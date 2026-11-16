#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 2 - 1 km monthly normals 1991-2020 by PRISM-type local weighted regression.

For each target (grid cell, or station for leave-one-out) and calendar month, the K nearest
stations with an accepted normal (step 1) enter a weighted regression
    precipitation : ln N = b0 + b1 z
    temperature   : N    = b0 + b1 z + b2 TPI20
whose station weights combine (PRISM-style) distance (adaptive bisquare), elevation
difference, topographic facet (8 orientations on a ~5 km smoothed DEM, plus flat) and,
for temperature, topographic position (valley floor vs slope/ridge, TPI over 41 km).
Slopes are bounded to physically plausible ranges. Residual corrections (full IDW, or a short-range Gaussian
taper of range RESID_L) are tested by LOO and adopted only if they improve it; for the v4 network both degrade
LOO at every range, so the normals are the regression alone (residuals are station-local, not spatially coherent
at the network spacing). The same code predicts each station's normal
WITHOUT that station, which step 3 needs for an honest LOO of the anomaly method.

Usage:  python s2_normals_grid.py --variable precip|temp [--tune]
Output: <OUT_ROOT>/normals/normals_<var>.nc, station_normals_loo_<var>.csv, tuning_<var>.csv
"""
import argparse, itertools, json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
import po_config as C, po_common as P

DEFAULTS = dict(K=40, DZ0=300.0, FACET_W=1.0, TPI0=1e9, USE_TPI=1, RESID="none", RESID_K=7, RESID_L=0.0, EXTRAP_DZ=1e9, EXTRAP_LAPSE=0.0)
EXTRAP_LAPSE_GRID = [0.0, -0.0055, -0.0065, -0.0075, -0.0085, -0.0098]   # degC/m used above the highest of the K stations (0 = local lapse)
EXTRAP_GRID = [1e9, 1000.0, 600.0, 400.0, 250.0, 100.0, 0.0]   # m above the highest of the K stations
RESID_L_GRID = [500.0, 1000.0, 2000.0, 3000.0, 5000.0, 8000.0]   # m, taper ranges tried in stage-2 tuning
SLOPE_BOUNDS = {"precip": (-0.3e-3, 0.8e-3),      # d lnP / dz per metre  (-30 % .. +120 % per km, e^0.8)
                "temp": (-10.0e-3, 3.0e-3),        # degC per metre
                "wetfreq": (-1.0e-3, 2.0e-3)}      # d logit(wet-day frequency) / dz per metre (two-part precipitation)


def design(variable, z, tpi):
    return np.stack([np.ones_like(z), z] + ([tpi] if variable == "temp" else []), -1)


def _tpi_scale(prm):
    return 1.0 if int(prm.get("USE_TPI", 1)) else 0.0


def predict(variable, prm, txy, tz, ttpi, tfacet, tself, sxy, sz, stpi, sfacet, y, resid=None, fixed_slope=None):
    """Batched local weighted regression at targets. tself: index of the target station or -1.
    fixed_slope (precip): per-target elevation slope imposed instead of the monthly fit."""
    K = int(prm["K"]); tree = cKDTree(sxy); kq = min(K + 1, len(sxy))
    d, il = tree.query(txy, k=kq)
    order = np.argsort(il == tself[:, None], axis=1, kind="stable")[:, :K]          # drop self (LOO)
    il = np.take_along_axis(il, order, 1); d = np.take_along_axis(d, order, 1)
    h = d[:, -1:] * 1.01 + 1.0
    w = (1 - (d / h) ** 2) ** 2
    w *= 1.0 / (1.0 + (np.abs(sz[il] - tz[:, None]) / prm["DZ0"]) ** 2)
    same = (sfacet[il] == tfacet[:, None]) | (sfacet[il] == 8) | (tfacet[:, None] == 8)
    w *= np.where(same, 1.0, prm["FACET_W"])
    if variable == "temp":
        w *= 1.0 / (1.0 + (np.abs(stpi[il] - ttpi[:, None]) / prm["TPI0"]) ** 2)
    X = design(variable, sz[il], stpi[il] * _tpi_scale(prm)); Y = y[il]                           # (M,K,p), (M,K)
    # centre predictors on the target so b0 is the estimate at the target. Precipitation: the target elevation used for
    # the elevation term is capped at EXTRAP_DZ above the highest of the K stations (no unbounded vertical
    # extrapolation of a locally fitted ln-slope into terrain above the local network)
    tz_eff = np.minimum(tz, sz[il].max(1) + prm["EXTRAP_DZ"]) if variable in ("precip", "wetfreq") else tz
    Xc = X.copy(); Xc[..., 1] -= tz_eff[:, None]
    if variable == "temp": Xc[..., 2] -= ttpi[:, None] * _tpi_scale(prm)
    XtW = Xc.transpose(0, 2, 1) * w[:, None, :]
    A = XtW @ Xc; b = (XtW @ Y[..., None])[..., 0]
    p = A.shape[-1]; lam = np.zeros(p); lam[1:] = 1e-9 * np.maximum(A[:, 1:, 1:].diagonal(axis1=1, axis2=2).mean(0), 1e-12)
    beta = np.linalg.solve(A + np.diag(lam)[None], b[..., None])[..., 0]
    lo, hi = SLOPE_BOUNDS[variable]
    slope = np.clip(beta[:, 1], lo, hi) if fixed_slope is None else fixed_slope
    if variable == "temp":
        # re-solve intercept and TPI coefficient with the bounded lapse held fixed
        Yr = Y - slope[:, None] * Xc[..., 1]; X2 = Xc[..., [0, 2]]; XtW2 = X2.transpose(0, 2, 1) * w[:, None, :]
        b2 = np.linalg.solve(XtW2 @ X2 + np.diag([0, lam[2]])[None], (XtW2 @ Yr[..., None]))[..., 0]
        est = b2[:, 0]
        if prm.get("EXTRAP_LAPSE", 0.0) < 0:
            # above the highest of the K stations the local lapse (often flattened by valley inversions) is replaced by
            # EXTRAP_LAPSE (degC/m) for the extrapolated part; 0 = keep the local lapse
            zt = sz[il].max(1); up = np.maximum(tz - zt, 0.0)
            est = est - slope * up + np.minimum(slope, prm["EXTRAP_LAPSE"]) * up
    else:
        est = ((Y - slope[:, None] * Xc[..., 1]) * w).sum(1) / w.sum(1)
    if resid is not None and prm["RESID"] == "idw":
        k = int(prm["RESID_K"]); dr = d[:, :k]; ir = il[:, :k]
        wr = 1.0 / np.maximum(dr, 100.0) ** 2; est = est + (wr * resid[ir]).sum(1) / wr.sum(1)
    elif resid is not None and prm["RESID"] == "taper" and prm["RESID_L"] > 0:
        # short-range residual correction: Gaussian kernel of range RESID_L; normalised only when the summed
        # weight exceeds 1, so the correction honours stations locally and fades to the regression far away
        k = int(prm["RESID_K"]); dr = d[:, :k]; ir = il[:, :k]
        wr = np.exp(-(dr / prm["RESID_L"]) ** 2); est = est + (wr * resid[ir]).sum(1) / np.maximum(wr.sum(1), 1.0)
    return est, slope


def annual_slopes(prm, grid, st, obs):
    """Precipitation: local ln-elevation slope fitted for each month on the stations of that month,
    then averaged over the 12 months, at every grid cell and (leave-one-out) at every station.
    One slope for all months removes the month-to-month artefacts of independent monthly fits on
    small, changing station sets, while keeping every station (also those with incomplete months)."""
    N = len(grid["xy"]); gs = np.zeros(N); ss = np.zeros(len(obs))
    for m in range(12):
        sidx = np.flatnonzero(np.isfinite(obs[:, m])); pos = np.full(len(obs), -1); pos[sidx] = np.arange(len(sidx))
        S = dict(sxy=st["xy"][sidx], sz=st["elev"][sidx], stpi=st["tpi20"][sidx], sfacet=st["facet"][sidx], y=np.log(np.maximum(obs[sidx, m], 0.1)))
        gs += predict("precip", prm, grid["xy"], grid["elev"], grid["tpi20"].astype(float), grid["facet"], np.full(N, -1), **S)[1] / 12
        ss += predict("precip", prm, st["xy"], st["elev"], st["tpi20"], st["facet"], pos, **S)[1] / 12
    return gs, ss


def run_month(variable, prm, grid, st, y, ok, tune_only=False, slopes=None):
    """Returns grid estimates (N,), station LOO estimates (S,), grid slope (N,).
    slopes: (grid_slope, station_slope) from annual_slopes(), or None for monthly fits."""
    S = len(st["ids"]); sidx = np.flatnonzero(ok)
    sxy, sz, stpi, sfc, ys = st["xy"][sidx], st["elev"][sidx], st["tpi20"][sidx], st["facet"][sidx], y[sidx]
    # station LOO residuals of the regression part (used both for the correction and as LOO estimate base)
    gsl, ssl = slopes if slopes is not None else (None, None)
    loo_reg, _ = predict(variable, dict(prm, RESID="none"), sxy, sz, stpi, sfc, np.arange(len(sidx)), sxy, sz, stpi, sfc, ys,
                         fixed_slope=None if ssl is None else ssl[sidx])
    resid = ys - loo_reg
    # station LOO estimates: regression without the station + residual correction from other stations
    loo_est = np.full(S, np.nan)
    allxy = st["xy"]; self_all = np.full(S, -1); pos = np.full(S, -1); pos[sidx] = np.arange(len(sidx)); self_all = pos
    e_loo, _ = predict(variable, prm, allxy, st["elev"], st["tpi20"], st["facet"], self_all, sxy, sz, stpi, sfc, ys, resid, fixed_slope=ssl)
    loo_est[:] = e_loo
    if tune_only: return None, loo_est, None
    g_est, g_slope = predict(variable, prm, grid["xy"], grid["elev"], grid["tpi20"].astype(float), grid["facet"], np.full(len(grid["xy"]), -1),
                             sxy, sz, stpi, sfc, ys, resid, fixed_slope=gsl)
    return g_est, loo_est, g_slope


def to_space(variable, v): return np.log(np.maximum(v, 0.1)) if variable == "precip" else v
def from_space(variable, v): return np.exp(v) if variable == "precip" else v


def score(variable, obs, est):
    ok = np.isfinite(obs) & np.isfinite(est)
    if variable == "precip":
        ann_o = np.nansum(np.where(ok, obs, 0), 1); ann_e = np.nansum(np.where(ok, est, 0), 1); full = ok.all(1)
        e = (ann_e - ann_o)[full]; return dict(n=int(full.sum()), annual_mae_mm=np.mean(np.abs(e)), annual_bias_mm=np.mean(e),
                                             annual_rel_mae_pct=100 * np.mean(np.abs(e) / ann_o[full]), monthly_rel_mae_pct=100 * np.nanmean(np.abs(est - obs)[ok] / obs[ok]))
    e = (est - obs)[ok]; return dict(n=int(ok.sum()), monthly_mae_c=np.mean(np.abs(e)), monthly_rmse_c=np.sqrt(np.mean(e * e)), monthly_bias_c=np.mean(e))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--variable", choices=["precip", "temp"], required=True)
    ap.add_argument("--tune", action="store_true", help="grid-search the weighting parameters by LOO before producing the grid")
    ap.add_argument("--params", default=None, help="JSON file with fixed parameters (e.g. params_<var>.json of a tuned run)")
    ap.add_argument("--monthly-slopes", action="store_true", help="precip: use each month's own elevation slope (v1 behaviour) "
                    "instead of the 12-month mean slope")
    a = ap.parse_args(); v = a.variable; t0 = time.time()
    grid = P.load_grid(); st = P.attach_station_predictors(P.load_station_data(v), grid)
    nd = pd.read_csv(C.OUT_ROOT / "normals" / f"station_normals_{v}.csv", dtype={"station": str})
    pv = nd.pivot(index="station", columns="month", values="normal").reindex(st["ids"])
    obs = pv.to_numpy(float, copy=True)                                                    # (S,12)
    prm = dict(DEFAULTS)
    if a.params: prm.update(json.load(open(a.params)))
    annual = (v == "precip") and not a.monthly_slopes
    def mask_slopes(p_, m):
        if not annual: return np.isfinite(obs[:, m]), None
        key = tuple(sorted(p_.items()))
        if key not in cache: cache[key] = annual_slopes(p_, grid, st, obs)
        return np.isfinite(obs[:, m]), cache[key]
    cache = {}
    if a.tune:
        space = dict(K=[10, 15, 25, 40, 60], DZ0=[150.0, 300.0, 500.0, 1e9], FACET_W=[0.5, 1.0], TPI0=[100.0, 1e9] if v == "temp" else [1e9],
                     USE_TPI=[0, 1] if v == "temp" else [0], RESID=["none", "idw"])
        rows = []
        for combo in itertools.product(*space.values()):
            p_ = dict(DEFAULTS, **dict(zip(space, combo))); est = np.full_like(obs, np.nan)
            for m in range(12):
                ok, sl_ = mask_slopes(p_, m); _, le, _ = run_month(v, p_, grid, st, to_space(v, obs[:, m]), ok, tune_only=True, slopes=sl_)
                est[:, m] = from_space(v, le)
            rows.append(dict(p_, **score(v, obs, est)))
        tdf = pd.DataFrame(rows); key = "annual_mae_mm" if v == "precip" else "monthly_mae_c"
        best = tdf.sort_values(key).iloc[0]; prm.update({k: (int(best[k]) if k in ("K", "RESID_K", "USE_TPI") else best[k]) for k in DEFAULTS})
        # stage 2: short-range residual correction with the stage-1 regression parameters held fixed, adopted only
        # if it improves LOO (v4, 2026-10-10: LOO degraded monotonically with range for both variables -> 'none').
        rows2 = []
        for L in RESID_L_GRID:
            p_ = dict(prm, RESID="taper", RESID_L=L); est = np.full_like(obs, np.nan)
            for m in range(12):
                ok, sl_ = mask_slopes(p_, m); _, le, _ = run_month(v, p_, grid, st, to_space(v, obs[:, m]), ok, tune_only=True, slopes=sl_)
                est[:, m] = from_space(v, le)
            rows2.append(dict(p_, **score(v, obs, est)))
        t2 = pd.DataFrame(rows2).sort_values(key); b2 = t2.iloc[0]
        if b2[key] < best[key]: prm.update(RESID="taper", RESID_L=float(b2["RESID_L"]))   # adopted only if it improves LOO
        # stage 3 (precipitation): cap on vertical extrapolation above the highest of the K stations, LOO-best of
        # EXTRAP_GRID (1e9 = no cap); evaluated with the stage-1/2 parameters held fixed
        t3 = pd.DataFrame()
        if v == "precip":
            hasn = np.isfinite(obs).all(1); ii = np.flatnonzero(hasn); tr_ = cKDTree(st["xy"][ii]); kk = int(prm["K"])
            _, nb = tr_.query(st["xy"], k=kk + 1); nb = ii[nb]
            nb = np.where(nb == np.arange(len(nb))[:, None], -1, nb)          # drop self
            zmax_nb = np.array([st["elev"][r[r >= 0][:kk]].max() for r in nb])
            rows3 = []
            for dz in EXTRAP_GRID:
                p_ = dict(prm, EXTRAP_DZ=dz); est = np.full_like(obs, np.nan)
                for m in range(12):
                    ok, sl_ = mask_slopes(p_, m); _, le, _ = run_month(v, p_, grid, st, to_space(v, obs[:, m]), ok, tune_only=True, slopes=sl_)
                    est[:, m] = from_space(v, le)
                sc_ = score(v, obs, est); ext = st["elev"] - zmax_nb
                hi = (ext > 0) & hasn & np.isfinite(est).all(1); e_hi = np.nanmean(np.abs(np.nansum(est[hi], 1) - np.nansum(obs[hi], 1)) / np.nansum(obs[hi], 1)) * 100
                rows3.append(dict(p_, **sc_, above_network_n=int(hi.sum()), above_network_annual_rel_mae_pct=e_hi))
            t3 = pd.DataFrame(rows3).sort_values(key); prm.update(EXTRAP_DZ=float(t3.iloc[0]["EXTRAP_DZ"]))
            print("stage 3 (extrapolation cap):"); print(t3[["EXTRAP_DZ", key, "above_network_n", "above_network_annual_rel_mae_pct"]].to_string(index=False))
        if v == "temp":
            rows3 = []; hasn = np.isfinite(obs).all(1); ii = np.flatnonzero(hasn); kk = int(prm["K"])
            _, nb = cKDTree(st["xy"][ii]).query(st["xy"], k=kk + 1); nb = ii[nb]; nb = np.where(nb == np.arange(len(nb))[:, None], -1, nb)
            zmax_nb = np.array([st["elev"][r[r >= 0][:kk]].max() for r in nb]); hi = (st["elev"] > zmax_nb) & hasn
            for lp in EXTRAP_LAPSE_GRID:
                p_ = dict(prm, EXTRAP_LAPSE=lp); est = np.full_like(obs, np.nan)
                for m in range(12):
                    ok, sl_ = mask_slopes(p_, m); _, le, _ = run_month(v, p_, grid, st, to_space(v, obs[:, m]), ok, tune_only=True, slopes=sl_)
                    est[:, m] = from_space(v, le)
                e = est - obs
                rows3.append(dict(p_, **score(v, obs, est), above_network_n=int(hi.sum()),
                                  above_network_bias_c=float(np.nanmean(e[hi])), above_network_mae_c=float(np.nanmean(np.abs(e[hi])))))
            t3 = pd.DataFrame(rows3).sort_values(key); prm.update(EXTRAP_LAPSE=float(t3.iloc[0]["EXTRAP_LAPSE"]))
            print("stage 3 (extrapolation lapse):"); print(t3[["EXTRAP_LAPSE", key, "above_network_n", "above_network_bias_c", "above_network_mae_c"]].to_string(index=False))
        tdf = pd.concat([tdf, t2, t3], ignore_index=True).sort_values(key); tdf.to_csv(C.OUT_ROOT / "normals" / f"tuning_{v}.csv", index=False)
        print("stage 2 (taper range):"); print(t2[["RESID_L", key]].to_string(index=False))
        print("tuning (best 5):"); print(tdf.head(5).to_string(index=False))
    N = len(grid["xy"]); G = np.full((12, N), np.nan); SL = np.full((12, N), np.nan); LOO = np.full_like(obs, np.nan)
    for m in range(12):
        ok, sl_ = mask_slopes(prm, m)
        g, le, sl = run_month(v, prm, grid, st, to_space(v, obs[:, m]), ok, slopes=sl_)
        G[m], SL[m], LOO[:, m] = from_space(v, g), sl, from_space(v, le)
    sc = score(v, obs, LOO); print("LOO skill of the normals:", {k: round(float(x), 3) for k, x in sc.items()})
    (C.OUT_ROOT / "normals" / f"params_{v}.json").write_text(json.dumps({k: (float(x) if not isinstance(x, str) else x) for k, x in prm.items()}, indent=1))
    import xarray as xr
    ny, nx = int(grid["ny"]), int(grid["nx"]); iy, ix = grid["iy"], grid["ix"]
    def to2d(a):
        out = np.full((12, ny, nx), np.nan, np.float32); out[:, iy, ix] = a; return out
    units = "mm month-1" if v == "precip" else "degC"
    ds = xr.Dataset(dict(normal=(("month", "y", "x"), to2d(G), dict(units=units, long_name=f"1991-2020 monthly normal of {v}")),
                         elevation_slope=(("month", "y", "x"), to2d(SL * 1000), dict(units=("ln(mm) km-1" if v == "precip" else "degC km-1"),
                                          long_name="local regression slope on elevation"))),
                    coords=dict(month=np.arange(1, 13), y=grid["y"], x=grid["x"]),
                    attrs=dict(method="PRISM-type local weighted regression + IDW residual correction", crs=C.METRIC_CRS,
                               parameters=str(prm), loo_skill=str(sc),
                               slope_mode=("12-month mean of the monthly local ln-elevation slopes" if annual else "monthly fit"), period=f"{C.PERIOD[0]}..{C.PERIOD[1]}"))
    enc = {k: dict(zlib=True, complevel=4) for k in ds.data_vars}
    od = C.OUT_ROOT / "normals"; ds.to_netcdf(od / f"normals_{v}.nc", encoding=enc)
    out = pd.DataFrame(dict(station=np.repeat(st["ids"], 12), month=np.tile(np.arange(1, 13), len(st["ids"])),
                            normal_station=obs.ravel(), normal_grid_loo=LOO.ravel()))
    out.to_csv(od / f"station_normals_loo_{v}.csv", index=False)
    print(f"done in {time.time()-t0:.0f} s -> {od}")


if __name__ == "__main__":
    main()
