#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 1w - 1991-2020 monthly normals of WET-DAY FREQUENCY (two-part precipitation method).

The monthly precipitation normal is M = n * i: wet-day frequency times mean wet-day intensity. The two-part
method keeps M (steps 1-2) and adds a separate normal of the wet-day frequency p = n / days, so that daily
occurrence and daily amounts can be interpolated separately (see s3_run.py, PO_PRECIP_METHOD=twopart).

Station value: for each station and calendar month, p = wet days (>= OCC_THRESHOLD mm) / valid days, pooled
over the months of 1991-2020 with >= MONTH_COMPLETENESS valid days; accepted with >= MIN_OBS_YEARS such months.
Grid: the same PRISM-type local weighted regression as the precipitation normals (s2_normals_grid.predict),
applied to logit(p) on elevation, with distance and elevation-difference weights; the elevation slope is the
12-month mean of the monthly local slopes; an optional cap on vertical extrapolation above the highest of the K
stations. K, DZ0 and the cap are selected by leave-one-out of the station normals (MAE of the annual number
of wet days). The same code predicts every station WITHOUT that station (for an honest LOO in step 3).

Usage:  python s1w_wetfreq_normals.py [--tune] [--params file.json] [--exclude-stations file]
Output: <OUT_ROOT>/normals/wetfreq_precip.nc, station_wetfreq_loo_precip.csv, params_wetfreq_precip.json,
        tuning_wetfreq_precip.csv (with --tune)
"""
import argparse, itertools, json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import po_config as C, po_common as P
import s2_normals_grid as S2

DEFAULTS = dict(K=25, DZ0=300.0, FACET_W=1.0, TPI0=1e9, USE_TPI=0, RESID="none", RESID_K=7, RESID_L=0.0, EXTRAP_DZ=1e9, EXTRAP_LAPSE=0.0)
SPACE = dict(K=[10, 15, 25, 40], DZ0=[150.0, 300.0, 500.0, 1e9], EXTRAP_DZ=[0.0, 500.0, 1e9])
PMIN, PMAX = 0.005, 0.995
DIM = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])


def logit(p): p = np.clip(p, PMIN, PMAX); return np.log(p / (1 - p))
def expit(x): return 1.0 / (1.0 + np.exp(-x))


def station_frequency(st):
    """(S,12) wet-day frequency, (S,12) number of accepted months."""
    d = st["dates"]; sel = (d >= pd.Timestamp(C.PERIOD[0])) & (d <= pd.Timestamp(C.PERIOD[1]))
    v = st["values"][sel]; d = d[sel]; per = pd.PeriodIndex(d, freq="M")
    codes = np.asarray((per.year - per.year.min()) * 12 + per.month - 1); edges = np.r_[0, np.flatnonzero(np.diff(codes)) + 1]
    fin = np.isfinite(v); wet = fin & (v >= C.OCC_THRESHOLD)
    nval = np.add.reduceat(fin, edges, 0); nwet = np.add.reduceat(wet, edges, 0); ndays = np.diff(np.r_[edges, len(v)])[:, None]
    acc = nval >= C.MONTH_COMPLETENESS * ndays; cal = (np.asarray(per.month)[edges] - 1)
    S = v.shape[1]; p = np.full((S, 12), np.nan); nm = np.zeros((S, 12), int)
    for m in range(12):
        r = cal == m; a = acc[r]; nm[:, m] = a.sum(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            p[:, m] = np.where(a, nwet[r], 0).sum(0) / np.where(a, nval[r], 0).sum(0)
    p[nm < C.MIN_OBS_YEARS] = np.nan
    return p, nm


def slopes(prm, grid, st, y, stations_only=False):
    """12-month mean of the monthly local logit-elevation slopes at grid cells and (LOO) at stations."""
    N = len(grid["xy"]); gs = np.zeros(N); ss = np.zeros(len(y))
    for m in range(12):
        sidx = np.flatnonzero(np.isfinite(y[:, m])); pos = np.full(len(y), -1); pos[sidx] = np.arange(len(sidx))
        A = dict(sxy=st["xy"][sidx], sz=st["elev"][sidx], stpi=st["tpi20"][sidx], sfacet=st["facet"][sidx], y=y[sidx, m])
        ss += S2.predict("wetfreq", prm, st["xy"], st["elev"], st["tpi20"], st["facet"], pos, **A)[1] / 12
        if not stations_only:
            gs += S2.predict("wetfreq", prm, grid["xy"], grid["elev"], grid["tpi20"].astype(float), grid["facet"], np.full(N, -1), **A)[1] / 12
    return gs, ss


def loo(prm, grid, st, y, ss):
    est = np.full_like(y, np.nan)
    for m in range(12):
        ok = np.isfinite(y[:, m])
        _, le, _ = S2.run_month("wetfreq", prm, grid, st, y[:, m], ok, tune_only=True, slopes=(None, ss))
        est[:, m] = le
    return est


def score(p_obs, p_est):
    ok = np.isfinite(p_obs).all(1) & np.isfinite(p_est).all(1)
    no = (p_obs[ok] * DIM).sum(1); ne = (p_est[ok] * DIM).sum(1)
    return dict(n=int(ok.sum()), annual_wetdays_mae=float(np.mean(np.abs(ne - no))), annual_wetdays_bias=float(np.mean(ne - no)),
                annual_wetdays_rel_mae_pct=float(100 * np.mean(np.abs(ne - no) / no)))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tune", action="store_true"); ap.add_argument("--params", default=None)
    ap.add_argument("--exclude-stations", default=None)
    a = ap.parse_args(); t0 = time.time()
    excl = [l.strip() for l in open(a.exclude_stations)] if a.exclude_stations else None
    grid = P.load_grid(); st = P.attach_station_predictors(P.load_station_data("precip"), grid)
    p_st, nm = station_frequency(st)
    if excl:   # withheld stations keep their row (they are predicted, LOO-style) but do not inform any normal
        p_st[np.isin(st["ids"], excl)] = np.nan
    y = logit(p_st)
    od = C.OUT_ROOT / "normals"; od.mkdir(parents=True, exist_ok=True)
    prm = dict(DEFAULTS)
    pf = Path(a.params) if a.params else (C.DATA_ROOT / "normals" / "params_wetfreq_precip.json")
    if not a.tune and pf.exists(): prm.update(json.load(open(pf)))
    if a.tune:
        rows = []
        for combo in itertools.product(*SPACE.values()):
            p_ = dict(DEFAULTS, **dict(zip(SPACE, combo)))
            _, ss = slopes(p_, grid, st, y, stations_only=True)
            rows.append(dict(p_, **score(p_st, expit(loo(p_, grid, st, y, ss)))))
        T = pd.DataFrame(rows).sort_values("annual_wetdays_mae"); T.to_csv(od / "tuning_wetfreq_precip.csv", index=False)
        print(T[list(SPACE) + ["n", "annual_wetdays_mae", "annual_wetdays_bias", "annual_wetdays_rel_mae_pct"]].head(8).to_string(index=False))
        b = T.iloc[0]; prm.update(K=int(b.K), DZ0=float(b.DZ0), EXTRAP_DZ=float(b.EXTRAP_DZ))
    gs, ss = slopes(prm, grid, st, y)
    N = len(grid["xy"]); G = np.full((12, N), np.nan); L = expit(loo(prm, grid, st, y, ss))
    for m in range(12):
        ok = np.isfinite(y[:, m])
        g, _, _ = S2.run_month("wetfreq", prm, grid, st, y[:, m], ok, slopes=(gs, ss))
        G[m] = expit(g)
    sc = score(p_st, L); print("LOO skill of the wet-day frequency normals:", {k: round(float(x), 3) for k, x in sc.items()})
    ok = np.isfinite(p_st).all(1)
    b_obs = np.polyfit(st["elev"][ok] / 1e3, np.log((p_st[ok] * DIM).sum(1)), 1)[0]
    b_loo = np.polyfit(st["elev"][ok] / 1e3, np.log((L[ok] * DIM).sum(1)), 1)[0]
    print(f"district-wide gradient of annual wet days: stations {100*(np.exp(b_obs)-1):+.1f} %/km, LOO normals {100*(np.exp(b_loo)-1):+.1f} %/km")
    (od / "params_wetfreq_precip.json").write_text(json.dumps({k: (float(x) if not isinstance(x, str) else x) for k, x in prm.items()}, indent=1))
    import xarray as xr
    ny, nx = int(grid["ny"]), int(grid["nx"]); iy, ix = grid["iy"], grid["ix"]
    def to2d(arr):
        out = np.full((12, ny, nx), np.nan, np.float32); out[:, iy, ix] = arr; return out
    ds = xr.Dataset(dict(wetfreq=(("month", "y", "x"), to2d(G), dict(units="1", long_name=f"1991-2020 wet-day frequency (days >= {C.OCC_THRESHOLD} mm / days)")),
                         elevation_slope=(("month", "y", "x"), to2d(np.tile(gs * 1000, (12, 1))), dict(units="logit km-1", long_name="local slope of logit(frequency) on elevation"))),
                    coords=dict(month=np.arange(1, 13), y=grid["y"], x=grid["x"]),
                    attrs=dict(method="PRISM-type local weighted regression of logit(wet-day frequency) on elevation", crs=C.METRIC_CRS,
                               parameters=str(prm), loo_skill=str(sc), wet_threshold_mm=C.OCC_THRESHOLD, period=f"{C.PERIOD[0]}..{C.PERIOD[1]}"))
    ds.to_netcdf(od / "wetfreq_precip.nc", encoding={k: dict(zlib=True, complevel=4) for k in ds.data_vars})
    pd.DataFrame(dict(station=np.repeat(st["ids"], 12), month=np.tile(np.arange(1, 13), len(st["ids"])), wetfreq_station=p_st.ravel(),
                      wetfreq_grid_loo=L.ravel(), n_months=nm.ravel())).to_csv(od / "station_wetfreq_loo_precip.csv", index=False)
    print(f"stations with 12 frequency normals: {int(np.isfinite(p_st).all(1).sum())} of {len(p_st)}; done in {time.time()-t0:.0f} s -> {od}")


if __name__ == "__main__":
    main()
