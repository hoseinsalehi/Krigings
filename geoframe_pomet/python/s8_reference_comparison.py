#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Comparison of the v3 product with independent gridded reference datasets (reviewer point 3).

The 1 km product is aggregated to each reference grid (area mean of the 1 km cells whose centres fall
in a reference cell; >= 50 % of the expected cells required), so every comparison is made at the
reference's native resolution without re-interpolating the reference. For references finer than or
equal to 1 km (IH-GAR, 30") the comparison is cell-to-cell at the nearest reference cell.

Precipitation references          temperature references
  ARCIS31  (Pavan et al. 2019)      E-OBS v31.0e tg 0.1 deg (Cornes et al. 2018)
  IH-GAR   (Manara/Brunetti, UNIMI) EMO-1arcmin ta6 (Thiemig et al. 2022)
                                    ARCIS (Tmax+Tmin)/2 (Pavan et al. 2026)
  RAPD     (MeteoSwiss RapdD, successor of APGD, Isotta et al. 2014)
  EMO-1arcmin pr (Thiemig et al. 2022; day 06-06 UTC, labelled by its end -> shifted -1 day)
  E-OBS v31.0e rr 0.1 deg
Per reference cell, over the common days: means, bias, daily r, MAE, RMSE, KGE, seasonal biases,
precipitation wet-day frequency (>= 1 mm) and mean annual maximum (Rx1day), event totals/means
(Storm Alex 2-5 Oct 2020, summer 2003; temperature: summer 2003, 1-14 Feb 2012), and for the area-mean
series of each region the daily correlation at lags -1, 0, +1 day (day-definition check).

Usage:  python s8_reference_comparison.py --variable precip [--products ARCIS31,IH-GAR] [--run precip_anomaly]
Output: <OUT_ROOT>/reference_comparison/
          cells_<var>_<product>.parquet   per reference cell metrics (lon, lat, elevation, region)
          summary_<var>.csv               per product x region / elevation band: medians of cell metrics
                                          and area-mean statistics
          lags_<var>.csv                  area-mean daily r at lags -1/0/+1
          events_<var>.csv                event spatial r, bias, means
"""
import argparse, glob, os, sys, tempfile, time, zipfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import po_config as C, po_common as P

DL = Path(os.environ.get("PO_INPUT_ROOT", "data") + "")
PRODUCTS = {
    "precip": {
        "ARCIS31": dict(kind="decade", files=sorted(DL.glob("ARCIS31_GG_PR_*.nc*")), var="rr", crs="EPSG:4326", shift=[(1900, -1), (2016, 0)],   # end-day labels until 2015, start-day from 2016 (yearly lag test)
                        ref="Pavan et al. (2019), ARCIS"),
        "IH-GAR": dict(kind="yearly", pattern=str(DL / "precip_IH-GAR" / "{y}.nc"), var="total_precipitation", crs="EPSG:4326", shift=0,
                       ref="Manara, Brunetti et al., UNIMI doi:10.13130/RD_UNIMI/CSYYAP"),
        "RAPD": dict(kind="rapd_zip", files=sorted((DL / "transfer-01a08b46").glob("RapdD_*.zip")), var="RapdD", crs="EPSG:3035", shift=0,
                     ref="MeteoSwiss RapdD (successor of APGD, Isotta et al. 2014)"),
        "EMO": dict(kind="yearly", pattern=str(DL / "EMO_pr" / "EMO-1arcmin-pr_{y}_clipped_32632.nc"), var="pr", crs="EPSG:32632", shift=-1,
                    ref="EMO-1arcmin (Thiemig et al. 2022), reprojected to EPSG:32632 by the user"),
        "E-OBS": dict(kind="multi", files=sorted((DL / "E-OBS").glob("rr_ens_mean_0.1deg_reg_*_v31.0e*")), var="rr", crs="EPSG:4326", shift=0,
                      ref="E-OBS v31.0e ensemble mean 0.1 deg (Cornes et al. 2018)"),
    },
    "temp": {
        "E-OBS": dict(kind="multi", files=sorted((DL / "E-OBS").glob("tg_ens_mean_0.1deg_reg_*_v31.0e.nc")), var="tg", crs="EPSG:4326", shift=0,
                      ref="E-OBS v31.0e ensemble mean 0.1 deg (Cornes et al. 2018)"),
        "EMO": dict(kind="yearly", pattern=str(DL / "EMO_temp" / "EMO-temp_1arcmin_{y}.nc"), var="ta6", crs="EPSG:4326", shift=0,
                    ref="EMO-1arcmin daily mean of 6-hourly temperature (Thiemig et al. 2022)"),
        "ARCIS_T": dict(kind="decade", files=sorted(DL.glob("ARCIS1_GG_TX_*.nc*")), files2=sorted(DL.glob("ARCIS1_GG_TN_*.nc*")),
                        var="tx", var2="tn", crs="EPSG:4326", shift=0,
                        ref="ARCIS daily (Tmax+Tmin)/2 (Pavan et al. 2026); not a 24-h mean"),
    },
}
EVENTS = {"precip": {"storm_alex_2020": ("2020-10-02", "2020-10-05", "sum"), "summer_2003": ("2003-06-01", "2003-08-31", "sum")},
          "temp": {"summer_2003": ("2003-06-01", "2003-08-31", "mean"), "cold_feb_2012": ("2012-02-01", "2012-02-14", "mean")}}
SEASONS = {"DJF": (12, 1, 2), "MAM": (3, 4, 5), "JJA": (6, 7, 8), "SON": (9, 10, 11)}
BANDS = ((0, 500), (500, 1000), (1000, 1500), (1500, 2000), (2000, 9000))


# ----------------------------------------------------------------------------- reading
def _xy_names(ds):
    for xn, yn in (("lon", "lat"), ("longitude", "latitude"), ("x", "y"), ("X", "Y")):
        if xn in ds.coords and yn in ds.coords and ds[xn].ndim == 1: return xn, yn
    raise KeyError(f"no 1-D horizontal coordinates in {list(ds.coords)}")


def read_year(cfg, year, bbox):
    """Returns DatetimeIndex (shifted to the station day), array (T, ny, nx) float32, x axis, y axis."""
    import xarray as xr
    tmp = None
    if cfg["kind"] == "yearly":
        f = cfg["pattern"].format(y=year)
        if not os.path.exists(f): return None
        dss = [xr.open_dataset(f)]
    elif cfg["kind"] in ("decade", "multi"):
        if "files2" in cfg:                        # paired Tmax/Tmin files -> (Tmax+Tmin)/2 under cfg["var"]
            dss = []
            for f1, f2 in zip(cfg["files"], cfg["files2"]):
                d1, d2 = xr.open_dataset(f1), xr.open_dataset(f2)
                assert d1.sizes == d2.sizes and (d1["time"].values == d2["time"].values).all(), (f1, f2)
                d1[cfg["var"]] = 0.5 * (d1[cfg["var"]] + d2[cfg["var2"]]); dss.append(d1)
        else:
            dss = [xr.open_dataset(f) for f in cfg["files"]]
    elif cfg["kind"] == "rapd_zip":
        z = [f for f in cfg["files"] if any(int(a) <= year <= int(b) for a, b in [Path(f).stem.split("_")[-2:]])]
        if not z: return None
        tmp = tempfile.TemporaryDirectory(); zz = zipfile.ZipFile(z[0])
        names = [n for n in zz.namelist() if f"_{year}" in n and n.endswith(".nc")]
        if not names: return None
        for n in names: zz.extract(n, tmp.name)
        parts = [xr.open_dataset(str(Path(tmp.name) / n)).load() for n in sorted(names)]  # noqa
        dss = [xr.concat(parts, dim="time")]
    out = []
    for ds in dss:
        t0 = pd.DatetimeIndex(ds["time"].values).floor("D"); sh = cfg["shift"]
        if isinstance(sh, (list, tuple)):          # year-dependent labelling convention: [(from_year, shift_days), ...]
            k = np.zeros(len(t0), int)
            for fy, d in sorted(sh): k[t0.year >= fy] = d
            t = t0 + pd.to_timedelta(k, unit="D")
        else: t = t0 + pd.Timedelta(days=sh)
        sel = t.year == year
        if not sel.any(): continue
        xn, yn = _xy_names(ds); x = ds[xn].values; y = ds[yn].values
        ix = np.flatnonzero((x >= bbox[0]) & (x <= bbox[2])); iy = np.flatnonzero((y >= bbox[1]) & (y <= bbox[3]))
        a = ds[cfg["var"]].isel(time=np.flatnonzero(sel), **{xn: slice(ix.min(), ix.max() + 1), yn: slice(iy.min(), iy.max() + 1)})
        a = a.transpose("time", yn, xn).values.astype(np.float32)
        out.append((t[sel], a, x[ix.min():ix.max() + 1], y[iy.min():iy.max() + 1]))
    for ds in dss: ds.close()
    if tmp is not None: tmp.cleanup()
    if not out: return None
    t = pd.DatetimeIndex(np.concatenate([o[0] for o in out])); a = np.concatenate([o[1] for o in out]); x, y = out[0][2], out[0][3]
    a[(a < -900) | (a > 1e6)] = np.nan
    o = np.argsort(t, kind="stable"); t, a = t[o], a[o]; keep = ~t.duplicated(); return t[keep], a[keep], x, y


def nearest_index(axis, v):
    """Index of the nearest axis cell and validity (inside half a spacing) for an ascending/descending axis."""
    asc = axis[1] > axis[0]; ax = axis if asc else axis[::-1]
    j = np.clip(np.searchsorted(ax, v), 1, len(ax) - 1); jl = np.where(np.abs(v - ax[j - 1]) <= np.abs(v - ax[j]), j - 1, j)
    sp = np.abs(np.diff(ax)).mean(); ok = np.abs(v - ax[jl]) <= 0.5 * sp * 1.0001
    return (jl if asc else len(ax) - 1 - jl), ok, sp


# ----------------------------------------------------------------------------- mapping 1 km -> reference cells
def build_mapping(cfg, grid, x, y):
    from pyproj import Transformer
    from scipy import sparse
    tr = Transformer.from_crs(C.METRIC_CRS, cfg["crs"], always_xy=True)
    px, py = tr.transform(grid["xy"][:, 0], grid["xy"][:, 1])
    jx, okx, sx = nearest_index(x, px); jy, oky, sy = nearest_index(y, py); ok = okx & oky
    cid = jy * len(x) + jx; cid = np.where(ok, cid, -1)
    # expected number of 1 km cells per reference cell (area in km2)
    if cfg["crs"] == "EPSG:4326":
        lat = y[np.clip(jy, 0, len(y) - 1)]; area = sx * 111.32 * np.cos(np.radians(lat)) * sy * 110.57
    else: area = np.full(len(px), sx * sy / 1e6)
    u, inv, cnt = np.unique(cid[ok], return_inverse=True, return_counts=True)
    exp = pd.Series(area[ok]).groupby(inv).mean().to_numpy()
    good = cnt >= np.maximum(1, 0.5 * np.minimum(exp, 1e9)) if exp.min() >= 1 else np.ones(len(u), bool)
    keep_cells = u[good]; pos = -np.ones(len(u), int); pos[good] = np.arange(good.sum())
    m = pos[inv]; src = np.flatnonzero(ok)[m >= 0]; dst = m[m >= 0]
    w = 1.0 / np.bincount(dst)[dst]
    A = sparse.csr_matrix((w, (dst, src)), shape=(len(keep_cells), len(px)))
    iyc, ixc = np.divmod(keep_cells, len(x))
    meta = pd.DataFrame(dict(ref_iy=iyc, ref_ix=ixc, ref_x=x[ixc], ref_y=y[iyc], n_1km=np.bincount(dst),
                             elevation=(A @ grid["elev"]), lon=(A @ grid["lon"]), lat=(A @ grid["lat"])))
    meta["region"] = "Rest_PRD"
    for b, (lo0, lo1, la0, la1) in C.BOXES.items():
        meta.loc[(meta.lon >= lo0) & (meta.lon <= lo1) & (meta.lat >= la0) & (meta.lat <= la1), "region"] = b
    meta["elev_band"] = pd.cut(meta.elevation, [b[0] for b in BANDS] + [9000], right=False,
                               labels=[f"{a}-{b}" if b < 9000 else f">={a}" for a, b in BANDS]).astype(object).fillna("unknown").astype(str)
    return A, iyc, ixc, meta, dict(spacing_x=sx, spacing_y=sy, mean_1km_per_cell=float(np.bincount(dst).mean()))


# ----------------------------------------------------------------------------- one product
def compare_product(v, name, run, years):
    import xarray as xr
    cfg = PRODUCTS[v][name]; grid = P.load_grid(); vn = "precipitation" if v == "precip" else "temperature"
    gy, gx = grid["iy"], grid["ix"]
    from pyproj import Transformer
    tr = Transformer.from_crs("EPSG:4326", cfg["crs"], always_xy=True)
    bx, by = tr.transform([6.0, 6.0, 14.5, 14.5], [43.5, 47.2, 43.5, 47.2]); bbox = (min(bx), min(by), max(bx), max(by))
    A = None; acc = {}; series = []; ev = {k: None for k in EVENTS[v]}; t0 = time.time(); yrs_used = []
    for yr in years:
        f = C.OUT_ROOT / run / f"{v}_{yr}.nc"
        if not f.exists(): continue
        r = read_year(cfg, yr, bbox)
        if r is None: continue
        tr_, ra, x, y = r
        if A is None:
            A, iyc, ixc, meta, info = build_mapping(cfg, grid, x, y); M = len(meta); x0, y0 = x, y
        elif len(x) != len(x0) or len(y) != len(y0) or not (np.allclose(x, x0) and np.allclose(y, y0)):
            # the yearly files of some products (EMO clipped by year) have slightly different extents:
            # re-index this year onto the grid of the first year (missing cells -> NaN)
            jx, okx, _ = nearest_index(x, x0); jy, oky, _ = nearest_index(y, y0)
            al = np.full((ra.shape[0], len(y0), len(x0)), np.nan, np.float32)
            al[:, np.ix_(oky, okx)[0], np.ix_(oky, okx)[1]] = ra[:, jy[oky]][:, :, jx[okx]]; ra = al
        if not acc:
            for k in ("n", "sx", "sy", "sxx", "syy", "sxy", "sad", "wx", "wy", "rxx", "rxy", "nyr"): acc[k] = np.zeros(M)
            for s in SEASONS: acc[f"n_{s}"] = np.zeros(M); acc[f"sx_{s}"] = np.zeros(M); acc[f"sy_{s}"] = np.zeros(M)
        with xr.open_dataset(f) as ds:
            tv = pd.DatetimeIndex(ds["time"].values).floor("D"); xv = ds[vn].values[:, gy, gx].astype(np.float32)
        common = tv.intersection(tr_)
        if len(common) == 0: continue
        X = np.asarray((A @ xv[tv.get_indexer(common)].T).T)                                     # (T, M) v3 aggregated
        Y = ra[tr_.get_indexer(common)][:, iyc, ixc]                                             # (T, M) reference
        ok = np.isfinite(X) & np.isfinite(Y); Xz = np.where(ok, X, 0.0); Yz = np.where(ok, Y, 0.0)
        acc["n"] += ok.sum(0); acc["sx"] += Xz.sum(0); acc["sy"] += Yz.sum(0); acc["sxx"] += (Xz * Xz).sum(0)
        acc["syy"] += (Yz * Yz).sum(0); acc["sxy"] += (Xz * Yz).sum(0); acc["sad"] += np.abs(Xz - Yz).sum(0)
        mon = common.month
        for s, mm in SEASONS.items():
            sel = np.isin(mon, mm); acc[f"n_{s}"] += ok[sel].sum(0); acc[f"sx_{s}"] += Xz[sel].sum(0); acc[f"sy_{s}"] += Yz[sel].sum(0)
        if v == "precip":
            acc["wx"] += ((X >= 1) & ok).sum(0); acc["wy"] += ((Y >= 1) & ok).sum(0)
            full = ok.sum(0) >= 0.9 * len(pd.date_range(f"{yr}-01-01", f"{yr}-12-31"))          # per cell: >= 90 % of the year
            with np.errstate(all="ignore"):
                acc["rxx"] += np.where(full, np.nanmax(np.where(ok, X, np.nan), 0), 0); acc["rxy"] += np.where(full, np.nanmax(np.where(ok, Y, np.nan), 0), 0)
            acc["nyr"] += full
        allok = ok.all(1) if ok.shape[1] else ok.any(1)
        for rg in ["All_PRD"] + list(C.BOXES) + ["Rest_PRD"]:
            msk = np.ones(M, bool) if rg == "All_PRD" else (meta.region.to_numpy() == rg)
            if not msk.any(): continue
            with np.errstate(all="ignore"):
                series.append(pd.DataFrame(dict(date=common, region=rg, v3=np.nanmean(np.where(ok[:, msk], X[:, msk], np.nan), 1),
                                                ref=np.nanmean(np.where(ok[:, msk], Y[:, msk], np.nan), 1))))
        for k, (a_, b_, how) in EVENTS[v].items():
            sel = (common >= a_) & (common <= b_)
            if sel.any():
                agg = (lambda z: np.nansum(z, 0)) if how == "sum" else (lambda z: np.nanmean(z, 0))
                e = (agg(np.where(ok[sel], X[sel], np.nan)), agg(np.where(ok[sel], Y[sel], np.nan)), int(sel.sum()))
                ev[k] = e if ev[k] is None else (ev[k][0] + e[0], ev[k][1] + e[1], ev[k][2] + e[2]) if how == "sum" else e
        yrs_used.append(yr)
        print(f"  {v} {name} {yr}: {len(common)} days, {M} cells ({time.time()-t0:.0f} s)", flush=True)
    if A is None: return None
    n = np.maximum(acc["n"], 1); mx = acc["sx"] / n; my = acc["sy"] / n
    vx = acc["sxx"] / n - mx ** 2; vy = acc["syy"] / n - my ** 2; cxy = acc["sxy"] / n - mx * my
    with np.errstate(all="ignore"):
        r = cxy / np.sqrt(vx * vy); alpha = np.sqrt(vx / vy); beta = mx / my
        kge = 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)
        rmse = np.sqrt(np.maximum(acc["sxx"] / n - 2 * acc["sxy"] / n + acc["syy"] / n, 0))
    cell = meta.copy(); cell["n_days"] = acc["n"].astype(int)
    if v == "precip":
        cell["mean_annual_v3"] = mx * 365.25; cell["mean_annual_ref"] = my * 365.25
        cell["bias_pct"] = 100 * (mx - my) / my; cell["wetfreq_v3"] = acc["wx"] / n; cell["wetfreq_ref"] = acc["wy"] / n
        with np.errstate(all="ignore"):
            cell["rx1day_v3"] = acc["rxx"] / acc["nyr"]; cell["rx1day_ref"] = acc["rxy"] / acc["nyr"]
    else:
        cell["mean_v3"] = mx; cell["mean_ref"] = my; cell["bias"] = mx - my
    cell["r_daily"] = r; cell["mae_daily"] = acc["sad"] / n; cell["rmse_daily"] = rmse; cell["kge_daily"] = kge
    for s in SEASONS:
        ns = np.maximum(acc[f"n_{s}"], 1)
        cell[f"bias_{s}"] = (acc[f"sx_{s}"] - acc[f"sy_{s}"]) / ns if v == "temp" else 100 * (acc[f"sx_{s}"] - acc[f"sy_{s}"]) / np.where(acc[f"sy_{s}"] > 0, acc[f"sy_{s}"], np.nan)
    evrows = []
    for k, e in ev.items():
        if e is None: continue
        cell[f"{k}_v3"] = e[0]; cell[f"{k}_ref"] = e[1]
        okk = np.isfinite(e[0]) & np.isfinite(e[1])
        evrows.append(dict(variable=v, product=name, event=k, n_cells=int(okk.sum()), spatial_r=float(np.corrcoef(e[0][okk], e[1][okk])[0, 1]),
                           mean_v3=float(np.mean(e[0][okk])), mean_ref=float(np.mean(e[1][okk])),
                           max_v3=float(np.max(e[0][okk])), max_ref=float(np.max(e[1][okk]))))
    S = pd.concat(series, ignore_index=True); lagrows = []
    for rg, g in S.groupby("region"):
        g = g.set_index("date").sort_index()
        for lag in (-1, 0, 1):
            ref = g["ref"].shift(lag); okk = g["v3"].notna() & ref.notna()
            lagrows.append(dict(variable=v, product=name, region=rg, lag_days=lag, r=float(np.corrcoef(g["v3"][okk], ref[okk])[0, 1]), n=int(okk.sum())))
    cell = cell[cell.n_days >= 0.5 * cell.n_days.max()].reset_index(drop=True)      # drop cells without a common record
    info.update(product=name, variable=v, reference=cfg["ref"], years=f"{min(yrs_used)}-{max(yrs_used)}", n_cells=len(cell))
    return cell, pd.DataFrame(evrows), pd.DataFrame(lagrows), S, info


def summarise(v, name, cell, S, info):
    rows = []
    groups = [("all", "All_PRD", np.ones(len(cell), bool))] + [("region", rg, (cell.region == rg).to_numpy()) for rg in list(C.BOXES) + ["Rest_PRD"]] \
             + [("elev_band", b, (cell.elev_band == b).to_numpy()) for b in sorted(set(cell.elev_band.astype(str))) if b != "unknown"]
    for grp, lev, m in groups:
        c = cell[m]
        if len(c) == 0: continue
        d = dict(variable=v, product=name, group=grp, level=lev, n_cells=len(c), years=info["years"])
        if v == "precip":
            d.update(mean_annual_v3=c.mean_annual_v3.mean(), mean_annual_ref=c.mean_annual_ref.mean(),
                     bias_pct_area=100 * (c.mean_annual_v3.mean() - c.mean_annual_ref.mean()) / c.mean_annual_ref.mean(),
                     bias_pct_cell_median=c.bias_pct.median(), wetfreq_v3=c.wetfreq_v3.mean(), wetfreq_ref=c.wetfreq_ref.mean(),
                     rx1day_v3=c.rx1day_v3.mean(), rx1day_ref=c.rx1day_ref.mean(),
                     spatial_r_climatology=np.corrcoef(c.mean_annual_v3, c.mean_annual_ref)[0, 1] if len(c) > 2 else np.nan)
        else:
            d.update(mean_v3=c.mean_v3.mean(), mean_ref=c.mean_ref.mean(), bias_area=c.mean_v3.mean() - c.mean_ref.mean(),
                     bias_cell_median=c.bias.median(), bias_cell_p10=c.bias.quantile(.1), bias_cell_p90=c.bias.quantile(.9),
                     spatial_r_climatology=np.corrcoef(c.mean_v3, c.mean_ref)[0, 1] if len(c) > 2 else np.nan)
        for k in ("r_daily", "mae_daily", "rmse_daily", "kge_daily"): d[f"{k}_cell_median"] = c[k].median()
        for s in SEASONS: d[f"bias_{s}_cell_median"] = c[f"bias_{s}"].median()
        rows.append(d)
    return rows


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--variable", choices=["precip", "temp"], required=True)
    ap.add_argument("--products", default=None, help="comma list (default: all available)"); ap.add_argument("--run", default=None)
    ap.add_argument("--years", default="1991-2020"); ap.add_argument("--workers", type=int, default=1,
                    help="parallel products (fork); 1 = sequential (macOS: forking after netCDF/xarray import can crash)")
    a = ap.parse_args(); v = a.variable; run = a.run or f"{v}_anomaly"
    y0, y1 = (int(s) for s in a.years.split("-")) if "-" in a.years else (int(a.years),) * 2
    names = a.products.split(",") if a.products else list(PRODUCTS[v])
    od = C.OUT_ROOT / "reference_comparison"; od.mkdir(parents=True, exist_ok=True)
    import multiprocessing as mp
    job = lambda nm: compare_product(v, nm, run, range(y0, y1 + 1))
    global _JOB; _JOB = job
    if a.workers <= 1: res = [_run(nm) for nm in names]
    else:
        with mp.get_context("fork").Pool(min(a.workers, len(names))) as pool: res = pool.map(_run, names)
    summ, evs, lags, infos = [], [], [], []
    for nm, r in zip(names, res):
        if r is None: print(f"{nm}: no overlapping data"); continue
        cell, ev, lag, S, info = r
        cell.to_parquet(od / f"cells_{v}_{nm}.parquet", index=False); S.to_parquet(od / f"series_{v}_{nm}.parquet", index=False)
        summ += summarise(v, nm, cell, S, info); evs.append(ev); lags.append(lag); infos.append(info)
    pd.DataFrame(summ).to_csv(od / f"summary_{v}.csv", index=False); pd.concat(evs).to_csv(od / f"events_{v}.csv", index=False)
    pd.concat(lags).to_csv(od / f"lags_{v}.csv", index=False); pd.DataFrame(infos).to_csv(od / f"products_{v}.csv", index=False)
    sm = pd.DataFrame(summ); key = ["bias_pct_area", "spatial_r_climatology", "r_daily_cell_median", "kge_daily_cell_median"] if v == "precip" \
        else ["bias_area", "spatial_r_climatology", "r_daily_cell_median", "mae_daily_cell_median"]
    print(sm[sm.group.isin(["all", "region"])][["product", "level", "n_cells", "years"] + key].round(3).to_string(index=False))
    L = pd.concat(lags); print(L[L.region == "All_PRD"].pivot_table(index="product", columns="lag_days", values="r").round(3).to_string())
    print(pd.concat(evs).round(2).to_string(index=False)); print(f"-> {od}")


def _run(nm):
    return _JOB(nm)


if __name__ == "__main__":
    main()
