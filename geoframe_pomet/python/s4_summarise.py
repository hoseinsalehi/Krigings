#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 4 - reviewer-ready summaries of one or more step-3 runs.

  python s4_summarise.py --runs precip_direct precip_anomaly [--climatology]

Writes to <OUT_ROOT>/summary/:
  loo_metrics_<run>.csv          overall, elevation band, season, decade, diagnostic box, intensity class
  loo_comparison_<var>.csv        the same key metrics side by side for all runs of a variable
  support_by_period_<run>.csv     distributions (p10/p50/p90) of D1, D7, D50, R_pos_half, R_abs90, N_eff per decade
  cell_climatology_<run>.nc       (with --climatology) 30-yr per-cell mean, mean kriging std, mean D7 by decade,
                                  % positive-lapse days by season (temperature), % oscillatory days
"""
import argparse, glob, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import po_config as C, po_common as P


def kge(p, o):
    if len(o) < 3 or o.std() == 0: return np.nan
    r = np.corrcoef(p, o)[0, 1]; return 1 - np.sqrt((r - 1) ** 2 + (p.std() / o.std() - 1) ** 2 + (p.mean() / o.mean() - 1) ** 2)


def metrics(g, variable):
    ok = np.isfinite(g.pred) & np.isfinite(g.obs); p = g.pred[ok].to_numpy(float); o = g.obs[ok].to_numpy(float); e = p - o
    d = dict(n=int(ok.sum()), bias=e.mean(), mae=np.abs(e).mean(), rmse=np.sqrt((e * e).mean()),
             r=np.corrcoef(p, o)[0, 1] if len(o) > 2 else np.nan, kge=kge(p, o))
    s = g.kriging_std[ok].to_numpy(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        pos = s > 0                                  # estimates with zero kriging std (all-dry shortcut) are excluded
        d["cover95_pct"] = 100 * np.mean(np.abs(e[pos]) <= 1.96 * s[pos]) if np.any(pos) else np.nan
    if variable == "precip":
        w = o >= C.WET_THRESHOLD
        d.update(wet_n=int(w.sum()), wet_mae=np.abs(e[w]).mean() if w.any() else np.nan,
                 wet_pbias_pct=100 * e[w].sum() / o[w].sum() if w.any() else np.nan,
                 pod_wet=np.mean(p[w] >= C.WET_THRESHOLD) if w.any() else np.nan,
                 far_wet=np.mean(o[p >= C.WET_THRESHOLD] < C.WET_THRESHOLD) if (p >= C.WET_THRESHOLD).any() else np.nan)
    return pd.Series(d)


def load_loo(run_dir, variable, st):
    files = sorted(glob.glob(str(run_dir / f"loo_{variable}_*.parquet")))
    if not files: return None
    d = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    meta = pd.DataFrame(dict(station=st["ids"], elev=st["elev"], lon=np.nan, lat=np.nan))
    from pyproj import Transformer
    meta["lon"], meta["lat"] = Transformer.from_crs(C.METRIC_CRS, "EPSG:4326", always_xy=True).transform(st["xy"][:, 0], st["xy"][:, 1])
    d = d.merge(meta, on="station", how="left")
    dt = pd.to_datetime(d.date); d["year"] = dt.dt.year
    d["season"] = dt.dt.month.map({12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"})
    d["decade"] = pd.cut(d.year, [1990, 2000, 2010, 2020], labels=["1991-2000", "2001-2010", "2011-2020"]).astype(str)
    d["elev_band"] = pd.cut(d.elev, [-np.inf, 500, 1000, 1500, 2000, np.inf], labels=["<500", "500-1000", "1000-1500", "1500-2000", ">=2000"]).astype(str)
    d["box"] = "Rest_PRD"
    for b, (x0, x1, y0, y1) in C.BOXES.items():
        d.loc[(d.lon >= x0) & (d.lon <= x1) & (d.lat >= y0) & (d.lat <= y1), "box"] = b
    if variable == "precip":
        d["intensity"] = pd.cut(d.obs, [-0.01, 0.1, 1, 10, 20, 50, np.inf], labels=["dry", "0.1-1", "1-10", "10-20", "20-50", ">=50"]).astype(str)
    return d


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--runs", nargs="+", required=True); ap.add_argument("--climatology", action="store_true")
    a = ap.parse_args(); od = C.OUT_ROOT / "summary"; od.mkdir(parents=True, exist_ok=True)
    grid = P.load_grid(); cache = {}; comp = {}
    for run in a.runs:
        rd = C.OUT_ROOT / run; cfg = json.loads((rd / "run_config.json").read_text()); v = cfg["variable"]
        if v not in cache: cache[v] = P.attach_station_predictors(P.load_station_data(v), grid)
        d = load_loo(rd, v, cache[v])
        if d is not None:
            parts = [metrics(d, v).to_frame().T.assign(group="all", level="all")]
            for col in ["elev_band", "season", "decade", "box"] + (["intensity"] if v == "precip" else []):
                for lev, g in d.groupby(col): parts.append(metrics(g, v).to_frame().T.assign(group=col, level=lev))
            t = pd.concat(parts, ignore_index=True); t.insert(0, "run", run)
            t.to_csv(od / f"loo_metrics_{run}.csv", index=False); comp.setdefault(v, []).append(t)
            print(f"== {run}"); print(t[t.group.isin(["all", "box", "elev_band"])][["group", "level", "n", "bias", "mae", "rmse", "r", "kge", "cover95_pct"]].round(3).to_string(index=False))
        ds_files = sorted(glob.glob(str(rd / f"daily_summary_{v}_*.csv")))
        if ds_files:
            s = pd.concat([pd.read_csv(f) for f in ds_files]); s["decade"] = pd.cut(pd.to_datetime(s.date).dt.year, [1990, 2000, 2010, 2020], labels=["1991-2000", "2001-2010", "2011-2020"])
            cols = [c for c in s.columns if c.endswith("_p50") or c.startswith("pct_") or c in ("n_active",)]
            s.groupby("decade", observed=True)[cols].median().round(3).to_csv(od / f"support_by_period_{run}.csv")
        if a.climatology: climatology(rd, v, run, grid, od)
    for v, ts in comp.items():
        t = pd.concat(ts); key = ["n", "bias", "mae", "rmse", "r", "kge", "cover95_pct"] + (["wet_mae", "wet_pbias_pct"] if v == "precip" else [])
        w = t.pivot_table(index=["group", "level"], columns="run", values=key); w.to_csv(od / f"loo_comparison_{v}.csv")
        print(f"\n== comparison {v} (MAE)"); print(w["mae"].round(3).to_string())


def climatology(rd, v, run, grid, od):
    import netCDF4, xarray as xr
    files = sorted(glob.glob(str(rd / f"{v}_[12][0-9][0-9][0-9].nc")))
    if not files: return
    iy, ix = grid["iy"], grid["ix"]; N = len(iy); vn = "precipitation" if v == "precip" else "temperature"
    tot = np.zeros(N); ks = np.zeros(N); nday = np.zeros(N); years = 0
    d7 = {}; osc = np.zeros(N);
    # support radii (reviewer: effective resolution): histograms (0.5 km bins = storage precision) per region x decade
    SUP = ("D1", "D7", "R_pos_2_3", "R_pos_half", "R_pos_1_3", "R_abs90"); edges = np.arange(0, 400.5, 0.5)
    regions = {"All_PRD": np.ones(N, bool)}
    for b, (lo0, lo1, la0, la1) in C.BOXES.items(): regions[b] = (grid["lon"] >= lo0) & (grid["lon"] <= lo1) & (grid["lat"] >= la0) & (grid["lat"] <= la1)
    regions["Rest_PRD"] = ~np.any([m for k_, m in regions.items() if k_ != "All_PRD"], 0)
    for lo, hi in ((0, 500), (500, 1000), (1000, 1500), (1500, 2000), (2000, 9000)): regions[f"elev_{lo}-{hi}"] = (grid["elev"] >= lo) & (grid["elev"] < hi)
    H = {}; supmean = {}
    posl = {s: np.zeros(N) for s in ("DJF", "MAM", "JJA", "SON")}; nseas = {s: 0 for s in posl}
    for f in files:
        with netCDF4.Dataset(f) as ds:
            x = ds[vn][:].filled(np.nan)[:, iy, ix]; k = ds["kriging_std"][:].filled(np.nan)[:, iy, ix]; t = pd.to_datetime(ds["time"][:], unit="D", origin="1991-01-01")
        tot += np.nansum(x, 0); ks += np.nansum(k, 0); nday += np.isfinite(x).sum(0); years += 1
        dfile = Path(f).with_name(f"diagnostics_{Path(f).name}")
        if dfile.exists():
            with netCDF4.Dataset(dfile) as dg:
                dec = "1991-2000" if t[0].year <= 2000 else "2001-2010" if t[0].year <= 2010 else "2011-2020"
                a7 = dg["D7"][:].filled(np.nan)[:, iy, ix]; d7.setdefault(dec, []).append(np.nanmean(a7, 0))
                for sv in SUP:
                    if sv not in dg.variables: continue
                    a_ = a7 if sv == "D7" else dg[sv][:].filled(np.nan)[:, iy, ix]
                    supmean.setdefault((sv, dec), []).append(np.nanmean(a_, 0))
                    for rg, m in regions.items():
                        vals = a_[:, m].ravel(); vals = vals[np.isfinite(vals)]
                        h = np.histogram(vals, edges)[0]; H[(sv, rg, dec)] = H.get((sv, rg, dec), 0) + h; H[(sv, rg, "1991-2020")] = H.get((sv, rg, "1991-2020"), 0) + h
                dg["qc_flags"].set_auto_maskandscale(False); q = dg["qc_flags"][:][:, iy, ix]
                osc += ((q & 8) > 0).sum(0)
                if v == "temp":
                    for s_, mm in (("DJF", (12, 1, 2)), ("MAM", (3, 4, 5)), ("JJA", (6, 7, 8)), ("SON", (9, 10, 11))):
                        sel = t.month.isin(mm); posl[s_] += ((q[sel] & 4) > 0).sum(0); nseas[s_] += int(sel.sum())
    ny, nx = int(grid["ny"]), int(grid["nx"])
    def to2d(a):
        o = np.full((ny, nx), np.nan, np.float32); o[iy, ix] = a; return o
    out = {}
    if v == "precip": out["mean_annual_precipitation"] = (("y", "x"), to2d(tot / years), dict(units="mm yr-1"))
    else: out["mean_temperature"] = (("y", "x"), to2d(tot / np.maximum(nday, 1)), dict(units="degC"))
    out["mean_kriging_std"] = (("y", "x"), to2d(ks / np.maximum(nday, 1)), dict(units="mm day-1" if v == "precip" else "degC"))
    out["pct_oscillatory_days"] = (("y", "x"), to2d(100 * osc / np.maximum(nday, 1)), dict(units="%"))
    for dec, lst in d7.items(): out[f"mean_D7_{dec.replace('-', '_')}"] = (("y", "x"), to2d(np.nanmean(lst, 0)), dict(units="km"))
    for (sv, dec), lst in supmean.items():
        if sv != "D7": out[f"mean_{sv}_{dec.replace('-', '_')}"] = (("y", "x"), to2d(np.nanmean(lst, 0)), dict(units="km"))
    if H:
        rows = []; mid = 0.5 * (edges[1:] + edges[:-1])
        for (sv, rg, dec), h in H.items():
            c = np.cumsum(h); n = c[-1]
            if n == 0: continue
            qv = {f"p{p_}": mid[np.searchsorted(c, p_ / 100 * n)] for p_ in (10, 25, 50, 75, 90)}
            rows.append(dict(variable=v, support=sv, region=rg, period=dec, n_cell_days=int(n), mean=float((h * mid).sum() / n), **qv))
        tb = pd.DataFrame(rows).sort_values(["support", "region", "period"]); tb.to_csv(od / f"support_radii_by_region_period_{run}.csv", index=False)
        print(tb[(tb.period != "1991-2020") & tb.region.isin(["All_PRD", "Aosta_box", "Valtellina_box"]) & tb.support.isin(["D7", "R_pos_2_3", "R_pos_half", "R_pos_1_3"])]
              .pivot_table(index=["region", "period"], columns="support", values="p50").round(1).to_string())
    if v == "temp":
        for s_ in posl: out[f"pct_positive_lapse_{s_}"] = (("y", "x"), to2d(100 * posl[s_] / max(nseas[s_], 1)), dict(units="% of days"))
    xr.Dataset(out, coords=dict(y=grid["y"], x=grid["x"]), attrs=dict(run=run, years=f"{len(files)} files", crs=C.METRIC_CRS)).to_netcdf(
        od / f"cell_climatology_{run}.nc", encoding={k: dict(zlib=True, complevel=5) for k in out})
    print(f"cell climatology written for {run}")


if __name__ == "__main__":
    main()
