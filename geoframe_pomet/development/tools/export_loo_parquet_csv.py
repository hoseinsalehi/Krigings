#!/usr/bin/env python3
"""Export the LOO parquet tables of a Python run (loo_<var>_<Y>.parquet) to the CSV layout written by PoKriging
(loo_<var>_<Y>.csv.gz, columns date,station,obs,pred,kriging_std; float32 values printed with 9 significant digits
= exact float32 round trip). Row order is kept (the station halves of the calibration depend on first appearance).
Optionally (--recompute stations.csv out.csv) recompute the s9 table with the s9 code on the same data, with station
elevations from the converted stations CSV instead of the shapefile (test that the inputs are identical).
With Java LOO CSVs in <run_dir> (no parquet) only the recomputation is done.
Usage: export_loo_parquet_csv.py <run_dir> <var> <out_dir> [--recompute stations.csv out.csv]
"""
import glob, sys
from pathlib import Path
import numpy as np, pandas as pd
run, v, od = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3]); od.mkdir(parents=True, exist_ok=True)
fs = sorted(glob.glob(str(run / f"loo_{v}_*.parquet")))
cols = ["date", "station", "obs", "pred", "kriging_std"]
csv_in = not fs                                   # Java LOO tables (loo_<var>_<Y>.csv): only --recompute, read as float32
if csv_in: fs = sorted(glob.glob(str(run / f"loo_{v}_*.csv")))
read = (lambda f, columns: pd.read_csv(f, usecols=columns, dtype={"station": str, "date": str, "obs": np.float32, "pred": np.float32,
        "kriging_std": np.float32})[columns]) if csv_in else pd.read_parquet
for f in ([] if csv_in else fs):
    d = pd.read_parquet(f, columns=cols)
    d.to_csv(od / (Path(f).stem + ".csv.gz"), index=False, float_format="%.9g", compression=dict(method="gzip", compresslevel=1))
print("exported", len(fs), "files")
if "--recompute" in sys.argv:
    i = sys.argv.index("--recompute"); stf, outf = sys.argv[i + 1], sys.argv[i + 2]
    SEAS = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}
    # ---- verbatim s9_uncertainty_calibration.py:main (only the elevation source differs) ----
    L = pd.concat([read(f, columns=cols) for f in fs], ignore_index=True).dropna()
    L["station"] = L.station.astype(str); L["season"] = pd.to_datetime(L.date).dt.month.map(SEAS)
    S = pd.read_csv(stf, dtype={"id": str}); el = pd.Series(S["elevation"].to_numpy(float), index=S["id"]); L["elev"] = L.station.map(el)
    print("stations without elevation:", L.loc[L.elev.isna(), "station"].nunique())
    if v == "precip":
        L["cls"] = pd.cut(L.pred, [-1, 0.0, 1, 5, 20, 50, 1e9], labels=["0", "0-1", "1-5", "5-20", "20-50", ">50"]).astype(str)
    else:
        L["cls"] = pd.cut(L.elev, [-1e3, 500, 1000, 1500, 2000, 9e3], labels=["<500", "500-1000", "1000-1500", "1500-2000", ">=2000"]).astype(str)
    L["e"] = L.obs - L.pred; L["half"] = pd.factorize(L.station)[0] % 2
    rows = []
    for (s, c), g in list(L.groupby(["season", "cls"])) + [(("all", "all"), L)]:
        pos = g.kriging_std > 0; z = (g.e.abs() / g.kriging_std.where(pos)).dropna()
        k = float(np.quantile(z, 0.95) / 1.96) if len(z) > 50 else np.nan
        cov_cal = []
        for h in (0, 1):
            tr = g[pos & (g.half == h)]; te = g[pos & (g.half != h)]
            if len(tr) > 50 and len(te) > 50:
                kh = np.quantile(tr.e.abs() / tr.kriging_std, 0.95) / 1.96; cov_cal.append(np.mean(te.e.abs() <= 1.96 * kh * te.kriging_std))
        q = g.e.quantile([.05, .25, .5, .75, .95]).to_numpy()
        rows.append(dict(variable=v, season=s, cls=c, n=len(g), share_zero_std=float((~pos).mean()),
                         coverage95_raw=float(np.mean(g.e[pos].abs() <= 1.96 * g.kriging_std[pos])) if pos.any() else np.nan, scale_k=k,
                         coverage95_cal_crossval=float(np.mean(cov_cal)) if cov_cal else np.nan,
                         e_q05=q[0], e_q25=q[1], e_q50=q[2], e_q75=q[3], e_q95=q[4]))
    pd.DataFrame(rows).to_csv(outf, index=False); print("dtypes e/z:", L.e.dtype, type(np.quantile(L.e.abs().to_numpy()[:1000], .95)))
