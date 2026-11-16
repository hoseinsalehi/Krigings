#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 5 - repeated 50 % network thinning (reviewer comment: robustness to missing data).

For each realization r a fixed random half of the stations is withheld for the whole test
period (default: the 2011-2020 representative days of the 360-day set, the data-rich
decade). For the anomaly method the normals are REBUILT from the retained stations only
(steps 1-2 with the tuned parameters), so the climatological layer gets no information
from the withheld half. Each method then predicts:
  (a) every withheld station-day   (error against observations), and
  (b) a fixed sample of grid cells (difference from the full-network field).

  python s5_thinning.py --variable precip --dates-csv selected_days.csv --realizations 30 --workers 12

Output: <OUT_ROOT>/thinning_<var>/summary_withheld.csv, summary_field.csv, per-realization files.
"""
import argparse, json, os, shutil, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import po_config as C, po_common as P

HERE = Path(__file__).resolve().parent


def sh(args, env):
    r = subprocess.run([sys.executable] + args, env=env, capture_output=True, text=True)
    if r.returncode: raise RuntimeError(f"{' '.join(args)}\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    return r.stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variable", choices=["precip", "temp"], required=True)
    ap.add_argument("--dates-csv", required=True); ap.add_argument("--start-year", type=int, default=2011)
    ap.add_argument("--realizations", type=int, default=30); ap.add_argument("--retain", type=float, default=0.5)
    ap.add_argument("--modes", default="direct,anomaly"); ap.add_argument("--grid-sample", type=int, default=5000)
    ap.add_argument("--workers", type=int, default=12); ap.add_argument("--seed", type=int, default=20261008)
    a = ap.parse_args(); v = a.variable; base = C.OUT_ROOT; root = base / f"thinning_{v}"; root.mkdir(parents=True, exist_ok=True)
    dc = pd.read_csv(a.dates_csv); col = "Date" if "Date" in dc else "date"; dc = dc[pd.to_datetime(dc[col]).dt.year >= a.start_year]
    dfile = root / "dates.csv"; dc[[col]].rename(columns={col: "Date"}).to_csv(dfile, index=False)
    ids = P.load_station_data(v)["ids"]; modes = a.modes.split(",")
    params = base / "normals" / f"params_{v}.json"
    env0 = dict(os.environ, PO_DATA_ROOT=str(C.DATA_ROOT))      # realizations keep reading the harmonised, cleaned station data
    # full-network reference fields at the grid sample
    for m in modes:
        tag = f"full_{m}"
        if not list((root / tag).glob("grid_sample_*.parquet")):
            sh([str(HERE / "s3_run.py"), "--variable", v, "--mode", m, "--dates-csv", str(dfile), "--grid-sample", str(a.grid_sample),
                "--no-loo", "--no-diagnostics", "--workers", str(a.workers), "--tag", f"thinning_{v}/{tag}"], env0)
    rows_w, rows_f = [], []
    for r in range(1, a.realizations + 1):
        rd = root / f"r{r:02d}"; rd.mkdir(exist_ok=True)
        rng = np.random.default_rng(a.seed + 1009 * r); keep = rng.random(len(ids)) < a.retain
        (rd / "withheld.txt").write_text("\n".join(ids[~keep]))
        env = dict(env0, PO_OUT_ROOT=str(rd))
        shutil.copy(base / "grid_predictors.npz", rd / "grid_predictors.npz")
        if "anomaly" in modes and not (rd / "normals" / f"normals_{v}.nc").exists():
            sh([str(HERE / "s1_station_normals.py"), "--variable", v, "--exclude-stations", str(rd / "withheld.txt")], env)
            sh([str(HERE / "s2_normals_grid.py"), "--variable", v] + (["--params", str(params)] if params.exists() else []), env)
        if "anomaly" in modes and v == "precip" and C.PRECIP_METHOD == "twopart" and not (rd / "normals" / "wetfreq_precip.nc").exists():
            sh([str(HERE / "s1w_wetfreq_normals.py"), "--exclude-stations", str(rd / "withheld.txt")], env)   # two-part: frequency normals of the retained network
        for m in modes:
            tag = f"thin_{m}"
            if not (rd / tag / "run_config.json").exists() or not list((rd / tag).glob("loo_*.parquet")):
                sh([str(HERE / "s3_run.py"), "--variable", v, "--mode", m, "--dates-csv", str(dfile), "--grid-sample", str(a.grid_sample),
                    "--exclude-stations", str(rd / "withheld.txt"), "--eval-excluded", "--no-diagnostics", "--workers", str(a.workers), "--tag", tag], env)
            w = pd.concat([pd.read_parquet(f) for f in (rd / tag).glob("loo_*.parquet")]); e = w.pred - w.obs; ok = e.notna()
            row = dict(realization=r, mode=m, n=int(ok.sum()), bias=e[ok].mean(), mae=e[ok].abs().mean(), rmse=np.sqrt((e[ok] ** 2).mean()),
                       r=np.corrcoef(w.pred[ok], w.obs[ok])[0, 1])
            g = w[ok].groupby("station").agg(o=("obs", "mean"), p=("pred", "mean"), n=("obs", "size")); g = g[g.n >= 0.5 * g.n.max()]
            if v == "precip":
                wet = ok & (w.obs >= C.WET_THRESHOLD); row.update(wet_mae=e[wet].abs().mean(), wet_pbias_pct=100 * e[wet].sum() / w.obs[wet].sum())
                g = g[g.o > 0.5]; rel = 100 * (g.p - g.o) / g.o          # bias of the station's period-mean precipitation
                row.update(stmean_abs_relbias_median_pct=rel.abs().median(), stmean_abs_relbias_p90_pct=rel.abs().quantile(.9), stmean_relbias_mean_pct=rel.mean())
            else:
                b_ = g.p - g.o; row.update(stmean_abs_bias_median_c=b_.abs().median(), stmean_abs_bias_p90_c=b_.abs().quantile(.9), stmean_bias_mean_c=b_.mean())
            rows_w.append(row)
            th = pd.concat([pd.read_parquet(f) for f in (rd / tag).glob("grid_sample_*.parquet")])
            fu = pd.concat([pd.read_parquet(f) for f in (root / f"full_{m}").glob("grid_sample_*.parquet")])
            j = th.merge(fu, on=["date", "cell"], suffixes=("_thin", "_full"))
            sc = 0.1 if v == "precip" else 0.01; d = (j.value_thin.astype(float) - j.value_full.astype(float)) * sc
            rows_f.append(dict(realization=r, mode=m, n=len(d), mean_diff=d.mean(), mad=d.abs().mean(), rmsd=np.sqrt((d ** 2).mean()),
                               p95_absdiff=d.abs().quantile(0.95), rel_mad_pct=100 * d.abs().sum() / (j.value_full.astype(float) * sc).abs().sum()))
        print(f"realization {r}: " + " | ".join(f"{x['mode']} MAE {x['mae']:.3f}" for x in rows_w[-len(modes):]), flush=True)
    W = pd.DataFrame(rows_w); F = pd.DataFrame(rows_f)
    W.to_csv(root / "withheld_by_realization.csv", index=False); F.to_csv(root / "field_by_realization.csv", index=False)
    def summ(T):
        num = [c for c in T.columns if c not in ("realization", "mode")]
        return T.groupby("mode")[num].agg(["mean", lambda x: x.quantile(.05), "median", lambda x: x.quantile(.95)]).rename(columns={"<lambda_0>": "p05", "<lambda_1>": "p95"}, level=1)
    summ(W).to_csv(root / "summary_withheld.csv"); summ(F).to_csv(root / "summary_field.csv")
    print(W.groupby("mode")[["mae", "rmse", "bias"]].describe().round(3).to_string())
    print(F.groupby("mode")[["mad", "rmsd", "rel_mad_pct"]].mean().round(3).to_string())


if __name__ == "__main__":
    main()
