#!/usr/bin/env python3
"""Step 5b - temporal-missingness test (reviewer comment 2.3: 'randomly deleting 50 % of the data').
Same days and grid sample as s5_thinning (2011-2020 test days). In each realization, every valid station-day of the
test days is withheld with probability 0.5 (independently per station and day). Stations are NOT removed from the
network: the climatological normals are those of the full network (production normals). Each method predicts the
withheld station-days from the remaining observations of the same day; the target normal at a withheld station is
its leave-one-out grid normal (conservative). Output: <OUT_ROOT>/thinning_days_<var>/...
  python s5b_missing_days.py --variable precip --realizations 30 --workers 14"""
import argparse, os, subprocess, sys
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
    ap = argparse.ArgumentParser(); ap.add_argument("--variable", choices=["precip", "temp"], required=True)
    ap.add_argument("--realizations", type=int, default=30); ap.add_argument("--withhold", type=float, default=0.5)
    ap.add_argument("--modes", default="direct,anomaly"); ap.add_argument("--grid-sample", type=int, default=5000)
    ap.add_argument("--workers", type=int, default=14); ap.add_argument("--seed", type=int, default=20261020)
    ap.add_argument("--start", type=int, default=1)
    a = ap.parse_args(); v = a.variable; base = C.OUT_ROOT
    src = base / f"thinning_{v}"; root = base / f"thinning_days_{v}"; root.mkdir(exist_ok=True)
    dfile = src / "dates.csv"; dates = pd.to_datetime(pd.read_csv(dfile)["Date"])
    st = P.load_station_data(v); ids = st["ids"]; ti = st["dates"].get_indexer(dates); assert (ti >= 0).all()
    V = st["values"][ti]                                       # (D, S) values on the test days (cleaned, harmonised)
    env0 = dict(os.environ, PO_DATA_ROOT=str(C.DATA_ROOT)); modes = a.modes.split(",")
    rows_w, rows_f = [], []
    for r in range(a.start, a.realizations + 1):
        rd = root / f"r{r:02d}"; rd.mkdir(exist_ok=True); xf = rd / "withheld_days.parquet"
        if not xf.exists():
            rng = np.random.default_rng(a.seed + 1009 * r); W = np.isfinite(V) & (rng.random(V.shape) < a.withhold)
            di, si = np.nonzero(W); pd.DataFrame(dict(date=dates.values[di], station=ids[si])).to_parquet(xf)
        for m in modes:
            tag = f"thinning_days_{v}/r{r:02d}/{m}"; od = base / tag
            if not list(od.glob("loo_*.parquet")):
                sh([str(HERE / "s3_run_days.py"), "--variable", v, "--mode", m, "--dates-csv", str(dfile), "--grid-sample", str(a.grid_sample),
                    "--exclude-station-days", str(xf), "--eval-excluded", "--no-diagnostics", "--workers", str(a.workers), "--tag", tag], env0)
            w = pd.concat([pd.read_parquet(f) for f in od.glob("loo_*.parquet")]); e = w.pred - w.obs; ok = e.notna()
            row = dict(realization=r, mode=m, n=int(ok.sum()), bias=e[ok].mean(), mae=e[ok].abs().mean(), rmse=np.sqrt((e[ok] ** 2).mean()),
                       r=np.corrcoef(w.pred[ok], w.obs[ok])[0, 1])
            g = w[ok].groupby("station").agg(o=("obs", "mean"), p=("pred", "mean"), n=("obs", "size")); g = g[g.n >= 20]
            if v == "precip":
                wo, wp = (w.obs[ok] >= 1), (w.pred[ok] >= 1)
                row.update(wetfreq_obs=100 * wo.mean(), wetfreq_pred=100 * wp.mean(), int_obs=w.obs[ok][wo].mean(), int_pred=w.pred[ok][wp].mean(),
                           pod=(wo & wp).sum() / wo.sum(), far=(wp & ~wo).sum() / wp.sum())
                g = g[g.o > 0.5]; rel = 100 * (g.p - g.o) / g.o
                row.update(stmean_abs_relbias_median_pct=rel.abs().median(), stmean_relbias_mean_pct=rel.mean(), n_st=len(g))
            else:
                b_ = g.p - g.o; row.update(stmean_abs_bias_median_c=b_.abs().median(), stmean_bias_mean_c=b_.mean(), n_st=len(g))
            rows_w.append(row)
            th = pd.concat([pd.read_parquet(f) for f in od.glob("grid_sample_*.parquet")])
            fu = pd.concat([pd.read_parquet(f) for f in (src / f"full_{m}").glob("grid_sample_*.parquet")])
            j = th.merge(fu, on=["date", "cell"], suffixes=("_thin", "_full"))
            sc = 0.1 if v == "precip" else 0.01; d = (j.value_thin.astype(float) - j.value_full.astype(float)) * sc
            rows_f.append(dict(realization=r, mode=m, n=len(d), mean_diff=d.mean(), mad=d.abs().mean(), p95_absdiff=d.abs().quantile(0.95),
                               rel_mad_pct=100 * d.abs().sum() / (j.value_full.astype(float) * sc).abs().sum()))
        print(f"realization {r}: " + " | ".join(f"{x['mode']} MAE {x['mae']:.3f} stbias {x.get('stmean_abs_relbias_median_pct', x.get('stmean_abs_bias_median_c')):.3f}" for x in rows_w[-len(modes):]), flush=True)
        pd.DataFrame(rows_w).to_csv(root / f"withheld_by_realization_from{a.start}.csv", index=False)
        pd.DataFrame(rows_f).to_csv(root / f"field_by_realization_from{a.start}.csv", index=False)

if __name__ == "__main__":
    main()
