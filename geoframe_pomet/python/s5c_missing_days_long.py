#!/usr/bin/env python3
"""Step 5c - temporal-missingness test on a long window (all days of 2017-2020, 1461 days), precipitation.
Same design as s5b (50 % of valid station-days withheld, full-network normals), but long enough for the
production occurrence calibration (per station and calendar month, >= CALIB_MIN_DAYS days) and for stable
station-mean (climatological) bias. Station LOO only (no grid). Baseline = full-network LOO on the same days.
  python s5c_missing_days_long.py --realizations 10 --workers 14"""
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
    ap = argparse.ArgumentParser(); ap.add_argument("--realizations", type=int, default=10)
    ap.add_argument("--withhold", type=float, default=0.5); ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--seed", type=int, default=20261021); ap.add_argument("--y0", type=int, default=2017); ap.add_argument("--y1", type=int, default=2020)
    a = ap.parse_args(); v = "precip"; base = C.OUT_ROOT; root = base / f"thinning_days_long_{v}"; root.mkdir(exist_ok=True)
    dates = pd.date_range(f"{a.y0}-01-01", f"{a.y1}-12-31", freq="D"); dfile = root / "dates.csv"
    pd.DataFrame(dict(Date=dates.strftime("%Y-%m-%d"))).to_csv(dfile, index=False)
    st = P.load_station_data(v); ids = st["ids"]; ti = st["dates"].get_indexer(dates); assert (ti >= 0).all()
    V = st["values"][ti]; env0 = dict(os.environ, PO_DATA_ROOT=str(C.DATA_ROOT)); modes = ["direct", "anomaly"]
    for m in modes:
        tag = f"thinning_days_long_{v}/baseline_{m}"
        if not list((base / tag).glob("loo_*.parquet")):
            sh([str(HERE / "s3_run_days.py"), "--variable", v, "--mode", m, "--dates-csv", str(dfile), "--loo-only", "--no-diagnostics",
                "--workers", str(a.workers), "--tag", tag], env0); print("baseline", m, "done", flush=True)
    for r in range(1, a.realizations + 1):
        rd = root / f"r{r:02d}"; rd.mkdir(exist_ok=True); xf = rd / "withheld_days.parquet"
        if not xf.exists():
            rng = np.random.default_rng(a.seed + 1009 * r); W = np.isfinite(V) & (rng.random(V.shape) < a.withhold)
            di, si = np.nonzero(W); pd.DataFrame(dict(date=dates.values[di], station=ids[si])).to_parquet(xf)
        for m in modes:
            tag = f"thinning_days_long_{v}/r{r:02d}/{m}"
            if not list((base / tag).glob("loo_*.parquet")):
                sh([str(HERE / "s3_run_days.py"), "--variable", v, "--mode", m, "--dates-csv", str(dfile), "--loo-only",
                    "--exclude-station-days", str(xf), "--eval-excluded", "--no-diagnostics", "--workers", str(a.workers), "--tag", tag], env0)
        print(f"realization {r} done", flush=True)

if __name__ == "__main__":
    main()
