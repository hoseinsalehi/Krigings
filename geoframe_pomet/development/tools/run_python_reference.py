#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python reference outputs for the GEOframe/OMS3 port (equivalence tests).

Runs the UNMODIFIED PoGrid_v3_code engine (s3_run.process_day / po_core.estimate_targets) for a few dates on a
clipped target set, plus station leave-one-out on the same dates. Nothing in PoGrid_v3_code is edited; the driver
only changes, in memory, three things:

 1. targets: whole production chunks of 2048 grid cells (production order = row-major valid cells of
    grid_predictors.npz, CHUNK = 2048). Because the PSO random stream is seeded per chunk with
    seed_for(day_tag, s0) and consumed in target order, a chunk-aligned subset reproduces the production seeds:
    the chunk offset s0 of the subset run is mapped back to the production offset 2048*k (grid calls only;
    LOO calls are untouched, they are a single chunk starting at 0 exactly as in production).
 2. encoders: s3_run.enc / enc_tp / value_enc are replaced by float64 pass-through (value_enc scale 2^-40,
    i.e. rounding below 1e-12), so the reference holds unquantised values. Production-quantised values are
    re-derived in compare_with_production.py.
 3. optional --salt N (N > 0): every seed_for() call gets an extra part 'saltN' -> an independent PSO random
    stream, used to measure the PSO seed sensitivity (the tolerance a non-bitwise Java port can be held to).

Environment (set by the caller, see run_reference.sh):
  PO_STATION_SET=v4  PO_OUT_ROOT=<root with normals/ + grid_predictors.npz>  PO_DATA_ROOT=~/Downloads/PoGrid_v5
  PO_PRECIP_METHOD=ratio|twopart   PYTHONDONTWRITEBYTECODE=1  NUMBA_CACHE_DIR=<scratch>

Usage:
  python run_python_reference.py --method precip_twopart --out <dir> [--salt 0] [--workers 10]
methods: precip_twopart precip_ratio precip_direct temp_anomaly temp_anomaly_nolapse temp_direct
"""
import argparse, json, os, sys, time, platform
from pathlib import Path
from types import SimpleNamespace

CODE = Path(os.environ.get("PO_CODE", str(Path.home() / "Downloads/PoGrid_v3_code")))
sys.path.insert(0, str(CODE))
import numpy as np, pandas as pd
import multiprocessing as mp

DATES = ["2020-01-05", "2020-03-02", "2020-04-11", "2020-06-06", "2020-07-04",
         "2020-07-15", "2020-08-08", "2020-08-13", "2020-10-02", "2020-12-05"]
CHUNKS = [2, 7, 9, 20]          # production chunks: Valtellina, Gran San Bernardo/Mont Blanc, Aosta, Po plain
METHODS = {  # variable, mode, anom_lapse, PO_PRECIP_METHOD required
    "precip_twopart": ("precip", "anomaly", "local50", "twopart"),
    "precip_ratio": ("precip", "anomaly", "local50", "ratio"),
    "precip_direct": ("precip", "direct", "local50", None),
    "temp_anomaly": ("temp", "anomaly", "local50", None),
    "temp_anomaly_nolapse": ("temp", "anomaly", "none", None),
    "temp_direct": ("temp", "direct", "local50", None),
}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--method", choices=list(METHODS), required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--salt", type=int, default=0)
    ap.add_argument("--workers", type=int, default=10); ap.add_argument("--no-loo", action="store_true")
    a = ap.parse_args()
    var, mode, anom_lapse, pm = METHODS[a.method]
    if pm is not None and os.environ.get("PO_PRECIP_METHOD", "ratio") != pm:
        sys.exit(f"{a.method} needs PO_PRECIP_METHOD={pm}")
    import po_config as C, po_core as K, s3_run as S3
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------- in-memory patches (see docstring)
    _seed0, _est0 = K.seed_for, K.estimate_targets
    STATE = {"map": None}
    chunk_map = {i * C.CHUNK: k * C.CHUNK for i, k in enumerate(CHUNKS)}

    def seed_for(*parts):
        if STATE["map"] is not None and len(parts) == 2 and isinstance(parts[1], (int, np.integer)):
            parts = (parts[0], STATE["map"][int(parts[1])])
        if a.salt: parts = tuple(parts) + (f"salt{a.salt}",)
        return _seed0(*parts)

    def estimate_targets(txy, *args, **kw):
        tg = S3.G.get("targets")
        STATE["map"] = chunk_map if (tg is not None and txy is tg["xy"]) else None
        try: return _est0(txy, *args, **kw)
        finally: STATE["map"] = None
    K.seed_for = seed_for; K.estimate_targets = estimate_targets
    S3.enc = lambda x, name: np.asarray(x, float)
    S3.enc_tp = lambda x, name: np.asarray(x, float)
    S3.value_enc = lambda variable: ("f8", 2.0 ** -40, np.nan)

    # ---------------------------------------------------------------- shared state as in s3_run.main
    args = SimpleNamespace(variable=var, mode=mode, anom_lapse=anom_lapse, exclude_stations=None, eval_excluded=False,
                           no_diagnostics=False, no_loo=a.no_loo, loo_only=True, grid_sample=0, points=None)
    t_setup = time.time()
    shared, grid, _, sel = S3.build_shared(args)
    N = len(grid["xy"]); cells = np.concatenate([np.arange(k * C.CHUNK, min(N, (k + 1) * C.CHUNK)) for k in CHUNKS])
    shared["targets"] = dict(xy=grid["xy"][cells], elev=grid["elev"][cells], cell=cells, elev_filled=grid["elev_filled"][cells])
    dates = shared["dates"]; want = set(pd.to_datetime(DATES))
    day_idx = np.array([i for i in np.flatnonzero(sel) if dates[i] in want])
    t_setup = time.time() - t_setup
    tp = bool(shared.get("twopart"))
    print(f"{a.method}: kernels={K.KERNELS} lapse={shared['lapse_mode']} twopart={tp} days={len(day_idx)} "
          f"targets={len(cells)} stations={int(shared['usable'].sum())} salt={a.salt} workers={a.workers}", flush=True)

    t0 = time.time(); res = []
    with mp.get_context("fork").Pool(a.workers, initializer=S3.init, initargs=(shared,)) as pool:
        S3.G.update(shared)
        for r in pool.imap_unordered(S3.process_day, day_idx, chunksize=1): res.append(r)
    wall = time.time() - t0
    res.sort(key=lambda r: r["date"])

    # ---------------------------------------------------------------- grid table
    rows = []
    for r in res:
        if "grid" not in r: continue
        g = r.pop("grid"); d = pd.DataFrame({k: np.asarray(x) for k, x in g.items()})
        if not tp:
            d["value"] = d["value"] * 2.0 ** -40; d["kstd"] = d["kstd"] * 2.0 ** -40
        d.insert(0, "cell", cells); d.insert(0, "date", r["date"]); rows.append(d)
    G_ = pd.concat(rows, ignore_index=True)
    geo = pd.DataFrame(dict(cell=cells, iy=grid["iy"][cells], ix=grid["ix"][cells], x=grid["xy"][cells, 0], y=grid["xy"][cells, 1],
                            elevation=grid["elev"][cells], elev_filled=grid["elev_filled"][cells].astype(int),
                            chunk=cells // C.CHUNK))
    loo = pd.concat([r.pop("loo") for r in res if "loo" in r], ignore_index=True) if not a.no_loo else None
    info = {}
    if tp:
        # occurrence calibration exactly as s3_run.finalize_twopart for a run of < CALIB_MIN_DAYS days per month
        # (pooled over the run's days); pi and amount quantised as production stores them before calibration
        def _q(x, name):
            dt, sc, fill = S3.TP[name]; info_ = np.iinfo(np.dtype(dt)); qq = np.rint(np.asarray(x, float) / sc)
            return np.where(np.isfinite(qq), np.clip(qq, info_.min + (dt[0] == "i"), info_.max - (dt[0] == "u")), fill).astype(dt)
        PIf = G_.pivot(index="date", columns="cell", values="pi").reindex(columns=cells)
        PI = np.stack([S3.dec_tp(_q(PIf.values[i], "pi"), "pi") for i in range(len(PIf))]).astype(np.float32)
        mon = pd.DatetimeIndex(PIf.index).month.to_numpy() - 1
        PGF = shared["PGF"][:, cells]
        monthly = np.bincount(mon, minlength=12).min() >= C.CALIB_MIN_DAYS
        assert not monthly
        k = np.nansum(np.where(np.isfinite(PI), PGF[mon], 0.0), 0); t = S3._thresholds(PI, k); WET = PI >= t[None]
        info["grid_calibration"] = "pooled"
        lut = {d: i for i, d in enumerate(PIf.index)}; ci = {c: j for j, c in enumerate(cells)}
        ii = G_.date.map(lut).to_numpy(); jj = G_.cell.map(ci).to_numpy()
        G_["pi_q"] = PI[ii, jj]; G_["wet"] = WET[ii, jj].astype(np.uint8); G_["threshold"] = t[jj]
        G_["amount_q"] = S3.dec_tp(_q(G_.amount, "amount"), "amount"); G_["samount_q"] = S3.dec_tp(_q(G_.samount, "samount"), "samount")
        G_["p_target"] = PGF[mon[ii], jj]
        val = np.where(G_.wet == 1, G_.amount_q, 0.0); G_["value"] = np.where(np.isfinite(G_.pi_q), val, np.nan)
        G_["kstd"] = S3._mix_std(G_.pi_q.to_numpy(float), G_.amount_q.to_numpy(float), G_.samount_q.to_numpy(float))
        if loo is not None:
            L, info["loo_calibration"] = S3._calibrate(loo, "station")
            L["pred"] = np.where(L.wet, L.amount, 0.0).astype(np.float32)
            L["kriging_std"] = S3._mix_std(L.pi.to_numpy(), L.amount.to_numpy(), L.samount.to_numpy()).astype(np.float32)
            L["wet_pred"] = L.wet.astype(np.uint8); loo = L.drop(columns=["wet"])
    G_ = G_.merge(geo[["cell", "iy", "ix"]], on="cell", how="left")
    G_.to_csv(out / "grid.csv.gz", index=False, float_format="%.10g")
    geo.to_csv(out.parent / "targets_cells.csv", index=False) if a.salt == 0 else None
    if loo is not None: loo.to_csv(out / "loo.csv.gz", index=False, float_format="%.8g")
    summ = pd.DataFrame(res); summ.to_csv(out / "daily_summary.csv", index=False)
    import numba, scipy
    info.update(method=a.method, variable=var, mode=mode, anom_lapse=anom_lapse, lapse_mode=shared["lapse_mode"],
                precip_method=os.environ.get("PO_PRECIP_METHOD", "ratio"), kernels=K.KERNELS, salt=a.salt, dates=DATES,
                chunks=CHUNKS, n_targets=int(len(cells)), n_stations_usable=int(shared["usable"].sum()),
                workers=a.workers, blas_threads=os.environ.get("OPENBLAS_NUM_THREADS"), setup_s=round(t_setup, 2),
                wall_s=round(wall, 2), worker_s_per_day_mean=round(float(summ.seconds.mean()), 3),
                worker_s_per_day_max=round(float(summ.seconds.max()), 3), python=platform.python_version(),
                numpy=np.__version__, scipy=scipy.__version__, numba=numba.__version__, pandas=pd.__version__,
                machine=f"{platform.system()} {platform.machine()}, {os.cpu_count()} logical CPUs",
                out_root=str(C.OUT_ROOT), data_root=str(C.DATA_ROOT), obs_file=str(C.DATA[var]["obs"]),
                method_text=shared["method_text"])
    (out / "run_info.json").write_text(json.dumps(info, indent=1, default=str))
    print(json.dumps({k: info[k] for k in ("method", "salt", "wall_s", "worker_s_per_day_mean", "setup_s")}), flush=True)


if __name__ == "__main__":
    main()
