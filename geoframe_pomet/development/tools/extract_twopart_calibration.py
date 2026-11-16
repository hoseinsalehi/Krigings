#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Production (1991-2020, monthly) occurrence-calibration thresholds of the two-part method for the reference cells,
recomputed from ~/Downloads/PoGrid_v5/precip_anomaly/_twopart/grid_<Y>.npz exactly as s3_run.finalize_twopart does
(read-only). Then applies them to the reference pass-1 pi/amount and checks the production final precip_2020.nc.
Writes reference/twopart_production_calibration/{thresholds_cells.csv, grid_final_production.csv.gz, loo_final_production.csv.gz,
check.json}.
"""
import json, os, sys
from pathlib import Path
CODE = Path.home() / "Downloads/PoGrid_v3_code"; sys.path.insert(0, str(CODE)); sys.dont_write_bytecode = True
import numpy as np, pandas as pd, xarray as xr

B = Path.home() / "Documents/Po_precip/GEOframe_PoKriging/reference"; V5 = Path.home() / "Downloads/PoGrid_v5"
OUT = B / "twopart_production_calibration"; OUT.mkdir(exist_ok=True)


def main():
    os.environ.update(PO_STATION_SET="v4", PO_PRECIP_METHOD="twopart", PO_OUT_ROOT=str(V5))
    import s3_run as S3, po_config as C
    cells = pd.read_csv(B / "daily/targets_cells.csv"); cid = cells.cell.to_numpy()
    dates = json.load(open(B / "daily/precip_twopart/run_info.json"))["dates"]
    PI, days = [], []
    for y in range(1991, 2021):
        z = np.load(V5 / f"precip_anomaly/_twopart/grid_{y}.npz")
        PI.append(S3.dec_tp(z["pi"][:, cid], "pi").astype(np.float32)); days += list(z["days"])
    PI = np.concatenate(PI); days = pd.DatetimeIndex(days); mon = days.month.to_numpy() - 1
    with xr.open_dataset(V5 / "normals/wetfreq_precip.nc") as ds:
        PGF = ds["wetfreq"].values[:, cells.iy.to_numpy(), cells.ix.to_numpy()].astype(float)
    assert np.bincount(mon, minlength=12).min() >= C.CALIB_MIN_DAYS
    T = np.full((12, len(cid)), np.nan); K = np.full((12, len(cid)), np.nan)
    for m in range(12):
        r = np.flatnonzero(mon == m); k = PGF[m] * np.isfinite(PI[r]).sum(0)
        T[m] = S3._thresholds(PI[r], k); K[m] = np.clip(np.rint(k), 0, np.isfinite(PI[r]).sum(0))
    th = pd.DataFrame(dict(cell=np.repeat(cid, 12), month=np.tile(np.arange(1, 13), len(cid)),
                           threshold=T.T.ravel(), wet_days_target=K.T.ravel().astype(int), p_normal=PGF.T.ravel()))
    th.to_csv(OUT / "thresholds_cells.csv", index=False, float_format="%.10g")
    # apply to the reference pass-1 fields (quantised as production) and compare with production final values
    g = pd.read_csv(B / "daily/precip_twopart/grid.csv.gz")
    g["month"] = pd.to_datetime(g.date).dt.month
    g = g.merge(th[["cell", "month", "threshold"]], on=["cell", "month"], suffixes=("_pooled10d", ""))
    wet = g.pi_q.astype(np.float32) >= g.threshold.astype(np.float32);  # float32 as in production (ties = wet)
    val = np.where(wet, g.amount_q, 0.0)
    with xr.open_dataset(V5 / "precip_anomaly/precip_2020.nc") as ds:
        P = ds["precipitation"].sel(time=pd.to_datetime(dates)).values
    di = {d: i for i, d in enumerate(dates)}; ci = dict(zip(cells.cell, zip(cells.iy, cells.ix)))
    iy = g.cell.map(lambda c: ci[c][0]).to_numpy(); ix = g.cell.map(lambda c: ci[c][1]).to_numpy()
    prod = P[g.date.map(di).to_numpy(), iy, ix]
    d = np.abs(np.rint(val / 0.1) * 0.1 - prod)
    res = dict(n=int(len(d)), max_abs_diff_mm=float(np.nanmax(d)), frac_exact=float((d < 1e-6).mean()),
               wet_fraction=float(wet.mean()), calibration="monthly 1991-2020 (production)", days_used=int(len(days)))
    pd.DataFrame(dict(date=g.date, cell=g.cell, wet=wet.astype(np.uint8), value=val, value_production_file=prod)) \
        .to_csv(OUT / "grid_final_production.csv.gz", index=False, float_format="%.6g")
    L = pd.read_parquet(V5 / "precip_anomaly/loo_precip_2020.parquet"); L = L[L.date.astype(str).isin(dates)]
    L.to_csv(OUT / "loo_final_production.csv.gz", index=False, float_format="%.8g")
    res["loo_rows"] = int(len(L)); res["loo_wet_pred_frac"] = float(L.wet_pred.mean())
    (OUT / "check.json").write_text(json.dumps(res, indent=1)); print(res)


if __name__ == "__main__":
    main()
