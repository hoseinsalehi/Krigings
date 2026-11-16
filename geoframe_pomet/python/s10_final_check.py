#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 10 - final completeness and plausibility check of every deliverable. Writes
<OUT_ROOT>/summary/final_check.csv and final_check.md; exits 1 if any check FAILs."""
import glob, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd, xarray as xr
import po_config as C, po_common as P

O = C.OUT_ROOT; rows = []
def chk(item, ok, detail="", warn=False):
    stt = "PASS" if ok else ("WARN" if warn else "FAIL")
    rows.append(dict(item=item, status=stt, detail=detail)); print(f"{stt}  {item}  {detail}", flush=True)

grid = P.load_grid(); mask = np.zeros((int(grid["ny"]), int(grid["nx"])), bool); mask[grid["iy"], grid["ix"]] = True
for v, vn in (("precip", "precipitation"), ("temp", "temperature")):
    rd = O / f"{v}_anomaly"; years = range(1991, 2021)
    for kind in (f"{v}_", f"diagnostics_{v}_", f"loo_{v}_", f"daily_summary_{v}_"):
        n = len([y for y in years if list(rd.glob(f"{kind}{y}.*"))]); chk(f"{v}: {kind}<year> files", n == 30, f"{n}/30")
    chk(f"{v}: no .part files", not list(rd.glob("*.part")))
    st = P.load_station_data(v); obs = pd.DataFrame(st["values"], index=st["dates"])
    nan_in, bad_len, incomplete, worst, viol, jumps, hot_cd, cold_cd, tot_cd = 0, [], [], 0.0, 0, 0, 0, 0, 0
    for y in years:
        f = rd / f"{v}_{y}.nc"
        if not f.exists(): continue
        with xr.open_dataset(f) as ds:
            a = ds[vn].values; t = pd.to_datetime(ds.time.values)
            if ds.attrs.get("complete") != "yes": incomplete.append(y)
            if len(t) != len(pd.date_range(f"{y}-01-01", f"{y}-12-31")): bad_len.append(y)
        nan_in += int(np.isnan(a[:, mask]).sum()); o = obs.reindex(t)
        if v == "precip":
            r = np.nanmax(a[:, mask], 1) / np.maximum(o.max(1).to_numpy(), 1.0); worst = max(worst, float(np.nanmax(r))); viol += int((np.nanmin(a[:, mask], 1) < 0).sum())
        else:
            mx = np.nanmax(a[:, mask], 1) - o.max(1).to_numpy()
            zc = grid["dem2d"][mask] if "dem2d" in grid else None
            ztop = np.array([st["elev"][np.isfinite(r_)].max() if np.isfinite(r_).any() else np.nan for r_ in o.to_numpy()])
            inside = zc[None, :] <= (ztop[:, None] + 300.0)                         # cells not above the highest active station + 300 m
            mn = o.min(1).to_numpy() - np.nanmin(np.where(inside, a[:, mask], np.nan), 1)
            hot_cd += int((a[:, mask] - o.max(1).to_numpy()[:, None] > 3).sum()); tot_cd += a[:, mask].size
            cold_cd += int(((o.min(1).to_numpy()[:, None] - np.where(inside, a[:, mask], np.nan)) > 10).sum())
            worst = max(worst, float(np.nanmax(mx))); viol += int((mn > 10).sum()); jumps += int((np.abs(np.diff(a[:, mask], axis=0)) > 15).sum())
        with xr.open_dataset(rd / f"diagnostics_{v}_{y}.nc") as dg: miss = [k for k in ("D1", "D7", "R_pos_2_3", "R_pos_half", "R_pos_1_3", "R_abs90", "qc_flags") if k not in dg]
        if miss: chk(f"{v} {y}: diagnostics fields", False, str(miss))
    chk(f"{v}: complete attribute", not incomplete, str(incomplete)); chk(f"{v}: days per year", not bad_len, str(bad_len))
    chk(f"{v}: no NaN inside the domain", nan_in == 0, f"{nan_in} NaN cell-days")
    if v == "precip":
        chk("precip: no negative values", viol == 0, f"{viol} days"); chk("precip: daily grid max <= 3 x station max", worst <= 3.0, f"worst ratio {worst:.2f}")
    else:
        chk("temp: daily grid max <= station max + 3 C", worst <= 3.0, f"worst excess {worst:.2f} C; {hot_cd} cell-days ({1e6*hot_cd/tot_cd:.2f} per million)", warn=hot_cd < 1e-5 * tot_cd)
        chk("temp: daily grid min >= station min - 10 C (cells up to highest active station + 300 m)", viol == 0, f"{viol} days; {cold_cd} cell-days", warn=cold_cd < 1e-5 * tot_cd); chk("temp: >15 C day-to-day jumps < 5000 cell-days", jumps < 5000, f"{jumps}")
    for f in (C.DATA_ROOT / C.SHIFT_FILE.format(variable=v), C.DATA_ROOT / "qc" / f"station_day_flags_{v}.parquet",
              O / "normals" / f"normals_{v}.nc", O / "summary" / f"loo_metrics_{v}_anomaly.csv", O / "summary" / f"support_radii_by_region_period_{v}_anomaly.csv",
              O / "summary" / f"cell_climatology_{v}_anomaly.nc", O / "summary" / f"uncertainty_calibration_{v}.csv",
              O / f"thinning_{v}" / "summary_withheld.csv", O / f"point_test_{v}" / "summary.csv", O / "reference_comparison" / f"summary_{v}.csv",
              O / f"loo360_{v}_direct" / f"daily_summary_{v}_2020.csv", O / f"loo360_{v}_anomaly" / f"daily_summary_{v}_2020.csv"):
        chk(f"{v}: {f.relative_to(O) if str(f).startswith(str(O)) else f}", f.exists())
    m = O / "summary" / f"loo_metrics_{v}_anomaly.csv"
    if m.exists():
        t = pd.read_csv(m); a_ = t[t.group == "all"].iloc[0]; chk(f"{v}: LOO skill (all)", True, f"MAE {a_.mae:.3f} r {a_.r:.3f} KGE {a_.kge:.3f} bias {a_.bias:.3f}")
R = pd.DataFrame(rows); (O / "summary").mkdir(exist_ok=True); R.to_csv(O / "summary" / "final_check.csv", index=False)
(O / "summary" / "final_check.md").write_text("| item | status | detail |\n|---|---|---|\n" + "\n".join(f"| {r.item} | {r.status} | {r.detail} |" for r in R.itertuples()))
print(f"\n{(R.status == 'PASS').sum()} PASS, {(R.status == 'WARN').sum()} WARN, {(R.status == 'FAIL').sum()} FAIL"); sys.exit(int((R.status == "FAIL").any()))
