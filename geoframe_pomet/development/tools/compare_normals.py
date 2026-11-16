#!/usr/bin/env python3
"""Java PoNormals vs Python normals (reference_run/normals_rebuild/normals): max relative differences per file/column."""
import sys, json, numpy as np, pandas as pd, xarray as xr
from pathlib import Path
B = Path(__file__).resolve().parents[1]; P = B / "reference_run/normals_rebuild/normals"; J = Path(sys.argv[1] if len(sys.argv) > 1 else B / "output/normals_java")
res = {}
def cmp(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float); nanm = int((np.isnan(a) != np.isnan(b)).sum()); f = np.isfinite(a) & np.isfinite(b)
    rel = np.abs(a[f] - b[f]) / np.maximum(np.abs(b[f]), 1e-12)
    return {"n": int(f.sum()), "nan_mismatch": nanm, "max_rel": float(rel.max()) if f.any() else 0.0, "frac_rel_le_1e-9": float(np.mean(rel <= 1e-9)) if f.any() else 1.0,
            "frac_bitwise": float(np.mean(a[f] == b[f])) if f.any() else 1.0}
for f, cols in [("station_normals_precip.csv", ["normal", "normal_obs_only", "n_obs_years", "n_total_years"]), ("station_normals_temp.csv", ["normal", "normal_obs_only", "n_obs_years", "n_total_years"]),
                ("station_normals_loo_precip.csv", ["normal_station", "normal_grid_loo"]), ("station_normals_loo_temp.csv", ["normal_station", "normal_grid_loo"]),
                ("station_wetfreq_loo_precip.csv", ["wetfreq_station", "wetfreq_grid_loo", "n_months"])]:
    a = pd.read_csv(J / f, dtype={"station": str}); b = pd.read_csv(P / f, dtype={"station": str}); m = b.merge(a, on=["station", "month"], suffixes=("_py", "_j"))
    res[f] = {"rows": [len(b), len(a), len(m)], **{c: cmp(m[c + "_j"], m[c + "_py"]) for c in cols}}
for f, v in [("normals_precip.nc", "normal"), ("normals_precip.nc", "elevation_slope"), ("normals_temp.nc", "normal"), ("normals_temp.nc", "elevation_slope"), ("wetfreq_precip.nc", "wetfreq"), ("wetfreq_precip.nc", "elevation_slope")]:
    with xr.open_dataset(J / f) as a, xr.open_dataset(P / f) as b: res[f + ":" + v] = cmp(a[v].values.ravel(), b[v].values.ravel())
print(json.dumps(res, indent=1))
