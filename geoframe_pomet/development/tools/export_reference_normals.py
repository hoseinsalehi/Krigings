#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Normals references for the Java equivalence tests (from the rebuilt normals, bit-identical to PoGrid_v5/normals).
Writes reference/normals/:
  grid_normals_cells.csv    cell,month, precip_M_mm_month, precip_lnslope_per_km, wetfreq_p, wetfreq_logitslope_per_km,
                            temp_normal_c, temp_lapse_c_per_km   (reference target cells, 12 months)
  station_*.csv             station normals (s1), LOO normals (s2), wet-day frequency + LOO (s1w), copied
  station_monthly_<var>_{observed,filled}.csv.gz   s1 monthly matrices (station x 360 months), for the gap-fill test
  params_*.json             locked normals parameters
"""
import shutil
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

B = Path.home() / "Documents/Po_precip/GEOframe_PoKriging"; R = B / "reference_run/normals_rebuild/normals"; O = B / "reference/normals"
O.mkdir(parents=True, exist_ok=True)
cells = pd.read_csv(B / "reference/daily/targets_cells.csv"); iy, ix = cells.iy.to_numpy(), cells.ix.to_numpy()
cols = {}
for f, v, name in (("normals_precip.nc", "normal", "precip_M_mm_month"), ("normals_precip.nc", "elevation_slope", "precip_lnslope_per_km"),
                   ("wetfreq_precip.nc", "wetfreq", "wetfreq_p"), ("wetfreq_precip.nc", "elevation_slope", "wetfreq_logitslope_per_km"),
                   ("normals_temp.nc", "normal", "temp_normal_c"), ("normals_temp.nc", "elevation_slope", "temp_lapse_c_per_km")):
    with xr.open_dataset(R / f) as ds: cols[name] = ds[v].values[:, iy, ix].T.ravel()      # (cells*12) cell-major
T = pd.DataFrame(dict(cell=np.repeat(cells.cell.to_numpy(), 12), month=np.tile(np.arange(1, 13), len(cells)), **cols))
T.to_csv(O / "grid_normals_cells.csv", index=False, float_format="%.9g")
for f in ("station_normals_precip.csv", "station_normals_temp.csv", "station_normals_loo_precip.csv", "station_normals_loo_temp.csv",
          "station_wetfreq_loo_precip.csv", "params_precip.json", "params_temp.json"):
    shutil.copy2(R / f, O / f)
shutil.copy2(Path.home() / "Downloads/PoGrid_v5/normals/params_wetfreq_precip.json", O / "params_wetfreq_precip.json")
for v in ("precip", "temp"):
    z = np.load(R / f"station_monthly_{v}.npz")
    for k in ("observed", "filled"):
        pd.DataFrame(z[k], index=pd.Index(z["ids"], name="station"), columns=z["months"]).to_csv(O / f"station_monthly_{v}_{k}.csv.gz", float_format="%.7g")
print(T.describe().T[["count", "mean", "min", "max"]].to_string())
