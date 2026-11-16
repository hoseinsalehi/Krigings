#!/usr/bin/env python3
"""Check that the numba kernels give bit-identical results to the numpy kernels on real data.
   python test_numba.py            (~30 s; prints True for every quantity)"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.spatial import cKDTree
import po_common as P, po_core as K, po_numba as N

grid = P.load_grid()
for var in ("precip", "temp"):
    st = P.attach_station_predictors(P.load_station_data(var), grid); sd = P.station_distance_km(st["xy"]); ok = True
    for date in ("2020-10-03", "2020-01-15", "1993-07-20"):
        di = int(np.flatnonzero(st["dates"] == date)[0]); y = st["values"][di].astype(float); act = np.flatnonzero(np.isfinite(y))
        idx = np.random.default_rng(0).choice(len(grid["xy"]), 2048, replace=False)
        _, il = cKDTree(st["xy"][act]).query(grid["xy"][idx], k=50); n_ = act[il]
        Y = y[n_]; D = sd[n_[:, :, None], n_[:, None, :]]
        A = K.empirical_batch_np(Y, D); B = N.empirical_batch(Y, D)
        pa = K.pso_spherical_np(*A[:4], 7); pb = N.pso_spherical(*B[:4], 7)
        ok &= all(np.array_equal(a, b, equal_nan=True) for a, b in zip(A + pa, B + pb))
    print(var, "numba == numpy (bit-identical):", ok)
