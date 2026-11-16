#!/usr/bin/env python3
"""Compare a Java uncertainty table (PoUncertainty) with the Python s9 table: row keys, max abs/rel difference per
column, number of bit-identical values. Usage: compare_uncertainty.py <java.csv> <python.csv>"""
import sys, numpy as np, pandas as pd
a, b = pd.read_csv(sys.argv[1], dtype={"cls": str}), pd.read_csv(sys.argv[2], dtype={"cls": str})
keys_equal = (a.season + "|" + a.cls).tolist() == (b.season + "|" + b.cls).tolist()
print(f"rows java {len(a)} python {len(b)}; season/cls order identical: {keys_equal}")
for c in [c for c in b.columns if c not in ("variable", "season", "cls")]:
    x, y = a[c].to_numpy(float), b[c].to_numpy(float); d = np.abs(x - y); nanm = np.isnan(x) == np.isnan(y)
    rel = np.nanmax(d / np.maximum(np.abs(y), 1e-300)) if np.isfinite(d).any() else 0
    print(f"  {c:24s} max|d| {np.nanmax(d) if np.isfinite(d).any() else 0:.3e}  max rel {rel:.3e}  identical {int(np.sum((x == y) | (np.isnan(x) & np.isnan(y))))}/{len(x)}  NaN pattern equal {nanm.all()}")
