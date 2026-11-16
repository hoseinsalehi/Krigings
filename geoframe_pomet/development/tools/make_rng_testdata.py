#!/usr/bin/env python3
"""Reference draws of numpy.random.default_rng for the Java RNG unit test (tests/data/rng_numpy.csv) and
po_core.seed_for values (tests/data/seed_for.csv)."""
import hashlib, numpy as np, sys
from pathlib import Path
out = Path(sys.argv[1] if len(sys.argv) > 1 else "tests/data"); out.mkdir(parents=True, exist_ok=True)
def seed_for(*parts):
    h = hashlib.blake2b(digest_size=8)
    for p in parts: h.update(str(p).encode())
    return int.from_bytes(h.digest(), "little") & 0xFFFFFFFF
seeds = [0, 1, 42, 12345, 2**31 - 1, 4294967295, seed_for("precip|anomaly|2020-10-02", 4096), seed_for("temp|direct|2020-01-05", "global")]
with open(out / "rng_numpy.csv", "w") as f:
    f.write("seed,kind,i,value\n")
    for s in seeds:
        g = np.random.default_rng(s); u = g.random(1000)
        g = np.random.default_rng(s); z = g.standard_normal(1000)
        g = np.random.default_rng(s); mixed = np.concatenate([g.random(30), g.normal(0, .06, 30), g.random(40)])   # draw order of the PSO
        for i, v in enumerate(u): f.write(f"{s},random,{i},{v.hex()}\n")
        for i, v in enumerate(z): f.write(f"{s},normal,{i},{v.hex()}\n")
        for i, v in enumerate(mixed): f.write(f"{s},mixed,{i},{v.hex()}\n")
tags = [("precip|anomaly|2020-10-02", 0), ("precip|anomaly|2020-10-02", 4096), ("temp|direct|2020-01-05", "global"),
        ("precip|twopart|2020-07-04|occ", 18432), ("precip|twopart|2020-07-04|amt", "global"), ("x", "y", "salt1")]
with open(out / "seed_for.csv", "w") as f:
    f.write("concatenated_parts,seed\n")
    for t in tags: f.write("".join(map(str, t)) + "," + str(seed_for(*t)) + "\n")
print("written", out)
