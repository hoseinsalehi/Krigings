#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 6 - point test between Aosta and Gran San Bernardo (reviewer example).

Locates the precipitation/temperature stations closest to Aosta (45.737 N, 7.320 E, ~580 m)
and to the Gran San Bernardo pass (45.869 N, 7.171 E, ~2470 m), then produces 1991-2020
daily series at test points between them under four networks:
  full | without Aosta | without Gran San Bernardo | without both
for each method. Points cost almost nothing (no grid), so the full 30 years are used.

  python s6_point_test.py --variable precip --modes direct,anomaly --workers 12
Output: <OUT_ROOT>/point_test_<var>/stations_used.csv, points.csv, series_<mode>_<scenario>.csv,
        summary.csv (annual mean, p99, annual maxima, and differences w.r.t. the full network).
"""
import argparse, os, shutil, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from pyproj import Transformer
import po_config as C, po_common as P

HERE = Path(__file__).resolve().parent
SITES = {"Aosta": (7.320, 45.737), "Gran_San_Bernardo": (7.171, 45.869)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--variable", choices=["precip", "temp"], required=True)
    ap.add_argument("--modes", default="direct,anomaly"); ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--station-ids", default=None, help="override: 'AOSTA_ID,GSB_ID'")
    ap.add_argument("--years", default=None, help="e.g. 2020 for a quick check (default 1991-2020)")
    a = ap.parse_args(); v = a.variable; root = C.OUT_ROOT / f"point_test_{v}"; root.mkdir(parents=True, exist_ok=True)
    st = P.load_station_data(v); tr = Transformer.from_crs("EPSG:4326", C.METRIC_CRS, always_xy=True)
    used = []
    for k, (lon, lat) in SITES.items():
        x, y = tr.transform(lon, lat); d = np.hypot(st["xy"][:, 0] - x, st["xy"][:, 1] - y) / 1000; i = int(np.argmin(d))
        used.append(dict(site=k, station=st["ids"][i], distance_km=d[i], elevation=st["elev"][i], valid_days_1991_2020=int(np.isfinite(st["values"][:, i]).sum())))
    if a.station_ids:
        for u, s in zip(used, a.station_ids.split(",")): u["station"] = s
    U = pd.DataFrame(used); U.to_csv(root / "stations_used.csv", index=False); print(U.to_string(index=False))
    # test points: 1/4, 1/2, 3/4 of the way from Aosta to Gran San Bernardo
    (lo0, la0), (lo1, la1) = SITES["Aosta"], SITES["Gran_San_Bernardo"]
    pts = pd.DataFrame([dict(name=f"AO_GSB_{int(f*100)}", lon=lo0 + f * (lo1 - lo0), lat=la0 + f * (la1 - la0)) for f in (0.25, 0.5, 0.75, 1.0)])
    pts.to_csv(root / "points.csv", index=False)
    ao, gs = U.station.tolist()
    scen = {"full": [], "no_Aosta": [ao], "no_GSB": [gs], "no_both": [ao, gs]}
    out = []
    for m in a.modes.split(","):
        for s, ex in scen.items():
            exf = root / f"exclude_{s}.txt"; exf.write_text("\n".join(ex))
            # without a station the 1991-2020 normals must also be rebuilt without it (anomaly mode), otherwise the
            # withheld station still informs the climatological layer: each scenario gets its own normals
            sub = root / f"scenario_{s}"; sub.mkdir(exist_ok=True); env = dict(os.environ, PO_OUT_ROOT=str(sub), PO_DATA_ROOT=str(C.DATA_ROOT))
            if not (sub / "grid_predictors.npz").exists(): shutil.copy(C.OUT_ROOT / "grid_predictors.npz", sub / "grid_predictors.npz")
            if not (sub / "normals" / f"normals_{v}.nc").exists():
                prm = C.OUT_ROOT / "normals" / f"params_{v}.json"
                for script, extra in (("s1_station_normals.py", (["--exclude-stations", str(exf)] if ex else [])), ("s2_normals_grid.py", ["--params", str(prm)])):
                    r = subprocess.run([sys.executable, str(HERE / script), "--variable", v] + extra, capture_output=True, text=True, env=env)
                    if r.returncode: raise RuntimeError(r.stdout[-1500:] + r.stderr[-1500:])
            if v == "precip" and C.PRECIP_METHOD == "twopart" and not (sub / "normals" / "wetfreq_precip.nc").exists():
                r = subprocess.run([sys.executable, str(HERE / "s1w_wetfreq_normals.py")] + (["--exclude-stations", str(exf)] if ex else []),
                                   capture_output=True, text=True, env=env)
                if r.returncode: raise RuntimeError(r.stdout[-1500:] + r.stderr[-1500:])
            tag = f"{m}_{s}"
            cmd = [sys.executable, str(HERE / "s3_run.py"), "--variable", v, "--mode", m, "--loo-only", "--no-loo", "--points", str(root / "points.csv"),
                   "--workers", str(a.workers), "--tag", tag] + (["--exclude-stations", str(exf)] if ex else []) + (["--years", a.years] if a.years else [])
            r = subprocess.run(cmd, capture_output=True, text=True, env=env)
            if r.returncode: raise RuntimeError(r.stdout[-1500:] + r.stderr[-1500:])
            ser = pd.concat([pd.read_csv(f) for f in sorted((sub / tag).glob("points_*.csv"))]); ser["mode"], ser["scenario"] = m, s
            ser.to_csv(root / f"series_{m}_{s}.csv", index=False); out.append(ser); print(f"{m} {s}: {len(ser)} rows", flush=True)
    A = pd.concat(out); A["year"] = pd.to_datetime(A.date).dt.year
    agg = (lambda g: pd.Series(dict(annual_total=g.groupby("year").value.sum().mean(), p99=g.value.quantile(.99), mean_annual_max=g.groupby("year").value.max().mean()))) \
        if v == "precip" else (lambda g: pd.Series(dict(mean=g.value.mean(), p01=g.value.quantile(.01), p99=g.value.quantile(.99))))
    S = A.groupby(["mode", "name", "scenario"]).apply(agg).reset_index()
    full = S[S.scenario == "full"].drop(columns="scenario").set_index(["mode", "name"])
    for c in [c for c in S.columns if c not in ("mode", "name", "scenario")]:
        S[f"{c}_minus_full"] = S[c].to_numpy() - full.loc[list(zip(S["mode"], S["name"])), c].to_numpy()
    S.to_csv(root / "summary.csv", index=False); print(S.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
