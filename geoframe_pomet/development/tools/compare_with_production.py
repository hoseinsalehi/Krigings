#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Check the Python reference (reference/daily/*) against the production files of 2020 and measure the PSO seed
sensitivity (reference/seed_sensitivity/*). Read-only on PoGrid_v4 / PoGrid_v5.

Production files: v2.0 ratio + temperature = ~/Downloads/PoGrid_v4/{precip_anomaly,temp_anomaly};
two-part = ~/Downloads/PoGrid_v5/precip_anomaly (pass-1 fields in _twopart/grid_2020.npz, loo_2020.parquet).
Writes reference/check_vs_production.csv and reference/seed_sensitivity.csv.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

B = Path.home() / "Documents/Po_precip/GEOframe_PoKriging/reference"
V4, V5 = Path.home() / "Downloads/PoGrid_v4", Path.home() / "Downloads/PoGrid_v5"
rows = []


def stat(name, field, a, b, quantum=None):
    a = np.asarray(a, float); b = np.asarray(b, float); both = np.isfinite(a) & np.isfinite(b)
    d = np.abs(a - b)[both]
    r = dict(check=name, field=field, n=int(both.size), n_both_finite=int(both.sum()),
             nan_pattern_equal=bool(np.array_equal(np.isfinite(a), np.isfinite(b))),
             max_abs_diff=float(d.max()) if d.size else np.nan, mean_abs_diff=float(d.mean()) if d.size else np.nan,
             p99_abs_diff=float(np.percentile(d, 99)) if d.size else np.nan,
             frac_exact=float((d == 0).mean()) if d.size else np.nan)
    if quantum: r["frac_within_1_quantum"] = float((d <= quantum * 1.0001).mean()) if d.size else np.nan
    rows.append(r); return r


def q(x, sc): return np.rint(np.asarray(x, float) / sc) * sc


def nc_at(path, var, dates, iy, ix):
    with xr.open_dataset(path) as ds:
        da = ds[var].sel(time=pd.to_datetime(dates)).values            # (nd, ny, nx), decoded
    return da[:, iy, ix]


def main():
    dates = pd.read_json(B / "daily/precip_ratio/run_info.json", typ="series")["dates"]
    cells = pd.read_csv(B / "daily/targets_cells.csv"); iy, ix, cid = cells.iy.to_numpy(), cells.ix.to_numpy(), cells.cell.to_numpy()
    # ---------------------------------------------------------------- ratio (v2.0) and temperature anomaly: grid
    for m, var, nc, sc in (("precip_ratio", "precipitation", V4 / "precip_anomaly/precip_2020.nc", 0.1),
                           ("temp_anomaly", "temperature", V4 / "temp_anomaly/temp_2020.nc", 0.01)):
        g = pd.read_csv(B / f"daily/{m}/grid.csv.gz")
        for fld, pv in (("value", var), ("kstd", "kriging_std")):
            ref = g.pivot(index="date", columns="cell", values=fld).reindex(index=dates, columns=cid).to_numpy()
            prod = nc_at(nc, pv, dates, iy, ix)
            stat(f"{m} grid vs production (quantised {sc})", fld, q(ref, sc), prod, sc)
        diag = nc.parent / f"diagnostics_{nc.stem.split('_')[0]}_2020.nc"
        ref = g.pivot(index="date", columns="cell", values="range").reindex(index=dates, columns=cid).to_numpy()
        stat(f"{m} grid vs production (quantised 0.5 km)", "variogram range", q(ref, 0.5), nc_at(diag, "range", dates, iy, ix), 0.5)
        # LOO
        L = pd.read_csv(B / f"daily/{m}/loo.csv.gz", dtype={"station": str})
        Pq = pd.read_parquet(nc.parent / f"loo_{nc.stem.split('_')[0]}_2020.parquet"); Pq["station"] = Pq.station.astype(str)
        Pq = Pq[Pq.date.astype(str).isin(dates)]
        J = L.merge(Pq, on=["date", "station"], suffixes=("_ref", "_prod"), how="outer", indicator=True)
        rows.append(dict(check=f"{m} LOO vs production", field="row match", n=len(J), n_both_finite=int((J._merge == "both").sum())))
        for fld in ("pred", "kriging_std"): stat(f"{m} LOO vs production", fld, J[f"{fld}_ref"], J[f"{fld}_prod"])
    # ---------------------------------------------------------------- two-part: pass-1 grid fields, LOO, final values
    g = pd.read_csv(B / "daily/precip_twopart/grid.csv.gz")
    Z = np.load(V5 / "precip_anomaly/_twopart/grid_2020.npz"); zd = list(Z["days"]); di = [zd.index(d) for d in dates]
    for fld, sc in (("pi", 1e-4), ("amount", 0.1), ("samount", 0.1)):
        ref = g.pivot(index="date", columns="cell", values=fld).reindex(index=dates, columns=cid).to_numpy()
        code = Z[fld][di][:, cid]; fill = {"pi": 65535, "amount": -32768, "samount": -32768}[fld]
        prod = np.where(code == fill, np.nan, code.astype(float) * sc)
        stat(f"precip_twopart pass-1 grid vs production (quantised {sc})", fld, q(ref, sc), prod, sc)
    ref = g.pivot(index="date", columns="cell", values="range").reindex(index=dates, columns=cid).to_numpy()
    stat("precip_twopart grid vs production (quantised 0.5 km)", "occurrence variogram range", q(ref, 0.5),
         nc_at(V5 / "precip_anomaly/diagnostics_precip_2020.nc", "range", dates, iy, ix), 0.5)
    # final value: reference = pooled 10-day calibration; production = monthly 1991-2020 calibration (not expected equal)
    ref = g.pivot(index="date", columns="cell", values="value").reindex(index=dates, columns=cid).to_numpy()
    prod = nc_at(V5 / "precip_anomaly/precip_2020.nc", "precipitation", dates, iy, ix)
    stat("precip_twopart FINAL grid vs production (different calibration sample: 10-day pooled vs 30-yr monthly)", "value", q(ref, 0.1), prod, 0.1)
    rw, pw = ref >= 0.05, prod >= 0.05; ok = np.isfinite(ref) & np.isfinite(prod)
    rows.append(dict(check="precip_twopart FINAL wet/dry agreement (10-day pooled vs 30-yr monthly calibration)", field="wet flag",
                     n=int(ok.sum()), frac_exact=float((rw == pw)[ok].mean()), ref_wet_frac=float(rw[ok].mean()), prod_wet_frac=float(pw[ok].mean())))
    L = pd.read_csv(B / "daily/precip_twopart/loo.csv.gz", dtype={"station": str})
    Pq = pd.read_parquet(V5 / "precip_anomaly/_twopart/loo_2020.parquet"); Pq["station"] = Pq.station.astype(str)
    Pq = Pq[Pq.date.astype(str).isin(dates)]
    J = L.merge(Pq, on=["date", "station"], suffixes=("_ref", "_prod"), how="outer", indicator=True)
    rows.append(dict(check="precip_twopart pass-1 LOO vs production", field="row match", n=len(J), n_both_finite=int((J._merge == "both").sum())))
    for fld in ("pi", "amount", "samount"): stat("precip_twopart pass-1 LOO vs production", fld, J[f"{fld}_ref"], J[f"{fld}_prod"])
    out = pd.DataFrame(rows); out.to_csv(B / "check_vs_production.csv", index=False, float_format="%.6g")
    with pd.option_context("display.width", 250, "display.max_columns", 20): print(out.to_string(index=False))
    # ---------------------------------------------------------------- PSO seed sensitivity
    srows = []
    for m, flds in (("precip_ratio", ("value", "kstd", "range")), ("temp_anomaly", ("value", "kstd", "range")),
                    ("precip_twopart", ("pi", "amount", "value", "range"))):
        g0 = pd.read_csv(B / f"daily/{m}/grid.csv.gz"); L0 = pd.read_csv(B / f"daily/{m}/loo.csv.gz", dtype={"station": str})
        for s in (1, 2):
            d = B / f"seed_sensitivity/{m}_salt{s}"
            g1 = pd.read_csv(d / "grid.csv.gz"); L1 = pd.read_csv(d / "loo.csv.gz", dtype={"station": str})
            J = g0.merge(g1, on=["date", "cell"], suffixes=("_0", "_1"))
            for f in flds:
                a, b = J[f"{f}_0"].to_numpy(float), J[f"{f}_1"].to_numpy(float); ok = np.isfinite(a) & np.isfinite(b); dd = np.abs(a - b)[ok]
                srows.append(dict(method=m, salt=s, target="grid", field=f, n=int(ok.sum()), frac_identical=float((dd == 0).mean()),
                                  median_abs_diff=float(np.median(dd)), p95_abs_diff=float(np.percentile(dd, 95)),
                                  p99_abs_diff=float(np.percentile(dd, 99)), max_abs_diff=float(dd.max()),
                                  mean_value=float(np.nanmean(np.abs(a)))))
            if m == "precip_twopart":
                J = g0.merge(g1, on=["date", "cell"], suffixes=("_0", "_1"))
                srows.append(dict(method=m, salt=s, target="grid", field="wet flag", n=len(J), frac_identical=float((J.wet_0 == J.wet_1).mean())))
            lf = "pred"
            J = L0.merge(L1, on=["date", "station"], suffixes=("_0", "_1")); a, b = J[f"{lf}_0"].to_numpy(float), J[f"{lf}_1"].to_numpy(float)
            ok = np.isfinite(a) & np.isfinite(b); dd = np.abs(a - b)[ok]; o = J.obs_0.to_numpy(float)[ok]
            srows.append(dict(method=m, salt=s, target="loo", field=lf, n=int(ok.sum()), frac_identical=float((dd == 0).mean()),
                              median_abs_diff=float(np.median(dd)), p95_abs_diff=float(np.percentile(dd, 95)), p99_abs_diff=float(np.percentile(dd, 99)),
                              max_abs_diff=float(dd.max()), mean_value=float(np.mean(np.abs(a[ok]))),
                              loo_mae_salt0=float(np.mean(np.abs(a[ok] - o))), loo_mae_salt=float(np.mean(np.abs(b[ok] - o)))))
    S = pd.DataFrame(srows); S.to_csv(B / "seed_sensitivity.csv", index=False, float_format="%.6g")
    with pd.option_context("display.width", 250, "display.max_columns", 20): print(S.to_string(index=False))


if __name__ == "__main__":
    main()
