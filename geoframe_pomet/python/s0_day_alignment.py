#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 0 - detection and harmonisation of the daily labelling convention of every station.

The merged station files mix networks whose daily values are labelled by different days: for many
stations the accumulation (precipitation) or mean (temperature) of day D is stored under D+1 (value
labelled by the end of its interval). This is detected per station and year from the correlation with
independent daily gridded references at lags -1/0/+1 day, using
    precipitation: EMO-1arcmin (06-06 UTC, labelled by the start day)  - daily values
    temperature:   EMO-1arcmin ta6 AND E-OBS v31.0e tg (both must agree) - day-to-day changes
                   (differencing removes the seasonal cycle and persistence, which hide a 1-day shift)
A station-year is 'later' (shift -1 day), 'earlier' (+1) or 'aligned' when that lag beats the others by
GAIN; ambiguous years take the class of the nearest classified year of the same station. Stations
without a reference classification are tested against the harmonised mean of their 5 nearest classified
neighbours (within 25 km). The harmonised convention is the calendar day 00-24 (the 'aligned' group).

Usage:  python s0_day_alignment.py --variable precip|temp
Output: <OUT_ROOT>/harmonisation/day_shift_<var>.csv              station, year, shift_days, source
        <OUT_ROOT>/harmonisation/day_alignment_station_summary_<var>.csv
        <OUT_ROOT>/harmonisation/neighbour_pair_check_<var>.csv   before/after check (no reference used)
Every later step applies the shift automatically when reading the raw CSV (po_common.load_station_data).
"""
import argparse, os, sys
from pathlib import Path
os.environ["PO_RAW"] = "1"; os.environ["PO_QC"] = "0"; os.environ["PO_SHIFT"] = "0"     # detect on raw, unshifted data
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
import po_config as C, po_common as P
import s8_reference_comparison as R

GAIN = {"precip": 0.15, "temp": 0.20}


def corr(A, B, minn):
    ok = np.isfinite(A) & np.isfinite(B); n = ok.sum(0); A = np.where(ok, A, 0); B = np.where(ok, B, 0); k = np.maximum(n, 1)
    ma = A.sum(0) / k; mb = B.sum(0) / k; c = (A * B).sum(0) / k - ma * mb; va = (A * A).sum(0) / k - ma ** 2; vb = (B * B).sum(0) / k - mb ** 2
    with np.errstate(all="ignore"): return np.where(n >= minn, c / np.sqrt(va * vb), np.nan)


def lag_r(O, F, v):
    """r at lag 0, station-later (station day t ~ reference day t-1), station-earlier."""
    if v == "temp": O = np.diff(O, axis=0); F = np.diff(F, axis=0)
    m = 150 if v == "temp" else 200
    return corr(O, F, m), corr(O[1:], F[:-1], m), corr(O[:-1], F[1:], m)


def classify(r0, rl, re, g):
    return np.select([rl - r0 > g, re - r0 > g, r0 - np.maximum(rl, re) > g], ["later", "earlier", "aligned"], "ambiguous")


def reference_classes(v, st, product):
    from pyproj import Transformer
    cfg = R.PRODUCTS[v][product]; tr = Transformer.from_crs(C.METRIC_CRS, cfg["crs"], always_xy=True)
    px, py = tr.transform(st["xy"][:, 0], st["xy"][:, 1]); tb = Transformer.from_crs("EPSG:4326", cfg["crs"], always_xy=True)
    bx, by = tb.transform([6.0, 6.0, 14.5, 14.5], [43.5, 47.2, 43.5, 47.2]); bbox = (min(bx), min(by), max(bx), max(by))
    obs = pd.DataFrame(st["values"], index=st["dates"]); rows = []
    for y in range(st["dates"].year.min(), st["dates"].year.max() + 1):
        r = R.read_year(cfg, y, bbox)
        if r is None: continue
        t, a, x, yy = r; jx, okx, _ = R.nearest_index(x, px); jy, oky, _ = R.nearest_index(yy, py)
        F = a[:, jy, jx].astype(float); F[:, ~(okx & oky)] = np.nan; O = obs.reindex(t).to_numpy().astype(float)
        r0, rl, re = lag_r(O, F, v)
        rows.append(pd.DataFrame(dict(station=st["ids"], year=y, r0=r0, rl=rl, re=re)))
        print(f"  {product} {y}", flush=True)
    D = pd.concat(rows).dropna(subset=["r0", "rl", "re"]); D["cls"] = classify(D.r0, D.rl, D.re, GAIN[v]); return D


def fill_years(cls_by_year, years):
    """Ambiguous / missing years take the class of the nearest classified year."""
    known = {y: c for y, c in cls_by_year.items() if c in ("later", "earlier", "aligned")}
    if not known: return {}
    ky = np.array(sorted(known)); return {y: known[int(ky[np.argmin(np.abs(ky - y))])] for y in years}


def apply_shift_matrix(values, dates, ids, table):
    """Same mapping as po_common.apply_day_shift, local copy for the checks."""
    return P.apply_day_shift(values, dates, ids, table)


def pair_check(v, values, xy, radius=10000.0):
    X = np.diff(values, axis=0) if v == "temp" else values
    pairs = cKDTree(xy).query_pairs(radius, output_type="ndarray"); out = []
    for i, j in pairs:
        a, b = X[:, i], X[:, j]; ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 500: continue
        r0 = np.corrcoef(a[ok], b[ok])[0, 1]
        ok1 = np.isfinite(a[1:]) & np.isfinite(b[:-1]); ok2 = np.isfinite(a[:-1]) & np.isfinite(b[1:])
        r1 = np.corrcoef(a[1:][ok1], b[:-1][ok1])[0, 1]; r2 = np.corrcoef(a[:-1][ok2], b[1:][ok2])[0, 1]
        out.append((r0, max(r1, r2)))
    o = np.array(out); return dict(pairs=len(o), median_r_same_day=float(np.median(o[:, 0])), pct_pairs_offset=float(100 * np.mean(o[:, 1] > o[:, 0])))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--variable", choices=["precip", "temp"], required=True)
    ap.add_argument("--reuse-classes", action="store_true", help="reuse station_year_classes_<var>.parquet (skip reading references)")
    a = ap.parse_args(); v = a.variable; od = C.DATA_ROOT / "harmonisation"; od.mkdir(parents=True, exist_ok=True)
    st = P.attach_station_predictors(P.load_station_data(v), P.load_grid()); ids = st["ids"]; years = list(range(st["dates"].year.min(), st["dates"].year.max() + 1))
    if a.reuse_classes and (od / f"station_year_classes_{v}.parquet").exists():
        D = pd.read_parquet(od / f"station_year_classes_{v}.parquet"); D["station"] = D["station"].astype(str)
    elif v == "precip":
        D = reference_classes(v, st, "EMO"); D["source"] = "EMO"
    else:
        D1 = reference_classes(v, st, "EMO"); D2 = reference_classes(v, st, "E-OBS")
        D = D1.merge(D2[["station", "year", "cls"]], on=["station", "year"], how="outer", suffixes=("", "_eobs"))
        D["cls"] = np.where(D.cls == D.cls_eobs, D.cls, "ambiguous"); D["source"] = "EMO+E-OBS"
    if not a.reuse_classes: D.to_parquet(od / f"station_year_classes_{v}.parquet", index=False)
    rows = []; summ = []
    by = {s: dict(zip(g.year, g.cls)) for s, g in D.groupby("station")}
    unresolved = []
    for s in ids:
        f = fill_years(by.get(s, {}), years)
        if not f: unresolved.append(s); continue
        for y, c in f.items(): rows.append(dict(station=s, year=y, shift_days={"later": -1, "earlier": 1, "aligned": 0}[c], source="reference"))
        cs = pd.Series(list(by.get(s, {}).values())); summ.append(dict(station=s, n_years_classified=int(cs.isin(["later", "earlier", "aligned"]).sum()),
                                                                    n_later=int((cs == "later").sum()), n_aligned=int((cs == "aligned").sum()), n_earlier=int((cs == "earlier").sum()), source="reference"))
    T = pd.DataFrame(rows)
    # stations without a reference classification: iterative consistency with harmonised neighbours.
    # The lag with the highest correlation to the mean of the 5 nearest stations (within 40 km) is taken when
    # it beats lag 0 by >= 0.05; repeated until no station changes (neighbours may themselves be unresolved).
    ref_cls = set(T.station.unique()); open_ = [s for s in ids if s not in ref_cls]; cur = {s: 0 for s in open_}
    tree = cKDTree(st["xy"]); X = st["values"].astype(float)
    # initial guess: neighbours restricted to reference-classified stations
    Hk = P.apply_day_shift(X, st["dates"], ids, T); kid = np.flatnonzero(np.isin(ids, list(ref_cls))); tk = cKDTree(st["xy"][kid])
    for s in open_:
        i = int(np.flatnonzero(ids == s)[0]); d, nb = tk.query(st["xy"][i], k=5); nb = kid[nb[d <= 40000]]
        if len(nb) < 2: continue
        with np.errstate(all="ignore"): ref = np.nanmean(Hk[:, nb], 1)[:, None]
        r0, rl, re = (float(z[0]) for z in lag_r(X[:, [i]], ref, v))
        if np.isfinite(r0):
            best = max([(r0, 0), (rl, -1), (re, 1)]); cur[s] = best[1] if best[0] - r0 >= 0.05 else 0
    for it in range(5):
        Tall = pd.concat([T, pd.DataFrame([dict(station=s, year=y, shift_days=cur[s]) for s in open_ for y in years])], ignore_index=True) if open_ else T
        Hv = P.apply_day_shift(X, st["dates"], ids, Tall); changed = 0
        for s in open_:
            i = int(np.flatnonzero(ids == s)[0]); d, nb = tree.query(st["xy"][i], k=6); nb = nb[(d <= 40000) & (nb != i)][:5]
            if len(nb) < 2: continue
            with np.errstate(all="ignore"): ref = np.nanmean(Hv[:, nb], 1)[:, None]
            r0, rl, re = (float(z[0]) for z in lag_r(X[:, [i]], ref, v))
            if not np.isfinite(r0): continue
            best = max([(r0, 0), (rl, -1), (re, 1)]); sh = best[1] if best[0] - r0 >= 0.05 else 0
            if sh != cur[s]: cur[s] = sh; changed += 1
            if it == 4 or True: summ_n = (s, r0, rl, re)
        print(f"  neighbour iteration {it + 1}: {changed} changes", flush=True)
        if changed == 0: break
    for s in open_:
        T = pd.concat([T, pd.DataFrame(dict(station=s, year=years, shift_days=cur[s], source="neighbours"))], ignore_index=True)
        summ.append(dict(station=s, source=f"neighbours (shift {cur[s]})"))
    unresolved = open_
    T.to_csv(od / f"day_shift_{v}.csv", index=False); S = pd.DataFrame(summ); S.to_csv(od / f"day_alignment_station_summary_{v}.csv", index=False)
    stn = T.groupby("station").shift_days.agg(lambda x: "/".join(map(str, sorted(set(x)))))
    print(f"{v}: stations {len(ids)} | shift pattern: {stn.value_counts().to_dict()} | resolved by neighbours: {len(unresolved)}")
    # independent check on the station data alone: neighbouring pairs before / after
    before = pair_check(v, st["values"].astype(float), st["xy"]); after = pair_check(v, P.apply_day_shift(st["values"].astype(float), st["dates"], ids, T), st["xy"])
    chk = pd.DataFrame([dict(state="raw", **before), dict(state="harmonised", **after)]); chk.to_csv(od / f"neighbour_pair_check_{v}.csv", index=False)
    print(chk.round(3).to_string(index=False)); print(f"-> {od}")


if __name__ == "__main__":
    main()
