#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Station-day quality control of the daily observations and writing of cleaned CSV files.

Tests (all on the raw CSV; flagged station-days are set to the missing code in the cleaned CSV):
  common   R  physical range / invalid token (precip < 0, temp outside TEMP_RANGE, non-numeric)
           S  stuck sensor: the same value on >= STUCK_DAYS consecutive days
              (temp: any value; precip: values >= 0.5 mm)
  temp     T  spatial consistency on the leave-one-out kriging error of a first run (pred - obs), after
              removing the station's median error of the same calendar month (its systematic
              representativeness bias); |debiased error| > QC_TEMP_ABS_C. Opposite-sign, smaller errors
              within QC_CONTAM_KM of a flagged station the same day are its neighbour contamination and
              are not flagged.
  precip   H  isolated high value: obs >= 30 mm and obs > 5 * (max of >= 3 neighbours within 25 km) + 20 mm
              and obs > 4 * LOO prediction + 20 mm
           A  accumulated total after missing days: day after >= 1 missing day, obs >= 30 mm and
              obs > 3 * neighbour max + 10 mm
           Z  false zero: obs == 0 while >= 3 neighbours within 15 km all report >= 10 mm and the LOO
              prediction is >= 15 mm
The neighbours are the 10 nearest stations (excluding the station itself) reporting that day.

Usage:  python s7_station_qc.py --variable temp   [--run temp_anomaly]
        python s7_station_qc.py --variable precip [--run precip_anomaly]
Output: <OUT_ROOT>/qc/station_day_flags_<var>.parquet  (date, station, obs, test, info)
        <OUT_ROOT>/qc/station_qc_summary_<var>.csv     (per station counts by test)
        <OUT_ROOT>/qc/<raw name>_QC.csv                  cleaned CSV, same layout and date format
Then set DATA[<var>]["obs"] in po_config.py to the cleaned CSV (see README) and rerun s1, s2, s3.
"""
import argparse, glob, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
import po_config as C, po_common as P

STUCK_DAYS = {"temp": 5, "precip": 4}
K_NB = 10


def load_raw(v):
    cfg = C.DATA[v]; raw = pd.read_csv(cfg["raw"], dtype=str, keep_default_na=False)        # always the raw file
    dates = pd.to_datetime(raw["Date"], format=cfg["date_format"], errors="coerce")
    raw = raw.loc[dates.notna()].reset_index(drop=True); dates = pd.DatetimeIndex(dates[dates.notna()])
    cols = [c for c in raw.columns if c != "Date"]
    num = np.array(raw[cols].apply(pd.to_numeric, errors="coerce").to_numpy(np.float64), copy=True)
    tok = (~np.isfinite(num)) & (raw[cols].to_numpy() != "")
    sh = P.load_day_shift(v)
    if sh is not None:                              # harmonise day labelling (s0); vacated days -> missing
        cols_s = [str(c) for c in cols]
        strs = P.apply_day_shift(raw[cols].to_numpy(object), dates, cols_s, sh, fill=str(int(C.MISSING_VALUE)))
        raw = raw.copy(); raw[cols] = strs
        num = np.array(pd.DataFrame(strs).apply(pd.to_numeric, errors="coerce").to_numpy(np.float64), copy=True)
        tok = (~np.isfinite(num)) & (strs != "")
    num[np.isclose(num, C.MISSING_VALUE)] = np.nan
    return raw, dates, cols, num, tok


def station_xy(v, cols):
    import geopandas as gpd
    g = gpd.read_file(C.DATA[v]["stations"]).to_crs(C.METRIC_CRS)
    lut = dict(zip(g[C.STATION_ID_FIELD].astype(str), np.c_[g.geometry.x, g.geometry.y]))
    return np.array([lut.get(str(c), (np.nan, np.nan)) for c in cols], float)


def stuck_mask(x, n, min_val=None):
    """True where x is part of a run of >= n identical consecutive finite values."""
    T, S = x.shape; out = np.zeros_like(x, bool)
    same = np.isfinite(x[1:]) & (x[1:] == x[:-1])
    if min_val is not None: same &= x[1:] >= min_val
    run = np.zeros(S, int)
    for t in range(1, T):
        run = np.where(same[t - 1], run + 1, 0)
        hit = run >= n - 1
        if hit.any():
            for s in np.flatnonzero(hit): out[t - run[s]:t + 1, s] = True
    return out


def loo_matrix(run, v, dates, cols):
    fs = sorted(glob.glob(str(C.OUT_ROOT / run / f"loo_{v}_*.parquet")))
    if not fs: sys.exit(f"no LOO files in {C.OUT_ROOT / run} (a first production run is needed)")
    L = pd.concat([pd.read_parquet(f, columns=["date", "station", "pred"]) for f in fs], ignore_index=True)
    ti = dates.get_indexer(pd.to_datetime(L["date"])); cp = {str(c): j for j, c in enumerate(cols)}
    si = L["station"].astype(str).map(cp).to_numpy(float); ok = (ti >= 0) & np.isfinite(si)
    M = np.full((len(dates), len(cols)), np.nan); M[ti[ok], si[ok].astype(int)] = L["pred"].to_numpy(float)[ok]
    return M


def neighbour_stats(x, xy, radius_m):
    """Per station-day: max and count of the K_NB nearest other stations within radius reporting."""
    good = np.isfinite(xy).all(1); tree = cKDTree(xy[good]); gi = np.flatnonzero(good)
    d, il = tree.query(xy[good], k=K_NB + 1); il = gi[il[:, 1:]]; d = d[:, 1:]
    nbmax = np.full_like(x, np.nan); nbmin = np.full_like(x, np.nan); nbn = np.zeros(x.shape, np.int16)
    V = x[:, il]                                                     # (T, Sgood, K)
    V = np.where((d <= radius_m)[None], V, np.nan)
    with np.errstate(all="ignore"):
        nbmax[:, gi] = np.nanmax(V, 2); nbmin[:, gi] = np.nanmin(V, 2)
    nbn[:, gi] = np.isfinite(V).sum(2)
    return nbmax, nbmin, nbn


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--variable", choices=["temp", "precip"], required=True)
    ap.add_argument("--run", default=None, help="run folder whose loo_*.parquet are used (default <var>_anomaly)")
    ap.add_argument("--fresh", action="store_true", help="ignore flags of earlier passes (default: cumulative)")
    a = ap.parse_args(); v = a.variable; run = a.run or f"{v}_anomaly"
    raw, dates, cols, x, tok = load_raw(v); xy = station_xy(v, cols); T, S = x.shape
    F = {}                                                            # test -> bool mask
    if v == "precip": F["R"] = tok | (np.isfinite(x) & (x < 0))
    else: F["R"] = tok | (np.isfinite(x) & ((x < C.TEMP_RANGE[0]) | (x > C.TEMP_RANGE[1])))
    xv = np.where(F["R"], np.nan, x)
    F["S"] = stuck_mask(xv, STUCK_DAYS[v], 0.5 if v == "precip" else None)
    xv = np.where(F["S"], np.nan, xv)
    pred = loo_matrix(run, v, dates, cols)
    if v == "temp":
        err = pred - xv
        mon = dates.month.to_numpy()
        med = np.full_like(err, np.nan)
        for m in range(1, 13):
            r = mon == m
            with np.errstate(all="ignore"): med[r] = np.nanmedian(err[r], 0)
        eb = err - med; cand = np.isfinite(eb) & (np.abs(eb) > C.QC_TEMP_ABS_C)
        Tm = np.zeros_like(cand); R = C.QC_CONTAM_KM * 1000.0
        for t in np.flatnonzero(cand.any(1)):
            ss = np.flatnonzero(cand[t]); ss = ss[np.argsort(-np.abs(eb[t, ss]))]; acc = []
            for s in ss:
                if any(np.hypot(*(xy[s] - xy[q])) < R and np.sign(eb[t, s]) != np.sign(eb[t, q]) and abs(eb[t, s]) < abs(eb[t, q]) for q in acc):
                    continue
                acc.append(s); Tm[t, s] = True
        F["T"] = Tm; info = eb
    else:
        nbmax25, _, n25 = neighbour_stats(xv, xy, 25000.0)
        _, nbmin15, n15 = neighbour_stats(xv, xy, 15000.0)
        with np.errstate(invalid="ignore"):
            F["H"] = (xv >= 30) & (n25 >= 3) & (xv > 5 * nbmax25 + 20) & (xv > 4 * np.nan_to_num(pred, nan=0.0) + 20)
            prev_missing = np.r_[np.zeros((1, S), bool), ~np.isfinite(x[:-1])]
            F["A"] = prev_missing & (xv >= 30) & (n25 >= 3) & (xv > 3 * nbmax25 + 10) & ~F["H"]
            F["Z"] = (xv == 0) & (n15 >= 3) & (nbmin15 >= 10) & (pred >= 15)
        info = xv - pred
    # cumulative QC: station-days flagged by an earlier pass have no LOO prediction in a run that excluded
    # them, so they cannot be re-detected; they are kept (test letter of the earlier pass)
    od = C.DATA_ROOT / "qc"; prev_f = od / C.QC_FLAGS_FILE.format(variable=v)
    if prev_f.exists() and not a.fresh:
        pv = pd.read_parquet(prev_f); ti = dates.get_indexer(pd.to_datetime(pv["date"])); cp = {str(c): j for j, c in enumerate(cols)}
        si = pv["station"].astype(str).map(cp).to_numpy(float); ok = (ti >= 0) & np.isfinite(si)
        for k in pv["test"].unique():
            m = np.zeros((T, S), bool); sel = ok & (pv["test"].to_numpy() == k); m[ti[sel], si[sel].astype(int)] = True
            F[k] = F.get(k, np.zeros((T, S), bool)) | m
        n_prev = int(ok.sum())
    else: n_prev = 0
    any_ = np.zeros((T, S), bool); rows = []
    for k, m in F.items():
        m = m & ~any_; any_ |= m; ti, si = np.nonzero(m)
        rows.append(pd.DataFrame(dict(date=dates[ti], station=np.asarray(cols)[si], obs=x[ti, si], test=k, info=info[ti, si])))
    FL = pd.concat(rows, ignore_index=True).sort_values(["date", "station"])
    od.mkdir(parents=True, exist_ok=True)
    FL.to_parquet(od / C.QC_FLAGS_FILE.format(variable=v), index=False)
    nobs = np.isfinite(x).sum(0)
    summ = FL.pivot_table(index="station", columns="test", values="obs", aggfunc="size", fill_value=0)
    summ["n_flagged"] = summ.sum(1); summ["n_obs"] = pd.Series(nobs, index=cols).reindex(summ.index)
    summ["pct_flagged"] = 100 * summ["n_flagged"] / summ["n_obs"]; summ = summ.sort_values("n_flagged", ascending=False)
    summ.to_csv(od / f"station_qc_summary_{v}.csv")
    # cleaned CSV: same layout, untouched strings except the flagged cells
    clean = raw.copy(); arr = clean[cols].to_numpy(object); arr[any_] = str(int(C.MISSING_VALUE)); clean[cols] = arr
    out = od / (Path(C.DATA[v]["raw"]).stem + "_QC.csv"); tmp = out.with_suffix(".part"); clean.to_csv(tmp, index=False); tmp.replace(out)
    n = int(np.isfinite(x).sum() + tok.sum())
    print(f"{v}: reported station-days {n} | flagged {len(FL)} ({100*len(FL)/n:.3f} %) | by test "
          + ", ".join(f"{k}={int((FL.test==k).sum())}" for k in F) + f" | stations {FL.station.nunique()} | from earlier passes {n_prev}, new {len(FL) - n_prev}")
    print(summ.head(10).round(2).to_string()); print(f"cleaned CSV -> {out}")


if __name__ == "__main__":
    main()
