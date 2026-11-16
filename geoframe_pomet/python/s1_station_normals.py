#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 1 - station monthly normals 1991-2020 with missing-data reconstruction.

For every station and month of 1991-2020 a monthly value is formed when at least
MONTH_COMPLETENESS of the days are observed (precipitation: daily mean x days in month;
temperature: monthly mean). Missing station-months are reconstructed from the best
correlated neighbours (<= FILL_RADIUS_KM, >= FILL_MIN_OVERLAP_MONTHS common months,
r >= FILL_MIN_R on calendar-month standardised series): precipitation by the neighbour's
value times the station/neighbour ratio of that calendar month, temperature by the
neighbour's value plus the station-neighbour difference; up to FILL_TOP_K donors weighted
by r^2. A normal is accepted when the calendar month has >= MIN_OBS_YEARS observed and
>= MIN_TOTAL_YEARS observed + reconstructed years. This removes the bias of normals
computed only from the years a station happened to report.

Usage:  python s1_station_normals.py --variable precip|temp
Output: <OUT_ROOT>/normals/station_normals_<var>.csv  (+ station_monthly_<var>.npz)
"""
import argparse, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
import po_config as C, po_common as P


def monthly_matrix(dates, values, variable):
    per = pd.PeriodIndex(dates, freq="M"); start = pd.Period(C.PERIOD[0], "M"); end = pd.Period(C.PERIOD[1], "M")
    months = pd.period_range(start, end, freq="M"); sel = (per >= start) & (per <= end)
    v = values[sel]; per = per[sel]
    codes = (per.year - start.year) * 12 + (per.month - 1); edges = np.r_[0, np.flatnonzero(np.diff(codes)) + 1]
    fin = np.isfinite(v); cnt = np.add.reduceat(fin, edges, axis=0); sm = np.add.reduceat(np.where(fin, v, 0.0), edges, axis=0)
    ndays = np.diff(np.r_[edges, len(v)])[:, None]
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = sm / cnt
    M = np.where(cnt >= C.MONTH_COMPLETENESS * ndays, mean * (ndays if variable == "precip" else 1.0), np.nan)
    return M.T.astype(float), months          # (S, 360)


def fill_gaps(M, xy, variable):
    S, NT = M.shape; cal = np.arange(NT) % 12
    mu = np.stack([np.nanmean(M[:, cal == m], 1) for m in range(12)], 1)          # (S,12)
    sd = np.stack([np.nanstd(M[:, cal == m], 1) for m in range(12)], 1)
    Z = (M - mu[:, cal]) / np.where(sd[:, cal] > 0, sd[:, cal], np.nan)
    tree = cKDTree(xy); F = M.copy(); n_filled = np.zeros(S, int); donors_used = np.zeros(S, int)
    k = min(C.FILL_MAX_CANDIDATES + 1, S)
    for s in range(S):
        miss = ~np.isfinite(M[s])
        if not miss.any(): continue
        d, idx = tree.query(xy[s], k=k); cand = idx[(d > 0) & (d <= C.FILL_RADIUS_KM * 1000)]
        if len(cand) == 0: continue
        zs = Z[s]; zc = Z[cand]; both = np.isfinite(zs)[None] & np.isfinite(zc); nov = both.sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            a = np.where(both, zs[None], 0); b = np.where(both, zc, 0)
            am = a.sum(1) / nov; bm = b.sum(1) / nov
            cov = (np.where(both, (zs[None] - am[:, None]) * (zc - bm[:, None]), 0)).sum(1)
            r = cov / np.sqrt((np.where(both, (zs[None] - am[:, None]) ** 2, 0)).sum(1) * (np.where(both, (zc - bm[:, None]) ** 2, 0)).sum(1))
        good = (nov >= C.FILL_MIN_OVERLAP_MONTHS) & (r >= C.FILL_MIN_R[variable])
        if not good.any(): continue
        o = np.argsort(-r[good]); cand = cand[good][o]; rr = r[good][o]; donors_used[s] = len(cand)
        # per-calendar-month ratio/difference between station and donor on common months
        both = np.isfinite(M[s])[None] & np.isfinite(M[cand])
        rel = np.full((len(cand), 12), np.nan)
        for m in range(12):
            bm_ = both[:, cal == m]; xs = np.where(bm_, M[s, cal == m][None], np.nan); xn = np.where(bm_, M[cand][:, cal == m], np.nan)
            with np.errstate(invalid="ignore", divide="ignore"):
                if variable == "precip":
                    rel[:, m] = np.nansum(xs, 1) / np.nansum(xn, 1)
                else:
                    rel[:, m] = np.nanmean(xs - xn, 1)
            rel[bm_.sum(1) < 5, m] = np.nan
        est = (M[cand] * rel[:, cal]) if variable == "precip" else (M[cand] + rel[:, cal])   # (c, NT)
        est = np.where(np.isfinite(est), est, np.nan)
        avail = np.isfinite(est); rank = np.cumsum(avail, 0); use = avail & (rank <= C.FILL_TOP_K)
        w = np.where(use, rr[:, None] ** 2, 0.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            val = (np.where(use, est, 0) * w).sum(0) / w.sum(0)
        put = miss & np.isfinite(val)
        if variable == "precip": val = np.maximum(val, 0.0)
        F[s, put] = val[put]; n_filled[s] = int(put.sum())
    return F, n_filled, donors_used


def normals(M, F):
    S, NT = M.shape; cal = np.arange(NT) % 12; rows = []
    for m in range(12):
        Mo, Fo = M[:, cal == m], F[:, cal == m]
        nobs = np.isfinite(Mo).sum(1); ntot = np.isfinite(Fo).sum(1)
        n_obs_only = np.nanmean(Mo, 1); n_fill = np.nanmean(Fo, 1)
        ok = (nobs >= C.MIN_OBS_YEARS) & (ntot >= C.MIN_TOTAL_YEARS)
        rows.append(pd.DataFrame(dict(month=m + 1, n_obs_years=nobs, n_total_years=ntot, normal_obs_only=n_obs_only,
                                      normal=np.where(ok, n_fill, np.nan), ok=ok)))
    return rows


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--variable", choices=["precip", "temp"], required=True)
    ap.add_argument("--exclude-stations", default=None, help="text file with one station ID per line to drop")
    a = ap.parse_args(); t0 = time.time()
    excl = [l.strip() for l in open(a.exclude_stations)] if a.exclude_stations else None
    st = P.load_station_data(a.variable, exclude=excl)
    M, months = monthly_matrix(st["dates"], st["values"], a.variable)
    F, n_filled, donors = fill_gaps(M, st["xy"], a.variable)
    parts = normals(M, F)
    out = pd.concat([p.assign(station=st["ids"]) for p in parts], ignore_index=True)
    dim = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    if a.variable == "precip": out["normal_daily_mm"] = out["normal"] / dim[out["month"] - 1]
    out = out[["station", "month", "normal", "normal_obs_only", "n_obs_years", "n_total_years", "ok"] + (["normal_daily_mm"] if a.variable == "precip" else [])]
    od = C.OUT_ROOT / "normals"; od.mkdir(parents=True, exist_ok=True)
    out.to_csv(od / f"station_normals_{a.variable}.csv", index=False)
    np.savez_compressed(od / f"station_monthly_{a.variable}.npz", ids=st["ids"], observed=M.astype(np.float32), filled=F.astype(np.float32),
                        n_filled=n_filled, donors=donors, months=np.array([str(p) for p in months]))
    ok = out["ok"]
    rel = (out.normal - out.normal_obs_only) / (out.normal_obs_only.abs() if a.variable == "precip" else 1.0)
    print(f"stations={len(st['ids'])}  monthly values observed={np.isfinite(M).mean()*100:.1f}%  after reconstruction={np.isfinite(F).mean()*100:.1f}%")
    print(f"station-months with accepted normal: {ok.sum()} of {len(ok)} ({100*ok.mean():.1f}%); stations with all 12 normals: {(out.groupby('station').ok.sum()==12).sum()}")
    print(f"effect of reconstruction on accepted normals ({'relative' if a.variable=='precip' else 'degC'}): "
          f"median {np.nanmedian(rel[ok]):+.4f}, p05 {np.nanpercentile(rel[ok],5):+.4f}, p95 {np.nanpercentile(rel[ok],95):+.4f}")
    print(f"done in {time.time()-t0:.0f} s -> {od}")


if __name__ == "__main__":
    main()
