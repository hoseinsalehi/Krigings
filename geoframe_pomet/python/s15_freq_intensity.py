#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 15 - daily frequency / intensity test of the elevation dependence (supervisor's question, 2026-10-11).

Paired-station design: pairs closer than 20 km with an elevation difference >= 300 m, compared on their COMMON
valid days 1991-2020, so that the same events are compared. For each pair (higher h, lower l):
    ln(P_h/P_l) = ln(F_h/F_l) + ln(I_h/I_l)      F: wet-day frequency (>= 1 mm), I: mean wet-day amount
plus P99 of wet days and the same-day ratio on days when both are wet. Gradients are % per km of elevation,
median over pairs, with 90 % station-block bootstrap intervals, by elevation band of the higher station.
Sources, all sampled on the same station-days:
  obs                 station observations (QC'd, day-harmonised; from the leave-one-out files)
  <run>_loo / _grid   this run (OUT_ROOT): leave-one-out prediction and the gridded product at the station cell
  cmp_loo / cmp_grid  comparison run (PO_COMPARE_ROOT, default ~/Downloads/PoGrid_v4 = dataset v2.0, ratio method)
  v1_grid             submitted dataset v1.0 (GEOframe ordinary kriging; ~/Downloads/precipitation_data)
Output: <OUT_ROOT>/freq_intensity/{pairs.csv, gradients_ci.csv, stations_by_band.csv, freq_intensity.png/.pdf}
"""
import glob, os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd, xarray as xr
from scipy.spatial import cKDTree
import po_config as C

HOME = Path.home(); CMP = Path(os.environ.get("PO_COMPARE_ROOT", HOME / "Downloads/PoGrid_v4")); V1 = HOME / "Downloads/precipitation_data"
OD = C.OUT_ROOT / "freq_intensity"; OD.mkdir(parents=True, exist_ok=True)
WET, MAXD, MINDZ, MINCOMMON, NB = 1.0, 20.0, 300.0, 1500, 300
BANDS = [(0, 1000, "< 1.0"), (1000, 1500, "1.0–1.5"), (1500, 2000, "1.5–2.0"), (2000, 5000, "> 2.0"), (1500, 5000, "> 1.5 (all)")]
MET = ["lnP", "lnF", "lnI", "p99_ln", "ev_ln_median"]


def loo_table(root):
    L = pd.concat([pd.read_parquet(f, columns=["date", "station", "obs", "pred"]) for f in sorted(glob.glob(str(root / "precip_anomaly/loo_precip_*.parquet")))])
    L = L.dropna(subset=["obs", "pred"]); L["station"] = L.station.astype(str); L["date"] = pd.to_datetime(L.date); return L


def main():
    st = pd.read_csv(CMP / "paper_tables/stations_precip.csv", dtype={"id": str}).set_index("id")
    Lr = loo_table(C.OUT_ROOT); Lc = loo_table(CMP)
    O = Lr.pivot_table(index="date", columns="station", values="obs")
    ids = [s for s in O.columns if s in st.index and O[s].notna().sum() >= 3000]; O = O[ids]
    R = Lr.pivot_table(index="date", columns="station", values="pred").reindex(index=O.index, columns=ids)
    Cm = Lc.pivot_table(index="date", columns="station", values="pred").reindex(index=O.index, columns=ids)
    xy = st.loc[ids, ["x", "y"]].to_numpy(); z = st.loc[ids, "elevation"].to_numpy()
    with xr.open_dataset(V1 / "precipitation_1991.nc") as d: gx, gy = d.x.values, d.y.values
    ix = np.abs(gx[None, :] - xy[:, :1]).argmin(1); iy = np.abs(gy[None, :] - xy[:, 1:]).argmin(1)
    def grid(files, var, nodata=None):
        out = []
        for f in files:
            with xr.open_dataset(f) as d:
                a = d[var].values[:, iy, ix].astype(float); t = pd.to_datetime(d.time.values).normalize()
            if nodata is not None: a[a < nodata] = np.nan
            out.append(pd.DataFrame(a, index=t, columns=ids))
        return pd.concat(out).reindex(O.index)
    yrs = range(1991, 2021)
    SRC = dict(obs=O, run_loo=R.clip(lower=0), run_grid=grid([C.OUT_ROOT / f"precip_anomaly/precip_{y}.nc" for y in yrs], "precipitation").clip(lower=0),
               cmp_loo=Cm.clip(lower=0), cmp_grid=grid([CMP / f"precip_anomaly/precip_{y}.nc" for y in yrs], "precipitation").clip(lower=0),
               v1_grid=grid([V1 / f"precipitation_{y}.nc" for y in yrs], "precipitation", -9000).clip(lower=0))
    A = {k: v.to_numpy() for k, v in SRC.items()}
    mask = np.logical_and.reduce([np.isfinite(x) for x in A.values()])

    def stats(a, b):
        wa, wb = a >= WET, b >= WET; both = wa & wb
        return dict(lnP=np.log(a.mean() / b.mean()), lnF=np.log(wa.mean() / wb.mean()), lnI=np.log(a[wa].mean() / b[wb].mean()),
                    p99_ln=np.log(np.quantile(a[wa], .99) / np.quantile(b[wb], .99)), ev_ln_median=np.median(np.log(a[both] / b[both])))
    rows = []
    for i, j in sorted(cKDTree(xy / 1e3).query_pairs(MAXD)):
        h, l = (i, j) if z[i] > z[j] else (j, i); dz = z[h] - z[l]
        if dz < MINDZ: continue
        ok = mask[:, h] & mask[:, l]
        if ok.sum() < MINCOMMON or any(((A[k][ok, h] >= WET).sum() < 100 or (A[k][ok, l] >= WET).sum() < 100) for k in A): continue
        r = dict(high=ids[h], low=ids[l], z_high=z[h], dz_km=dz / 1e3)
        for k in A: r.update({f"{k}_{m}": x for m, x in stats(A[k][ok, h], A[k][ok, l]).items()})
        rows.append(r)
    D = pd.DataFrame(rows); D.to_csv(OD / "pairs.csv", index=False)
    rng = np.random.default_rng(1); g = lambda d, c: 100 * (np.exp(np.nanmedian(d[c] / d.dz_km)) - 1); out = []
    for lo, hi, lab in BANDS:
        d = D[(D.z_high >= lo) & (D.z_high < hi)].reset_index(drop=True); stn = np.unique(np.r_[d.high, d.low])
        boots = [d[d.high.isin(p) & d.low.isin(p)] for p in (set(rng.choice(stn, len(stn), replace=True)) for _ in range(NB))]
        for k in A:
            for m in MET:
                c = f"{k}_{m}"; b = [g(x, c) for x in boots]
                out.append(dict(band=lab, source=k, metric=m, median=g(d, c), ci05=np.nanquantile(b, .05), ci95=np.nanquantile(b, .95), n_pairs=len(d)))
    T = pd.DataFrame(out); T.to_csv(OD / "gradients_ci.csv", index=False)
    # absolute station-level biases by elevation band
    rows = []
    for k_, s_ in enumerate(ids):
        ok = mask[:, k_]; r = dict(station=s_, z=z[k_])
        for k in A:
            v = A[k][ok, k_]; w = v >= WET
            r.update({f"{k}_P": 365.25 * v.mean(), f"{k}_F": w.mean(), f"{k}_I": v[w].mean() if w.any() else np.nan, f"{k}_P99": np.quantile(v[w], .99) if w.sum() > 50 else np.nan})
        rows.append(r)
    S = pd.DataFrame(rows); S["band"] = pd.cut(S.z, [0, 500, 1000, 1500, 2000, 4000])
    B = S.groupby("band", observed=True).apply(lambda q: pd.Series({"n": len(q), **{f"{k}_{x}_bias_pct": 100 * (np.nanmedian(q[f"{k}_{x}"] / q[f"obs_{x}"]) - 1)
                                                                            for k in A if k != "obs" for x in ("P", "F", "I", "P99")}}))
    B.to_csv(OD / "stations_by_band.csv")
    pd.set_option("display.width", 250)
    W = T[T.band == "> 1.5 (all)"].pivot(index="metric", columns="source", values="median").round(1)
    print(f"pairs {len(D)}; gradients above 1.5 km (% per km):"); print(W.to_string())
    print(B.round(1).T.to_string())
    fig_plot(T)


def fig_plot(T):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    SRC = [("obs", "Stations (observed)", "#000000"), ("v1_grid", "v1.0 (ordinary kriging)", "#E69F00"), ("cmp_grid", "v2.0 ratio method, grid", "#0072B2"),
           ("cmp_loo", "v2.0 ratio method, station withheld", "#56B4E9"), ("run_grid", "two-part method, grid", "#009E73"), ("run_loo", "two-part method, station withheld", "#8FD9B6")]
    TITLE = dict(lnP="Annual total", lnF="Wet-day frequency (≥ 1 mm)", lnI="Mean wet-day intensity", p99_ln="P99 of wet days", ev_ln_median="Same-day ratio (both wet)")
    bands = [b[2] for b in BANDS[:4]]
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300})
    fig, AX = plt.subplots(2, 3, figsize=(175 / 25.4, 115 / 25.4), sharey=True); x = np.arange(len(bands)); w = 0.13
    for ax, m in zip(AX.flat, MET):
        for n, (s, l, col) in enumerate(SRC):
            q = T[(T.source == s) & (T.metric == m)].set_index("band").reindex(bands)
            ax.bar(x + (n - 2.5) * w, q["median"], w * 0.95, color=col, label=l)
            ax.errorbar(x + (n - 2.5) * w, q["median"], yerr=[q["median"] - q.ci05, q.ci95 - q["median"]], fmt="none", ecolor="#555555", elinewidth=0.5)
        ax.axhline(0, color="k", lw=0.6); ax.set_title(TITLE[m], fontsize=7.5); ax.set_xticks(x); ax.set_xticklabels(bands, fontsize=6.5); ax.grid(axis="y", lw=0.3, alpha=0.5)
    for a in AX[:, 0]: a.set_ylabel("Change per km of elevation (%)")
    L = AX[1, 2]; L.axis("off"); h, l = AX[0, 0].get_legend_handles_labels(); L.legend(h, l, frameon=False, loc="upper left", fontsize=6.5)
    L.text(0, 0.25, "Station pairs < 20 km apart, Δz ≥ 300 m,\ncommon days 1991–2020. Bars: median\nover pairs; whiskers: 90 % station\nbootstrap interval.", transform=L.transAxes, fontsize=6.5, va="top")
    fig.supxlabel("Elevation of the higher station of the pair (km)", fontsize=7.5); fig.tight_layout()
    for e in ("png", "pdf"): fig.savefig(OD / f"freq_intensity.{e}", bbox_inches="tight")


if __name__ == "__main__":
    main()
