"""Figs. 2-9 of the revised manuscript (dataset v2.0). Usage: python figs_02_09.py [F02 F03 ...]
Inputs: ../data/stations_<v>.csv, ../data/loo_station[_year]_<v>.csv (precipitation rebuilt from the final two-part
product by build_precip_tables.py; temperature from s14_paper_tables.py), <V5>/harmonisation, <V5>/normals,
<V5>/summary/cell_climatology_*.nc. v1.0 reference values in Fig. 6: station-wise means reported in the submitted
manuscript (Table R5 of the reply to Reviewer 1)."""
import sys
import numpy as np, pandas as pd, xarray as xr, matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import cmcrameri.cm as cmc
from figlib import *

GRP = [("inside district", "Inside the District", C["blue"]), ("buffer (20 km)", "20 km buffer", C["orange"]), ("MeteoSwiss", "MeteoSwiss", C["vermil"])]
st = {v: pd.read_csv(TAB / f"stations_{v}.csv", dtype={"id": str}) for v in ("precip", "temp")}
YRS = [c for c in st["precip"].columns if c.startswith("days_")]
NDAYS = np.array([366 if int(c[5:]) % 4 == 0 else 365 for c in YRS])


def loo(v): return pd.read_csv(TAB / f"loo_station_{v}.csv", dtype={"station": str})


def F02():
    style(); fig = plt.figure(figsize=(DOUBLE, 205 * MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 0.95, 0.75], hspace=0.55, wspace=0.22)
    for i, v in enumerate(("precip", "temp")):
        s = st[v]; ax = fig.add_subplot(gs[0, i]); basemap(ax, ticks=True, pad=5)
        for g, lab, col in GRP:
            q = s[s.group == g]; ax.scatter(q.x / 1e3, q.y / 1e3, s=2.2, color=col, lw=0, zorder=4, label=lab)
        if i == 1: ax.set_ylabel(""); ax.set_yticklabels([])
        ax.set_title(f"{VLAB[v]} ({len(s)} stations)", pad=3); panel_label(ax, "ab"[i], x=-0.02 if i else -0.17, y=1.01)
        if i == 0:
            h = [Line2D([], [], marker="o", ls="", ms=4, color=c) for _, _, c in GRP]
            fig.legend(h, [l for _, l, _ in GRP], loc="upper center", bbox_to_anchor=(0.5, 0.635), ncol=3, frameon=False, handletextpad=0.2, columnspacing=1.5)
        ax = fig.add_subplot(gs[1, i]); ss = s.sort_values("elevation"); A = 100 * ss[YRS].to_numpy(float) / NDAYS
        im = ax.imshow(A, aspect="auto", cmap=cmc.batlow, vmin=0, vmax=100, interpolation="nearest")
        ax.set_xticks(range(4, len(YRS), 5)); ax.set_xticklabels([c[5:] for c in YRS][4::5])
        zt = [0, 250, 500, 1000, 2000]; pos = [int(np.searchsorted(ss.elevation.to_numpy(), z)) for z in zt]
        ax.set_yticks(pos); ax.set_yticklabels(zt); ax.set_ylabel("Stations by elevation (m)" if i == 0 else "")
        ax.set_title(VLAB[v], pad=3); panel_label(ax, "cd"[i], x=-0.17 if i == 0 else -0.12)
    cax = fig.add_axes([0.93, 0.36, 0.012, 0.22]); cb = fig.colorbar(im, cax=cax); cb.set_label("Valid days per year (%)")
    ax = fig.add_subplot(gs[2, 0]); bins = np.arange(0, 31, 1)
    for v, col in (("precip", C["blue"]), ("temp", C["orange"])):
        rl = st[v][YRS].sum(1) / 365.25
        ax.hist(rl, bins=bins, histtype="step", lw=1.3, color=col, label=VLAB[v])
        ax.axvline(rl.median(), color=col, lw=0.8, ls="--")
        print(v, "median record (yr)", round(rl.median(), 2), "< 10 yr:", int((rl < 10).sum()), "n", len(rl))
    ax.set_xlabel("Record length (years of valid data)"); ax.set_ylabel("Stations"); ax.set_xlim(0, 30.5); panel_label(ax, "e", x=-0.17)
    ax = fig.add_subplot(gs[2, 1])
    for v, col in (("precip", C["blue"]), ("temp", C["orange"])):
        n = st[v][YRS].sum(0).to_numpy() / NDAYS; ax.plot([int(c[5:]) for c in YRS], n, color=col, marker="o", ms=2.5, label=VLAB[v])
        print(v, "stations/day 1991", round(n[0]), "2020", round(n[-1]))
    ax.set_xlabel("Year"); ax.set_ylabel("Stations reporting per day"); ax.set_ylim(0, None); panel_label(ax, "f", x=-0.2)
    h, l = ax.get_legend_handles_labels(); fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 0.045), ncol=2, frameon=False)
    for a in fig.axes[-2:]: a.grid(lw=0.3, alpha=0.5)
    save(fig, "Fig02_station_network")


def F03():
    style(); fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 68 * MM))
    bins = np.arange(0, 4001, 400); mid = (bins[:-1] + bins[1:]) / 2; area = np.histogram(G["elev"], bins)[0] / len(G["elev"]) * 100
    ax = axes[0]; ax2 = ax.twinx(); w = 150
    for k, (v, col) in enumerate((("precip", C["blue"]), ("temp", C["orange"]))):
        s = st[v]; h = np.histogram(s.elevation, bins)[0]; ax.bar(mid + (k - 0.5) * w, h, width=w, color=col, label=f"{VLAB[v]} stations")
        av = 100 * s[YRS].sum(1) / NDAYS.sum(); b = pd.cut(s.elevation, bins)
        m = av.groupby(b, observed=False).mean().to_numpy().copy(); n_ = av.groupby(b, observed=False).size().to_numpy(); m[n_ < 10] = np.nan
        axes[1].plot(mid, m, marker="o", ms=3, color=col, label=VLAB[v]); print(v, "stations per bin", dict(zip(mid.astype(int), n_)))
    ax2.plot(mid, area, color=C["black"], lw=1, ls="--", label="District area"); ax2.set_ylabel("Share of District area (%)"); ax2.spines["right"].set_visible(True)
    ax.set_xlabel("Elevation (m)"); ax.set_ylabel("Number of stations"); panel_label(ax, "a", x=-0.2)
    axes[1].set_xlabel("Elevation (m)"); axes[1].set_ylabel("Mean data availability (%)"); axes[1].set_ylim(0, 100); axes[1].grid(lw=0.3, alpha=0.5); panel_label(axes[1], "b", x=-0.2)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0.1, 1, 1), w_pad=3)
    fig.legend(h1 + h2, l1 + l2, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
    save(fig, "Fig03_station_elevation")


def F04():
    style(); fig = plt.figure(figsize=(DOUBLE, 118 * MM)); gs = fig.add_gridspec(2, 2, width_ratios=[1.6, 1], hspace=0.35, wspace=0.25)
    CC = {"unchanged": ("Unchanged", C["grey"]), "-1": ("Shifted by −1 day", C["blue"]), "+1": ("Shifted by +1 day", C["green"]), "mixed": ("Changed over time", C["vermil"])}
    for i, v in enumerate(("precip", "temp")):
        d = pd.read_csv(V5 / "harmonisation" / f"day_shift_{v}.csv", dtype={"station": str})
        cls = d.groupby("station").shift_days.agg(lambda s: "unchanged" if (s == 0).all() else ("-1" if (s == -1).all() else ("+1" if (s == 1).all() else "mixed")))
        m = st[v].set_index("id").join(cls.rename("cls")); ax = fig.add_subplot(gs[i, 0]); basemap(ax, grey=False, pad=5)
        for c, (lab, col) in CC.items():
            q = m[m.cls == c]
            if len(q): ax.scatter(q.x / 1e3, q.y / 1e3, s=2.2 if c == "unchanged" else 3.5, color=col, lw=0, zorder=4 if c == "unchanged" else 5)
        print(v, "classes", m.cls.value_counts().to_dict())
        ax.set_title(VLAB[v], pad=2); panel_label(ax, "ac"[i], x=0.0, y=0.98)
        pc = pd.read_csv(V5 / "harmonisation" / f"neighbour_pair_check_{v}.csv").set_index("state")
        ax = fig.add_subplot(gs[i, 1]); vals = [pc.loc["raw", "pct_pairs_offset"], pc.loc["harmonised", "pct_pairs_offset"]]
        ax.bar([0, 1], vals, color=[C["grey"], C["blue"]], width=0.6); ax.set_xticks([0, 1]); ax.set_xticklabels(["Original labels", "Harmonised"])
        ax.set_ylabel("Pairs offset by one day (%)"); ax.set_ylim(0, 8); ax.grid(axis="y", lw=0.3, alpha=0.5); panel_label(ax, "bd"[i], x=-0.3)
        for x, val in zip([0, 1], vals): ax.text(x, val + 0.15, f"{val:.2f}" if val < 1 else f"{val:.1f}", ha="center", va="bottom", fontsize=8)
        print(v, pc.to_dict())
    h = [Line2D([], [], marker="o", ls="", ms=4, color=col) for _, col in CC.values()]
    fig.legend(h, [l for l, _ in CC.values()], loc="upper center", bbox_to_anchor=(0.36, 0.06), ncol=4, frameon=False, handletextpad=0.2)
    save(fig, "Fig04_day_harmonisation")


def F05():
    style(); fig, axes = plt.subplots(2, 2, figsize=(DOUBLE, 128 * MM))
    nP = xr.open_dataset(V5 / "normals" / "normals_precip.nc").normal.values.sum(0); nT = xr.open_dataset(V5 / "normals" / "normals_temp.nc").normal.values.mean(0)
    cP = xr.open_dataset(V5 / "summary" / "cell_climatology_precip_anomaly.nc"); cT = xr.open_dataset(V5 / "summary" / "cell_climatology_temp_anomaly.nc")
    pan = [(nP, "Annual precipitation normal (mm yr$^{-1}$)", dict(cmap=cmc.oslo_r, norm=LogNorm(500, 2600)), None),
           (nT, "Annual mean temperature normal (°C)", dict(cmap=cmc.batlow, vmin=-10, vmax=15), "both"),
           (cP.mean_R_pos_half_2011_2020.values, "Precipitation, ½-weight radius (km)", dict(cmap=cmc.lajolla, vmin=0, vmax=30), "max"),
           (cT.mean_R_pos_half_2011_2020.values, "Temperature, ½-weight radius (km)", dict(cmap=cmc.lajolla, vmin=0, vmax=30), "max")]
    for k, (a, lab, kw, ext) in enumerate(pan):
        ax = axes.flat[k]; basemap(ax, grey=False, pad=3); im = raster(ax, a, **kw)
        cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02, orientation="horizontal", extend=ext or "neither", aspect=35); cb.set_label(lab)
        if k == 0: cb.set_ticks([500, 700, 1000, 1500, 2000, 2500]); cb.set_ticklabels(["500", "700", "1000", "1500", "2000", "2500"]); cb.ax.xaxis.set_minor_formatter(plt.NullFormatter())
        panel_label(ax, "abcd"[k], x=0.0, y=0.97)
        print(lab, "district mean", np.nanmean(np.where(MASK, a, np.nan)))
    fig.tight_layout(h_pad=0.5, w_pad=0.5); save(fig, "Fig05_normals_support")


V1 = {  # station-wise means of the submitted dataset (v1.0), Table R5 of the reply to Reviewer 1
    ("precip", "all"): dict(kge=0.84, cc=0.90, mae=1.28, rmse=3.79, bias=1.1),
    ("precip", "wet"): dict(kge=0.82, cc=0.88, mae=3.05, rmse=6.02, bias=-4.3),
    ("temp", "all"): dict(kge=0.88, cc=0.98, mae=1.14, rmse=1.50, bias=0.02)}


def F06():
    style(); fig, axes = plt.subplots(2, 5, figsize=(DOUBLE, 100 * MM))
    P, Tm = loo("precip"), loo("temp")
    M = [("kge", "KGE"), ("cc", "r"), ("mae", "MAE"), ("rmse", "RMSE"), ("bias", "Bias")]
    bp = dict(widths=0.55, showfliers=False, whis=(5, 95), patch_artist=True, medianprops=dict(color="k", lw=1), boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    for j, (m, lab) in enumerate(M):
        for i, (v, sets) in enumerate((("precip", [("all", "All", C["blue"]), ("wet", "Wet", C["sky"])]), ("temp", [("all", "All", C["orange"])]))):
            ax = axes[i, j]; d = P if v == "precip" else Tm
            data = [d[d.days == s][m].dropna() for s, _, _ in sets]
            b = ax.boxplot(data, **bp)
            for patch, (_, _, col) in zip(b["boxes"], sets): patch.set_facecolor(col); patch.set_alpha(0.8)
            for x, (s, _, _) in enumerate(sets, 1):
                ax.scatter(x, data[x - 1].mean(), marker="o", s=14, facecolor="white", edgecolor="k", lw=0.7, zorder=5)
                ax.scatter(x + 0.38, V1[(v, s)][m], marker="D", s=14, color=C["vermil"], zorder=5)
            ax.set_xticks(range(1, len(sets) + 1)); ax.set_xticklabels([l for _, l, _ in sets]); ax.set_xlim(0.4, len(sets) + 0.7)
            unit = {"mae": " (mm)" if v == "precip" else " (°C)", "rmse": " (mm)" if v == "precip" else " (°C)", "bias": " (%)" if v == "precip" else " (°C)"}.get(m, "")
            ax.set_title(lab + unit, pad=3); ax.grid(axis="y", lw=0.3, alpha=0.5)
            if j == 0: ax.set_ylabel(VLAB[v])
            if m == "bias": ax.axhline(0, color=C["grey"], lw=0.6)
            panel_label(ax, "abcdefghij"[i * 5 + j], x=-0.35, y=1.02)
    print("P n", P.groupby("days").size().to_dict(), "T n", len(Tm))
    h = [plt.Rectangle((0, 0), 1, 1, color=C["grey"], alpha=0.6), Line2D([], [], marker="o", ls="", mfc="white", mec="k", ms=4), Line2D([], [], marker="D", ls="", color=C["vermil"], ms=4)]
    fig.tight_layout(rect=(0, 0.07, 1, 1), w_pad=0.8, h_pad=1.2)
    fig.legend(h, ["v2.0, station-wise values (box: quartiles; whiskers: 5–95 %)", "v2.0, mean", "v1.0 (submitted), mean"], loc="lower center", ncol=3, frameon=False)
    save(fig, "Fig06_loo_distributions")


def F07():
    style(); fig, axes = plt.subplots(2, 2, figsize=(DOUBLE, 128 * MM)); P, Tm = loo("precip"), loo("temp"); P = P[P.days == "all"]
    pan = [(P, "kge", "Precipitation KGE", cmc.batlow, (0.5, 1.0)), (P, "mae", "Precipitation MAE (mm)", cmc.imola, (0.4, 2.5)),
           (Tm, "kge", "Temperature KGE", cmc.batlow, (0.7, 1.0)), (Tm, "mae", "Temperature MAE (°C)", cmc.imola, (0.4, 1.8))]
    for k, (d, m, lab, cm, (lo, hi)) in enumerate(pan):
        ax = axes.flat[k]; basemap(ax, pad=3); o = d.sort_values(m, ascending=(m != "kge"))
        sc = ax.scatter(o.x / 1e3, o.y / 1e3, c=o[m], s=3.5, cmap=cm, vmin=lo, vmax=hi, lw=0, zorder=4)
        cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02, orientation="horizontal", extend="both", aspect=35); cb.set_label(lab); panel_label(ax, "abcd"[k], x=0.0, y=0.97)
    fig.tight_layout(h_pad=0.5, w_pad=0.5); save(fig, "Fig07_loo_maps")


def F08():
    style(); fig, axes = plt.subplots(2, 3, figsize=(DOUBLE, 105 * MM)); P, Tm = loo("precip"), loo("temp")
    bands = [0, 500, 1000, 1500, 2000, 4000]; blab = ["<500", "500–\n1000", "1000–\n1500", "1500–\n2000", ">2000"]
    bp = dict(showfliers=False, whis=(5, 95), patch_artist=True, medianprops=dict(color="k", lw=1), boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    for i, (v, d) in enumerate((("precip", P), ("temp", Tm))):
        sets = [("all", C["blue"]), ("wet", C["sky"])] if v == "precip" else [("all", C["orange"])]
        for j, (m, lab) in enumerate((("kge", "KGE"), ("mae", "MAE"), ("bias", "Bias"))):
            ax = axes[i, j]; w = 0.36 if len(sets) == 2 else 0.55
            for k, (s, col) in enumerate(sets):
                q = d[d.days == s]; b = pd.cut(q.elevation, bands, right=False)
                data = [q.loc[b == c, m].dropna() for c in b.cat.categories]
                if i == 0 and j == 0 and s == "all": print(v, "stations per band", [len(x) for x in data])
                if i == 1 and j == 0: print(v, "stations per band", [len(x) for x in data])
                off = (k - (len(sets) - 1) / 2) * w
                bb = ax.boxplot(data, positions=np.arange(len(data)) + off, widths=w * 0.9, **bp)
                for p in bb["boxes"]: p.set_facecolor(col); p.set_alpha(0.85)
            ax.set_xticks(range(len(blab))); ax.set_xticklabels(blab if i == 1 else [""] * len(blab))
            unit = {"mae": " (mm)" if v == "precip" else " (°C)", "bias": " (%)" if v == "precip" else " (°C)", "kge": ""}[m]
            ax.set_ylabel(f"{lab}{unit}"); ax.grid(axis="y", lw=0.3, alpha=0.5)
            if m == "bias": ax.axhline(0, color=C["grey"], lw=0.6)
            if i == 1: ax.set_xlabel("Station elevation (m)")
            ax.set_title(f"{VLAB[v]}, {lab}", pad=3); panel_label(ax, "abcdef"[i * 3 + j], x=-0.28)
    h = [plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.85) for c in (C["blue"], C["sky"], C["orange"])]
    fig.tight_layout(rect=(0, 0.06, 1, 1), w_pad=1.0)
    fig.legend(h, ["Precipitation, all days", "Precipitation, wet days (≥ 0.1 mm)", "Temperature"], loc="lower center", ncol=3, frameon=False)
    save(fig, "Fig08_loo_elevation")


def F09():
    style(); fig, axes = plt.subplots(3, 1, figsize=(DOUBLE, 130 * MM))
    sets = [("precip", "all", "Precipitation, all days"), ("precip", "wet", "Precipitation, wet days (≥ 0.1 mm)"), ("temp", "all", "Temperature")]
    M = [("kge", "KGE", True), ("cc", "r", True), ("bias", "Bias", None), ("mae", "MAE", False), ("rmse", "RMSE", False)]
    for k, (v, days, title) in enumerate(sets):
        d = pd.read_csv(TAB / f"loo_station_year_{v}.csv"); d = d[d.days == days]; g = d.groupby("year")[[m for m, _, _ in M]].mean()
        if k == 0: print("station-years per year (P all)", d.groupby("year").size().iloc[[0, -1]].to_dict())
        Z = np.zeros((len(M), len(g)))
        for r, (m, _, hb) in enumerate(M):
            x = g[m].to_numpy(); x = -np.abs(x) if hb is None else (x if hb else -x); Z[r] = (x - x.min()) / (np.ptp(x) + 1e-12)
        ax = axes[k]; ax.imshow(Z, aspect="auto", cmap=cmc.batlow, vmin=0, vmax=1)
        for r, (m, lab, _) in enumerate(M):
            for c, val in enumerate(g[m].to_numpy()):
                ax.text(c, r, f"{val:.2f}" if abs(val) < 10 else f"{val:.0f}", ha="center", va="center", fontsize=5.6, color="w" if Z[r, c] < 0.45 else "k")
        unit = ("mm", "%") if v == "precip" else ("°C", "°C")
        ax.set_yticks(range(len(M))); ax.set_yticklabels(["KGE", "r", f"Bias ({unit[1]})", f"MAE ({unit[0]})", f"RMSE ({unit[0]})"])
        ax.set_xticks(range(0, len(g), 2)); ax.set_xticklabels(g.index[::2]); ax.set_title(title, pad=3, loc="left")
        for s in ax.spines.values(): s.set_visible(False)
        ax.tick_params(length=0); panel_label(ax, "abc"[k], x=-0.1)
    fig.tight_layout(h_pad=0.8); save(fig, "Fig09_loo_annual")


FIGS = dict(F02=F02, F03=F03, F04=F04, F05=F05, F06=F06, F07=F07, F08=F08, F09=F09)
if __name__ == "__main__":
    w = sys.argv[1:] or list(FIGS)
    for k in w: FIGS[k]()
