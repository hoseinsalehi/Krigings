"""Reviewer 1, points 1b and 3 - physiographic structure of the 1-km normals vs the ISAC-CNR 30-arcsec climatologies
(precipitation: Crespi et al. 2018; temperature: Brunetti et al. 2014; both 1961-1990, annual).
Reference period differs from ours (1991-2020): compare pattern and elevation dependence, not absolute level.
(a, c) cell-by-cell scatter (density); (b, d) ratio / difference per 250-m elevation bin (median and 10-90 %).
Source: <V4>/isac_comparison/{cells.parquet, gradients.csv} (s12_isac_normals_comparison.py)."""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from figlib import V4, DOUBLE, C, style, panel_label, save

d = pd.read_parquet(V4 / "isac_comparison" / "cells.parquet"); gr = pd.read_csv(V4 / "isac_comparison" / "gradients.csv").set_index(["variable", "group"])
style()
fig, axes = plt.subplots(2, 2, figsize=(DOUBLE, 125 / 25.4))
bins = np.arange(0, 4001, 250); mid = (bins[:-1] + bins[1:]) / 2
for i, (v, o, r, lab, unit, log) in enumerate((("precip", "P_ours_1991_2020", "P_crespi_1961_1990", "precipitation", "mm yr$^{-1}$", True),
                                               ("temp", "T_ours_1991_2020", "T_brunetti_1961_1990", "temperature", "°C", False))):
    q = d[[o, r, "elev"]].dropna(); ax = axes[i, 0]
    hb = ax.hexbin(q[r], q[o], gridsize=70, mincnt=1, cmap="viridis", norm=LogNorm(), xscale="log" if log else "linear", yscale="log" if log else "linear", lw=0, rasterized=True)
    # rasterized: the vector hexbin collection was dropped/clipped by PDF renderers (empty panel a in round 1)
    lo, hi = np.nanpercentile(np.r_[q[o], q[r]], [0.1, 99.9]); ax.plot([lo, hi], [lo, hi], color=C["vermil"], lw=0.8)
    rr = np.corrcoef(np.log(q[o]) if log else q[o], np.log(q[r]) if log else q[r])[0, 1]
    print(v, "r", round(rr, 3), "n", len(q))
    ax.set_xlabel(f"{'Crespi et al. (2018)' if log else 'Brunetti et al. (2014)'} ({unit})"); ax.set_ylabel(f"GEOframe-PoMet ({unit})")
    if log:
        from matplotlib.ticker import FixedLocator, NullFormatter
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(FixedLocator([500, 1000, 2000, 3000])); axis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}")); axis.set_minor_formatter(NullFormatter())
    cb = fig.colorbar(hb, ax=ax, pad=0.02); cb.set_label("Cells")
    ax = axes[i, 1]; e = 100 * (q[o] / q[r] - 1) if log else q[o] - q[r]; k = np.digitize(q["elev"], bins) - 1
    st = pd.DataFrame(dict(k=k, e=e)).groupby("k").e.quantile([0.1, 0.5, 0.9]).unstack().reindex(range(len(mid)))
    cnt = pd.Series(k).value_counts().reindex(range(len(mid))).fillna(0).to_numpy(); ok = cnt >= 50
    ax.fill_betweenx(mid[ok], st[0.1][ok], st[0.9][ok], color=C["sky"], alpha=0.35, lw=0, label="10–90 %")
    ax.plot(st[0.5][ok], mid[ok], color=C["blue"], marker="o", ms=2.5, label="Median")
    ax.axvline(0, color=C["grey"], lw=0.6); ax.set_ylabel("Cell elevation (m)")
    ax.set_xlabel("Ratio GEOframe-PoMet / reference − 1 (%)" if log else "Difference GEOframe-PoMet − reference (°C)")
    g = gr.loc[(v, "All_PRD")]
    un = g.unit.replace("degC", "°C")
    print(v, "gradient ours", g.gradient_ours, "ref", g.gradient_ref, un)
    hL, lL = ax.get_legend_handles_labels(); ax.grid(lw=0.3, alpha=0.5)
for k, ax in enumerate(axes.flat): panel_label(ax, "abcd"[k], x=-0.3, y=1.03)
fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.legend(hL, lL, loc="lower center", ncol=2, frameon=False); save(fig, "Fig10_normals_vs_crespi_brunetti")
