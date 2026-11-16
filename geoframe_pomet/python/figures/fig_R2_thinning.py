"""Reviewer 1, point 2 - robustness to missing data: 50 % random station deletion, 30 realisations.
Each realisation withholds half of the stations; the withheld stations are predicted from the other half with the
anomaly (climatology-aided) and the direct mode. Lines join the two modes of the same realisation.
Source: <V4>/thinning_<var>/withheld_by_realization.csv (s5_thinning.py)."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from figlib import V4, DOUBLE, MODE_COL, C, style, panel_label, save

PANELS = [("precip", "mae", "Precipitation\ndaily MAE (mm)"),
          ("precip", "stmean_abs_relbias_median_pct", "Precipitation\nstation-mean |bias| (%)"),
          ("temp", "mae", "Temperature\ndaily MAE (°C)"),
          ("temp", "stmean_abs_bias_median_c", "Temperature\nstation-mean |bias| (°C)")]
style()
fig, axes = plt.subplots(1, 4, figsize=(DOUBLE, 62 / 25.4))
for k, (v, col, lab) in enumerate(PANELS):
    d = pd.read_csv(V4 / f"thinning_{v}" / "withheld_by_realization.csv").pivot(index="realization", columns="mode", values=col)
    ax = axes[k]
    for _, r in d.iterrows(): ax.plot([0, 1], [r["anomaly"], r["direct"]], color=C["grey"], lw=0.4, alpha=0.6)
    for x, m in enumerate(("anomaly", "direct")):
        ax.scatter(np.full(len(d), x), d[m], s=8, color=MODE_COL[m], zorder=3)
        ax.plot([x - 0.18, x + 0.18], [d[m].median()] * 2, color=MODE_COL[m], lw=2, zorder=4)
    ax.set_xticks([]); ax.set_xlim(-0.5, 1.5); ax.set_ylabel(lab)
    worse = (d["direct"] > d["anomaly"]).sum()
    print(v, col, "medians", d.median().round(3).to_dict(), "direct worse", worse, "/", len(d))
    panel_label(ax, "abcd"[k], x=-0.4, y=1.02)
from matplotlib.lines import Line2D
h = [Line2D([], [], marker="o", ls="-", lw=2, ms=4, color=MODE_COL[m]) for m in ("anomaly", "direct")]
fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.legend(h, ["GEOframe-PoMet", "Direct interpolation of daily values (without the climatological layer)"], loc="lower center", ncol=2, frameon=False)
save(fig, "Fig07_missing_data")
