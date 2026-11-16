"""Reviewer 1, point 1a - effective resolution.
Distance (km) from each 1-km cell to the farthest kriging station whose weight is >= 2/3, 1/2 and 1/3 of the
largest positive weight of that cell-day (R_pos_2_3, R_pos_half, R_pos_1_3), per region and decade. Points = median over all cell-days,
bars = interquartile range. Source: <V4>/summary/support_radii_by_region_period_<var>_anomaly.csv (s4_summarise.py)."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from figlib import V4, DOUBLE, PERIOD_COL, REGION_LAB, style, panel_label, save

SUPPORTS = [("R_pos_2_3", "Weight ≥ 2/3 of maximum"), ("R_pos_half", "Weight ≥ 1/2 of maximum"), ("R_pos_1_3", "Weight ≥ 1/3 of maximum")]
REGIONS = ["All_PRD", "Aosta_box", "Valtellina_box", "Rest_PRD"]; PERIODS = ["1991-2000", "2001-2010", "2011-2020"]
style()
fig, axes = plt.subplots(2, 3, figsize=(DOUBLE, 120 / 25.4), sharey="row")
for i, (v, vlab) in enumerate((("precip", "Precipitation"), ("temp", "Temperature"))):
    d = pd.read_csv(V4 / "summary" / f"support_radii_by_region_period_{v}_anomaly.csv")
    for j, (s, slab) in enumerate(SUPPORTS):
        ax = axes[i, j]
        for k, per in enumerate(PERIODS):
            q = d[(d.support == s) & (d.period == per)].set_index("region").reindex(REGIONS)
            x = np.arange(len(REGIONS)) + (k - 1) * 0.22
            ax.errorbar(x, q.p50, yerr=[q.p50 - q.p25, q.p75 - q.p50], fmt="o", ms=3.5, lw=0.8, capsize=0,
                        color=PERIOD_COL[per], label=per.replace("-", "–"))
        ax.set_xticks(range(len(REGIONS))); ax.set_xticklabels([REGION_LAB[r].replace(" of ", " of\n").replace("Whole ", "Whole\n").replace("Aosta ", "Aosta\n") for r in REGIONS])
        ax.set_ylim(0, 1.05 * d[d.support.isin([x for x, _ in SUPPORTS]) & d.region.isin(REGIONS) & d.period.isin(PERIODS)].p75.max()); ax.grid(axis="y", lw=0.3, alpha=0.5)
        if i == 0: ax.set_title(slab)
        if j == 0: ax.set_ylabel(f"{vlab}, distance (km)")
        panel_label(ax, "abcdef"[i * 3 + j])
h, l = axes[0, 0].get_legend_handles_labels()
fig.tight_layout(rect=(0, 0.06, 1, 1)); fig.legend(h, l, frameon=False, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0))
save(fig, "Fig06_effective_resolution")
