"""Point check at MeteoSwiss Col du Grand St-Bernard (GSB, 2472 m; MeteoSwiss Open Data, rka150d0 00-00 UTC,
tre200d0 24-h mean), 1991-2020.
  v3 (without GSB): GSB was not in the earlier network -> independent estimate
  v4 LOO          : GSB withheld in the leave-one-out run of the final product
  v4 (with GSB)   : nearest 1-km cell of the final product (0.44 km from the station)
(a) monthly mean precipitation; (b) annual totals (years with >= 330 valid days); (c) monthly mean temperature
difference product minus station. References that assimilate GSB (EMO, APGD, E-OBS) are dashed.
Sources: <V3>/gsb_check/gsb_daily.parquet, <V4>/gsb_check/gsb_daily.parquet, <V4>/{precip,temp}_anomaly/loo_*.parquet."""
import glob
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from figlib import V3, V4, DOUBLE, C, style, panel_label, save

def own_col(df, var):
    c = [k for k in df.columns if k.startswith(f"{var}_PoGrid") and k.endswith("_cell")]; return df[c[0]]
def loo(var):
    fs = sorted(glob.glob(str(V4 / f"{var}_anomaly" / f"loo_{var}_*.parquet")))
    d = pd.concat([pd.read_parquet(f, columns=["date", "station", "pred"]) for f in fs])
    d = d[d.station.astype(str) == "CH_GSB"]                     # filter in pandas (parquet filters are ignored by some engines)
    d.index = pd.to_datetime(d.date); s = d.pred.astype(float)
    return s[~s.index.duplicated(keep="first")].sort_index()

g4 = pd.read_parquet(V4 / "gsb_check" / "gsb_daily.parquet"); g3 = pd.read_parquet(V3 / "gsb_check" / "gsb_daily.parquet")
P = {"Station (MeteoSwiss)": (g4["obs_rka150d0"], C["black"], "-", 2.0),
     "Leave-one-out (GSB withheld)": (loo("precip").reindex(g4.index), C["vermil"], "-", 1.2),
     "GEOframe-PoMet (GSB included)": (own_col(g4, "P"), C["blue"], "-", 1.2), "Manara et al. (2026)": (g4.get("P_IH-GAR_cell"), C["green"], "-", 0.8),
     "ARCIS": (g4.get("P_ARCIS31_cell"), C["purple"], "-", 0.8), "APGD": (g4.get("P_RAPD_cell"), C["sky"], "--", 0.8),
     "EMO-1arcmin": (g4.get("P_EMO_cell"), C["grey"], "--", 0.8), "E-OBS": (g4.get("P_E-OBS_cell"), C["orange"], "--", 0.8)}
T = {"Leave-one-out (GSB withheld)": (loo("temp").reindex(g4.index), C["vermil"], "-"),
     "GEOframe-PoMet (GSB included)": (own_col(g4, "T"), C["blue"], "-"), "ARCIS (Tmax+Tmin)/2": (g4.get("T_ARCIS_T_cell"), C["purple"], "-"),
     "EMO-1arcmin": (g4.get("T_EMO_cell"), C["grey"], "--"), "E-OBS": (g4.get("T_E-OBS_cell"), C["orange"], "--")}
ok = g4["obs_rka150d0"].notna(); cnt = ok.groupby(g4.index.year).sum(); full = cnt[cnt >= 330].index
style()
fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 105 / 25.4))
for lab, (s, col, ls, lw) in P.items():
    if s is None: continue
    s = s.where(ok); m = s.groupby([s.index.year, s.index.month]).sum(min_count=25).groupby(level=1).mean()
    axes[0].plot(m.index, m.values, color=col, ls=ls, lw=lw, label=lab)
    n = s.notna().groupby(s.index.year).sum(); a = s.groupby(s.index.year).sum().where(n >= 330).reindex(full)
    axes[1].plot(a.index, a.values, color=col, ls=ls, lw=lw)
axes[0].set_xticks(range(1, 13)); axes[0].set_xticklabels("JFMAMJJASOND"); axes[0].set_ylabel("Precipitation (mm month$^{-1}$)")
axes[1].set_ylabel("Annual precipitation (mm)"); axes[1].set_xlabel("Year"); axes[0].set_xlabel("Month")
to = g4["obs_tre200d0"]
for lab, (s, col, ls) in T.items():
    if s is None: continue
    e = (s - to).dropna(); m = e.groupby(e.index.month).mean(); axes[2].plot(m.index, m.values, color=col, ls=ls, label=lab)
axes[2].axhline(0, color=C["black"], lw=0.6); axes[2].set_xticks(range(1, 13)); axes[2].set_xticklabels("JFMAMJJASOND")
axes[2].set_ylabel("Product minus station (°C)"); axes[2].set_xlabel("Month")
for k, ax in enumerate(axes): ax.grid(lw=0.3, alpha=0.5); panel_label(ax, "abc"[k], x=-0.35, y=1.02)
fig.tight_layout(rect=(0, 0.38, 1, 1), w_pad=2.5)
axes[0].legend(frameon=False, loc="upper center", ncol=3, bbox_to_anchor=(1.2, -0.25), title="Precipitation (a, b)")
axes[2].legend(frameon=False, loc="upper center", ncol=1, bbox_to_anchor=(0.5, -0.25), title="Temperature (c)")
save(fig, "FigS6_gran_san_bernardo")
