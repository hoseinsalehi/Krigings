"""Figures 8 (precipitation comparison) and 9 (temperature comparison) of the revised manuscript.
Inputs (written by fig08_09_compute.py): ../data/fig08_09_fields.npz (1-km fields of every product, masked to
~/Downloads/Po_District_clean.shp), ../data/fig08_station_climatology.csv, ../data/fig08_storm_alex_stations.csv,
../data/fig09_station_values.csv, ../stats_fig08.csv, ../stats_fig09.csv.
Output: ../Fig08_precip_comparison.{pdf,png}, ../Fig09_temp_comparison.{pdf,png} (300 dpi)."""
from pathlib import Path
import string
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, ListedColormap

HERE = Path(__file__).resolve().parent; OUT = HERE.parent; D = OUT / "data"
BND = gpd.read_file(Path.home() / "Downloads/Po_District_clean.shp").to_crs(3035)
Z = np.load(D / "fig08_09_fields.npz"); EXT = Z["EXT"]; CMASK = Z["CMASK"]
MM = 1 / 25.4
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman"], "mathtext.fontset": "stix",
                     "font.size": 9, "savefig.dpi": 300, "pdf.fonttype": 42, "axes.unicode_minus": True})

bx, by = [], []
for g in BND.geometry:
    for p in getattr(g, "geoms", [g]):
        for ring in [p.exterior] + list(p.interiors):
            x, y = ring.xy; bx += list(x) + [np.nan]; by += list(y) + [np.nan]
bx, by = np.array(bx) / 1e3, np.array(by) / 1e3


def labels():
    a = list(string.ascii_lowercase); return a + ["a" + c for c in a]


def frame(ax):
    ax.plot(bx, by, color="k", lw=0.35, zorder=4)
    ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3]); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)


def grey(ax, shade=0.90):
    ax.imshow(np.where(CMASK, 1.0, np.nan), extent=EXT, origin="upper", cmap=ListedColormap([(shade,) * 3]),
              interpolation="nearest", zorder=0)


def corner(ax, txt):
    ax.text(0.995, 0.985, txt, transform=ax.transAxes, ha="right", va="top", fontsize=7.5, linespacing=1.05, zorder=6)


def plot_grid(rows, cols, fname, stat_txt, station_size):
    """rows: list of (key, row label); cols: list of (key, title, norm, cmap, ticks, unit)."""
    nr, nc = len(rows), len(cols)
    w = 170 * MM; cw = w * 0.94 / nc; ch = cw * (EXT[3] - EXT[2]) / (EXT[1] - EXT[0])
    top, bot_cb, left = 0.32, 0.62, 0.06 * w
    h = top + nr * ch * 1.04 + bot_cb
    fig = plt.figure(figsize=(w, h)); L = labels(); n = 0; ims = {}
    for r, (rk, rl) in enumerate(rows):
        for c, (ck, title, norm, cmap, ticks, unit) in enumerate(cols):
            x0 = (left + c * cw) / w; y0 = (h - top - (r + 1) * ch * 1.04) / h
            ax = fig.add_axes([x0, y0, cw / w * 0.985, ch / h])
            if rk == "obs":
                grey(ax); s = station(ck)
                o = np.argsort(s.v.values) if ck in ("Alex",) else np.arange(len(s))
                im = ax.scatter(s.x.values[o] / 1e3, s.y.values[o] / 1e3, c=s.v.values[o], s=station_size[ck], cmap=cmap,
                                norm=norm, linewidths=0, zorder=3)
            else:
                k = f"{rk}_{ck}"
                if k not in Z.files:
                    ax.axis("off"); n += 0; continue          # no data (APGD has no 2020): panel left empty
                im = ax.imshow(Z[k], extent=EXT, origin="upper", cmap=cmap, norm=norm, interpolation="nearest", zorder=1)
                t = stat_txt(rk, ck)
                if t: corner(ax, t)
            ims[c] = im; frame(ax)
            ax.text(0.01, 0.99, f"({L[n]})", transform=ax.transAxes, ha="left", va="top", fontsize=8.5, zorder=6); n += 1
            if r == 0: ax.set_title(title, fontsize=9, pad=3)
            if c == 0: ax.text(-0.03, 0.5, rl, transform=ax.transAxes, rotation=90, ha="right", va="center", fontsize=9)
    for c, (ck, title, norm, cmap, ticks, unit) in enumerate(cols):
        x0 = (left + c * cw) / w; cax = fig.add_axes([x0 + 0.1 * cw / w, (bot_cb - 0.20) / h, 0.78 * cw / w, 0.075 / h])
        cb = fig.colorbar(ims[c], cax=cax, orientation="horizontal", ticks=ticks, extend="both")
        cb.ax.set_xticklabels([f"{t:g}".replace("-", "\u2212") for t in ticks]); cb.ax.minorticks_off()
        cb.ax.tick_params(labelsize=8, length=2, pad=1); cb.set_label(unit, fontsize=8.5, labelpad=1)
    for ext in ("png", "pdf"): fig.savefig(OUT / f"{fname}.{ext}", dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig); print("->", OUT / f"{fname}.png")


def m(v):
    v = round(float(v), 1)
    return "0.0" if v == 0 else f"{v:+.1f}".replace("-", "\u2212")


# ------------------------------------------------------------------ Fig. 8 precipitation
S8 = pd.read_csv(OUT / "stats_fig08.csv")
SC = pd.read_csv(D / "fig08_station_climatology.csv"); SA = pd.read_csv(D / "fig08_storm_alex_stations.csv")
PN = {"tp": "GEOframe-PoMet", "IH-GAR": "Manara et al. (2026)", "RAPD": "APGD", "ARCIS31": "ARCIS", "EMO": "EMO-1arcmin", "E-OBS": "E-OBS"}
CN = {"P": "mean_annual_precipitation_mm", "F": "wet_day_frequency_pct", "I": "wet_day_intensity_mm_per_wet_day",
      "Alex": "storm_alex_total_2_5_Oct_2020_mm"}


def st8(rk, ck):
    q = S8[(S8["product"] == PN[rk]) & (S8.column == CN[ck])]
    if q.empty or q.n_stations.iloc[0] == 0: return ""
    q = q.iloc[0]
    if ck == "Alex": return f"{m(q.total_bias_pct)} %\n{q.mae_mm:.1f} mm"
    return f"{m(q.median_bias_pct)} %\n{q.median_abs_error_pct:.1f} %"


def station8(ck):
    if ck == "Alex": return pd.DataFrame(dict(x=SA.x, y=SA.y, v=SA.obs))
    return pd.DataFrame(dict(x=SC.x, y=SC.y, v=SC[ck]))


cols8 = [("P", "Mean annual precipitation", LogNorm(400, 3000), "YlGnBu", [500, 1000, 2000], "mm yr$^{-1}$"),
         ("F", "Wet-day frequency (\u2265 1 mm)", Normalize(15, 40), "PuBuGn", [15, 20, 25, 30, 35, 40], "% of days"),
         ("I", "Mean wet-day intensity", Normalize(6, 20), "YlOrRd", [6, 10, 14, 18], "mm per wet day"),
         ("Alex", "Storm Alex, 2\u20135 Oct 2020", LogNorm(10, 600), "BuPu", [10, 30, 100, 300], "mm")]
rows8 = [("obs", "Stations"), ("tp", "GEOframe-PoMet"), ("IH-GAR", "Manara et al. (2026)"), ("RAPD", "APGD"), ("ARCIS31", "ARCIS"),
         ("EMO", "EMO-1arcmin"), ("E-OBS", "E-OBS")]
station = station8
plot_grid(rows8, cols8, "Fig08_precip_comparison", st8, {"P": 1.8, "F": 1.8, "I": 1.8, "Alex": 2.0})

# ------------------------------------------------------------------ Fig. 9 temperature
S9 = pd.read_csv(OUT / "stats_fig09.csv"); TS = pd.read_csv(D / "fig09_station_values.csv")
PN9 = {"v2": "GEOframe-PoMet", "ARCIS_T": "ARCIS", "EMO": "EMO-1arcmin", "E-OBS": "E-OBS"}
CN9 = {"Tann": "mean_annual_T_1991_2020", "Tjja": "mean_T_JJA_2003", "Tcold": "daily_mean_T_1991-02-07", "Thot": "daily_mean_T_2019-06-27"}


def st9(rk, ck):
    q = S9[(S9["product"] == PN9[rk]) & (S9.column == CN9[ck])].iloc[0]
    return f"{m(q.median_bias_degC)} \u00b0C\n{q.median_abs_error_degC:.1f} \u00b0C"


def station9(ck):
    s = TS[np.isfinite(TS[ck])]; return pd.DataFrame(dict(x=s.x, y=s.y, v=s[ck]))


cols9 = [("Tann", "Mean annual, 1991\u20132020", Normalize(-2, 16), "RdYlBu_r", [0, 4, 8, 12, 16], "\u00b0C"),
         ("Tjja", "Summer (JJA) 2003", Normalize(8, 28), "RdYlBu_r", [8, 12, 16, 20, 24, 28], "\u00b0C"),
         ("Tcold", "7 February 1991", Normalize(-22, -2), "RdYlBu_r", [-22, -17, -12, -7, -2], "\u00b0C"),
         ("Thot", "27 June 2019", Normalize(14, 32), "RdYlBu_r", [14, 18, 22, 26, 30], "\u00b0C")]
rows9 = [("obs", "Stations"), ("v2", "GEOframe-PoMet"), ("ARCIS_T", "ARCIS"), ("EMO", "EMO-1arcmin"), ("E-OBS", "E-OBS")]
station = station9
plot_grid(rows9, cols9, "Fig09_temp_comparison", st9, {"Tann": 6, "Tjja": 3, "Tcold": 6, "Thot": 1.5})
