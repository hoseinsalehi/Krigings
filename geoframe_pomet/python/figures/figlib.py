"""Shared paths, style and map helpers for the figures of the revised manuscript (dataset v2.0).
Environment: PO_V4 (product root, default ~/Downloads/PoGrid_v5), FIG_OUT (default ../ relative to this file).
District boundary: ~/Downloads/Po_District_clean.shp (EPSG:3035, sliver holes removed)."""
import os
from pathlib import Path
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HOME = Path.home()
V5 = Path(os.environ.get("PO_V4", "output"))
V3 = HOME / "Downloads/PoGrid_v3"
HERE = Path(__file__).resolve().parent
OUT = Path(os.environ.get("FIG_OUT", HERE.parent)); OUT.mkdir(parents=True, exist_ok=True)
TAB = HERE.parent / "data"           # station-level tables (precip LOO rebuilt from the final product)
BND = gpd.read_file(HOME / "Downloads/Po_District_clean.shp").to_crs(3035)
G = np.load(V5 / "grid_predictors.npz")
EXT = [G["x"].min() / 1e3 - 0.5, G["x"].max() / 1e3 + 0.5, G["y"].min() / 1e3 - 0.5, G["y"].max() / 1e3 + 0.5]
ORIG = "upper" if G["y"][0] > G["y"][-1] else "lower"
MM = 1 / 25.4
SINGLE, DOUBLE = 83 * MM, 170 * MM
C = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73", yellow="#F0E442",
         blue="#0072B2", vermil="#D55E00", purple="#CC79A7", grey="#888888")
PERIOD_COL = {"1991-2000": C["sky"], "2001-2010": C["blue"], "2011-2020": C["black"]}
REGION_LAB = {"All_PRD": "Whole District", "Aosta_box": "Aosta Valley", "Valtellina_box": "Valtellina", "Rest_PRD": "Rest of District"}
BOXES = {"Aosta_box": (6.75, 7.95, 45.42, 46.02), "Valtellina_box": (9.15, 10.55, 45.90, 46.65)}   # lon0, lon1, lat0, lat1 (po_config.py)
VLAB = {"precip": "Precipitation", "temp": "Temperature"}


def style():
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman"], "mathtext.fontset": "stix",
                         "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
                         "legend.fontsize": 8, "legend.title_fontsize": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "lines.linewidth": 1.0,
                         "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42, "axes.unicode_minus": True})


def panel_label(ax, s, x=-0.12, y=1.02, **kw):
    ax.text(x, y, f"({s})", transform=ax.transAxes, fontweight="bold", fontsize=10, va="bottom", ha="left", **kw)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight", pad_inches=0.03, dpi=300)
    plt.close(fig); print("->", OUT / f"{name}.png")


def to2d(vals):
    a = np.full(G["dem2d"].shape, np.nan, np.float32); a[G["iy"], G["ix"]] = vals; return a


MASK = np.isfinite(to2d(np.ones(len(G["elev"]))))


def bnd_xy():
    xs, ys = [], []
    for geom in BND.geometry:
        for poly in (geom.geoms if hasattr(geom, "geoms") else [geom]):
            for ring in [poly.exterior] + list(poly.interiors):
                x, y = ring.xy; xs += list(x) + [np.nan]; ys += list(y) + [np.nan]
    return np.array(xs) / 1e3, np.array(ys) / 1e3


def basemap(ax, grey=True, ticks=False, pad=10):
    if grey:
        ax.imshow(to2d(G["elev"]), extent=EXT, origin=ORIG, cmap="Greys", vmin=-1500, vmax=5000, zorder=0, interpolation="nearest")
    bx, by = bnd_xy(); ax.plot(bx, by, color="k", lw=0.5, zorder=3)
    ax.set_aspect("equal"); ax.set_xlim(EXT[0] - pad, EXT[1] + pad); ax.set_ylim(EXT[2] - pad, EXT[3] + pad)
    if ticks:
        km_axes(ax)
    else:
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)


def km_axes(ax, xlabel=True, ylabel=True):
    for s in ax.spines.values(): s.set_visible(True)
    ax.set_xticks(np.arange(4100, 4600, 100)); ax.set_yticks(np.arange(2300, 2650, 100))
    if xlabel: ax.set_xlabel("Easting (km)")
    if ylabel: ax.set_ylabel("Northing (km)")
    ax.tick_params(labelsize=7.5)


def raster(ax, a, **kw):
    return ax.imshow(np.where(MASK, a, np.nan), extent=EXT, origin=ORIG, interpolation="nearest", zorder=1, **kw)


def box_xy(name, n=60):
    from pyproj import Transformer
    lo0, lo1, la0, la1 = BOXES[name]; t = Transformer.from_crs(4326, 3035, always_xy=True)
    lon = np.r_[np.linspace(lo0, lo1, n), np.full(n, lo1), np.linspace(lo1, lo0, n), np.full(n, lo0)]
    lat = np.r_[np.full(n, la0), np.linspace(la0, la1, n), np.full(n, la1), np.linspace(la1, la0, n)]
    X, Y = t.transform(lon, lat); return X / 1e3, Y / 1e3

# aliases for the reviewer-response scripts (fig_R*.py)
V4 = V5
MODE_COL = {"anomaly": C["blue"], "direct": C["orange"]}
