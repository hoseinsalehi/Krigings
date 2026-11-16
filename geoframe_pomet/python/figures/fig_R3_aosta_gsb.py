"""Reviewer 1, point 2 - Aosta (~580 m) - Gran San Bernardo (2472 m) point test.
(a) terrain and stations of the transect area; (b) mean annual precipitation and (c) mean annual temperature at
1/4, 1/2, 3/4 of the transect and at the Gran San Bernardo point, with all stations, without Aosta, without
Gran San Bernardo, without both (anomaly mode; normals rebuilt for every scenario). Open symbols: station
observations (Aosta: network station; GSB: MeteoSwiss rka150d0 / tre200d0, years with >= 330 valid days).
Sources: <V4>/point_test_<var>/{summary,points,stations_used}.csv, <V4>/normals/station_normals_<var>.csv,
<V4>/gsb_check/gsb_daily.parquet, <V4>/grid_predictors.npz, station shapefiles of the v4 network."""
import numpy as np, pandas as pd, geopandas as gpd, matplotlib.pyplot as plt
from pyproj import Transformer
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from figlib import V4, DOUBLE, C, style, panel_label, save

V4IN = V4.parent / "PoGrid_v4_input"
SCEN = [("full", "All stations", C["black"], "-"), ("no_Aosta", "Without Aosta", C["sky"], "--"),
        ("no_GSB", "Without GSB", C["vermil"], "-"), ("no_both", "Without both", C["orange"], ":")]
FRAC = {"AO_GSB_25": 0.25, "AO_GSB_50": 0.50, "AO_GSB_75": 0.75, "AO_GSB_100": 1.00}
tr = Transformer.from_crs(4326, 3035, always_xy=True)
style()
fig = plt.figure(figsize=(DOUBLE, 70 / 25.4)); gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1], wspace=0.45)

# (a) map
ax = fig.add_subplot(gs[0]); g = np.load(V4 / "grid_predictors.npz"); x, y, dem = g["x"], g["y"], g["dem2d"]
pts = pd.read_csv(V4 / "point_test_precip" / "points.csv"); su = pd.read_csv(V4 / "point_test_precip" / "stations_used.csv")
lon = np.r_[7.320, pts.lon.values]; lat = np.r_[45.737, pts.lat.values]; X, Y = tr.transform(lon, lat)
x0, x1, y0, y1 = X.min() - 15e3, X.max() + 15e3, Y.min() - 12e3, Y.max() + 12e3
ix = (x >= x0) & (x <= x1); iy = (y >= y0) & (y <= y1)
im = ax.imshow(dem[np.ix_(iy, ix)], extent=[x[ix].min() / 1e3, x[ix].max() / 1e3, y[iy].min() / 1e3, y[iy].max() / 1e3],
               origin="upper" if y[0] > y[-1] else "lower", cmap="Greys", vmin=300, vmax=4000)
st = gpd.read_file(V4IN / "STATIONS_plusCH_precip.shp").to_crs(3035)
st = st[(st.geometry.x.between(x0, x1)) & (st.geometry.y.between(y0, y1))]
ax.scatter(st.geometry.x / 1e3, st.geometry.y / 1e3, s=6, color=C["blue"], lw=0, label="Precipitation stations")
ax.plot(X / 1e3, Y / 1e3, color=C["vermil"], lw=1); ax.scatter(X[1:] / 1e3, Y[1:] / 1e3, s=14, color=C["vermil"], zorder=4, label="Test points")
for s_, (xx, yy) in zip(["Aosta", "GSB"], [(X[0], Y[0]), (X[-1], Y[-1])]):
    ax.annotate(s_, (xx / 1e3, yy / 1e3), xytext=(4, 3), textcoords="offset points", fontsize=7)
ax.set_xlabel("Easting (km)"); ax.set_ylabel("Northing (km)"); ax.set_aspect("equal")
cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.05, pad=0.25); cb.set_label("Elevation (m)")
panel_label(ax, "a", x=-0.3); hA, lA = ax.get_legend_handles_labels()

# observed references
gd = pd.read_parquet(V4 / "gsb_check" / "gsb_daily.parquet")
def annual_obs(s):
    n = s.groupby(s.index.year).count(); a = s.groupby(s.index.year).sum(); return a[n >= 330].mean()
def mean_obs(s):
    n = s.groupby(s.index.year).count(); a = s.groupby(s.index.year).mean(); return a[n >= 330].mean()
obs = {"precip": {}, "temp": {}}
obs["precip"][1.0] = annual_obs(gd["obs_rka150d0"]); obs["temp"][1.0] = mean_obs(gd["obs_tre200d0"])
for v in ("precip", "temp"):
    u = pd.read_csv(V4 / f"point_test_{v}" / "stations_used.csv"); aid = str(u.loc[u.site == "Aosta", "station"].iloc[0])
    n = pd.read_csv(V4 / "normals" / f"station_normals_{v}.csv", dtype={"station": str}); n = n[n.station == aid]
    obs[v][0.0] = n.normal.sum() if v == "precip" else n.normal.mean()

for k, (v, lab) in enumerate((("precip", "Mean annual precipitation (mm)"), ("temp", "Mean annual temperature (°C)"))):
    ax = fig.add_subplot(gs[k + 1]); s = pd.read_csv(V4 / f"point_test_{v}" / "summary.csv"); s = s[s["mode"] == "anomaly"]
    for sc, sl, col, ls in SCEN:
        q = s[s.scenario == sc].assign(f=lambda d: d["name"].map(FRAC)).sort_values("f")
        ax.plot(q.f, q["annual_total" if v == "precip" else "mean"], ls=ls, color=col, marker="o", ms=3, label=sl)
    for f, val in obs[v].items():
        ax.scatter([f], [val], s=30, facecolor="none", edgecolor=C["green"], lw=1.2, zorder=5,
                   label="Station observation" if f == 0.0 else None)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0]); ax.set_xticklabels(["Aosta", "1/4", "1/2", "3/4", "GSB"])
    ax.set_xlabel("Position along transect"); ax.set_ylabel(lab); ax.grid(lw=0.3, alpha=0.5); panel_label(ax, "bc"[k], x=-0.25)
    if k == 0: hB, lB = ax.get_legend_handles_labels()
fig.legend(hA + hB, lA + lB, loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=4, frameon=False)
save(fig, "FigS4_aosta_gsb_point_test")
