"""Fig. 1 - study area: hillshaded relief (Dem_Po_250m.tif), clean District boundary only, sea and lakes in the same
blue (sea = map frame minus European land polygons, Europe_coastline_poly.shp; lakes > 10 km2 from HydroLAKES),
evaluation sub-regions (Aosta Valley, Valtellina), location inset (Italy and neighbours) inside the frame."""
import numpy as np, geopandas as gpd, rasterio, matplotlib.pyplot as plt
from matplotlib.colors import LightSource, LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.patheffects as pe
from shapely.geometry import box as sbox
from figlib import *

style()
WATER = "#7fb2de"
XL, YL = (4040, 4560), (2270, 2632)
R = np.load(HERE.parent / "data" / "relief_250m_frame.npz")   # Dem_Po_250m.tif inside its coverage, Copernicus GLO-90 (data/cop90) outside, 250 m, EPSG:3035
z = R["z"].astype(float); ext = list(R["ext"])
frame = sbox(XL[0] * 1e3 - 2e4, YL[0] * 1e3 - 2e4, XL[1] * 1e3 + 2e4, YL[1] * 1e3 + 2e4)
land = gpd.read_file(HOME / "Downloads/Europe_coastline_shapefile/Europe_coastline_poly.shp", bbox=frame.bounds).to_crs(3035)
sea = gpd.GeoSeries([frame.difference(land.union_all())], crs=3035)
lakes = gpd.read_file(HOME / "Downloads/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp", bbox=(5.5, 43.0, 14.0, 47.5))
lakes = lakes[lakes.Lake_area > 10].to_crs(3035)
km = lambda g: g.scale(1e-3, 1e-3, origin=(0, 0))

terrain = LinearSegmentedColormap.from_list("relief", [(0.00, "#1fb51f"), (0.05, "#4cc236"), (0.15, "#a6d65a"), (0.25, "#ece27e"), (0.40, "#e3ad5e"), (0.58, "#bf7a3c"), (0.75, "#9a6440"), (0.88, "#d9c7b5"), (1.00, "#ffffff")])
fig = plt.figure(figsize=(DOUBLE, 122 * MM)); ax = fig.add_axes([0.07, 0.15, 0.80, 0.83])
ls = LightSource(azdeg=315, altdeg=40)
zz = np.nan_to_num(z, nan=0.0)
base = terrain(np.clip(zz / 4200, 0, 1))[..., :3]
hs = ls.hillshade(zz, vert_exag=3.0, dx=250, dy=250)          # 0-1, NW illumination
flat = ls.hillshade(np.zeros((3, 3)), vert_exag=3.0, dx=250, dy=250)[1, 1]
shade = np.clip(1.0 + 0.95 * (hs - flat), 0.35, 1.25)          # multiplicative, plains ~unchanged
rgb = np.clip(base * shade[..., None], 0, 1)
rgb[~np.isfinite(z)] = 1.0
# light white veil outside the District so the District relief stands out
from rasterio import features
from affine import Affine
ny, nx = z.shape; tf = Affine((ext[1] - ext[0]) * 1e3 / nx, 0, ext[0] * 1e3, 0, -(ext[3] - ext[2]) * 1e3 / ny, ext[3] * 1e3)
inside = features.rasterize([(g, 1) for g in BND.to_crs(3035).geometry], out_shape=z.shape, transform=tf, fill=0).astype(bool)
rgb[~inside] = 0.50 * rgb[~inside] + 0.50
ax.imshow(rgb, extent=ext, zorder=0, interpolation="bilinear")
km(sea).plot(ax=ax, color=WATER, lw=0, zorder=2)
km(lakes.geometry).plot(ax=ax, color=WATER, lw=0, zorder=2)
bx, by = bnd_xy(); ax.plot(bx, by, color="k", lw=1.1, zorder=5)
for nm in ("Aosta_box", "Valtellina_box"):
    X, Y = box_xy(nm); ax.plot(X, Y, color=C["vermil"], lw=1.1, zorder=6)
from pyproj import Transformer
tr = Transformer.from_crs(4326, 3035, always_xy=True)
for lon, lat, s in ((7.35, 46.05, "Aosta Valley"), (9.85, 45.80, "Valtellina")):
    x, y = tr.transform(lon, lat); ax.text(x / 1e3, y / 1e3, s, color=C["vermil"], fontsize=8, ha="center", va="top" if s == "Valtellina" else "bottom", fontweight="bold", zorder=7, path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])
km_axes(ax); ax.set_xlim(*XL); ax.set_ylim(*YL); ax.set_aspect("equal")
for s in ax.spines.values(): s.set_visible(True)
# colour bar (right, outside)
sm = plt.cm.ScalarMappable(cmap=terrain, norm=plt.Normalize(0, 4200))
cax = fig.add_axes([0.885, 0.25, 0.017, 0.62]); cb = fig.colorbar(sm, cax=cax, extend="max"); cb.set_label("Elevation (m a.s.l.)")
# legend (below, outside)
h = [Line2D([], [], color="k", lw=1.1), Patch(facecolor=WATER, lw=0), Line2D([], [], color=C["vermil"], lw=1.1)]
fig.legend(h, ["Po River District", "Sea and lakes", "Evaluation sub-regions"], loc="lower center", bbox_to_anchor=(0.47, 0.035), ncol=3, frameon=False)
ax.apply_aspect(); fig.canvas.draw()
# inset in the empty upper-right corner (outside the District)
x0, y0 = ax.transData.transform((4447, 2527)); x1, y1 = ax.transData.transform((4558, 2630))
inv = fig.transFigure.inverted(); (fx0, fy0), (fx1, fy1) = inv.transform((x0, y0)), inv.transform((x1, y1))
ia = fig.add_axes([fx0, fy0, fx1 - fx0, fy1 - fy0])
W = gpd.read_file(HOME / "Downloads/world/TM_WORLD_BORDERS-0.3.shp")
W = W[W.REGION == 150].to_crs(3035)
ia.set_facecolor(WATER)
km(W[W.ISO3 != "ITA"].geometry).plot(ax=ia, color="0.90", edgecolor="0.55", lw=0.25)
km(W[W.ISO3 == "ITA"].geometry).plot(ax=ia, color="0.70", edgecolor="0.45", lw=0.25)
km(BND.geometry).plot(ax=ia, color=C["vermil"], lw=0)
ia.set_xlim(3750, 5150); ia.set_ylim(2850 - 1400 * (y1 - y0) / (x1 - x0), 2850); ia.set_xticks([]); ia.set_yticks([]); ia.set_xlabel(""); ia.set_ylabel(""); ia.set_aspect("auto")
for s in ia.spines.values(): s.set_visible(True); s.set_linewidth(0.7)
assert BND.union_all().intersection(sbox(4447e3, 2527e3, 4558e3, 2630e3)).area == 0
save(fig, "Fig01_study_area")
