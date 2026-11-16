"""
s12_isac_normals_comparison.py -- compare the 1-km 1991-2020 normals (s2) with the ISAC-CNR 30-arcsec 1961-1990
annual climatologies: precipitation (Crespi et al. 2018, doi:10.1002/joc.5217) and mean temperature
(Brunetti et al. 2014, doi:10.1002/joc.3764). Source: https://www.isac.cnr.it/climstor/CLIMATE_DATA/ (CC BY-NC-ND 4.0;
direct-access int32 binary, tenths of mm / tenths of degC, 1801x1801, 5.0042E..20.0E, 50.0N..35.0N).

Because the reference period differs (1961-1990 vs 1991-2020), the comparison targets the SPATIAL PATTERN and the
ELEVATION GRADIENT (the PRISM-type physiographic structure questioned by the reviewer), not the absolute level:
  - spatial correlation of the annual fields (precip in log space as well), per region and elevation band
  - area-mean difference (precip ratio %, temperature difference; temperature includes ~1991-2020 warming)
  - elevation gradients: OLS of ln P (%/km) and T (degC/km) vs cell elevation, per region
Each 1-km cell takes the mean of the 30-arcsec cells whose centres fall within +-0.5 km (EPSG:3035).
Output: OUT_ROOT/isac_comparison/{cells.parquet, summary.csv, gradients.csv}
"""
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd, xarray as xr
import po_config as C

ISAC = Path(os.environ.get("PO_INPUT_ROOT", "data") + "/ISAC_climatologies")
N, X0, Y0, D = 1801, 5.00416666666667, 49.99583333333333, 0.00833333333333
BOXES = {"Aosta_box": (6.75, 45.42, 7.95, 46.02), "Valtellina_box": (9.15, 45.90, 10.55, 46.65)}
BANDS = ((0, 500), (500, 1000), (1000, 1500), (1500, 2000), (2000, 9000))
OD = C.OUT_ROOT / "isac_comparison"; OD.mkdir(parents=True, exist_ok=True)


def read_isac(name):
    a = np.fromfile(ISAC / name, "<i4").reshape(N, N).astype("float64"); a[a == -9999] = np.nan; return a / 10.0


def to_grid(a, g):
    """mean of 30" cells whose centres lie within the 1-km cell (EPSG:3035)"""
    from pyproj import Transformer
    lon = X0 + np.arange(N) * D; lat = Y0 - np.arange(N) * D
    x0, x1 = g["x"].min() - 1000, g["x"].max() + 1000; y0, y1 = g["y"].min() - 1000, g["y"].max() + 1000
    LO, LA = np.meshgrid(lon, lat); ok = np.isfinite(a)
    X, Y = Transformer.from_crs(4326, 3035, always_xy=True).transform(LO[ok], LA[ok]); v = a[ok]
    sel = (X >= x0) & (X <= x1) & (Y >= y0) & (Y <= y1); X, Y, v = X[sel], Y[sel], v[sel]
    dx = g["x"][1] - g["x"][0]; dy = g["y"][1] - g["y"][0]
    ix = np.rint((X - g["x"][0]) / dx).astype(int); iy = np.rint((Y - g["y"][0]) / dy).astype(int)
    ny, nx = int(g["ny"]), int(g["nx"]); m = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    s = np.zeros((ny, nx)); c = np.zeros((ny, nx))
    np.add.at(s, (iy[m], ix[m]), v[m]); np.add.at(c, (iy[m], ix[m]), 1)
    with np.errstate(invalid="ignore"): out = s / c
    return out[g["iy"], g["ix"]]


def grad(z, v, log):
    ok = np.isfinite(z) & np.isfinite(v) & ((v > 0) if log else True)
    if ok.sum() < 30 or np.ptp(z[ok]) < 300: return np.nan
    b = np.polyfit(z[ok] / 1000.0, np.log(v[ok]) if log else v[ok], 1)[0]
    return 100 * (np.exp(b) - 1) if log else b


def main():
    g = np.load(C.OUT_ROOT / "grid_predictors.npz"); z = g["elev"]; lon, lat = g["lon"], g["lat"]
    nP = xr.open_dataset(C.OUT_ROOT / "normals" / "normals_precip.nc").normal.values
    nT = xr.open_dataset(C.OUT_ROOT / "normals" / "normals_temp.nc").normal.values
    ours_P = nP[:, g["iy"], g["ix"]].sum(0); ours_T = nT[:, g["iy"], g["ix"]].mean(0)
    ref_P = to_grid(read_isac("CLINO_PCP_ITA_Y"), g); ref_T = to_grid(read_isac("CLINO_TMM_ITA_Y"), g)
    cells = pd.DataFrame(dict(lon=lon, lat=lat, elev=z, P_ours_1991_2020=ours_P, P_crespi_1961_1990=ref_P,
                              T_ours_1991_2020=ours_T, T_brunetti_1961_1990=ref_T))
    cells.to_parquet(OD / "cells.parquet", index=False)
    groups = {"All_PRD": np.ones(len(z), bool)}
    for k, (a, b, c_, d_) in BOXES.items(): groups[k] = (lon >= a) & (lon <= c_) & (lat >= b) & (lat <= d_)
    groups["Rest_PRD"] = ~(groups["Aosta_box"] | groups["Valtellina_box"])
    for lo_, hi_ in BANDS: groups[f"elev_{lo_}-{hi_}"] = (z >= lo_) & (z < hi_)
    rows, grows = [], []
    for gname, m in groups.items():
        for var, o, r in (("precip", ours_P, ref_P), ("temp", ours_T, ref_T)):
            ok = m & np.isfinite(o) & np.isfinite(r)
            if ok.sum() < 30: continue
            row = dict(variable=var, group=gname, n_cells=int(ok.sum()), mean_ours=o[ok].mean(), mean_ref=r[ok].mean(),
                       spatial_r=np.corrcoef(o[ok], r[ok])[0, 1])
            if var == "precip":
                row.update(diff_pct=100 * (o[ok].mean() / r[ok].mean() - 1), spatial_r_log=np.corrcoef(np.log(o[ok]), np.log(r[ok]))[0, 1],
                           cell_ratio_median_pct=100 * (np.median(o[ok] / r[ok]) - 1),
                           cell_ratio_p10_pct=100 * (np.percentile(o[ok] / r[ok], 10) - 1), cell_ratio_p90_pct=100 * (np.percentile(o[ok] / r[ok], 90) - 1))
            else:
                e = o[ok] - r[ok]; row.update(diff=e.mean(), diff_p10=np.percentile(e, 10), diff_p90=np.percentile(e, 90),
                                              anomaly_r=np.corrcoef(o[ok] - o[ok].mean(), r[ok] - r[ok].mean())[0, 1])
            rows.append(row)
            if not gname.startswith("elev_"):
                grows.append(dict(variable=var, group=gname, unit="% per km" if var == "precip" else "degC per km",
                                  gradient_ours=grad(z[ok], o[ok], var == "precip"), gradient_ref=grad(z[ok], r[ok], var == "precip")))
    S = pd.DataFrame(rows); G = pd.DataFrame(grows)
    S.to_csv(OD / "summary.csv", index=False); G.to_csv(OD / "gradients.csv", index=False)
    pd.set_option("display.width", 250)
    print(S.round(3).to_string(index=False)); print(G.round(3).to_string(index=False)); print("->", OD)


if __name__ == "__main__":
    main()
