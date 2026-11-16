# -*- coding: utf-8 -*-
"""Shared I/O: product grid + DEM-derived predictors, stations, daily observations."""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
import sys
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).resolve().parent))
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy import ndimage

import po_config as C


# ----------------------------------------------------------------------------- grid
def load_grid(cache=None):
    """Valid-cell grid of the product with elevation and topographic predictors.

    Returns dict with ny, nx, x (nx), y (ny), iy, ix (valid cells), xy (N,2), elev (N),
    elev_filled (N bool), slope_deg, northness, eastness, tpi5, tpi20 (N), facet (N int8:
    0-7 = N, NE, ... NW facing on the ~5 km smoothed DEM, 8 = flat), lon, lat (N).
    """
    cache = Path(cache) if cache else C.OUT_ROOT / "grid_predictors.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=False)
        return {k: z[k] for k in z.files}
    import xarray as xr, geopandas as gpd
    with xr.open_dataset(C.GRID_MASK_NC) as ds:
        x = ds["x"].values.astype(float); y = ds["y"].values.astype(float)
        mask = np.isfinite(ds["precipitation"].isel(time=0).values)
    ny, nx = mask.shape
    pts = gpd.read_file(C.DEM_POINTS)
    ix = np.rint((pts["x"].values - x[0]) / (x[1] - x[0])).astype(int)
    iy = np.rint((pts["y"].values - y[0]) / (y[1] - y[0])).astype(int)
    ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    dem = np.full((ny, nx), np.nan); dem[iy[ok], ix[ok]] = pts["elevation"].values[ok].astype(float)
    # fill DEM holes inside the mask (and everywhere, for gradient/TPI computation) by nearest value
    have = np.isfinite(dem)
    _, (fy, fx) = ndimage.distance_transform_edt(~have, return_indices=True)
    dem_f = dem[fy, fx]
    elev_filled_2d = mask & ~have
    # slope/aspect on the 1 km DEM, facets on a ~5 km smoothed DEM (PRISM-style facet scale)
    gy, gx = np.gradient(dem_f, y[1] - y[0], x[1] - x[0])          # dz/dy (y axis direction), dz/dx
    north_comp = gy if (y[1] - y[0]) > 0 else -gy                    # dz/dNorth
    slope = np.degrees(np.arctan(np.hypot(gx, north_comp)))
    aspect = (np.degrees(np.arctan2(-gx, -north_comp)) + 360.0) % 360.0   # direction the slope faces
    sm = ndimage.gaussian_filter(dem_f, sigma=2.5)
    sgy, sgx = np.gradient(sm, y[1] - y[0], x[1] - x[0]); snc = sgy if (y[1] - y[0]) > 0 else -sgy
    sslope = np.hypot(sgx, snc); sasp = (np.degrees(np.arctan2(-sgx, -snc)) + 360.0) % 360.0
    facet = np.where(sslope < 0.02, 8, (((sasp + 22.5) % 360) // 45).astype(int)).astype(np.int8)
    def tpi(r_km):   # topographic position index: cell elevation minus mean of a (2r+1) km window
        return dem_f - ndimage.uniform_filter(dem_f, size=2 * int(r_km) + 1, mode="nearest")
    tpi5 = tpi(5); tpi20 = tpi(20)
    vy, vx = np.nonzero(mask)
    xx, yy = x[vx], y[vy]
    from pyproj import Transformer
    lon, lat = Transformer.from_crs(C.METRIC_CRS, "EPSG:4326", always_xy=True).transform(xx, yy)
    out = dict(ny=np.int64(ny), nx=np.int64(nx), x=x, y=y, iy=vy.astype(np.int32), ix=vx.astype(np.int32),
               xy=np.c_[xx, yy], elev=dem_f[vy, vx], elev_filled=elev_filled_2d[vy, vx],
               slope_deg=slope[vy, vx].astype(np.float32), northness=np.cos(np.radians(aspect[vy, vx])).astype(np.float32),
               eastness=np.sin(np.radians(aspect[vy, vx])).astype(np.float32), tpi5=tpi5[vy, vx].astype(np.float32),
               tpi20=tpi20[vy, vx].astype(np.float32), facet=facet[vy, vx], lon=lon, lat=lat,
               dem2d=dem_f.astype(np.float32), tpi20_2d=tpi20.astype(np.float32), facet2d=facet)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, **out)
    return out


def sample_predictors_at(grid, xy):
    """Topographic predictors at arbitrary points (nearest 1 km cell of the full 2-D arrays)."""
    x, y = grid["x"], grid["y"]
    ix = np.clip(np.rint((xy[:, 0] - x[0]) / (x[1] - x[0])).astype(int), 0, int(grid["nx"]) - 1)
    iy = np.clip(np.rint((xy[:, 1] - y[0]) / (y[1] - y[0])).astype(int), 0, int(grid["ny"]) - 1)
    return dict(dem=grid["dem2d"][iy, ix].astype(float), tpi20=grid["tpi20_2d"][iy, ix].astype(float),
                facet=grid["facet2d"][iy, ix])


# ----------------------------------------------------------------------------- stations / observations
def parse_dates(raw, fmt):
    s = raw.astype(str).str.strip()
    nonempty = raw.notna() & ~s.isin(["", "nan", "NaT", "None"])
    dt = pd.to_datetime(s, format=fmt, errors="coerce")
    if not (dt.notna() | ~nonempty).all():
        raise ValueError(f"Date column does not parse with format {fmt!r}; first values {s.head(3).tolist()}")
    return dt


def apply_day_shift(values, dates, ids, table, fill=np.nan):
    """Move every value of station s observed in year Y by shift_days (table: station, year, shift_days).
    A value labelled D with shift -1 becomes D-1. Vacated days get `fill`; two values landing on the same
    day (at a change of convention) are both discarded. Works on float or object arrays (T, S)."""
    if table is None or len(table) == 0: return values
    out = values.copy(); yrs = np.asarray(dates.year); T = len(dates)
    tb = table.assign(station=table["station"].astype(str))
    for s, g in tb.groupby("station"):
        j = np.flatnonzero(np.asarray(ids).astype(str) == s)
        if len(j) == 0: continue
        j = int(j[0]); k = pd.Series(g["shift_days"].to_numpy(int), index=g["year"].to_numpy(int)).reindex(yrs).fillna(0).to_numpy(int)
        if not k.any(): continue
        src = np.arange(T); dst = src + k; ok = (dst >= 0) & (dst < T)
        col = np.full(T, fill, dtype=values.dtype); cnt = np.bincount(dst[ok], minlength=T)
        col[dst[ok]] = values[src[ok], j]; col[cnt > 1] = fill; out[:, j] = col
    return out


def load_day_shift(variable):
    f = C.DATA_ROOT / C.SHIFT_FILE.format(variable=variable)
    if not f.exists() or os.environ.get("PO_SHIFT", "1") == "0": return None
    return pd.read_csv(f, dtype={"station": str})


def load_station_data(variable, exclude=None):
    """Daily observations aligned to station metadata.

    Returns dict: dates (DatetimeIndex, strictly daily), ids (S), values (T,S float32, NaN missing),
    xy (S,2), elev (S; shapefile, DEM where missing), tpi20, facet (S).
    """
    import geopandas as gpd
    cfg = C.DATA[variable]
    g = gpd.read_file(cfg["stations"]).to_crs(C.METRIC_CRS)
    sid = g[C.STATION_ID_FIELD].astype(str).to_numpy()
    sxy = np.c_[g.geometry.x.to_numpy(float), g.geometry.y.to_numpy(float)]
    selev = pd.to_numeric(g[C.STATION_ELEV_FIELD], errors="coerce").to_numpy(float)
    df = pd.read_csv(cfg["obs"], low_memory=False)
    df["Date"] = parse_dates(df["Date"], cfg["date_format"])
    df = df.loc[df["Date"].notna()].reset_index(drop=True)
    if not (df["Date"].is_monotonic_increasing and df["Date"].is_unique and (df["Date"].diff().dropna().dt.days == 1).all()):
        raise ValueError("Observation dates are not a strictly daily increasing sequence")
    lut = {s: i for i, s in enumerate(sid)}
    excl = set(map(str, exclude or []))
    cols = [c for c in df.columns if c != "Date" and str(c) in lut and str(c) not in excl]
    v = df[cols].apply(pd.to_numeric, errors="coerce").to_numpy(np.float32, copy=True)
    v[np.isclose(v, C.MISSING_VALUE)] = np.nan
    if variable == "precip": v[v < 0] = np.nan
    else: v[(v < C.TEMP_RANGE[0]) | (v > C.TEMP_RANGE[1])] = np.nan
    idx = np.array([lut[str(c)] for c in cols])
    if cfg["obs"] == cfg["raw"]:                    # raw file: harmonise the day labelling (the cleaned CSV already is)
        sh = load_day_shift(variable)
        if sh is not None:
            v = apply_day_shift(v, pd.DatetimeIndex(df["Date"]), [str(c) for c in cols], sh)
            print(f"[shift] day labelling harmonised for {int((sh.shift_days != 0).groupby(sh.station).any().sum())} {variable} stations", file=sys.stderr)
    qf = C.DATA_ROOT / "qc" / C.QC_FLAGS_FILE.format(variable=variable)
    if qf.exists() and os.environ.get("PO_QC", "1") != "0" and cfg["obs"] == cfg["raw"]:   # cleaned CSV already excludes them
        fl = pd.read_parquet(qf, columns=["date", "station"])
        ti = pd.DatetimeIndex(df["Date"]).get_indexer(pd.to_datetime(fl["date"]))
        cpos = {str(c): j for j, c in enumerate(cols)}; si = fl["station"].astype(str).map(cpos).to_numpy(float)
        ok = (ti >= 0) & np.isfinite(si); v[ti[ok], si[ok].astype(int)] = np.nan
        print(f"[qc] {int(ok.sum())} {variable} station-days removed ({qf.name})", file=sys.stderr)
    return dict(dates=pd.DatetimeIndex(df["Date"]), ids=np.asarray([str(c) for c in cols]), values=v,
                xy=sxy[idx], elev=selev[idx])


def attach_station_predictors(st, grid):
    p = sample_predictors_at(grid, st["xy"])
    elev = st["elev"].copy(); bad = ~np.isfinite(elev); elev[bad] = p["dem"][bad]
    st.update(elev=elev, elev_from_dem=bad, tpi20=p["tpi20"], facet=p["facet"])
    return st


def station_distance_km(xy):
    d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1)) / 1000.0
    return d.astype(np.float32)
