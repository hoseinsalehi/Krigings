#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert the current Po inputs into the GEOframe PoKriging CSV format (see spec/CSV_AND_SIM_SCHEMA.md).

Source = exactly what the final Python runs read (PO_STATION_SET=v4, DATA_ROOT=~/Downloads/PoGrid_v5):
the cleaned, day-shift-harmonised CSVs <DATA_ROOT>/qc/Final_<var>_plusCH_QC.csv and the station shapefiles
~/Downloads/PoGrid_v4_input/STATIONS_plusCH_<var>.shp (EPSG:3035). Station metadata and values are taken through
po_common.load_station_data + attach_station_predictors, i.e. the engine's own view of the data:
  * only stations present in both the shapefile and the CSV, in CSV column order (= engine station index order);
  * MISSING_VALUE (-9999) -> missing; precipitation < 0 -> missing; temperature outside [-50, 50] degC -> missing;
  * station elevation from the shapefile, filled from the 1 km DEM where missing (column elev_from_dem = 1).
QC and day-shift harmonisation are NOT redone here (they stay in Python and are already in the *_QC.csv).

Values are written as the shortest decimal string that round-trips the engine's float32 value
(Java: Float.parseFloat(s) gives the identical float32; widen to double afterwards).

Outputs (in --out):
  stations_<var>.csv   id,x,y,elevation,elev_from_dem,tpi20,facet
  data_<var>.csv       date (YYYY-MM-DD), one column per station id, missing = --missing (default empty)
  dem_po_1km.asc / dem_po_1km.tif   product grid (494 x 330, 1 km, EPSG:3035), NODATA outside the product mask
  grid_cells.csv       id,x,y,elevation,iy,ix for every valid product cell in engine order (row-major from the north)
  points_reference.csv id,x,y,elevation of the reference target cells (if --reference-cells is given)

Usage:
  PO_STATION_SET=v4 PO_OUT_ROOT=~/Downloads/PoGrid_v5 python convert_po_inputs.py --out <dir> [--variables precip temp]
        [--start 1991-01-01 --end 2020-12-31] [--missing ""] [--reference-cells reference/targets_cells.csv]
"""
import argparse, os, sys, time
from pathlib import Path
os.environ.setdefault("PO_STATION_SET", "v4")
os.environ.setdefault("PO_OUT_ROOT", str(Path.home() / "Downloads/PoGrid_v5"))
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1"); sys.dont_write_bytecode = True
CODE = Path(os.environ.get("PO_CODE", str(Path.home() / "Downloads/PoGrid_v3_code")))
sys.path.insert(0, str(CODE))
import numpy as np, pandas as pd


def f32str(a, missing):
    """Shortest round-trip decimal strings of float32 values; NaN -> missing."""
    a = np.asarray(a, np.float32); out = np.full(a.shape, missing, dtype=object); ok = np.isfinite(a)
    out[ok] = [np.format_float_positional(x, unique=True, trim="-") for x in a[ok]]
    return out


def write_grid(grid, out):
    ny, nx = int(grid["ny"]), int(grid["nx"]); x, y = grid["x"], grid["y"]
    dx = float(x[1] - x[0]); dy = float(y[1] - y[0]); assert dx == 1000.0 and dy == -1000.0, (dx, dy)
    dem = np.full((ny, nx), -9999.0, np.float32); dem[grid["iy"], grid["ix"]] = grid["elev"]
    hdr = (f"ncols {nx}\nnrows {ny}\nxllcorner {x[0] - dx / 2:.1f}\nyllcorner {y[-1] + dy / 2:.1f}\n"
           f"cellsize {dx:.1f}\nNODATA_value -9999\n")
    with open(out / "dem_po_1km.asc", "w") as f:
        f.write(hdr)
        for r in dem: f.write(" ".join(np.format_float_positional(v, unique=True, trim="-") for v in r) + "\n")
    (out / "dem_po_1km.prj").write_text('PROJCS["ETRS89-extended / LAEA Europe",AUTHORITY["EPSG","3035"]]\n')
    try:
        import rasterio
        from rasterio.transform import from_origin
        with rasterio.open(out / "dem_po_1km.tif", "w", driver="GTiff", height=ny, width=nx, count=1, dtype="float32",
                           crs="EPSG:3035", transform=from_origin(x[0] - dx / 2, y[0] - dy / 2, dx, -dy), nodata=-9999.0,
                           compress="deflate") as ds:
            ds.write(dem, 1)
    except ImportError:
        print("rasterio not available: GeoTIFF skipped", file=sys.stderr)
    cells = pd.DataFrame(dict(id=np.arange(len(grid["xy"])), x=grid["xy"][:, 0], y=grid["xy"][:, 1], elevation=grid["elev"],
                              iy=grid["iy"], ix=grid["ix"], elev_filled=grid["elev_filled"].astype(int)))
    cells.to_csv(out / "grid_cells.csv", index=False)
    return cells


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True); ap.add_argument("--variables", nargs="+", default=["precip", "temp"])
    ap.add_argument("--start", default="1991-01-01"); ap.add_argument("--end", default="2020-12-31")
    ap.add_argument("--missing", default=""); ap.add_argument("--reference-cells", default=None)
    ap.add_argument("--no-grid", action="store_true")
    a = ap.parse_args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    import po_config as C, po_common as P
    grid = P.load_grid()
    if not a.no_grid:
        cells = write_grid(grid, out); print(f"grid: {len(cells)} valid cells -> dem_po_1km.asc/.tif, grid_cells.csv")
        if a.reference_cells:
            rc = pd.read_csv(a.reference_cells)
            cells.set_index("id").loc[rc["cell"].to_numpy(), ["x", "y", "elevation"]].reset_index() \
                .to_csv(out / "points_reference.csv", index=False)
    for v in a.variables:
        t0 = time.time()
        st = P.attach_station_predictors(P.load_station_data(v), grid)
        print(f"{v}: source {C.DATA[v]['obs']}")
        S = pd.DataFrame(dict(id=st["ids"], x=st["xy"][:, 0], y=st["xy"][:, 1], elevation=st["elev"],
                              elev_from_dem=st["elev_from_dem"].astype(int), tpi20=st["tpi20"], facet=st["facet"].astype(int)))
        assert S.id.is_unique
        S.to_csv(out / f"stations_{v}.csv", index=False)   # shortest round-trip repr
        d = st["dates"]; sel = (d >= pd.Timestamp(a.start)) & (d <= pd.Timestamp(a.end))
        vals = st["values"][sel]; dd = d[sel]
        with open(out / f"data_{v}.csv", "w") as f:
            f.write("date," + ",".join(st["ids"]) + "\n")
            for i in range(len(dd)):
                f.write(dd[i].strftime("%Y-%m-%d") + "," + ",".join(f32str(vals[i], a.missing)) + "\n")
        n_val = int(np.isfinite(vals).sum())
        print(f"{v}: {len(S)} stations, {len(dd)} days {dd[0]:%Y-%m-%d}..{dd[-1]:%Y-%m-%d}, {n_val} valid station-days, "
              f"{int(S.elev_from_dem.sum())} elevations from DEM, data file {os.path.getsize(out / f'data_{v}.csv') / 1e6:.1f} MB, "
              f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
