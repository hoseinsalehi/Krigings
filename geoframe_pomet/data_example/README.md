# Input data

The station observations used for GEOframe-PoMet belong to the data providers (regional environmental protection
agencies and meteorological services, EEAR-Clim, MeteoSwiss) and are not distributed in this repository.
Put the inputs in `geoframe_pomet/data/` (ignored by git). Their format is given in
`geoframe_pomet/spec/CSV_AND_SIM_SCHEMA.md`:

| file | content |
|---|---|
| `stations_<var>.csv` | `id,x,y,elevation` (ETRS89-LAEA, EPSG:3035, metres) |
| `data_<var>.csv` | `date,<id1>,<id2>,...` (YYYY-MM-DD; empty field = missing) |
| `cells/po_cells_1km.shp` | target cells with an elevation field |
| `dem_po_1km.asc` / `.tif` | 1-km DEM of the grid |

`geoframe_pomet/tools/convert_po_inputs.py` builds these files from the station tables of the Python pipeline.
Ask the authors for access to the prepared inputs.
