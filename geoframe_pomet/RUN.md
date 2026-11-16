# Quick start

1. **Java 11.** For example Temurin 11. If it is not the default `java`, set `JAVA11` to its `bin` folder.
2. **Data.** Put the input files you received in `geoframe_pomet/data/` (see `data_example/README.md`). You need:
   * `stations_precip.csv`, `data_precip.csv`
   * `stations_temp.csv`, `data_temp.csv`
   * `cells/po_cells_1km.shp`
   * `dem_po_1km.asc`, `grid_cells.csv`
3. **Build and run.** From `geoframe_pomet/` (macOS or Linux):

   ```bash
   ./build.sh                              # -> ../build/libs/Krigings-0.9.9-po1.jar, uses the jars in lib/
   ./run_oms.sh simulation/precip.sim      # daily precipitation, two-part method
   ./run_oms.sh simulation/temp.sim        # daily mean temperature, anomaly method
   ```

   On Windows, use the OMS Console with this folder as the project, or WSL.
4. **Test first.** Set `runFrom`/`runTo` in the USER CONFIGURATION block of the sim to a few days, e.g.
   2020-10-02 to 2020-10-05. On a laptop with 12 threads, one year of precipitation takes about 10 minutes.
   1991–2020 takes several hours.
5. **Outputs.** One NetCDF file per year, plus daily summaries and the run configuration, are written to `output/`.
   If `nccopy` (netCDF C utilities) is not installed, the files are written as NetCDF-3 instead of compressed NetCDF-4.
   The memory limit is `XMX` (default 8g).

Other sims:
* `simulation/advanced/` holds single components, point targets and leave-one-out.
* Most of them read precomputed normals from `output/normals_java`. Create these with `normals_precip.sim` and
  `normals_temp.sim`.
* `simulation/tests/` are development tests. Some of them need reference outputs that are not in the repository.
