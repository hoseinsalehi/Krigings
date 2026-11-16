# PoKriging — input files and .sim schema (GEOframe / OMS3)

Two OMS3 components, each runnable as a one-shot model with `oms3.CLI`:
* `PoNormals`  — station monthly normals (gap-filled), wet-day frequency normals, gridded normals + station LOO normals;
* `PoKriging`  — daily interpolation (grid or points) and/or station leave-one-out.
One normals .sim and one interpolation/LOO .sim per variable (examples in `sim_examples/`).
Package name proposed: `org.geoframe.interpolation.pokriging` (new repository; Krigings untouched).
Method details: `METHOD_SPEC.md`.

---------------------------------------------------------------------------------------------------------------

## 1. Input files

### 1.1 stations.csv (one file per variable)
```
id,x,y,elevation[,name,network,elev_from_dem,tpi20,facet,...]
2466000,4348255.127,2581754.33,950
```
| column | required | type | meaning |
|---|---|---|---|
| id | yes | string (unique, no commas) | station id = column header in data.csv |
| x, y | yes | double | coordinates in the CRS declared by `crs` (metric CRS required; Po: EPSG:3035) |
| elevation | yes | double, m | station elevation (missing not allowed; the converter fills it from the DEM and sets elev_from_dem = 1) |
| any other | no | – | carried through to outputs, ignored by the method (`tpi20`, `facet` are only needed if the inactive PRISM facet/TPI weights are ever switched on) |
* Header required; comma separator; `.` decimal; UTF-8; lines starting with `#` ignored.
* Station order in the engine = **column order of data.csv** (needed for equivalence: LOO chunk order, ties).
  Stations in stations.csv without a data column are ignored; data columns without a station row → error.

### 1.2 data.csv (wide, one file per variable)
```
date,2466000,2546000,2576000,...
1991-01-01,0,,12.4,...
```
* first column `date`, format `YYYY-MM-DD`, strictly daily and increasing without gaps (error otherwise);
* one column per station id; value in mm/day (precip) or °C (temp, daily mean);
* missing: the string given by `missingValue` (default empty field; `-9999` accepted as well when
  `missingValue = -9999`; both "" and the code are always treated as missing);
* values are parsed as **float32** (`Float.parseFloat`) and then widened to double (bitwise equal to the Python engine);
* validity after reading: precip < 0 → missing; temp outside [`tempMin`, `tempMax`] = [−50, 50] → missing;
* QC and day-labelling harmonisation must already be applied (Python s7/s0). The file is the clean input.

### 1.3 Targets
* **grid**: DEM as GeoTIFF or ESRI ASCII grid (`demFile`), cell centres = target locations, NODATA cells skipped;
  optional `maskFile` (same grid, NODATA/0 = outside). Target order = row-major from the first (northern) row,
  exactly as `grid_cells.csv` (needed for `rngMode = numpy`, METHOD_SPEC §A-2). The DEM must be the grid of the
  normals (PoNormals writes normals on the DEM grid).
* **points**: `points.csv` = `id,x,y,elevation[,...]`; same CRS. Normals at a point = normal of the nearest grid cell
  (as Python); elevation from the file (if missing → nearest DEM cell).
* elevation holes inside the mask: filled with the nearest valid DEM value (as Python) and flagged (qc value 64).

### 1.4 Optional tables
| file | used by | columns |
|---|---|---|
| `uncertaintyTable` | PoKriging | variable,season,cls,scale_k (s9 output `uncertainty_calibration_<var>.csv`) |
| `occurrenceThresholds` | PoKriging, twopart, `calibration = table` | cell,month,threshold[,wet_days_target,p_normal] |
| normals (from PoNormals) | PoKriging | see §3 |

---------------------------------------------------------------------------------------------------------------

## 2. .sim keys

Syntax: OMS3 Groovy DSL (`oms3.SimBuilder`); keys are `@In` fields of the component, written
`"<component>.<key>" <value>` inside `parameter { }`. Types: S = String, D = double, I = int, B = boolean.
Defaults = locked Po values; a .sim needs only the keys without a default.

### 2.1 Common keys (both components)
| key | type | default | values / meaning |
|---|---|---|---|
| variable | S | – | `precip` \| `temp` |
| stationsFile | S | – | path of stations.csv |
| dataFile | S | – | path of data.csv |
| missingValue | S | "" | missing code in data.csv (e.g. `-9999`) |
| crs | S | – | e.g. `EPSG:3035` (metric) |
| tempMin, tempMax | D | −50, 50 | temperature validity range |
| threads | I | number of cores − 2 | worker threads (one day per thread, as Python) |
| outputDir | S | – | output folder (created) |
| logLevel | S | INFO | |

### 2.2 PoNormals
| key | type | default | meaning (METHOD_SPEC section) |
|---|---|---|---|
| normalsStart, normalsEnd | S | 1991-01-01, 2020-12-31 | normal period (§1) |
| monthCompleteness | D | 0.80 | valid-day share for a station-month |
| minObsYears | I | 5 | observed years per calendar month |
| minTotalYears | I | 20 | observed + reconstructed years |
| fillRadiusKm | D | 100 | gap-fill donor radius |
| fillMaxCandidates | I | 30 | donors examined |
| fillMinOverlapMonths | I | 36 | common months |
| fillMinR | D | 0.70 (precip) / 0.85 (temp) | min correlation |
| fillTopK | I | 3 | donors used (weights r²) |
| demFile | S | – | target grid (§1.3) |
| regK | I | 15 (precip) / 25 (temp) | stations in the local regression (§2.2) |
| regDz0 | D | 500 (precip) / 300 (temp) | elevation-difference weight scale, m |
| regExtrapDz | D | 0 (precip) / 1e9 (temp) | cap of the target elevation above the highest of the K stations, m (0 = no upward extrapolation) |
| regExtrapLapse | D | 0 (precip) / −0.0075 (temp) | °C/m used above the highest of K stations (temp only; 0 = off) |
| slopeMin, slopeMax | D | −0.3e-3, 0.8e-3 (precip, ln per m) / −10e-3, 3e-3 (temp, °C/m) | slope bounds |
| slopeMode | S | annual (precip) / monthly (temp) | 12-month mean of bounded monthly slopes, or monthly |
| precipFloorMm | D | 0.1 | ln(max(M, floor)) |
| wetFrequency | B | true (precip) | also compute p normals (two-part) |
| occThreshold | D | 1.0 | wet day ≥ this (mm) for p |
| wfK, wfDz0, wfExtrapDz | I, D, D | 15, 1e9, 0 | regression parameters for logit(p) |
| wfSlopeMin, wfSlopeMax | D | −1.0e-3, 2.0e-3 | logit per m |
| wfPmin, wfPmax | D | 0.005, 0.995 | clip before logit |
| facetWeight, tpi0, useTpi | D, D, B | 1.0, 1e9, false | PRISM facet/TPI terms (inactive in the Po product) |
| writeLoo | B | true | station LOO normals (needed for PoKriging LOO) |
**Outputs** (`outputDir`): `normals_<var>.nc` (normal, elevation_slope; month × y × x), `wetfreq_precip.nc`
(wetfreq, elevation_slope), `station_normals_<var>.csv`, `station_normals_loo_<var>.csv`
(station,month,normal_station,normal_grid_loo), `station_wetfreq_loo_precip.csv`
(station,month,wetfreq_station,wetfreq_grid_loo,n_months), `station_monthly_<var>.csv.gz` (observed/filled),
`params_<var>.json`. Same names and columns as Python, so either engine can feed the other.

### 2.3 PoKriging
| key | type | default | meaning |
|---|---|---|---|
| mode | S | interpolation | `interpolation` \| `loo` \| `both` |
| target | S | grid | `grid` \| `points` (ignored for `mode = loo`) |
| demFile | S | – | grid targets (§1.3) |
| pointsFile | S | – | point targets |
| periodStart, periodEnd | S | 1991-01-01, 2020-12-31 | days to compute |
| datesFile | S | – | optional list of dates (column `date`), intersected with the period |
| normalsDir | S | – | PoNormals output (required unless `precipMethod`/temp method is `direct`) |
| precipMethod | S | twopart | `direct` \| `climatology` (ratio, v2.0) \| `twopart` (v5) |
| tempMethod | S | anomaly | `direct` \| `anomaly` |
| temperatureTrend | S | local | `local` (per-target daily lapse on the 50 neighbours) \| `none` \| `global` (§A-16: not in Python — to be decided) |
| lapseMinSpanM, lapseMinStations, lapseFallbackN | D, I, I | 300, 10, 100 | lapse fallback rule (§3.4) |
| anomLapseMaxExtrapM | D | 300 | elevation clip of the anomaly lapse (anomaly mode only) |
| nVario | I | 50 | stations defining the local variogram (and lapse) |
| nKrig | I | 7 | kriging neighbours |
| nLags, minPairsPerLag | I | 12, 5 | empirical variogram |
| variogramModel | S | spherical | only spherical in the reference |
| psoParticles, psoIterations | I | 24, 35 | |
| psoInertia, psoC1, psoC2 | D | 0.72, 1.45, 1.45 | |
| fitWeighting | S | pair-count | `pair-count` \| `equal` |
| rngMode | S | numpy | `numpy` = exact port of numpy PCG64/SeedSequence/ziggurat + BLAKE2b seeds + chunking (equivalence with v2/v5) \| `java` |
| chunk | I | 2048 | targets per random-stream chunk (numpy mode) |
| globalVariogramFallback | B | true | |
| zeroTol | D | 0.0 | precip zero shortcut: all nKrig neighbours ≤ zeroTol → 0 (direct, climatology) |
| minDailyNormalMm | D | 0.05 | floor of the daily climatology (divisor) |
| maxPrecipRatio | D | 50 | cap of station ratios |
| occThreshold | D | 1.0 | two-part wet day (mm) |
| pFloor | D | 0.01 | floor of p in i = M / (p·days) |
| calibration | S | auto | two-part occurrence: `auto` (monthly if every month has ≥ calibMinDays run days, else pooled) \| `monthly` \| `pooled` \| `table` |
| calibMinDays | I | 100 | |
| occurrenceThresholds | S | – | for `calibration = table` |
| writeTwopartPass1 | B | false | keep π / amount / s_amt per day (needed to calibrate later or in chunks of years) |
| uncertaintyTable | S | – | if set, writes `<variable>_uncertainty` = kriging_std × k(season, class) |
| diagnostics | B | true | write the diagnostics file (grid); used only with `outputVariables = default` |
| outputVariables | S | default | NetCDF grid variables: `default` (value, kriging_std [, uncertainty] + all diagnostics, as before) \| `all` \| comma list of names (README, table *Output variables*); a name not available for the method stops the run and lists the allowed names |
| minActiveStations | I | nVario + 1 = 51 | days with fewer active stations are skipped (fill value) |
| wetThresholdStats | D | 0.1 | wet threshold of the LOO summary metrics |
| outputFormat | S | netcdf | `netcdf` (grid, one file per year) \| `csv` (points, LOO; grid as long table only for tests) |
| quantise | B | true | NetCDF int16 encoding as Python (0.1 mm, 0.01 °C); false = float32 |
| equivalenceDump | S | – | optional: write unquantised value/kstd/diagnostics for a cell list (cell ids of grid_cells.csv) as CSV, for the automatic Python comparison |

**Outputs**
* grid, NetCDF: `<var>_<YYYY>.nc` (precipitation [mm day-1] | temperature [degC], kriging_std
  [, <var>_uncertainty], plus the method intermediates selected by `outputVariables`); `diagnostics_<var>_<YYYY>.nc`
  (METHOD_SPEC §3.10, selected diagnostics); CF time "days since 1991-01-01", grid_mapping crs EPSG:3035, attributes
  method/parameters (+ `output_variables` when not default). A file is not written if none of its variables is selected.
* points, CSV: `points_<var>_<YYYY>.csv` = id,date,value,kriging_std,D1,D7,qc_flags[,pi,amount].
* LOO, CSV: `loo_<var>_<YYYY>.csv` = date,station,obs,pred,kriging_std,D1,D7,qc_flags[,pi,amount,samount,p_target,wet_pred];
  `daily_summary_<var>_<YYYY>.csv` (n_active, loo_n, loo_bias, loo_mae, loo_rmse, loo_r, loo_wet_n,
  loo_wet_mae, loo_wet_pbias, seconds, …).
* `run_config.json`: all keys, version, timings.

---------------------------------------------------------------------------------------------------------------

## 3. Converter from the current Po data

`../tools/convert_po_inputs.py` (run 2026-10-19, output `../input_po/`, 87 MB):
* source = exactly the final runs' input: `~/Downloads/PoGrid_v5/qc/Final_<var>_plusCH_QC.csv` (= PoGrid_v4/qc,
  cleaned and day-shift harmonised) + `~/Downloads/PoGrid_v4_input/STATIONS_plusCH_<var>.shp`, read through
  `po_common.load_station_data` / `attach_station_predictors`;
* `stations_precip.csv` 1636 stations (1 elevation from DEM), `stations_temp.csv` 1570;
* `data_precip.csv` 1991-01-01…2020-12-31, 9 841 333 valid station-days (34.9 MB);
  `data_temp.csv` 8 926 076 valid station-days (50.5 MB); missing = empty field;
* `dem_po_1km.tif` / `.asc` (+ `.prj`): 494 × 330, 1 km, EPSG:3035, NODATA −9999 outside the 88 074-cell mask;
  `grid_cells.csv` (id = engine cell index, x, y, elevation, iy, ix, elev_filled); `points_reference.csv`
  (the 8192 reference cells as point targets);
* verified: values re-read as float32 are bit-identical to the engine arrays; station x/y/elevation, DEM cells and
  grid-cell coordinates identical to the engine.
The full 1991-2020 table is exported because the normals need it (station data only, < 100 MB).

---------------------------------------------------------------------------------------------------------------

## 4. Questions for the user before the Java core
1. `temperatureTrend = global`: implement as one daily OLS on all active stations (new, no Python reference), or
   expose only `local | none` (METHOD_SPEC §A-16)?
2. RNG: accept porting numpy's PCG64/SeedSequence/ziggurat for near-bitwise equivalence (recommended, §A-1), or a
   statistical equivalence test only?
3. Two-part production: run the 30-year calibration inside Java (two passes, π stored ≈ 2.2 GB for 30 years, as in
   Python), or keep `calibration = table` with thresholds from a first pass?
4. Grid output: one NetCDF per year (as Python, proposed) or a GEOframe-specific format; points: is point kriging
   enough for HRUs or is block kriging over the HRU area needed?

---------------------------------------------------------------------------------------------------------------

## 5. Status 2026-10-19 (implemented)
Decisions: `temperatureTrend = global` implemented as one daily OLS on all active stations (no Python reference); RNG =
exact numpy port (`rngMode = numpy`, default); two-part calibration two-pass inside Java plus `calibration = table`;
NetCDF one file per year (NetCDF-3 via netcdf-java, converted to NetCDF-4/zlib by nccopy); point kriging for HRUs (no block
kriging); precipMethod `hybrid` dropped. Implementation: `../Krigings` (component `pokriging.PoKriging`, `pokriging.PoNormals`),
sims in `../simulation`, results in `../EQUIVALENCE_REPORT.md`. Package names follow the GEOframe Krigings repository
(`pokriging`, `krigingsCore`, `utilities.random`) instead of the proposed `org.geoframe.interpolation.pokriging`.
