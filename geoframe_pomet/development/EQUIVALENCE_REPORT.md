# Java (GEOframe Krigings + pokriging) vs Python Po engine — equivalence and speed (2026-10-19)

Reference: `reference/` (unmodified `PoGrid_v3_code` engine, numba kernels; 10 dates of 2020; 8 192 grid cells = production
chunks 2, 7, 9, 20; station LOO on the same dates; normals rebuilt bit-identically to `PoGrid_v5/normals`).
Java: `Krigings/` build. The runs went **through `oms3.CLI -r`** (`simulation/tests/eq_*.sim`, Java 11 temurin, 10 threads,
`rngMode = numpy`, CSV output without quantisation). Comparison: `tools/compare_equivalence.py`, `tools/compare_normals.py`.
Raw results: `output/eq_oms/cmp_*.json`, `output/eq_oms/equivalence_summary.csv`, `output/normals_java/cmp_normals.json`.
Note: the Python reference grid CSV stores 10 significant digits and the LOO CSV 8 (float32), so "identical" means equal
to that precision.

## 1. Daily kriging with the Python normals (the port of the daily engine on its own)

"≤1e-6" = share of cell-days (or station-days) whose relative difference is ≤ 1e-6. Max columns are absolute differences
(mm or °C; range in km).

| method | cell-days | value ≤1e-6 [%] | value max | kriging-std max | range ≤1e-6 [%] | range max | LOO rows | LOO pred ≤1e-6 [%] | LOO pred max | LOO KGE py / java | LOO MAE py / java |
|---|---|---|---|---|---|---|---|---|---|---|---|
| precip_direct | 81920 | 100.0000 | 5e-08 | 5e-08 | 100.0000 | 5e-08 | 12091 | 100.0000 | 1e-05 | 0.923807 / 0.923807 | 3.825018 / 3.825018 |
| precip_ratio | 81920 | 99.9988 | 0.036 | 0.028 | 99.9981 | 0.06 | 12091 | 99.9917 | 0.73 | 0.866803 / 0.866807 | 3.912512 / 3.912452 |
| precip_twopart | 81920 | 100.0000 | 1.1e-13 | 0.052 | 99.9976 | 0.51 | 12091 | 100.0000 | 2e-05 | 0.516131 / 0.516131 | 7.469147 / 7.469147 |
| temp_direct | 81920 | 100.0000 | 5e-09 | 5e-10 | 100.0000 | 5e-08 | 11830 | 100.0000 | 1e-06 | 0.990289 / 0.990289 | 0.823834 / 0.823834 |
| temp_anomaly | 81920 | 100.0000 | 5e-09 | 5e-10 | 99.9988 | 0.0065 | 11830 | 100.0000 | 1e-06 | 0.989388 / 0.989388 | 0.817166 / 0.817166 |
| temp_anomaly_nolapse | 81920 | 100.0000 | 2.5e-06 | 6.1e-06 | 99.9988 | 0.027 | 11830 | 100.0000 | 1e-06 | 0.988340 / 0.988340 | 0.879517 / 0.879517 |

Two-part precipitation:
* **pooled 10-day calibration**: the wet/dry flag is identical in 81 920 / 81 920 cell-days and in all 12 091 LOO
  station-days. π agrees within 1e-6 in 99.998 % of cell-days.
* **production threshold table** (`calibration = table`, `reference/twopart_production_calibration/thresholds_cells.csv`):
  * versus the production final field `grid_final_production.csv.gz`, the wet flag is equal in
    100.0 % of the 81920 cell-days;
  * the final value differs by at most 1.1e-13 mm.

LOO wet-day statistics (wet = ≥ 1 mm, 10 dates):

| method | predicted wet frequency [%] py / java | predicted wet-day intensity [mm] py / java |
|---|---|---|
| precip_direct | 45.6455 / 45.6455 | 32.3681 / 32.3681 |
| precip_ratio | 45.4967 / 45.4967 | 31.4721 / 31.4722 |
| precip_twopart | 24.7870 / 24.7870 | 38.7378 / 38.7378 |

Residual differences:
* Direct precipitation, two-part and direct temperature agree to the precision of the reference files in every cell.
* The isolated larger differences (precip ratio: 1 cell-day 0.036 mm and 1 LOO station-day 0.73 mm; two-part: kriging-std
  0.05 mm in 2 cell-days) come from targets where the PSO trajectory diverges after a sub-ulp difference of its inputs. The
  OK solver uses LU in Java and LAPACK in numpy, and the detrending fit for the global fallback variogram is centred OLS in
  Java versus `np.polyfit` (SVD) in Python.
* This is the expected behaviour of a near-tie in the weakly identified variogram fit (METHOD_SPEC §A-1). It is far below
  the seed-sensitivity envelope: an independent seed gives p99 = 5.7 mm and 0.50 °C.
* The acceptance target of reference/README test 2 (≥ 99.9 % of cell-days within 1e-6) is met by every method.
* The new `temperatureTrend = global` option has no Python counterpart (METHOD_SPEC §A-16). It was run
  (`eq_temp_anomaly_globallapse.sim`) but not compared.

## 2. Normals (PoNormals vs Python s1/s2/s1w)

| file | column | n | NaN mismatch | max rel. diff | ≤1e-9 [%] | bitwise [%] |
|---|---|---|---|---|---|---|
| station_normals_precip.csv | normal | 17598 | 0 | 6.3e-16 | 100.00 | 71.4 |
| station_normals_precip.csv | normal_obs_only | 19306 | 0 | 4.7e-16 | 100.00 | 95.9 |
| station_normals_precip.csv | n_obs_years | 19632 | 0 | 0.0e+00 | 100.00 | 100.0 |
| station_normals_precip.csv | n_total_years | 19632 | 0 | 0.0e+00 | 100.00 | 100.0 |
| station_normals_temp.csv | normal | 15280 | 0 | 5.9e-14 | 100.00 | 58.3 |
| station_normals_temp.csv | normal_obs_only | 18736 | 0 | 9.7e-14 | 100.00 | 75.4 |
| station_normals_temp.csv | n_obs_years | 18840 | 0 | 0.0e+00 | 100.00 | 100.0 |
| station_normals_temp.csv | n_total_years | 18840 | 0 | 0.0e+00 | 100.00 | 100.0 |
| station_normals_loo_precip.csv | normal_station | 17598 | 0 | 6.3e-16 | 100.00 | 71.4 |
| station_normals_loo_precip.csv | normal_grid_loo | 19632 | 0 | 9.6e-15 | 100.00 | 41.5 |
| station_normals_loo_temp.csv | normal_station | 15280 | 0 | 5.9e-14 | 100.00 | 58.3 |
| station_normals_loo_temp.csv | normal_grid_loo | 18840 | 0 | 2.6e-12 | 100.00 | 33.2 |
| station_wetfreq_loo_precip.csv | wetfreq_station | 18205 | 0 | 0.0e+00 | 100.00 | 100.0 |
| station_wetfreq_loo_precip.csv | wetfreq_grid_loo | 19632 | 0 | 5.0e-15 | 100.00 | 69.7 |
| station_wetfreq_loo_precip.csv | n_months | 19632 | 0 | 0.0e+00 | 100.00 | 100.0 |
| normals_precip.nc | normal | 1056888 | 0 | 0.0e+00 | 100.00 | 100.0 |
| normals_precip.nc | elevation_slope | 1056888 | 0 | 0.0e+00 | 100.00 | 100.0 |
| normals_temp.nc | normal | 1056888 | 0 | 0.0e+00 | 100.00 | 100.0 |
| normals_temp.nc | elevation_slope | 1056888 | 0 | 0.0e+00 | 100.00 | 100.0 |
| wetfreq_precip.nc | wetfreq | 1056888 | 0 | 0.0e+00 | 100.00 | 100.0 |
| wetfreq_precip.nc | elevation_slope | 1056888 | 0 | 0.0e+00 | 100.00 | 100.0 |

* All gridded normals (`normals_*.nc`, `wetfreq_precip.nc`, float32) are **bitwise identical**.
* Station normals agree within 6e-14 relative (temperature LOO grid normals within 2.6e-12). The remaining differences come
  from the summation order of the gap-filling correlations (numpy pairwise over 360 months is reproduced, but not every
  numpy reduction).
* Accepted station-months: precip 17 598 of 19 632, temp 15 280 of 18 840, identical to Python.
* One subtlety reproduced: numpy `add.reduceat` on float32 is (first element) + pairwise(rest).

## 3. End to end (Java normals → Java daily kriging)

| method | value ≤1e-6 [%] | value max | range ≤1e-6 [%] | LOO KGE py / java | LOO MAE py / java |
|---|---|---|---|---|---|
| precip_twopart | 100.000 | 1.1e-13 | 99.996 | 0.516131 / 0.516142 | 7.469147 / 7.469256 |
| precip_ratio | 99.977 | 2 | 99.966 | 0.866803 / 0.866797 | 3.912512 / 3.912425 |
| temp_anomaly | 99.999 | 0.0016 | 99.999 | 0.989388 / 0.989388 | 0.817166 / 0.817167 |

The 1e-14 differences of the station normals occasionally change a PSO trajectory. This happened in 0.02–0.04 % of
cell-days for the precipitation ratio method (max 2 mm). LOO skill is unchanged at the 5th decimal.

## 4. Speed (same machine: macOS, 16 logical CPUs; 10 workers / threads)

| case | Python (numba) | Java 11 temurin x86_64 (Rosetta) | Java 25 zulu arm64 |
|---|---|---|---|
| 10 reference dates, temp anomaly (8 192 cells + ~1 180 LOO), wall of the daily loop | 0.69 s | 3.4–3.9 s | 1.73–1.83 s |
| 10 reference dates, precip ratio | 0.67–1.73 s | 2.4–2.6 s | 1.18–1.20 s |
| 10 reference dates, precip two-part | 1.18–1.22 s | 3.7–4.4 s | 2.29 s |
| steady state, temp anomaly, 90 days, 10 threads | ≈ 0.069 s/day (0.69 s per 10 days with 10 workers) | 0.112 s/day | **0.065 s/day** |
| steady state, 1 thread | 0.65 s/day/worker | 0.77 s/day | **0.46 s/day** |
| set-up (data, normals, distance matrix) | 0.9–1.2 s | 1.1–1.6 s | 0.6–0.8 s |

* The 10-day wall times are dominated by JIT warm-up and set-up. Over longer runs, the Java core on a native arm64 JVM is
  as fast as the numba engine at 10 threads and 1.4× faster single-threaded.
* The required Java 11 of the OMS console here is the x86_64 build running under Rosetta, which is 1.6–1.7× slower.
  Installing an arm64 (aarch64) Java 11, for example Temurin 11 aarch64, removes this penalty; the code needs no change.
* Optimisations done:
  * a 12-bin objective unrolled in numpy pairwise order;
  * array fills of the random stream;
  * per-block parallelism over (day × 2 048-target chunk);
  * an asynchronous ordered writer.
* The results are identical for any number of threads, because every chunk owns its random stream.
* Profile (JFR, temp anomaly): PSO 54 %, empirical variograms 13 %, CSV/data reading 10 %, k-d tree 5 %.


## Addendum 2026-10-19: output selection (`outputVariables`)
New key `outputVariables` (default | all | list) exposes the intermediates of each method (README, *Output variables*).
`krigingsCore.local.TargetResult` gained two arrays `resid`, `trend` (est = resid + trend), filled in `DayKriging`;
no computed value changed. Re-certification: all `simulation/tests/eq_*.sim` were run with the build before and after
the change; the 28 output files (grid CSV, LOO, daily summaries, NetCDF) are identical value by value, and
`tools/compare_equivalence.py` gives JSON results identical to the stored ones for the six methods (production
threshold-table case: wet flag equal in 100 % of cell-days, max |value diff| 1.1e-13 mm, as before). Whole district,
2020-10-02/03: `default` outputs of the new jar identical to the old jar for precip direct / climatology / two-part
and temperature anomaly.

## Uncertainty calibration (s9 → pokriging.PoUncertainty), 2026-10-19
Data: the production LOO of v5 two-part precipitation (`PoGrid_v5/precip_anomaly/loo_precip_1991..2020.parquet`,
9 841 333 station-days, 1 636 stations) and of the v4 anomaly temperature (`PoGrid_v4/temp_anomaly/loo_temp_*.parquet`,
8 926 076 station-days, 1 570 stations), exported to the PoKriging LOO CSV layout (`tools/export_loo_parquet_csv.py`,
float32 printed with 9 significant digits = exact round trip). Reference: the stored Python tables
`PoGrid_v5/summary/uncertainty_calibration_{precip,temp}.csv` (the v5 temp table is a copy of v4's); the s9 code
re-run on the same parquet files with elevations from `stations_temp.csv` reproduces them exactly (max |d| = 0).

| variable | rows | k identical | coverage95_raw | coverage95_cal_crossval | e_q05..e_q95 | n | share_zero_std | Java time |
|---|---|---|---|---|---|---|---|---|
| precip | 25/25, same order | 25/25 (max \|d\| 0) | 25/25 | 25/25 | 125/125 | 25/25 | 24/25, 1 value 4.1e-17 (1 ulp) | 20 s |
| temp | 21/21, same order | 21/21 (max \|d\| 0) | 21/21 | 21/21 | 104/105, 1 value 8.1e-17 (1 ulp) | 21/21 | 21/21 | 28 s |

Java LOO (OMS chain PoKriging mode = loo → PoUncertainty, 10 reference dates of 2020; `simulation/tests/test_uncertainty_chain_*.sim`)
vs the s9 code applied to the same Java CSVs: precip 21 rows, temp 21 rows, all values identical (max |d| = 0).
Feed-back: PoKriging with `uncertaintyTable` = the Java table gives `uncertainty` = kriging_std × k within 4.8e-7
(float32 rounding).

## Full 1991-2020 two-part precipitation run in GEOframe vs the Python v5 product (2026-10-20)
Run through `oms3.CLI -r` from the package (`simulation/precip.sim`, method twopart, mode both, calibration auto -> monthly,
12 threads, Intel Java 11 under Rosetta): 4.9 h wall including ~1.5 h of system sleep. Outputs compared with
`PoGrid_v5/precip_anomaly` (Python production):
* Occurrence thresholds (`occurrence_thresholds_1991-01-01_2020-12-31.csv`, 1,056,888 cell-months): identical (<= 1e-6) in
  99.998 %, max difference 0.008; wet_days_target identical in 100 %.
* Daily grids, every year 1991-2020 (32.1-32.2 M cell-days per year): identical (<= 0.05 mm, the int16 quantisation) in
  99.9998-99.9999 % of cell-days; wet/dry identical; yearly domain means identical to 4 decimals; 30-year totals per cell
  within 0.1 % everywhere (max 0.099 %). The few differing cell-days (max 45 mm in one cell-day of 1995) are targets where
  the PSO takes a different path after a tiny floating-point difference.
* Station LOO (9,841,333 station-days): predictions identical (<= 1e-4) in 99.998 %; KGE 0.8867 / 0.8867, r 0.9070 / 0.9070,
  MAE 1.2457 / 1.2457 mm, wet-day frequency 25.45 / 25.45 % (observed 25.64 %), wet-day intensity 11.92 / 11.92 mm
  (observed 11.67) for Java / Python.
* Uncertainty table from the Java LOO (`uncertainty_precip_1991-01-01_2020-12-31.csv`) vs the Python v5 table: k max
  difference 0.0035 (DJF, class 0-1), coverages max 0.005.
