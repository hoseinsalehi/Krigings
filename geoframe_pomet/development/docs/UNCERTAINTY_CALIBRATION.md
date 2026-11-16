# Calibration of the published uncertainty (s9) — method and GEOframe port

Source of the method: `PoGrid_v3_code/s9_uncertainty_calibration.py:main` (the only function of the script).
Use of the result: `PoGrid_v3_code/s13_zenodo_package.py:k_lookup` / `:calibrated_std` (Python product) and
`Krigings/src/main/java/pokriging/engine/KrigingEngine.java:readUncertainty` / `:uncertaintyK` (GEOframe PoKriging,
key `uncertaintyTable`). Port: `pokriging.uncertainty.UncertaintyCalibration` (OMS-independent) and the OMS component
`pokriging.PoUncertainty`.

## 1. Purpose

The ordinary-kriging standard deviation `kriging_std` (σ) is a model quantity. Its nominal ±1.96σ interval does not
contain 95 % of the leave-one-out (LOO) errors (temperature: 83 %, precipitation wet classes: 89–94 %). s9 estimates,
per class, a multiplicative factor k so that the interval ±1.96·k·σ contains 95 % of the LOO errors, checks the factor by
cross-validation over two disjoint halves of the stations, and lists empirical error quantiles as a distribution-free
alternative. The published uncertainty is `<var>_uncertainty = k(season, class) · σ` (1-σ equivalent).

## 2. Input data (s9:main, lines "fs = ..." to "dropna()")

* All files `<OUT_ROOT>/<run>/loo_<var>_*.parquet` (run = `precip_anomaly` = two-part v5 product,
  `temp_anomaly` = anomaly temperature v4/v5), sorted by file name (= by year), concatenated in that order.
* Columns used: `date, station, obs, pred, kriging_std` (all float32 in the files). No other LOO column is used
  (D1, D7, qc_flags, pi, amount, ... are ignored). For the two-part run `pred` is the final LOO estimate
  (0 on dry-predicted days, the amount otherwise) and `kriging_std` is the mixture standard deviation of the two parts
  (`TwoPart.mixStdF`).
* Rows with a missing value in any of the five columns are removed (`DataFrame.dropna()`), before anything else.
* Station elevation (temperature classes only): `po_common.load_station_data(var)["elev"]` = shapefile elevation,
  1-km DEM where missing — identical to the `elevation` column of the converted `stations_<var>.csv`
  (`tools/convert_po_inputs.py`).

## 3. Definitions

Error: e = obs − pred (float32 arithmetic, as in pandas on float32 columns).

Season of the date's calendar month m: DJF = {12, 1, 2}, MAM = {3, 4, 5}, JJA = {6, 7, 8}, SON = {9, 10, 11}
(December is pooled with January–February of the same calendar year; irrelevant here because all years are pooled).

Class (`pd.cut`, right-closed intervals (a, b], values outside all bins → class "nan"):

| variable | classified quantity | classes |
|---|---|---|
| precip | LOO estimate `pred` [mm d⁻¹] | "0" = (−1, 0], "0-1" = (0, 1], "1-5" = (1, 5], "5-20" = (5, 20], "20-50" = (20, 50], ">50" = (50, 10⁹] |
| temp | station elevation [m] | "<500" = (−1000, 500], "500-1000", "1000-1500", "1500-2000", ">=2000" = (2000, 9000] |

Note that the labels "<500" and ">=2000" are inexact: 500 m belongs to "<500", 2000 m to "1500-2000".
Precipitation is classified by the *estimate* (known at every cell), not by the observation, so that k can be
applied on the grid.

Station halves: h(station) = (index of first appearance of the station id in the concatenated, NaN-free table) mod 2
(`pd.factorize(L.station)[0] % 2`). The order is therefore "first year, then date, then the station order within the
file", not the sorted station ids. The halves split *stations* (all days of a station are in the same half), so the
test half is spatially independent of the training half at the station level.

Groups: every (season, class) combination that occurs, in lexicographic order of the strings (pandas `groupby`
sort; e.g. "0","0-1","1-5","20-50","5-20",">50"), followed by one pooled group (season = "all", cls = "all").

## 4. Statistics per group g (n = number of rows of g)

Let P = {rows of g with σ > 0} (σ = 0 occurs for precipitation when all kriging neighbours are dry; then the
estimate is exactly 0 with no stated uncertainty) and z = |e| / σ on P (float32).

| column | formula |
|---|---|
| `n` | n (all rows, including σ = 0) |
| `share_zero_std` | #(σ ≤ 0) / n |
| `coverage95_raw` | #{P : \|e\| ≤ 1.96·σ} / #P   (NaN if P is empty) |
| `scale_k` | k = Q₀.₉₅(z) / 1.96 if #P > 50, else NaN |
| `coverage95_cal_crossval` | for h ∈ {0, 1}: train = P ∩ half h, test = P ∩ half (1−h); if #train > 50 and #test > 50: k_h = Q₀.₉₅(z on train)/1.96, c_h = #{test : \|e\| ≤ 1.96·k_h·σ} / #test. Column = mean of the available c_h (NaN if none) |
| `e_q05 … e_q95` | Q_p(e) for p = 0.05, 0.25, 0.50, 0.75, 0.95 on all n rows (including σ = 0) |

Q_p is the linear-interpolation quantile (numpy default, Hyndman–Fan type 7): with the sorted sample x₍₀₎ ≤ … ≤ x₍N−1₎,
v = (N − 1)·p, j = ⌊v⌋, γ = v − j, d = x₍j+1₎ − x₍j₎ (float32), Q = x₍j₎ + d·γ if γ < 0.5, else x₍j+1₎ − d·(1 − γ)
(numpy `_lerp`, evaluated in float64). For the e quantiles pandas evaluates p as (100·p)/100.

Properties. Because 1.96·k·σ = Q₀.₉₅(|e|/σ)·σ, the in-sample coverage of the calibrated interval is 95 % by
construction; `coverage95_cal_crossval` is the honest out-of-sample estimate (station split). k < 1 (class "0",
JJA light rain) means the kriging σ is too wide there; k > 1 means too narrow. The calibration is unconditional
within a class: it does not correct heteroscedasticity inside a class, and it assumes the LOO error at a station
is representative of the error at a grid cell (the LOO estimate uses one station fewer than the grid estimate,
so k is slightly conservative where stations are dense).

## 5. Application on the grid

`<var>_uncertainty(cell, day) = σ(cell, day) · k(season(day), class(cell, day))`, class from the cell estimate
(precip) or the cell elevation (temp). Python (`s13_zenodo_package.py:calibrated_std`) uses the season×class rows
only (no "all" row; a NaN k gives NaN; lower precipitation edge −∞). Java (`KrigingEngine.uncertaintyK`) falls back to
the (all, all) k where the class k is NaN and gives NaN for an estimate ≤ −1. For the v4/v5 tables no class k is NaN,
so both give the same field.

## 6. GEOframe port

`pokriging.uncertainty.UncertaintyCalibration` reproduces §2–§4 on the LOO tables written by PoKriging
(`loo_<var>_<Y>.csv` or `.csv.gz`; columns date,station,obs,pred,kriging_std, any further column is ignored), read in
year order and in file row order, rows restricted to [periodStart, periodEnd]. Arithmetic follows the float32/float64
pattern of pandas/numpy described above (float32 e, z, 1.96·σ; float64 quantiles and k).

Outputs (period suffix = first and last LOO date used, as for `occurrence_thresholds_<start>_<end>.csv`):

* `uncertainty_<var>_<start>_<end>.csv` — exactly the s9 columns
  (variable,season,cls,n,share_zero_std,coverage95_raw,scale_k,coverage95_cal_crossval,e_q05,e_q25,e_q50,e_q75,e_q95);
  it is read back by PoKriging with `uncertaintyTable`.
* `uncertainty_summary_<var>_<start>_<end>.csv` — extension (not in s9): for each group and each nominal level
  p in `nominalLevels` (default 0.5,0.8,0.9,0.95): z_p = Φ⁻¹((1+p)/2) (1.96 for p = 0.95, the s9 convention), n, n_pos,
  coverage_raw (|e| ≤ z_p σ), scale_k = Q_p(z)/z_p, k_half0, k_half1, coverage_cal_crossval (same station-half scheme),
  coverage_empirical_crossval (interval [Q_(1−p)/2(e), Q_(1+p)/2(e)] of the training half, all rows, tested on the
  other half; the distribution-free alternative). The p = 0.95 rows equal the main table.

Chaining in OMS: PoKriging (mode = loo | both) → `@Out outLooDir` → PoUncertainty `looDir`. The table produced from a
long LOO run (e.g. 1991–2020) is then given to later grid runs with `uncertaintyTable` (it cannot be applied in the
same run: PoKriging reads `uncertaintyTable` at initialisation, before any LOO exists).
