# PoKriging — method specification of the Python reference engine (frozen 2026-10-19)

Scope: the CORE that the GEOframe/OMS3 Java component must reproduce — monthly normals (precipitation total M,
wet-day frequency p, temperature), daily interpolation (precipitation: direct / climatology(ratio, v2.0) /
two-part (v5, final); temperature: direct / anomaly with local or no daily lapse), grid or point targets, and
leave-one-out (LOO). QC (s7) and day-shift harmonisation (s0) stay in Python; they produce the clean CSVs.

Source frozen: `~/Downloads/PoGrid_v3_code` (files dated 2026-10-12), run with `PO_STATION_SET=v4`,
`PO_OUT_ROOT=~/Downloads/PoGrid_v5` (two-part) or `~/Downloads/PoGrid_v4` (ratio v2.0, temperature).
Reference citations are `file:function`. "§A-n" refers to the ambiguity list at the end.

Verified on 2026-10-19 (tools/ and reference/ of this folder):
* the normals rebuilt from the clean CSVs with the locked parameters are **bit-identical** to
  `PoGrid_v5/normals` (all grids and station tables, max |Δ| = 0);
* the daily engine re-run on 4 production chunks (8192 cells) × 10 dates of 2020 reproduces the production
  NetCDF values **exactly** after the production quantisation (precip ratio, temperature anomaly, two-part pass-1
  π/amount and variogram ranges: 100 % identical), and the production LOO to float32 precision (|Δ| ≤ 5e-6);
* applying the production 1991-2020 monthly occurrence thresholds to the re-run π reproduces the production
  two-part final field exactly (81 920 cell-days).

---------------------------------------------------------------------------------------------------------------

## 0. Conventions and data model

| item | value | source |
|---|---|---|
| CRS | EPSG:3035 (metres) | `po_config.METRIC_CRS` |
| product grid | 494 × 330 cells, 1 km, cell centres x = 4 055 500 … (+1000), y = 2 613 500 … (−1000, north→south); 88 074 valid cells | `po_common.load_grid` |
| target order | valid cells in row-major order from the north-west (`np.nonzero(mask)`), i.e. `grid_cells.csv` row order | `po_common.load_grid` |
| cell elevation | DEM point value; holes inside the mask filled with the nearest DEM value (`elev_filled` flag, qc bit 64) | `po_common.load_grid` |
| station elevation | shapefile `elevation`; where missing, DEM value of the nearest 1 km cell (`elev_from_dem`) | `po_common.attach_station_predictors` |
| station order | column order of the clean CSV restricted to stations present in the shapefile (= engine index) | `po_common.load_station_data` |
| values | parsed to **float32**, then used as float64; `-9999` → missing; precip < 0 → missing; temp outside [−50, 50] °C → missing | `po_common.load_station_data`, `po_config.MISSING_VALUE`, `TEMP_RANGE` |
| station–station distance (variogram binning) | float32 matrix in **km** | `po_common.station_distance_km` |
| target–station distance (kriging RHS, D1/D7/D50) | float64 metres from the k-d tree, /1000 → km | `po_core.estimate_targets` |
| station–station distance inside OK matrix | float64, recomputed from x/y, /1000 | `po_core.ok_solve` |
| period | 1991-01-01 … 2020-12-31 | `po_config.PERIOD` |
| days in month | actual days of the calendar month of the date (29 in leap Februaries) for the daily transforms; 28.25 only for the reported `normal_daily_mm` | `s3_run.build_shared` (`dim`), `s1_station_normals.main` (§A-7) |

Neighbour search: k nearest by Euclidean distance (scipy cKDTree, sorted ascending). The engine station sets
contain no co-located stations (checked: 0 pairs < 1 m among the 1636 precip / 1570 temp engine stations;
the shapefiles contain 9 such temperature pairs but one of each pair has no data column). Java: sort by
(distance, station index) for a deterministic tie rule (§A-15).

---------------------------------------------------------------------------------------------------------------

## 1. Station monthly values and station normals (step s1)

`s1_station_normals.monthly_matrix`, `fill_gaps`, `normals`

1. Month value for station s, month t of 1991-01 … 2020-12 (360 months): valid days n_v, days in month n_d.
   Accepted if n_v ≥ `MONTH_COMPLETENESS`·n_d = 0.80·n_d. Value: precip = (mean of valid days) × n_d [mm/month];
   temp = mean of valid days [°C].
2. Gap filling (only missing station-months):
   * calendar-month standardisation: Z = (M − μ_s,m)/σ_s,m (σ = population std, `np.nanstd`); σ = 0 → NaN;
   * candidates: the `FILL_MAX_CANDIDATES`+1 = 31 nearest stations (incl. self), keep 0 < d ≤ `FILL_RADIUS_KM` = 100 km;
   * Pearson r on common finite Z months; keep candidates with overlap ≥ `FILL_MIN_OVERLAP_MONTHS` = 36 and
     r ≥ `FILL_MIN_R` = 0.70 (precip) / 0.85 (temp); order by r descending;
   * per candidate and calendar month: precip ratio = Σ station / Σ donor over common months; temp difference =
     mean(station − donor) over common months; NaN when < 5 common months of that calendar month;
   * estimate = donor × ratio (precip) or donor + difference (temp); for every month use the first
     `FILL_TOP_K` = 3 donors (in r order) that have a finite estimate, weights r²; precip clipped at ≥ 0.
3. Normal per calendar month: mean of the filled series over the 30 years; accepted if observed years
   ≥ `MIN_OBS_YEARS` = 5 **and** observed + reconstructed years ≥ `MIN_TOTAL_YEARS` = 20; else NaN.
   Output `station_normals_<var>.csv` (station, month, normal, normal_obs_only, n_obs_years, n_total_years, ok
   [, normal_daily_mm = normal / {31,28.25,31,30,31,30,31,31,30,31,30,31}]).
   Reference: precip 17 598 of 19 632 station-months accepted (1375 stations with 12); temp 15 280 of 18 840 (1221).

## 1w. Station wet-day frequency (step s1w, two-part only)

`s1w_wetfreq_normals.station_frequency`

* wet day: P ≥ `OCC_THRESHOLD` = 1.0 mm (§A-8); per station-month accepted if valid days ≥ 0.80 n_d;
* p(s, m) = Σ wet days / Σ valid days over the accepted months of calendar month m (pooled ratio, not a mean of
  monthly ratios); NaN if fewer than `MIN_OBS_YEARS` = 5 accepted months. **No gap filling** (§A-14).
* Transform: logit(p) with p clipped to [`PMIN`, `PMAX`] = [0.005, 0.995]; back-transform expit.
  Reference: 1426 of 1636 stations with 12 frequency normals.

---------------------------------------------------------------------------------------------------------------

## 2. Gridded normals: PRISM-type local weighted regression (step s2, s1w)

`s2_normals_grid.predict`, `annual_slopes`, `run_month`; `s1w_wetfreq_normals.slopes`, `loo`

### 2.1 Locked parameters (`PoGrid_v5/normals/params_*.json`, selected by LOO with `--tune`)

| variable | space y | K | DZ0 [m] | FACET_W | TPI0 | USE_TPI | RESID | EXTRAP_DZ [m] | EXTRAP_LAPSE [°C/m] | slope mode | slope bounds (per m) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| precip M | ln(max(M, 0.1)) | 15 | 500 | 1.0 | 1e9 | 0 | none | **0** | – | annual mean | [−0.3e-3, +0.8e-3] |
| wet-day freq p | logit(p) | 15 | 1e9 | 1.0 | 1e9 | 0 | none | **0** | – | annual mean | [−1.0e-3, +2.0e-3] |
| temp | T | 25 | 300 | 1.0 | 1e9 | 0 | none | 1e9 (none) | **−0.0075** | monthly | [−10.0e-3, +3.0e-3] |

Precip slope bounds = −0.3 … +0.8 ln-units per km, i.e. ×e^−0.3 = −25.9 %/km … ×e^0.8 = +122.6 %/km
(the code comment says "−30 % .. +120 %", §A-4). With FACET_W = 1, TPI0 = 1e9 and USE_TPI = 0 the facet and TPI
terms are inactive (§A-6); RESID = none (residual correction rejected by LOO for v4).

### 2.2 Regression at one target (`predict`)

For a target (cell, or station in LOO) and calendar month m, with the station set S_m = stations with a finite
value in month m:
1. the K+1 nearest stations of S_m; drop the target's own station (stable sort) and keep K;
   d = distances [m] ascending.
2. weights: bandwidth h = d_K·1.01 + 1.0 m; w = (1 − (d/h)²)² · 1/(1 + (|z_s − z_t|/DZ0)²)
   · [facet term = 1] · [temp: 1/(1 + (|TPI_s − TPI_t|/TPI0)²) ≈ 1].
3. target elevation used for the elevation term: precip and wetfreq z_eff = min(z_t, max_K z_s + EXTRAP_DZ)
   = min(z_t, highest of the K stations) → **no upward extrapolation**; downward extrapolation is allowed.
   temp: z_eff = z_t.
4. design centred on the target: X = [1, z_s − z_eff (, TPI_s·0 − TPI_t·0)]; WLS normal equations
   A = Xᵀ W X, b = Xᵀ W y, ridge λ_0 = 0, λ_j = 1e-9·max(mean over targets of A_jj, 1e-12) for j ≥ 1 (λ computed
   per batch of targets — mean over the batch; numerically negligible, §A-6); β = solve(A + diag λ, b).
5. slope = clip(β_1, bounds) — or the imposed per-target slope (precip/wetfreq, §2.3).
6. estimate:
   * precip/wetfreq: est = Σ w (y_s − slope·(z_s − z_eff)) / Σ w;
   * temp: re-solve the intercept (and TPI coefficient, zero column) with the bounded slope held fixed:
     est = WLS intercept of (y_s − slope·(z_s − z_t)) (≡ weighted mean when the TPI column is zero);
     then **above the highest of the K stations**: up = max(z_t − max_K z_s, 0);
     est ← est − slope·up + min(slope, EXTRAP_LAPSE)·up, i.e. the steeper of the local lapse and −7.5 °C/km
     carries the extrapolated part.
7. back-transform: precip exp(est) [mm/month]; wetfreq expit(est); temp est [°C].

### 2.3 Slopes (`annual_slopes`, `s1w.slopes`)

* precip and wetfreq: for each month m fit (2.2 steps 1–5, bounds applied) on S_m, then slope = **mean of the 12
  bounded monthly slopes**; computed for every grid cell and, leave-one-out, for every station (station dropped
  from S_m when it belongs to it). The estimate (2.2 step 6) then uses this fixed slope for all months.
* temp: the bounded monthly slope of the same fit (`--monthly` behaviour; `annual` only for precip).
* Output `elevation_slope` in the NetCDF: slope ×1000 (per km).

### 2.4 Station LOO normals (no leakage at this stage)

`run_month`: every station is predicted with itself removed from S_m and from the slope fits
(`self_all = pos`). This gives `normal_grid_loo` (`station_normals_loo_<var>.csv`) and `wetfreq_grid_loo`
(`station_wetfreq_loo_precip.csv`), used as the LOO target climatology in step 3.
Remaining leakage: the station may have served as a gap-filling **donor** of its neighbours' normals (§A-11).
Reference LOO skill: precip annual MAE 89.96 mm (7.85 %), temp monthly MAE 0.605 °C, wet days annual MAE 4.73 d.

### 2.5 Tuning (optional, not needed for the locked product)

`s2_normals_grid.main --tune` grid: K ∈ {10, 15, 25, 40, 60}, DZ0 ∈ {150, 300, 500, 1e9}, FACET_W ∈ {0.5, 1},
TPI0 ∈ {100, 1e9} (temp), USE_TPI ∈ {0, 1} (temp), RESID ∈ {none, idw}; stage 2 taper RESID_L ∈ {500 … 8000 m};
stage 3 EXTRAP_DZ ∈ {1e9, 1000, 600, 400, 250, 100, 0} (precip) or EXTRAP_LAPSE ∈ {0, −0.0055 … −0.0098} (temp);
criterion annual MAE (precip) / monthly MAE (temp). s1w: K ∈ {10, 15, 25, 40}, DZ0 ∈ {150, 300, 500, 1e9},
EXTRAP_DZ ∈ {0, 500, 1e9}, criterion MAE of annual wet days. Recommendation: Java implements predict() with
these parameters as inputs; the grid search can stay in Python.

---------------------------------------------------------------------------------------------------------------

## 3. Daily interpolation (step s3)

`s3_run.transform`, `back`, `process_day`, `twopart_day`, `build_shared`; `po_core.estimate_targets`

### 3.1 Day set-up
* active stations: finite value (after the transform) and usable; **day skipped** (all targets = fill) if
  fewer than `N_VARIO`+1 = 51 active stations (`process_day`, `twopart_day`).
* station climatology NS(s, m): the accepted station normal (s1); where NaN, the station's LOO grid normal
  (s2 `normal_grid_loo`) — flag bit 128 in LOO output. Same rule for p_station (PSF) with `wetfreq_grid_loo`.
* target climatology: grid cell → gridded normal of that cell (NG, PGF); LOO station → its LOO grid normal
  (NL, PLF); point → gridded normal of the nearest grid cell (k-d tree on cell centres).

### 3.2 Value kriged at stations (transform) and back-transform

| method | y at station | back-transform at target | kriging std reported |
|---|---|---|---|
| precip direct | P | max(est, 0) | √kvar |
| precip climatology (ratio, v2.0) | min(P / c_s, 50) with c_s = max(NS/n_d, 0.05) | max(est, 0)·c_t, c_t = max(N_t/n_d, 0.05) | √kvar·c_t |
| temp direct | T (local daily lapse, §3.4) | est + trend | √kvar |
| temp anomaly | T − NS | est (+ trend if lapse) + N_t | √kvar |
| precip two-part | §3.9 | §3.9 | §3.9 |

Constants: `MIN_DAILY_NORMAL_MM` = 0.05 mm/day, `MAX_PRECIP_RATIO` = 50.
Temperature lapse mode: direct → `local50` always; anomaly → `local50` (production) or `none`
(`--anom-lapse none`). Precip: never a lapse.

### 3.3 Neighbourhoods (`estimate_targets`)
* query the min(51, A) nearest active stations; move the target's own station (LOO) to the end (stable sort) and
  keep the first `N_VARIO` = 50 (fewer if A < 51 — possible only for the two-part amount kriging);
* the kriging uses the first `N_KRIG` = 7 of these 50 (nearest, self excluded);
* reported D1, D7 (7th) and D50 (last kept) in km.
* targets are processed in chunks of `CHUNK` = 2048 consecutive targets (seed and random-stream unit, §3.6, §A-2).

### 3.4 Per-target daily lapse (temperature, `lapse_mode = local50`)
`po_core.lapse_batch` + `estimate_targets`
1. OLS y = b0 + b1·z on the 50 neighbours (finite pairs); r², span = max z − min z, n.
   Valid only if n ≥ 3 and Σ(z − z̄)² > 1e-8.
2. fallback if b1 not finite or span < `LAPSE_MIN_SPAN_M` = 300 m or n < `LAPSE_MIN_STATIONS` = 10: refit on the
   min(101, A) nearest active stations with the target's own station set to NaN (→ 100 stations for LOO,
   **101 for grid/points**, §A-9); flag bit 2. Still invalid → b0 = b1 = 0.
3. residuals y − (b0 + b1 z_s) at the 50 neighbours are kriged (variogram on residuals).
4. trend at target b0 + b1·z_e. Anomaly mode only: z_e = clip(z_t, min z₅₀ − 300, max z₅₀ + 300)
   (`ANOM_LAPSE_MAX_EXTRAP_M` = 300 m; z₅₀ = elevations of the 50 neighbours with finite residual, also when
   the fallback lapse was used, §A-9). Direct mode: z_e = z_t (no limit).
5. flags: bit 2 lapse fallback, bit 4 positive lapse (b1 > 0, "inversion").

### 3.5 Empirical variogram (`po_core.empirical_batch_np`, `po_numba._empirical`)
On the 50 neighbours' values (residuals if lapse):
* all pairs i < j with finite distance and semivariance and d > 0: γ_ij = ½(y_i − y_j)², d_ij from the float32 km
  matrix (§0);
* `N_LAGS` = 12 equal-width bins on [0, d_max·(1 + 1e-10)], d_max = largest valid pair distance (1.0 if none);
  bin index = number of upper edges e_k = k·step (k = 1…12, last edge = d_max·(1+1e-10) exactly) with d ≥ e_k,
  capped at 11;
* per bin: count, mean lag, mean γ; bins with count < `MIN_PAIRS_PER_LAG` = 5 are masked (count set 0);
* valid variogram: total pairs ≥ max(20, 3·5) = 20, ≥ 3 unmasked bins, max γ > 1e-12.
* invalid rows: the bins are replaced by a dummy (bins 0–2 unmasked, lag = 1, 2, 3, γ = 1, count = 1), PSO is run
  anyway (it consumes the random stream), and the parameters are then replaced by the day's **global variogram**
  if it exists (flag bit 32); if the global variogram is also invalid the dummy fit is kept (§A-10).

### 3.6 Spherical model and PSO fit (`po_core.pso_spherical_np` ≡ `po_numba.pso_spherical`)
* γ(h) = nug + c·(1.5 r − 0.5 r³), r = min(h/a, 1); γ(h < 1e-9) = 0 (no nugget on the diagonal / at a
  coincident station).
* per target, over unmasked bins: maxg = max(max γ̂, 1e-8); minh = max(min lag, 1e-4); maxh = max(max lag, 1.01 minh).
* bounds: nug ∈ [0, 1.5 maxg], c ∈ [1e-12, 3 maxg], a ∈ [max(0.25 minh, 1e-4), 2 maxh] (km).
* weights (`FIT_WEIGHTING` = pair-count): w_l = n_l / max(Σ n_l, 1); masked bins contribute 0 (h = 1, γ̂ = 0, w = 0).
* objective f = Σ_l w_l (nug + c·sph(h_l/max(a, 1e-12)) − γ̂_l)², summed in numpy pairwise order (numba
  replicates the 8-way unrolled pairwise sum, `po_numba._pairwise_sum12`).
* `PSO_PARTICLES` P = 24, `PSO_ITERATIONS` = 35, inertia 0.72, c1 = c2 = 1.45.
  init x = lo + U·(hi − lo), v = N(0, 0.06)·(hi − lo); personal best = x, global best = argmin (first index).
  each iteration: v ← 0.72 v + 1.45 r1 (pbest − x) + 1.45 r2 (gbest − x); x ← clip(x + v, lo, hi) (velocity not
  clipped); f(x); pbest replaced if f < pbest_f (strict); gbest ← pbest of argmin pbest_f (first index) if
  strictly smaller than current gbest_f.
* random stream: `numpy.random.default_rng(seed)` (SeedSequence → PCG64) for the whole chunk of B targets that
  need a variogram, drawn in this order: U (B,P,3) `random()`, N (B,P,3) `normal(0, 0.06)` (ziggurat),
  then for t = 0…34: r1 (B,P,3), r2 (B,P,3) `random()` (C order).
* seed = `po_core.seed_for(day_tag, s0)`: BLAKE2b (digest 8 bytes) of the UTF-8 strings str(part) concatenated,
  little-endian integer & 0xFFFFFFFF. day_tag = `"{variable}|{mode}|{YYYY-MM-DD}"` (variable precip|temp,
  mode direct|anomaly); two-part: `"precip|twopart|{date}|occ"` and `"…|amt"`; s0 = index of the first target
  of the chunk in the target list (0, 2048, 4096, …). Global variogram: `seed_for(day_tag, "global")`.
* range-bound flag (bit 16): a ≥ 0.975·hi_a.
* measured seed sensitivity (reference/seed_sensitivity.csv, 10 dates, 8192 cells): changing only the seed
  changes the fitted range by up to 275 km (precip ratio) / 308 km (temp) (p95 ≈ 125 km) and the daily values by p95/p99 = 1.32/5.7 mm
  (precip ratio), 0.18/0.50 °C (temp anomaly), 1.2/15.7 mm (two-part final) — see §A-1.

### 3.7 Global variogram (fallback only, `po_core.global_variogram`)
One empirical variogram + PSO on **all** active stations (incl. a LOO target, §A-10); y = the kriged quantity,
for lapse modes y minus a single linear fit of y on elevation over the active stations (`np.polyfit` deg 1).
Invalid → none.

### 3.8 Ordinary kriging (`po_core.ok_solve`)
* k = 7 nearest; A_(k+1)×(k+1): γ between neighbours (diagonal 0), last row/column 1, corner 0; RHS γ(d_target),
  1. Solve (LAPACK gesv batched; if any matrix of the batch is singular, solve one by one; a failed solve gives
  NaN estimate). est = Σ w_i y_i; kvar = max(Σ w_i γ_i0 + μ, 0).
* precip **zero shortcut** (only variable = precip, no lapse, i.e. direct and ratio; NOT two-part):
  if all 7 kriging neighbours have y ≤ `ZERO_TOL` = 0 (exactly zero rain) → est = 0, kvar = 0, no variogram
  fitted (and no random numbers drawn for that target), flag bit 1 (§A-8).
* weight diagnostics (`weight_diagnostics`): n_neg = #(w < −1e-10); N_eff = 1/Σ p_i², p = |w|/Σ|w|;
  R_pos_½, R_pos_⅔, R_pos_⅓ = largest distance among neighbours with w ≥ {½, ⅔, ⅓}·max positive w;
  R_abs90 = distance of the first neighbour (ascending) where cumulative p ≥ 0.9; Σ|w| > `OSC_SUM_ABS_W` = 2 →
  flag bit 8 (oscillatory).

### 3.9 Two-part precipitation (`s3_run.twopart_day`, `finalize_twopart`, `_thresholds`, `_calibrate`, `_mix_std`)
Pass 1 (per day, independent):
* station: p_s = PSF(s, m) (§3.1 fallback), i_s = NS(s, m) / (max(p_s, `P_FLOOR` = 0.01)·n_d);
  wet_s = P_s ≥ 1.0 mm.
* occurrence: y_occ = wet_s − p_s at all active stations; kriged with the full machinery (50/7, PSO, own seed tag
  `|occ`, global fallback) without zero shortcut and without lapse; π = clip(p_t + ô, 0, 1).
* amount: y_amt = min(P_s / max(i_s, 0.05), 50) on WET active stations only; if ≥ `N_KRIG`+1 = 8 wet stations:
  kriged (tag `|amt`, global fallback only if ≥ 10 wet), amount = max(est, 0)·i_t, s_amt = √kvar·i_t, with
  i_t = N_t / (max(p_t, 0.01)·n_d); else amount = median(y_amt over wet)·i_t (0 if none), s_amt = 0 (§A-21).
* production storage (quantisation applied BEFORE calibration): π u2 1e-4, amount i2 0.1 mm, s_amt i2 0.1 mm.

Pass 2 (calibration over the whole run, `finalize_twopart`):
* grid / points: per target and calendar month if every calendar month has ≥ `CALIB_MIN_DAYS` = 100 run days
  ("monthly", production 1991-2020), otherwise pooled over all run days ("pooled"); in each group
  k = rint(Σ_days with finite π of p_t) (half-to-even), clipped to [0, n]; threshold t = k-th largest π
  (t = +∞ if k = 0; t ≤ 0 → 1e-9); **wet ⇔ π ≥ t** (float32 comparison of the quantised π; ties are wet, so a group
  can have more than k wet days, §A-3).
* LOO and points: `_calibrate`: same grouping rule decided on the number of distinct dates per month; wet ⇔
  rank(π, descending, method "first") ≤ k and π > 0 (unquantised float32 π; tie rule differs from the grid, §A-3).
* final value = amount (quantised) if wet, else 0. std = √(π s_amt² + π(1 − π) amount²) (hurdle variance; uses
  π, not the wet flag).

### 3.10 Outputs of the daily stage (`s3_run.write_year`, `ENC`)
* `<var>_<Y>.nc`: value and kriging_std, int16, scale 0.1 mm (precip) / 0.01 °C (temp), fill −32768,
  dims (time, y, x), time = days since 1991-01-01, zlib 6, chunks (1, ny, nx).
* `diagnostics_<var>_<Y>.nc`: D1, D7 (u2 0.1 km), D50, range, R_pos_* , R_abs90 (u2 0.5 km), nugget_fraction
  (u1 0.01), log10_sill (i2 0.01), lapse (i2 0.01 °C/km), lapse_r2 (u1 0.01), n_negative_weights (u1),
  n_eff (u1 0.05), qc_flags (u1); two-part adds range_amount, nugget_fraction_amount, log10_sill_amount,
  R_pos_half_amount, wet_probability (u1 1/250). Encoding: rint(x/scale), clipped to [min+1 (signed), max−1
  (unsigned)], NaN → fill.
* `loo_<var>_<Y>.parquet`: date, station, obs, pred, kriging_std, D1, D7, qc_flags (+ two-part: pi, amount,
  p_target, wet_pred).
* qc flag values (bit-OR): 1 zero-shortcut, 2 lapse fallback, 4 positive lapse, 8 oscillatory, 16 range bound hit, 32 global
  variogram fallback, 64 elevation filled (grid), 128 station-normal fallback (LOO).

### 3.11 Calibrated uncertainty (post-processing; `s9_uncertainty_calibration.main`, `s13_zenodo_package.calibrated_std`)
k(season, class) = q95(|obs − pred| / kriging_std over LOO with std > 0) / 1.96, needs > 50 LOO values;
classes: precip by the estimate value, right-closed (−∞,0], (0,1], (1,5], (5,20], (20,50], (50,∞);
temp by elevation (−∞,500], (500,1000], (1000,1500], (1500,2000], (2000,∞) m; seasons DJF/MAM/JJA/SON.
Published `<variable>_uncertainty` = kriging_std × k. Java: apply a k table read from CSV (computing k may stay
in Python; it is a LOO post-process).

---------------------------------------------------------------------------------------------------------------

## 4. Leave-one-out rules (what the withheld station never sees)

1. its own daily value: removed from the 50/7 neighbourhood (`estimate_targets`, stable-sort drop), from the
   lapse fit (50 and fallback set: own value → NaN) and from the amount kriging;
2. its own normal: the target climatology is the LOO grid normal (§2.4) for M, p and T — the station's own
   normal is not used at the target;
3. its own normal still enters as **neighbour climatology for other targets** (correct: it is not the target);
4. known residual leakage (kept for equivalence, documented §A-10, §A-11): the day's global fallback variogram
   includes the station; the station may be a gap-filling donor of neighbouring station normals (s1);
5. LOO targets are all active stations of the day, in station-index order, as one chunk (≤ 2048), seed s0 = 0.
6. two-part LOO calibration uses the LOO frequency normal p_LOO (not the station's observed frequency).
7. evaluated variables: obs, pred, kriging_std, D1, D7, flags; precip daily metrics use wet = obs ≥ 0.1 mm
   (`WET_THRESHOLD`) for wet-day MAE/pbias.

---------------------------------------------------------------------------------------------------------------

## 5. Ambiguities and decisions needed (flagged)

| # | issue | impact | proposal |
|---|---|---|---|
| A-1 | **PSO randomness.** The fitted variogram is weakly identified: an independent seed changes ranges by up to 275–308 km and daily values by p99 5.7 mm / 0.50 °C (reference/seed_sensitivity.csv). LOO MAE over the 10 dates changes only slightly (precip ratio 3.913 → 3.906/3.916 mm; temp 0.8172 → 0.8162/0.8174 °C; two-part 7.47 → 7.42/7.33 mm with the 10-day calibration). | A Java port with its own RNG can only be "statistically equivalent"; per-cell tolerances would have to be ~6 mm / 0.5 °C, too loose to detect bugs. | **Port numpy's stream exactly**: BLAKE2b seed_for, SeedSequence(seed) → PCG64, `random()` (53-bit), `standard_normal` (ziggurat, 256 tables), same draw order. ~300 lines of Java; then the test can be near-bitwise (1e-9). Keep a `rngMode = numpy|java` switch. |
| A-2 | Random stream per **chunk of 2048 targets**: a cell's result depends on the other targets of its chunk (order, zero-shortcut mask). A point target at a cell centre does not reproduce the grid value. | Equivalence needs the same target list and chunking. | Java reproduces chunking in `numpy` mode (CHUNK = 2048, s0 = index of the chunk's first target). Optionally a future per-target seed (not equivalent to v2/v5). |
| A-3 | Two-part final values depend on the calibration sample: monthly over the whole run if every calendar month has ≥ 100 days, else pooled. 10-day run vs production: wet/dry agreement 82.1 %. Grid ties: π ≥ t (more than k wet days possible); LOO: rank "first" ≤ k and π > 0. Grid π quantised to 1e-4 before calibration, LOO π not. | Short Java runs cannot reproduce production final values. | Java: two-pass (pass-1 π stored), modes `auto|monthly|pooled|table`; `table` reads per cell×month thresholds (reference/twopart_production_calibration/thresholds_cells.csv, verified exact). Keep both tie rules for equivalence; compare integer π codes. |
| A-4 | Precip slope bounds stated as "−30 %..+120 %/km" are ln-slopes −0.3/+0.8 per km = −25.9 %/+122.6 %. | Documentation only. | Specify bounds in ln units per metre in the .sim. |
| A-5 | Slope mode differs: precip and wetfreq use the 12-month mean of bounded monthly slopes; temperature uses bounded monthly slopes. | Must be reproduced. | `.sim` key `slopeMode = annual|monthly` per variable. |
| A-6 | Facet, TPI and residual terms are inactive under the locked parameters; ridge λ uses a batch mean of A_jj (depends on batch). | Negligible numerically (λ ≈ 1e-9 relative). | Java implements K, DZ0, EXTRAP_DZ, EXTRAP_LAPSE, bounds; facet/TPI/RESID optional (not needed). Use λ per target = 1e-9·A_jj — accept ~1e-12 differences. |
| A-7 | Daily climatology divides the monthly normal by the actual days of the month (29 in leap Februaries) while the normal is a 30-year mean monthly total. | ±3.5 % in leap Februaries for ratio/two-part at the target (cancels at stations). | Keep for equivalence; note in paper if relevant. |
| A-8 | Three wet thresholds: zero shortcut uses P ≤ 0 (exact 0 mm), statistics use 0.1 mm, two-part occurrence 1.0 mm. | Must be reproduced. | Separate .sim keys `zeroTol`, `wetThresholdStats`, `occThreshold`. |
| A-9 | Lapse fallback uses 101 stations for grid/points but 100 for LOO; the ±300 m clip uses the 50-station range even after fallback; direct mode has no clip. | Reproduce. | As coded. |
| A-10 | Global fallback variogram includes the LOO station; if both local and global variograms are invalid the dummy-variogram PSO fit is used. | Small leakage; rare. | Keep for equivalence; flag bit 32 count reported. |
| A-11 | s1 gap filling may use the withheld station as a donor for neighbours' normals ("LOO without leakage" covers only the station's own normal). | Small optimistic bias of LOO. | Keep; optional strict mode later (re-run s1 without the station is too expensive per station). |
| A-12 | Mixed precision: variogram distances float32 km, kriging distances float64; station values float32. | Needed for bitwise equivalence. | Java mirrors it (float sdmat, double OK). |
| A-13 | Two-part amount fallback (< 8 wet stations): median ratio × i(x), std 0. | Rare. | As coded. |
| A-14 | M includes reconstructed months but p does not → i = M/(p·n_d) mixes samples. | Method-level, not a port issue. | Keep; note for the paper. |
| A-15 | k-nearest ties: cKDTree order. No ties in the engine set. | none now | Java: sort by (distance, index). |
| A-16 | `temperatureTrend = global` does not exist in Python (only `local50` and `none`; a global linear fit is used only to detrend the global fallback variogram). | No reference for "global". | Decide: (a) global = one daily OLS on all active stations (new, no Python reference), or (b) drop "global" and expose `local|none`. |
| A-17 | Uncertainty k classes: precip by the estimate itself (right-closed bins). | Reproduce. | As coded. |
| A-18 | Day skipped when < 51 active stations; fill values written. | Reproduce. | As coded. |
| A-19 | Python `process_day` grid value of the ratio method is `max(est,0)·c_t` where est includes the zero shortcut; kvar 0 there → kriging_std 0. | Reproduce. | As coded. |
