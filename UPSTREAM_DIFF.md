# Krigings working copy — changes with respect to the upstream code

Base: **GEOframe Kriging 0.9.9** (`Kriging-2023-01-11-0.9.9.jar` of the user's Brenta OMS project; no source was available,
so the jar was decompiled with CFR 0.152 into `../upstream/Kriging-0.9.9_from_Brenta/src_decompiled`).
GitHub `geoframecomponents/Krigings` master (`../upstream/Krigings-master`, June 2022) is used for the original comments,
for the files whose bytecode is unchanged in 0.9.9, and for its JUnit tests.
Package names and every public `@In`/`@Out` field of 0.9.9 are kept; all additions are new options whose defaults give the
0.9.9 behaviour.

## 1. Source reconstruction (45 files of 0.9.9)

Every class was compiled from GitHub master (with debug info), decompiled with the same CFR version, and compared with the
decompiled 0.9.9 class (comments and imports ignored).
* bytecode-equivalent → the **GitHub source** (with its comments) is used;
* the 6 variogram models whose only change is a new `computeGradient()` returning `null` → GitHub source + that method;
* everything else → the 0.9.9 decompilation, cleaned:
  * CFR banner removed and the GitHub licence header restored;
  * casts such as `(Object)((Object)this)` and `isNovalue((double)x)` simplified;
  * CFR's "removed try catching itself" notes dropped. They refer to try-with-resources on feature iterators, which behave
    the same as the original try/finally.

| file | origin |
|---|---|
| `curvefitter/KrigingParamValidator.java` | 0.9.9 decompiled, cleaned |
| `curvefitter/VariogramFitter.java` | 0.9.9 decompiled, cleaned |
| `curvefitter/VariogramFunction.java` | 0.9.9 decompiled, cleaned |
| `experimentalVariogram/ExperimentalVariogram.java` | 0.9.9 decompiled, cleaned |
| `experimentalVariogram/OmsVariogramOLD.java` | github source (bytecode-equivalent to 0.9.9) |
| `krigingsPointCase/GlobalParameterEvaluator.java` | 0.9.9 decompiled, cleaned |
| `krigingsPointCase/Krigings.java` | 0.9.9 decompiled, cleaned |
| `krigingsPointCase/MaxDistance.java` | github source (bytecode-equivalent to 0.9.9) |
| `krigingsPointCase/Model.java` | github source (bytecode-equivalent to 0.9.9) |
| `krigingsPointCase/NumberOfStations.java` | github source (bytecode-equivalent to 0.9.9) |
| `krigingsPointCase/ResidualsEvaluator.java` | 0.9.9 decompiled, cleaned |
| `krigingsPointCase/SimpleModelFactory.java` | github source (bytecode-equivalent to 0.9.9) |
| `krigingsPointCase/StationsSelection.java` | 0.9.9 decompiled, cleaned |
| `krigingsPointCase/Utility.java` | 0.9.9 decompiled, cleaned |
| `krigingsPointCase/VariogramParamsEvaluator.java` | 0.9.9 decompiled, cleaned |
| `krigingsRasterCase/Krigings.java` | github source (bytecode-equivalent to 0.9.9) |
| `krigingsRasterCase/MaxDistance.java` | github source (bytecode-equivalent to 0.9.9) |
| `krigingsRasterCase/Model.java` | github source (bytecode-equivalent to 0.9.9) |
| `krigingsRasterCase/NumberOfStations.java` | github source (bytecode-equivalent to 0.9.9) |
| `krigingsRasterCase/SimpleModelFactory.java` | github source (bytecode-equivalent to 0.9.9) |
| `krigingsRasterCase/StationsSelection.java` | github source (bytecode-equivalent to 0.9.9) |
| `leaveOneOut/Krigings.java` | github source (bytecode-equivalent to 0.9.9) |
| `leaveOneOut/MaxDistance.java` | github source (bytecode-equivalent to 0.9.9) |
| `leaveOneOut/Model.java` | github source (bytecode-equivalent to 0.9.9) |
| `leaveOneOut/NumberOfStations.java` | github source (bytecode-equivalent to 0.9.9) |
| `leaveOneOut/SimpleModelFactory.java` | github source (bytecode-equivalent to 0.9.9) |
| `leaveOneOut/StationsSelection.java` | github source (bytecode-equivalent to 0.9.9) |
| `linearSistemSolver/SimpleLinearSystemSolverFactory.java` | github source (bytecode-equivalent to 0.9.9) |
| `oldPointCase/Krigings.java` | 0.9.9 decompiled, cleaned |
| `pathGenerator/PathGenerator.java` | github source (bytecode-equivalent to 0.9.9) |
| `theoreticalVariogram/Bessel.java` | github source + 0.9.9 computeGradient() returning null |
| `theoreticalVariogram/Circular.java` | github source + 0.9.9 computeGradient() returning null |
| `theoreticalVariogram/Exponential.java` | 0.9.9 decompiled, cleaned |
| `theoreticalVariogram/Gaussian.java` | 0.9.9 decompiled, cleaned |
| `theoreticalVariogram/Hole.java` | github source + 0.9.9 computeGradient() returning null |
| `theoreticalVariogram/Linear.java` | 0.9.9 decompiled, cleaned |
| `theoreticalVariogram/Logarithmic.java` | 0.9.9 decompiled, cleaned |
| `theoreticalVariogram/Model.java` | 0.9.9 decompiled, cleaned |
| `theoreticalVariogram/Pentaspherical.java` | github source + 0.9.9 computeGradient() returning null |
| `theoreticalVariogram/Periodic.java` | github source + 0.9.9 computeGradient() returning null |
| `theoreticalVariogram/Power.java` | 0.9.9 decompiled, cleaned |
| `theoreticalVariogram/SimpleModelFactory.java` | github source (bytecode-equivalent to 0.9.9) |
| `theoreticalVariogram/Spherical.java` | 0.9.9 decompiled, cleaned |
| `theoreticalVariogram/Spline.java` | github source + 0.9.9 computeGradient() returning null |
| `theoreticalVariogram/TheoreticalVariogram.java` | github source (bytecode-equivalent to 0.9.9) |

Reconstruction fixes:
* `theoreticalVariogram/Bessel.bessi1`: CFR lost an assignment (`if (d < 3.75)` with `d` uninitialised). The GitHub source
  was used instead.
* `krigingsRasterCase/Krigings` and `experimentalVariogram/OmsVariogramOLD`: still import `com.vividsolutions.jts`. These two
  classes cannot even be loaded from the original 0.9.9 jar with the HortonMachine 0.10.4/0.10.5 libraries
  (`NoClassDefFoundError com/vividsolutions/jts/geom/Coordinate`). The imports were changed to `org.locationtech.jts`, the
  same fix N. Tubini applied to the point case in 2021.

## 2. Refactoring of the duplicated packages (behaviour unchanged)

New package **`krigingsCore`** holds the code that was copied three times in `krigingsPointCase`, `krigingsRasterCase` and
`leaveOneOut`:

| shared class | replaces | package wrappers kept (same names, same public API) |
|---|---|---|
| `krigingsCore.Model` | the 3 `Model` interfaces | `<pkg>.Model extends krigingsCore.Model` |
| `krigingsCore.NumberOfStations`, `MaxDistance`, `SimpleModelFactory` | the 3 copies | `<pkg>.NumberOfStations extends krigingsCore.NumberOfStations`, etc. |
| `krigingsCore.StationsSelection` | the 3 copies | point: no change; leaveOneOut: `idOut` through `include(id)`; raster: `rasterLegacy = true` |
| `krigingsCore.KrigingSystem` | private `covMatrixCalculating`, `knownTermsCalculation`, `variogram` of the 3 `Krigings` | the private methods delegate to it |

Upstream quirks deliberately kept, for identical results:
* `krigingsRasterCase.MaxDistance` never stores maxdist (it stays 0).
* The raster-only `StationsSelection` measures distances from `xStationInitialSet[n1]`. This class is not used by
  `krigingsRasterCase.Krigings`, which uses the point-case selection.
* In 0.9.9 the point case re-executes one `StationsSelection` object for every target. Its `n1` / `areAllEquals` state
  therefore accumulates. With `doLogarithmic` the values in `inData` are log-transformed in place on every re-execution.

## 3. Behaviour changes (explicit, small)

1. `krigingsPointCase.Krigings`:
   * **Fixed variogram parameters are now honoured.** This applies when `nugget`, `sill` and `range` (not all 0) are given
     together with `pSemivariogramType`. 0.9.9 accepted such parameters in `verifyInput` but then reset all three to 0 at
     every step, so it always fitted them.
   * The GitHub API is restored as additional inputs:
     * `inParamTS` gives the parameters per time step;
     * `doIncludezero` is the old spelling; setting either flag to false excludes zero values.
   * None of the Brenta sims set these keys, so their results are unchanged (section 5).
2. `krigingsCore.StationsSelection` carries the 0.9.9 `doLogarithmic` field (default false) for all packages.
3. `theoreticalVariogram.Spherical`: the new static method `semivariancePoOrder(h, nug, c, a)` serves the local PSO
   kriging. It evaluates the same formula in the operation order of the Python Po engine, using `StrictMath.pow` so that
   results do not depend on the JVM or CPU. `computeSemivariance()` is unchanged.

## 4. New code (all OMS-independent except the two components)

| package | content |
|---|---|
| `utilities.random` | `NumpyRandom` (bit-exact numpy `default_rng`: SeedSequence → PCG64 → `random()`; ziggurat `normal()`), `ZigguratTables` (numpy tables, BSD licence credited), `Blake2b` (`seed_for`), `Rng` (numpy or Java stream) |
| `krigingsCore.local` | `KdTree`, `Variogram` (Python binning), `SphericalPso` (24 × 35 PSO, numpy pairwise sums), `OrdinaryKriging` (+ weight diagnostics), `LocalLapse`, `DayKriging` (one day at many targets, chunked random streams, multithreaded), `NpMath` |
| `pokriging.io` | `stations.csv` + wide `data.csv` reader, ESRI-ASCII DEM / `grid_cells.csv` targets, normals readers (NetCDF/CSV), CF NetCDF writer (NetCDF-3 written by netcdf-java, converted to NetCDF-4/zlib by `nccopy`) |
| `pokriging.engine` | `KrigingConfig`, `KrigingEngine` (direct / climatology / two-part precipitation, direct / anomaly temperature, local / global / no lapse, grid / points, LOO without leakage in the same run, two-pass two-part calibration or threshold table, uncertainty table), `TwoPart`, `Cli` |
| `pokriging.normals` | `NormalsEngine` + `PrismRegression`: station normals with gap filling, PRISM-type gridded and LOO normals, wet-day frequency normals |
| `pokriging.PoKriging`, `pokriging.PoNormals` | thin OMS3 components (`@In` per setting, `@Execute`), generated from / mapped onto the config classes |

New option in `krigingsPointCase.Krigings` (default off): `doLocalVariogram`, plus `nVario`, `nKrig`, `localTrend`
(none / local / global), `doZeroShortcut`, `rngMode`, `seedTag` and `threads`, and the output `outKrigingVariance`.
With it, the existing per-time-step OMS chains (OMS time-series reader → Krigings → writer) use the local PSO variograms.
The existing Levenberg-Marquardt global/local fit of 0.9.9 remains the default.

## 5. Backward-compatibility tests

| test | result |
|---|---|
| Brenta `330.sim` (krigingsPointCase, temperature, exponential, global parameters), 2018-10-01 … 10-05 hourly, original 0.9.9 jar vs this build (same Brenta libs) | output CSV identical (128 lines; only the `Created` line differs) |
| Brenta `param_evaluator_precip.sim`, `param_evaluator_temp.sim` (GlobalParameterEvaluator over the whole hydrological year 2018/19, LM fit) | logs identical: precip linear nugget 0, sill 1.2892966988528525, range 5459.205645503307; temp exponential nugget 0, sill 238.89893330101154, range 220702.0965546152 |
| GitHub JUnit tests (6 classes, 8 tests) with this build | 6 pass. 2 point-case tests fail (`testKriging2`: ArrayIndexOutOfBounds with all-novalue data; `testKriging5`: expected 10.0, got 0.0) **exactly as with the original 0.9.9 jar**: these are 0.9.9 regressions with respect to GitHub master, not reconstruction errors. |
| outputs written by those tests: original 0.9.9 jar vs this build | 14 of 14 files identical. Versus GitHub master, 20 of 21 are identical; `theoretical_VGM.csv` differs because the 0.9.9 variogram models changed (e.g. `Power` uses dist² and `NaN`/`MAX_VALUE` guards). |
| refactor of the three packages (on GitHub master, before rebasing to 0.9.9) | 8/8 tests OK, 21/21 output files identical |
| new `TestLocalVariogram` | OK (finite estimates, variance ≥ 0, 1 vs 4 threads bitwise identical) |
| `pokrigingTests.RngTest` | 16 800 numpy draws and 6 `seed_for` values bitwise identical |

Scripts: `tools/test_upstream_regression.sh ../tests/github_baseline_outputs` (outputs of the unmodified GitHub build), `tools/source_origin.tsv`, `../tests/brenta_compat/` (sims and original outputs).

* 2026-10-19: `TargetResult.resid` / `TargetResult.trend` (components of `est`) for the `outputVariables` intermediates
  of PoKriging; `DayKriging` fills them. Results unchanged.
