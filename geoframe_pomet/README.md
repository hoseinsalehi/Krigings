# GEOframe-PoMet in GEOframe Krigings

This folder contains the OMS project and the Python pipeline of **GEOframe-PoMet**: 1-km daily precipitation and
daily mean air temperature for the Po River District, 1991–2020. The interpolation runs inside the GEOframe
`Krigings` component. The new code is in `src/main/java/pokriging` and `src/main/java/krigingsCore` at the
repository root. Package names and public `@In`/`@Out` fields of Krigings 0.9.9 are kept, and every addition is a
new option whose default gives the 0.9.9 behaviour. Every change to the upstream code is listed in
`../UPSTREAM_DIFF.md`. Existing sims of Krigings 0.9.9 give identical results: see `tests/brenta_compat`.

## Method

The method is a climatology-aided (anomaly) interpolation:

* **Monthly normals.** Normals of 1991–2020 are computed by PRISM-type local weighted regression on elevation:
  precipitation total, wet-day frequency (≥ 1 mm) and temperature. This is the OMS component `pokriging.PoNormals`.
* **Daily kriging.** Every cell and day gets a local spherical variogram fitted to the 50 nearest stations
  (particle swarm optimisation), then ordinary kriging with the 7 nearest stations. This is `pokriging.PoKriging`.
  * Temperature: the anomaly from the normal is kriged after removing a local daily lapse rate.
  * Precipitation (two-part scheme):
    * The occurrence anomaly (wet indicator minus wet-day frequency) is kriged and calibrated so that every cell
      keeps its normal number of wet days.
    * The amount ratio (observation / intensity normal) is kriged on the wet stations only.
* **Uncertainty.** A calibrated uncertainty is computed from the leave-one-out errors by season and class. This is
  `pokriging.PoUncertainty`.

The full specification is in `spec/METHOD_SPEC.md`. The input and sim keys are in `spec/CSV_AND_SIM_SCHEMA.md`.

## Layout

| path | content |
|---|---|
| `simulation/precip.sim`, `simulation/temp.sim` | main runs: cells → `PoNormals` → `PoKriging` (→ `PoUncertainty`); edit only the USER CONFIGURATION block |
| `simulation/advanced/` | single-component sims (normals, kriging per method, point targets, LOO) |
| `simulation/tests/` | equivalence and chain tests |
| `spec/` | method specification, CSV and sim schema, sim examples |
| `tests/` | RNG test data, equivalence date and cell lists, compatibility runs with an existing Krigings 0.9.9 project |
| `tools/convert_po_inputs.py` | builds the CSV inputs from the station tables of the Python pipeline |
| `python/` | Python pipeline: day-label harmonisation, quality control, normals, daily fields, evaluation, manuscript figures (see `python/README.md`) |
| `development/` | development record: equivalence between the Java and Python implementations, uncertainty notes, comparison tools. These files refer to internal folder names of the development runs. |
| `data_example/README.md` | input files expected in `data/` (station data are not distributed) |

## Build and run

See `RUN.md`. The jars needed to build and run are bundled in `lib/` (see `lib/README.md`); only Java 11 and the
input data are needed. `build.gradle` at the root lists the netcdf-java dependency as well; the tested build is
`build.sh`.
