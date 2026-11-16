# Files to copy into the package after the 1991-2020 run has finished (PoUncertainty, 2026-10-19)

DEV = ~/Documents/Po_precip/GEOframe_PoKriging, PKG = ~/Documents/Po_precip/GEOframe_PoKriging_share/PoKriging_GEOframe_0.9.9-po1.
Do NOT run PKG/build.sh while a run uses PKG/lib/Krigings-0.9.9-po1.jar (build.sh replaces that jar).

| DEV | PKG | change |
|---|---|---|
| Krigings/src/main/java/pokriging/uncertainty/UncertaintyCalibration.java | same path | new |
| Krigings/src/main/java/pokriging/PoUncertainty.java | same path | new OMS component |
| Krigings/src/main/java/pokriging/PoKriging.java | same path | + `@Out outLooDir` (3 lines; generator tail, regenerating is stable) |
| simulation/precip.sim, simulation/temp.sim | same paths | + USER CONFIGURATION keys makeUncertainty / uncertaintyTable / uncNominal, component `unc`, connect `krig.outLooDir -> unc.looDir`, `outputsAll` (adds `<var>_uncertainty` to an explicit output list when a table is given) |
| input_po/uncertainty_precip_v5.csv, input_po/uncertainty_temp_v5.csv | data/ | new: product k tables (= Python s9 v5/v4, made by PoUncertainty) |
| docs/UNCERTAINTY_CALIBRATION.md | docs/ | new |
| README.md section "Calibrated uncertainty: pokriging.PoUncertainty" + the 3 edited lines (Layout, Outputs, Remaining gaps) | merge into PKG/README.md | text |
| EQUIVALENCE_REPORT.md section "Uncertainty calibration" | append to PKG/docs/EQUIVALENCE_REPORT.md | text |
| tools/export_loo_parquet_csv.py, tools/compare_uncertainty.py | tools/ (optional) | Python test helpers |

Not for the package: simulation/tests/test_uncertainty_chain_{precip,temp}.sim (they use reference_run/ normals and tests/eq dates).
Then: `cd PKG && ./build.sh` (Java 11; commons-math3 is already in PKG/lib) and re-zip.
