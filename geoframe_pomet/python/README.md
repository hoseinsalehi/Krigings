# GEOframe-PoMet – Python pipeline

The Python implementation of GEOframe-PoMet. It contains the preprocessing and quality control, the normals, the
daily fields, the evaluation, and the figures and tables of the manuscript. The GEOframe components in
`../../src/main/java/pokriging` reproduce the interpolation (see `../development/EQUIVALENCE_REPORT.md`).

Python ≥ 3.10 with numpy, scipy, pandas, geopandas, pyproj, xarray, netCDF4, pyarrow and numba. Numba is optional:
without it, the numpy kernels are used.

## Configuration

All paths are read from environment variables in `po_config.py`. The defaults are relative to the working folder.

| variable | content |
|---|---|
| `PO_MERGED_GAUGES` | merged regional station archive |
| `PO_STATION_INPUT` | station tables including the MeteoSwiss stations (`s00_add_meteoswiss.py`) |
| `PO_GRID_MASK_NC`, `PO_DEM_POINTS` | 1-km grid mask and cell elevations |
| `PO_INPUT_ROOT` | external inputs (MeteoSwiss, District boundary, reference datasets) |
| `PO_OUT_ROOT`, `PO_DATA_ROOT` | output folder |
| `PO_STATION_SET`, `PO_PRECIP_METHOD` | network and precipitation method; the defaults (`v4`, `twopart`) give GEOframe-PoMet |

## Steps

The `run_all.sh` script runs the production chain:

| script | step |
|---|---|
| `s00_add_meteoswiss.py` | add the MeteoSwiss stations |
| `s0_day_alignment.py` | harmonise the day labels of all records |
| `s1_station_normals.py`, `s1w_wetfreq_normals.py`, `s2_normals_grid.py` | station and gridded monthly normals |
| `s3_run.py` | daily 1-km fields with leave-one-out |
| `s7_station_qc.py` | automatic station-day quality control (uses the LOO of a first run) |
| `s9_uncertainty_calibration.py` | calibration of the uncertainty |
| `s4_summarise.py`, `s10_final_check.py`, `s13_zenodo_package.py` | summaries, checks, dataset files |

These scripts produce the evaluation:

| script | evaluation |
|---|---|
| `s5_thinning.py`, `s5b_missing_days.py`, `s5c_missing_days_long.py` | missing-data experiments |
| `s6_point_test.py` | Aosta–GSB point test |
| `s8_reference_comparison.py` | comparison with the reference datasets |
| `s11_gsb_independent_check.py` | check at Col du Grand St-Bernard |
| `s12_isac_normals_comparison.py` | comparison with the reference climatologies |
| `s14_paper_tables.py`, `s15_freq_intensity.py` | tables and frequency/intensity analysis |

`figures/` contains the code of the manuscript figures. Run it with `PYTHONPATH=figures python figures/<script>.py`.

Some file names, environment values and comments still use development labels (for example `PO_STATION_SET=v4`);
they identify configurations of the development runs.
