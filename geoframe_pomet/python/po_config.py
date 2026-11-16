# -*- coding: utf-8 -*-
"""Paths and fixed method parameters for the GEOframe-PoMet Python pipeline (Po River District, 1 km, daily precipitation and temperature).

Edit only the PATHS block when moving machines. Method parameters are the values
selected by the LOO tests of 2026-10-07 (50 variogram / 7 kriging stations, fixed
spherical model, PSO 24x35, per-cell daily lapse for temperature).
"""
import os
from pathlib import Path

# ----------------------------------------------------------------------------- PATHS
MG = Path(os.environ.get("PO_MERGED_GAUGES", "data/Merged_Gauges"))
DATA = {
    "precip": dict(obs=MG / "Precip_New/Final_precip.csv",
                   stations=MG / "Final_Data/STATIONS_FINAL_CLEANED_V3.shp",
                   date_format="%d/%m/%y"),
    "temp":   dict(obs=MG / "Final_Data/Final_Merged_Temp_Timeseries_Cleaned_V2_Filtered.csv",
                   stations=MG / "Final_Data/Temp_new_all_Filtered.shp",
                   date_format="%d/%m/%y %H:%M"),
}
# PO_STATION_SET=v4 (default here, = GEOframe-PoMet network): network extended with the official MeteoSwiss SMN/NIME stations (s00_add_meteoswiss.py)
V4_INPUT = Path(os.environ.get("PO_STATION_INPUT", "data/station_input"))
if os.environ.get("PO_STATION_SET", "v4") == "v4":
    DATA["precip"].update(obs=V4_INPUT / "Final_precip_plusCH.csv", stations=V4_INPUT / "STATIONS_plusCH_precip.shp")
    DATA["temp"].update(obs=V4_INPUT / "Final_temp_plusCH.csv", stations=V4_INPUT / "STATIONS_plusCH_temp.shp")
GRID_MASK_NC = Path(os.environ.get("PO_GRID_MASK_NC", "data/grid_mask.nc"))   # 494 x 330, 88,074 valid cells
DEM_POINTS = Path(os.environ.get("PO_DEM_POINTS", "data/points_with_elevation.shp"))          # same grid, 88,028 cells
OUT_ROOT = Path(os.environ.get("PO_OUT_ROOT", "output"))
# station-data products (day-shift tables, QC flags, cleaned CSVs) live under DATA_ROOT, which stays the main
# output root when a sub-run (thinning realization) redirects OUT_ROOT
DATA_ROOT = Path(os.environ.get("PO_DATA_ROOT", str(OUT_ROOT)))
# day-labelling harmonisation (s0_day_alignment.py): per station-year shift in days, applied when reading the raw
# CSV (PO_SHIFT=0 disables); the cleaned *_QC.csv written by s7 is already harmonised
SHIFT_FILE = "harmonisation/day_shift_{variable}.csv"     # under OUT_ROOT
# raw observations are kept untouched; s7_station_qc.py writes cleaned copies to OUT_ROOT/qc/<name>_QC.csv,
# which every step then reads automatically (PO_RAW=1 forces the raw files)
for _v in DATA:
    DATA[_v]["raw"] = DATA[_v]["obs"]
    _q = DATA_ROOT / "qc" / (DATA[_v]["obs"].stem + "_QC.csv")
    if _q.exists() and os.environ.get("PO_RAW", "0") != "1": DATA[_v]["obs"] = _q

STATION_ID_FIELD = "ID"
STATION_ELEV_FIELD = "elevation"
METRIC_CRS = "EPSG:3035"
MISSING_VALUE = -9999.0
TEMP_RANGE = (-50.0, 50.0)
# station-day QC flags (written by s7_station_qc.py); applied in every step when the file exists, PO_QC=0 disables
QC_FLAGS_FILE = "station_day_flags_{variable}.parquet"     # under OUT_ROOT/qc/
QC_TEMP_ABS_C = 8.0         # |LOO error - station-month median LOO error| above which a temperature station-day is flagged
QC_CONTAM_KM = 30.0         # opposite-sign smaller errors within this distance of a flagged station are not flagged
# temperature, anomaly mode: the daily lapse on anomalies is not extrapolated more than this above/below the
# elevation range of the 50 local stations (above it the monthly normal alone carries the elevation signal)
ANOM_LAPSE_MAX_EXTRAP_M = 300.0
PERIOD = ("1991-01-01", "2020-12-31")

# ----------------------------------------------------------------------------- METHOD (locked by LOO)
N_VARIO = 50              # stations defining each local variogram (and the local lapse)
N_KRIG = 7                # kriging neighbours
N_LAGS = 12
MIN_PAIRS_PER_LAG = 5
PSO_PARTICLES = 24
PSO_ITERATIONS = 35
PSO_INERTIA = 0.72
PSO_C1 = 1.45
PSO_C2 = 1.45
FIT_WEIGHTING = "pair-count"
LAPSE_MIN_SPAN_M = 300.0  # local lapse falls back to LAPSE_FALLBACK_N stations below this elevation span
LAPSE_MIN_STATIONS = 10
LAPSE_FALLBACK_N = 100
WET_THRESHOLD = 0.1       # mm/day, for wet-day statistics only
ZERO_TOL = 0.0            # precipitation: all N_KRIG neighbours <= ZERO_TOL -> estimate 0
MIN_DAILY_NORMAL_MM = 0.05  # anomaly mode: floor on the daily climatological mean used as divisor
MAX_PRECIP_RATIO = 50.0   # anomaly mode: cap on station ratio P/Pclim (guards near-zero normals)
OSC_SUM_ABS_W = 2.0       # qc flag: oscillatory weights
CHUNK = 2048              # targets per vectorised chunk

# ----------------------------------------------------------------------------- NORMALS (step 1/2)
MONTH_COMPLETENESS = 0.80     # fraction of valid days required for a station-month value
MIN_OBS_YEARS = int(os.environ.get("PO_MIN_OBS_YEARS", 5))             # observed years required per calendar month for a station normal
MIN_TOTAL_YEARS = 20          # observed + reconstructed years required
FILL_RADIUS_KM = 100.0
FILL_MAX_CANDIDATES = 30
FILL_MIN_OVERLAP_MONTHS = int(os.environ.get("PO_FILL_MIN_OVERLAP", 36))   # 2026-10-10: 8 yr / 60 mo -> 5 yr / 36 mo (+173 precip, +131 temp normals; LOO neutral/better)
FILL_MIN_R = {"precip": 0.70, "temp": 0.85}
FILL_TOP_K = 3

# Diagnostic boxes (EPSG:4326, not administrative units) used in summaries
BOXES = {"Aosta_box": (6.75, 7.95, 45.42, 46.02), "Valtellina_box": (9.15, 10.55, 45.90, 46.65)}

# ----------------------------------------------------------------------------- TWO-PART PRECIPITATION (2026-10-12)
# PO_PRECIP_METHOD=twopart: precipitation '--mode anomaly' runs the two-part method (s3_run.py): daily occurrence from
# kriged occurrence anomalies (wet - p_station) calibrated to the wet-day frequency normal p(x) (s1w_wetfreq_normals.py),
# amounts from kriged intensity ratios P / i_station on wet stations only, i = M / (p * days). Default here: 'twopart' (GEOframe-PoMet).
PRECIP_METHOD = os.environ.get("PO_PRECIP_METHOD", "twopart")
OCC_THRESHOLD = 1.0           # mm/day defining a wet day for occurrence and frequency normals (ETCCDI R1mm)
P_FLOOR = 0.01                # floor on wet-day frequency when deriving the intensity normal i = M / (p * days)
CALIB_MIN_DAYS = 100          # occurrence calibration per cell and calendar month needs >= this many days, else pooled
