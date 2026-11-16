#!/usr/bin/env bash
# Two-part precipitation rebuild (occurrence x intensity), 2026-10-12.
# Temperature, station data (day-shift tables, QC) and the precipitation total normals M are taken unchanged from
# the v2.0 run in ~/Downloads/PoGrid_v4 (read-only: linked or copied, never written). Everything precipitation is
# recomputed into a NEW root ~/Downloads/PoGrid_v5, so the v2.0 results stay intact for the comparison.
# Safe to restart: finished years/steps are skipped.     Usage:  ./run_twopart.sh  2>&1 | tee ~/Downloads/run_twopart.log
set -euo pipefail
cd "$(dirname "$0")"
V4=${PO_INPUT_ROOT:-data}/PoGrid_v4
OUT=${PO_OUT_ROOT:-${PO_INPUT_ROOT:-data}/PoGrid_v5}
P360=${PO_INPUT_ROOT:-data}/variogram_precip_360days_v2/selected_days.csv
W=${WORKERS:-12}
log(){ echo "[$(date +%H:%M:%S)] $*"; }
export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES PO_STATION_SET=v4 PO_PRECIP_METHOD=twopart PO_OUT_ROOT=$OUT
unset PO_DATA_ROOT

free=$(df -g "$HOME/Downloads" | awk 'NR==2{print $4}')
[ "$free" -ge 20 ] || { echo "only ${free} GB free; at least 20 GB are needed (empty the Trash)"; exit 1; }

log "setup $OUT"
mkdir -p "$OUT/normals" "$OUT/summary" "$OUT/reference_comparison"
[ -f "$OUT/grid_predictors.npz" ] || cp "$V4/grid_predictors.npz" "$OUT/"
for d in harmonisation qc temp_anomaly loo360_temp_direct loo360_temp_anomaly loo360_precip_direct point_test_temp thinning_temp; do
  [ -e "$OUT/$d" ] || ln -s "$V4/$d" "$OUT/$d"          # read-only inputs from v2.0
done
[ -f "$OUT/normals/.copied" ] || { cp "$V4"/normals/*.* "$OUT/normals/"; cp "$V4"/normals/.final_* "$OUT/normals/" 2>/dev/null || true; touch "$OUT/normals/.copied"; }
[ -f "$OUT/summary/.copied" ] || { cp "$V4"/summary/*temp* "$OUT/summary/"; cp "$V4"/reference_comparison/*temp* "$OUT/reference_comparison/"; touch "$OUT/summary/.copied"; }

log "wet-day frequency normals (tuned by LOO)"
[ -f "$OUT/normals/wetfreq_precip.nc" ] || python s1w_wetfreq_normals.py --tune | tee "$OUT/normals/wetfreq_precip.log"

log "360-day leave-one-out (two-part)"
[ -f "$OUT/loo360_precip_anomaly/daily_summary_precip_2020.csv" ] || python s3_run.py --variable precip --mode anomaly --loo-only --dates-csv "$P360" --tag loo360_precip_anomaly --workers $W

log "production 1991-2020 (two-part)"
python s3_run.py --variable precip --mode anomaly --workers $W

log "summaries"
python s4_summarise.py --runs loo360_precip_direct loo360_precip_anomaly
python s4_summarise.py --runs precip_anomaly --climatology

log "uncertainty calibration"
python s9_uncertainty_calibration.py --variable precip

log "reference datasets"
python s8_reference_comparison.py --variable precip

log "Aosta - Gran San Bernardo point test"
[ -f "$OUT/point_test_precip/summary.csv" ] || python s6_point_test.py --variable precip --workers $W

log "50 % network thinning, 30 realizations"
[ -f "$OUT/thinning_precip/summary_withheld.csv" ] || python s5_thinning.py --variable precip --dates-csv "$P360" --realizations 30 --workers $W

log "independent check at MeteoSwiss Col du Grand St-Bernard"
python s11_gsb_independent_check.py

log "ISAC-CNR normals comparison"
python s12_isac_normals_comparison.py

log "frequency / intensity test (stations, v1.0, v2.0 ratio method, two-part)"
python s15_freq_intensity.py

log "final check"
python s10_final_check.py || true
du -sh "$OUT"; df -h "$HOME/Downloads" | tail -1
log "DONE"
