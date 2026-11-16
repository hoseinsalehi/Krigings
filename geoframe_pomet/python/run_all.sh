#!/usr/bin/env bash
# Full rebuild of PoGrid v3 after day-labelling harmonisation (s0 already run; tables in OUT_ROOT/harmonisation).
# Station QC is iterated on fast station-only LOO runs (about 1 min each) BEFORE the gridded production, so the
# production is run only once. Safe to restart: every step either skips finished work or is cheap.
set -euo pipefail
cd "$(dirname "$0")"
OUT=${PO_OUT_ROOT:-${PO_DATA_ROOT:-data}}
P360=${PO_INPUT_ROOT:-data}/variogram_precip_360days_v2/selected_days.csv
T360=${PO_INPUT_ROOT:-data}/variogram_temp_360days_v2/selected_days.csv
W=${WORKERS:-12}
log(){ echo "[$(date +%H:%M:%S)] $*"; }
export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES     # macOS: allow fork() in worker pools after Objective-C init

mkdir -p "$OUT"
[ -f "$OUT/grid_predictors.npz" ] || [ ! -f "${PO_DATA_ROOT:-data}/grid_predictors.npz" ] || cp "${PO_DATA_ROOT:-data}/grid_predictors.npz" "$OUT/"
for v in precip temp; do
  if [ -f "$OUT/harmonisation/day_shift_$v.csv" ]; then log "$v: day alignment already done"; continue; fi
  log "$v: day-labelling alignment (s0)"
  PO_SHIFT=0 PO_QC=0 python s0_day_alignment.py --variable $v
done

for v in precip temp; do
  if [ -f "$OUT/qc/.done_$v" ]; then log "$v: QC already final"; continue; fi
  for pass in 1 2 3; do
    log "$v QC pass $pass: normals"
    python s1_station_normals.py --variable $v
    if [ -f "$OUT/normals/params_$v.json" ]; then python s2_normals_grid.py --variable $v --params "$OUT/normals/params_$v.json"
    else python s2_normals_grid.py --variable $v --tune; fi
    log "$v QC pass $pass: station LOO, full period"
    rm -rf "$OUT/qcloo_${v}_p$pass"
    python s3_run.py --variable $v --mode anomaly --loo-only --tag qcloo_${v}_p$pass --workers $W
    log "$v QC pass $pass: tests"
    if [ $pass -eq 1 ]; then python s7_station_qc.py --variable $v --run qcloo_${v}_p$pass --fresh | tee "$OUT/qcloo_${v}_p$pass/qc.log"
    else python s7_station_qc.py --variable $v --run qcloo_${v}_p$pass | tee "$OUT/qcloo_${v}_p$pass/qc.log"; fi
    new=$(grep -o "new [0-9]*" "$OUT/qcloo_${v}_p$pass/qc.log" | awk '{print $2}')
    log "$v QC pass $pass: new flags = ${new:-all}"
    if [ $pass -gt 1 ] && [ "${new:-0}" -lt 100 ]; then break; fi
  done
  touch "$OUT/qc/.done_$v"
done

for v in precip temp; do
  [ -f "$OUT/normals/.final_$v" ] && { log "$v: final normals already built"; continue; }
  log "$v: final normals (tuned)"
  python s1_station_normals.py --variable $v
  python s2_normals_grid.py --variable $v --tune
  touch "$OUT/normals/.final_$v"
done

log "360-day method comparison"
for v in precip temp; do
  D=$P360; [ $v = temp ] && D=$T360
  for m in direct anomaly; do
    [ -f "$OUT/loo360_${v}_$m/daily_summary_${v}_2020.csv" ] || python s3_run.py --variable $v --mode $m --loo-only --dates-csv $D --tag loo360_${v}_$m --workers $W
  done
done

log "production 1991-2020"
python s3_run.py --variable precip --mode anomaly --workers $W
python s3_run.py --variable temp   --mode anomaly --workers $W

log "summaries"
if [ ! -f "$OUT/summary/cell_climatology_temp_anomaly.nc" ]; then
  python s4_summarise.py --runs loo360_precip_direct loo360_precip_anomaly loo360_temp_direct loo360_temp_anomaly
  python s4_summarise.py --runs precip_anomaly temp_anomaly --climatology
fi

log "uncertainty calibration"
for v in precip temp; do [ -f "$OUT/summary/uncertainty_calibration_$v.csv" ] || python s9_uncertainty_calibration.py --variable $v; done

log "reference datasets"
for v in precip temp; do [ -f "$OUT/reference_comparison/summary_$v.csv" ] || python s8_reference_comparison.py --variable $v; done

log "Aosta - Gran San Bernardo point test"
for v in precip temp; do [ -f "$OUT/point_test_$v/summary.csv" ] || python s6_point_test.py --variable $v --workers $W; done

log "50 % network thinning, 30 realizations"
[ -f "$OUT/thinning_precip/summary_withheld.csv" ] || python s5_thinning.py --variable precip --dates-csv $P360 --realizations 30 --workers $W
[ -f "$OUT/thinning_temp/summary_withheld.csv" ]   || python s5_thinning.py --variable temp   --dates-csv $T360 --realizations 30 --workers $W

log "independent check at MeteoSwiss Col du Grand St-Bernard"
[ -f "$OUT/gsb_check/gsb_metrics.csv" ] || python s11_gsb_independent_check.py

log "ISAC-CNR normals comparison (Crespi 2018 / Brunetti 2014)"
[ -f "$OUT/isac_comparison/summary.csv" ] || python s12_isac_normals_comparison.py

log "final check"
python s10_final_check.py
log "DONE"
