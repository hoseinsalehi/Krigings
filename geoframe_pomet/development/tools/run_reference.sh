#!/usr/bin/env bash
# Python reference run for the GEOframe port. Rebuilds the normals (bit-identical to PoGrid_v5/normals, checked)
# in a private root and runs the daily engine for 10 dates of 2020 on 4 production chunks (8192 cells) + LOO.
# Nothing is written to PoGrid_v3_code / PoGrid_v4 / PoGrid_v5 (bytecode and numba caches are redirected).
set -euo pipefail
B=$HOME/Documents/Po_precip/GEOframe_PoKriging
R=$B/reference_run/normals_rebuild          # PO_OUT_ROOT: normals/ rebuilt here, grid_predictors.npz -> symlink to v5
OUT=$B/reference
W=${WORKERS:-10}
export PO_STATION_SET=v4 PO_OUT_ROOT=$R PO_DATA_ROOT=$HOME/Downloads/PoGrid_v5 OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES
export PYTHONDONTWRITEBYTECODE=1 NUMBA_CACHE_DIR=$B/reference_run/.numba_cache
cd "$(dirname "$0")"
run(){ m=$1; s=$2; pm=$3; d=$OUT/daily/$m; [ "$s" -gt 0 ] && d=$OUT/seed_sensitivity/${m}_salt$s
  PO_PRECIP_METHOD=$pm python run_python_reference.py --method $m --salt $s --workers $W --out $d; }
for m in precip_ratio precip_direct temp_anomaly temp_anomaly_nolapse temp_direct; do run $m 0 ratio; done
run precip_twopart 0 twopart
for s in 1 2; do run precip_ratio $s ratio; run temp_anomaly $s ratio; run precip_twopart $s twopart; done
