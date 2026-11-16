#!/bin/bash
# Java (Krigings working copy, pokriging.engine.Cli) on the reference subset: 10 dates of 2020, 8192 cells (4 production
# chunks) + station LOO, Python normals (reference_run/normals_rebuild/normals), 10 threads, CSV output (unquantised).
set -euo pipefail
B="$(cd "$(dirname "$0")/.." && pwd)"; K=$B/Krigings
J=${JAVA11:-/Library/Java/JavaVirtualMachines/temurin-11.jdk/Contents/Home/bin}
G=${GEOET_LIB:-$HOME/Documents/GEPTOP/GEOET_CLEAN_TASK/OMS_Project_GEOET_clean_modified/lib}
CP="$K/build/classes/java/main:$G/*:$K/lib/flanagan.jar"; T=${THREADS:-10}; NORM=${NORMALS:-$B/reference_run/normals_rebuild/normals}
OUT=${OUT:-$B/output/eq}
common="mode=both target=grid demFile=$B/input_po/dem_po_1km.asc gridCellsFile=$B/input_po/grid_cells.csv targetCellsFile=$B/tests/eq/target_cells.csv datesFile=$B/tests/eq/dates_reference.csv normalsDir=$NORM outputFormat=csv quantise=false threads=$T"
P="stationsFile=$B/input_po/stations_precip.csv dataFile=$B/input_po/data_precip.csv variable=precip"
TT="stationsFile=$B/input_po/stations_temp.csv dataFile=$B/input_po/data_temp.csv variable=temp"
run(){ name=$1; shift; echo "== $name"; /usr/bin/time -p "$J/java" -Xmx8g -cp "$CP" pokriging.engine.Cli $common "$@" outputDir=$OUT/$name 2>&1 | grep -E "done:|real|Exception|Error" ; }
for m in ${METHODS:-precip_ratio precip_direct precip_twopart precip_twopart_table temp_anomaly temp_anomaly_nolapse temp_direct}; do
 case $m in
  precip_ratio) run $m $P precipMethod=climatology ;;
  precip_direct) run $m $P precipMethod=direct ;;
  precip_twopart) run $m $P precipMethod=twopart calibration=auto ;;
  precip_twopart_table) run $m $P precipMethod=twopart calibration=table occurrenceThresholds=$B/reference/twopart_production_calibration/thresholds_cells.csv ;;
  temp_anomaly) run $m $TT tempMethod=anomaly temperatureTrend=local ;;
  temp_anomaly_nolapse) run $m $TT tempMethod=anomaly temperatureTrend=none ;;
  temp_direct) run $m $TT tempMethod=direct temperatureTrend=local ;;
 esac
done
