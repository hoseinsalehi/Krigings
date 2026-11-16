#!/bin/bash
# Run a .sim through the OMS3 command line exactly like the OMS console: java -cp <oms-all + groovy + Krigings build + lib>
# oms3.CLI -r <sim>, with oms_prj = this folder. Java 11. Build first with ./build.sh.
set -euo pipefail
PRJ="$(cd "$(dirname "$0")" && pwd)"
J=${JAVA11:-/Library/Java/JavaVirtualMachines/temurin-11.jdk/Contents/Home/bin}
SIM="${1:?usage: run_oms.sh simulation/<file>.sim}"
L=${GEOET_LIB:?set GEOET_LIB to the lib folder of an OMS/GEOframe project}
CP="$L/oms-all.jar:$L/groovy-all-2.3.9.jar:$PRJ/../build/libs/Krigings-0.9.9-po1.jar:$L/*:$PRJ/../lib/flanagan.jar"
exec "$J/java" -Xmx${XMX:-8g} -cp "$CP" -Doms.prj="$PRJ" -Doms.home="$PRJ" oms3.CLI -r "$SIM"
