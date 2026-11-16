#!/bin/bash
# Run a .sim through the OMS3 command line exactly like the OMS console: java -cp <oms-all + groovy + Krigings build + lib>
# oms3.CLI -r <sim>, with oms_prj = this folder. Java 11. Build first with ./build.sh.
set -euo pipefail
B_LIB="$(cd "$(dirname "$0")" && pwd)/lib"      # jars bundled with this project (OMS, GeoTools, HortonMachine, netcdf-java)
PRJ="$(cd "$(dirname "$0")" && pwd)"
J=${JAVA11:-/Library/Java/JavaVirtualMachines/temurin-11.jdk/Contents/Home/bin}
[ -x "$J/java" ] || J="$(dirname "$(command -v java)")"   # fall back to the java on PATH (must be Java 11)
SIM="${1:?usage: run_oms.sh simulation/<file>.sim}"
L=${GEOET_LIB:-$B_LIB}
CP="$L/oms-all.jar:$L/groovy-all-2.3.9.jar:$PRJ/../build/libs/Krigings-0.9.9-po1.jar:$L/*:$PRJ/../lib/flanagan.jar"
exec "$J/java" -Xmx${XMX:-8g} -cp "$CP" -Doms.prj="$PRJ" -Doms.home="$PRJ" oms3.CLI -r "$SIM"
