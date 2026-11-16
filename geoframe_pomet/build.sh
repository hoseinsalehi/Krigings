#!/bin/bash
# Build the Krigings working copy (GEOframe Krigings 0.9.9 + krigingsCore + pokriging) with Java 11 and the jars of the
# GEOET OMS project (hm-* 0.10.4, geotools, commons-math3, netcdf-java, oms-all). Output: Krigings/build/classes/java/main
# and Krigings/build/libs/Krigings-0.9.9-po1.jar.
set -euo pipefail
B_LIB="$(cd "$(dirname "$0")" && pwd)/lib"      # jars bundled with this project (OMS, GeoTools, HortonMachine, netcdf-java)
B="$(cd "$(dirname "$0")" && pwd)"; K=$(cd "$B/.." && pwd)
J=${JAVA11:-/Library/Java/JavaVirtualMachines/temurin-11.jdk/Contents/Home/bin}
[ -x "$J/java" ] || J="$(dirname "$(command -v java)")"   # fall back to the java on PATH (must be Java 11)
G=${GEOET_LIB:-$B_LIB}
CP="${KRIG_CP:-$G/*}:$K/lib/flanagan.jar"
OUT=$K/build/classes/java/main; rm -rf "$OUT"; mkdir -p "$OUT" "$K/build/libs"
"$J/javac" -nowarn -encoding ISO-8859-1 -cp "$CP" -d "$OUT" $(find "$K/src/main/java" -name '*.java')
(cd "$OUT" && "$J/jar" cf "$K/build/libs/Krigings-0.9.9-po1.jar" .)
echo "compiled $(find "$OUT" -name '*.class' | wc -l | tr -d ' ') classes -> $K/build/libs/Krigings-0.9.9-po1.jar"
