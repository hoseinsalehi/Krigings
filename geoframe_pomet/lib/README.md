# Bundled runtime libraries

These are the third-party jars that `../build.sh` and `../run_oms.sh` need: OMS3 (`oms-all.jar`), Groovy 2.3.9,
GeoTools 25, HortonMachine 0.10.4, JAI/jai-ext, netcdf-java 5.5.3 and their dependencies. They were taken from the lib
folder of a GEOframe OMS project and reduced to the jars used to compile the component and to run the sims of
`../simulation` (grid, point and leave-one-out modes). Every jar keeps its own licence. To use a different set, point
`GEOET_LIB` to another folder.
