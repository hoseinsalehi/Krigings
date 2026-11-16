package pokriging.engine;

/** All settings of a PoKriging run (= the @In keys of the OMS component; defaults = locked Po values). */
public final class KrigingConfig {
    public String variable;                       // precip | temp
    public String mode = "interpolation";         // interpolation | loo | both
    public String target = "grid";                // grid | points
    public String stationsFile, dataFile, missingValue = "", crs = "EPSG:3035";
    public double tempMin = -50, tempMax = 50;
    public int threads = Math.max(1, Runtime.getRuntime().availableProcessors() - 2);
    public String outputDir;
    public String demFile, gridCellsFile, targetCellsFile, pointsFile;
    public pokriging.io.Grid targetGrid;          // set by PoKriging from inTargets (cells from a shapefile); replaces demFile
    public String periodStart = "1991-01-01", periodEnd = "2020-12-31", datesFile;
    public String normalsDir;
    public String precipMethod = "twopart";       // direct | climatology | twopart
    public String tempMethod = "anomaly";         // direct | anomaly
    public String temperatureTrend = "local";     // local | none | global
    public double lapseMinSpanM = 300, anomLapseMaxExtrapM = 300; public int lapseMinStations = 10, lapseFallbackN = 100;
    public int nVario = 50, nKrig = 7, nLags = 12, minPairsPerLag = 5;
    public String variogramModel = "spherical";
    public int psoParticles = 24, psoIterations = 35; public double psoInertia = 0.72, psoC1 = 1.45, psoC2 = 1.45;
    public String fitWeighting = "pair-count";
    public String rngMode = "numpy"; public int chunk = 2048;
    public boolean globalVariogramFallback = true;
    public double zeroTol = 0.0, minDailyNormalMm = 0.05, maxPrecipRatio = 50, occThreshold = 1.0, pFloor = 0.01;
    public String calibration = "auto";           // auto | monthly | pooled | table
    public int calibMinDays = 100;
    public String occurrenceThresholds;
    public boolean writeTwopartPass1 = false;
    public String uncertaintyTable;
    public boolean diagnostics = true;            // with outputVariables = default only
    public String outputVariables = "default";   // default | all | comma list of value/diagnostic variables (README)
    public int minActiveStations = -1;            // default nVario + 1
    public double wetThresholdStats = 0.1;
    public String outputFormat = "netcdf";        // netcdf | csv
    public boolean quantise = true;
    public String nccopy = "/opt/homebrew/bin/nccopy";   // NetCDF-3 -> NetCDF-4/zlib conversion; "" = keep NetCDF-3
    public String excludeStationsFile;
    public String logLevel = "INFO";
}
