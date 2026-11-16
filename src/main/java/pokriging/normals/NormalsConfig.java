package pokriging.normals;

/** Settings of a PoNormals run (= @In keys of pokriging.PoNormals); defaults = locked Po values (METHOD_SPEC §1-2). */
public final class NormalsConfig {
    public String variable;                                   // precip | temp
    public String stationsFile, dataFile, missingValue = "", demFile, gridCellsFile, outputDir, crs = "EPSG:3035";
    public pokriging.io.Grid targetGrid;                      // set by PoNormals from inTargets; replaces demFile
    public double tempMin = -50, tempMax = 50;
    public String normalsStart = "1991-01-01", normalsEnd = "2020-12-31";
    public double monthCompleteness = 0.80; public int minObsYears = 5, minTotalYears = 20;
    public double fillRadiusKm = 100; public int fillMaxCandidates = 30, fillMinOverlapMonths = 36, fillTopK = 3;
    public double fillMinR = Double.NaN;                      // NaN = 0.70 (precip) / 0.85 (temp)
    public int regK = -1; public double regDz0 = Double.NaN, regExtrapDz = Double.NaN, regExtrapLapse = Double.NaN;
    public double slopeMin = Double.NaN, slopeMax = Double.NaN; public String slopeMode = "";    // annual | monthly
    public double precipFloorMm = 0.1, facetWeight = 1.0, tpi0 = 1e9; public boolean useTpi = false;
    public boolean wetFrequency = true; public double occThreshold = 1.0;
    public int wfK = 15; public double wfDz0 = 1e9, wfExtrapDz = 0, wfSlopeMin = -1.0e-3, wfSlopeMax = 2.0e-3, wfPmin = 0.005, wfPmax = 0.995;
    public int threads = Math.max(1, Runtime.getRuntime().availableProcessors() - 2);
    public String nccopy = "/opt/homebrew/bin/nccopy";

    void resolveDefaults() {
        boolean p = "precip".equals(variable);
        if (Double.isNaN(fillMinR)) fillMinR = p ? 0.70 : 0.85;
        if (regK < 0) regK = p ? 15 : 25;
        if (Double.isNaN(regDz0)) regDz0 = p ? 500 : 300;
        if (Double.isNaN(regExtrapDz)) regExtrapDz = p ? 0 : 1e9;
        if (Double.isNaN(regExtrapLapse)) regExtrapLapse = p ? 0 : -0.0075;
        if (Double.isNaN(slopeMin)) slopeMin = p ? -0.3e-3 : -10.0e-3;
        if (Double.isNaN(slopeMax)) slopeMax = p ? 0.8e-3 : 3.0e-3;
        if (slopeMode == null || slopeMode.isEmpty()) slopeMode = p ? "annual" : "monthly";
    }
}
