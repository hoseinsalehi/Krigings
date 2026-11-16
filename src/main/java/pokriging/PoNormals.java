package pokriging;

import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;
import oms3.annotations.Unit;
import pokriging.normals.NormalsConfig;
import pokriging.normals.NormalsEngine;

/**
 * OMS3 component: monthly normals for PoKriging (station normals with gap filling, PRISM-type gridded normals, station
 * LOO normals, wet-day frequency normals). One-shot component. NaN numeric keys mean the variable-specific default.
 */
@Description("Monthly normals (station, gridded PRISM-type, LOO, wet-day frequency) for PoKriging")
@Author(name = "Po extension 2026 (port of the Python s1/s2/s1w steps)")
@Keywords("normals, PRISM, gap filling, wet-day frequency")
@Name("PoNormals")
@Status(Status.EXPERIMENTAL)
@License("General Public License Version 3 (GPLv3)")
public class PoNormals {
    @Description("precip | temp") @In public String variable;
    @Description("stations.csv: id,x,y,elevation[,tpi20,facet]") @In public String stationsFile;
    @Description("wide data.csv") @In public String dataFile;
    @Description("missing code (empty field always missing)") @In public String missingValue = "";
    @Description("ESRI ASCII DEM of the normals grid") @In public String demFile;
    @Description("optional grid_cells.csv (elev_filled, tpi20, facet)") @In public String gridCellsFile;
    @Description("output folder") @In public String outputDir;
    @Description("metric CRS") @In public String crs = "EPSG:3035";
    @Description("temperature validity, lower bound") @Unit("degC") @In public double tempMin = -50;
    @Description("temperature validity, upper bound") @Unit("degC") @In public double tempMax = 50;
    @Description("first day of the normal period") @In public String normalsStart = "1991-01-01";
    @Description("last day of the normal period") @In public String normalsEnd = "2020-12-31";
    @Description("valid-day share for a station-month") @In public double monthCompleteness = 0.80;
    @Description("observed years per calendar month") @In public int minObsYears = 5;
    @Description("observed + reconstructed years") @In public int minTotalYears = 20;
    @Description("gap-fill donor radius") @Unit("km") @In public double fillRadiusKm = 100;
    @Description("donors examined") @In public int fillMaxCandidates = 30;
    @Description("minimum common months") @In public int fillMinOverlapMonths = 36;
    @Description("donors used (weights r^2)") @In public int fillTopK = 3;
    @Description("minimum correlation (NaN = 0.70 precip / 0.85 temp)") @In public double fillMinR = Double.NaN;
    @Description("stations in the regression (-1 = 15 precip / 25 temp)") @In public int regK = -1;
    @Description("elevation-difference weight scale (NaN = 500 / 300)") @Unit("m") @In public double regDz0 = Double.NaN;
    @Description("cap of the target elevation above the highest of the K stations (NaN = 0 precip / 1e9 temp)") @Unit("m") @In public double regExtrapDz = Double.NaN;
    @Description("lapse used above the highest of the K stations (NaN = 0 precip / -0.0075 temp)") @Unit("degC m-1") @In public double regExtrapLapse = Double.NaN;
    @Description("slope lower bound (NaN = -0.3e-3 ln/m precip / -10e-3 degC/m temp)") @In public double slopeMin = Double.NaN;
    @Description("slope upper bound (NaN = 0.8e-3 / 3e-3)") @In public double slopeMax = Double.NaN;
    @Description("annual | monthly (empty = annual precip / monthly temp)") @In public String slopeMode = "";
    @Description("ln(max(M, floor))") @Unit("mm") @In public double precipFloorMm = 0.1;
    @Description("PRISM facet weight (1 = inactive)") @In public double facetWeight = 1.0;
    @Description("PRISM TPI weight scale (1e9 = inactive)") @Unit("m") @In public double tpi0 = 1e9;
    @Description("TPI predictor (temp)") @In public boolean useTpi = false;
    @Description("also wet-day frequency normals (precip, two-part)") @In public boolean wetFrequency = true;
    @Description("wet day threshold") @Unit("mm") @In public double occThreshold = 1.0;
    @Description("wet-frequency regression K") @In public int wfK = 15;
    @Description("wet-frequency elevation weight scale") @Unit("m") @In public double wfDz0 = 1e9;
    @Description("wet-frequency elevation cap") @Unit("m") @In public double wfExtrapDz = 0;
    @Description("wet-frequency slope lower bound (logit per m)") @In public double wfSlopeMin = -1.0e-3;
    @Description("wet-frequency slope upper bound (logit per m)") @In public double wfSlopeMax = 2.0e-3;
    @Description("clip before logit, lower") @In public double wfPmin = 0.005;
    @Description("clip before logit, upper") @In public double wfPmax = 0.995;
    @Description("threads") @In public int threads = Math.max(1, Runtime.getRuntime().availableProcessors() - 2);
    @Description("nccopy for NetCDF-4 conversion (empty = NetCDF-3)") @In public String nccopy = "/opt/homebrew/bin/nccopy";
    @Description("target cells (points = cell centres, or cell polygons) e.g. from OmsShapefileFeatureReader; when given, demFile and gridCellsFile are not used")
    @In
    public org.geotools.data.simple.SimpleFeatureCollection inTargets;

    @Description("elevation field of inTargets")
    @In
    public String fTargetZ = "elevation";

    @Description("optional 0/1 field of inTargets flagging cells whose elevation was gap-filled (empty = none)")
    @In
    public String fTargetElevFilled = "";

    @Description("write the target grid built from inTargets as <outputDir>/target_dem.asc")
    @In
    public boolean writeTargetDem = true;

    @Description("run summary") @Out public String outSummary;
    @Description("folder of the normals written (connect to PoKriging.normalsDir)") @Out public String outNormalsDir;

    @Execute
    public void execute() throws Exception {
        NormalsConfig c = new NormalsConfig();
        for (java.lang.reflect.Field f : NormalsConfig.class.getFields()) {
            java.lang.reflect.Field g;
            try { g = PoNormals.class.getField(f.getName()); } catch (NoSuchFieldException e) { continue; }
            f.set(c, g.get(this));
        }
        if (inTargets != null) {
            c.targetGrid = pokriging.io.Grid.fromFeatures(inTargets, fTargetZ, fTargetElevFilled);
            java.nio.file.Files.createDirectories(java.nio.file.Paths.get(outputDir));
            if (writeTargetDem) c.targetGrid.writeAscii(java.nio.file.Paths.get(outputDir, "target_dem.asc"));
            System.out.println("PoNormals: " + c.targetGrid.n + " target cells from inTargets (" + c.targetGrid.ny + " x " + c.targetGrid.nx + " grid)");
        }
        NormalsEngine e = new NormalsEngine(c);
        e.run();
        outSummary = e.info.toString();
        outNormalsDir = outputDir;
        System.out.println("PoNormals done: " + outSummary);
    }
}
