/*
 * GNU GPL v3 License (GEOframe Krigings, version 0.9.9 of 2023-01-11).
 */
/* Source reconstructed from Kriging-2023-01-11-0.9.9.jar (CFR 0.152 decompilation, cleaned 2026; see UPSTREAM_DIFF.md). */
package krigingsPointCase;

import experimentalVariogram.ExperimentalVariogram;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import krigingsPointCase.ResidualsEvaluator;
import krigingsPointCase.StationsSelection;
import krigingsPointCase.Utility;
import krigingsPointCase.VariogramParamsEvaluator;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.hortonmachine.gears.libs.monitor.IHMProgressMonitor;
import org.hortonmachine.gears.libs.monitor.LogProgressMonitor;
import org.opengis.feature.simple.SimpleFeatureType;

public class GlobalParameterEvaluator {
    @Description(value="The .shp of the measurement point, containing the position of the stations.")
    @In
    public SimpleFeatureCollection inStations = null;
    @Description(value="The field of the vector of stations, defining the id.")
    @In
    public String fStationsid = null;
    @Description(value="The field of the vector of stations, defining the elevation.")
    @In
    public String fStationsZ = null;
    @Description(value="The progress monitor.")
    @In
    public IHMProgressMonitor pm = new LogProgressMonitor();
    @Description(value="Include zeros in computations (default is true).")
    @In
    public boolean doIncludeZero = true;
    @Description(value="The type of theoretical semivariogram: exponential, gaussian, spherical, pentasphericallinear, circular, bessel, periodic, hole, logaritmic, power, spline")
    @In
    public String pSemivariogramType = null;
    @Description(value="Switch for detrended mode.")
    @In
    public boolean doDetrended;
    @Description(value="Specified cutoff")
    @In
    public double cutoffInput;
    @Description(value="Number of bins to consider in the anlysis")
    @In
    public int cutoffDivide;
    @Description(value="Distances input path")
    @In
    public String inHValuesPath;
    @In
    public String tStart = null;
    @In
    public String tEnd = null;
    @In
    public int tTimeStep = 60;
    @In
    public boolean doLogarithmic = false;
    @Out
    public double pAGlobal;
    @Out
    public double pSGlogal;
    @Out
    public double pNugGlogal;
    @Out
    public double pAGlobalDeTrended;
    @Out
    public double pSGlogalDeTrended;
    @Out
    public double pNugGlogalDeTrended;
    @Out
    public StationsSelection stations;
    @Description(value="In the case of kriging with neighbor, maxdist is the maximum distance within the algorithm has to consider the stations")
    @In
    public double maxdist;
    @In
    public String fileNoValue = "-9999";
    public String outGlobalVariogramType;
    public String outGlobalDetrendedVariogramType;

    @Execute
    public void execute() {
        this.verifyInput();
        this.stations = new StationsSelection();
        this.stations.inStations = this.inStations;
        this.stations.doIncludezero = this.doIncludeZero;
        this.stations.maxdist = this.maxdist;
        this.stations.fStationsid = this.fStationsid;
        this.stations.fStationsZ = this.fStationsZ;
        this.stations.doLogarithmic = this.doLogarithmic;
        try {
            this.createDefaulParams(this.stations);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void verifyInput() {
        int ff;
        if (this.doDetrended && (ff = ((SimpleFeatureType)this.inStations.getSchema()).indexOf(this.fStationsZ)) < 0) {
            throw new NullPointerException("check if the z field name is correct");
        }
    }

    private void createDefaulParams(StationsSelection stations) throws Exception {
        System.out.println("evaluate global params");
        double[] variance = new double[this.cutoffDivide];
        double[] distance = new double[this.cutoffDivide];
        double[] varianceDeTrended = new double[this.cutoffDivide];
        double[] distanceDeTrended = new double[this.cutoffDivide];
        OmsTimeSeriesIteratorReader readH = new OmsTimeSeriesIteratorReader();
        readH.file = this.inHValuesPath;
        readH.idfield = "ID";
        readH.fileNovalue = this.fileNoValue;
        readH.tStart = this.tStart;
        readH.tTimestep = this.tTimeStep;
        if (this.tEnd != null) {
            readH.tEnd = this.tEnd;
        }
        readH.initProcess();
        ExperimentalVariogram exp = Utility.getExperimentalVariogram(this.fStationsid, this.inStations, this.doIncludeZero, this.cutoffDivide, this.cutoffInput, 0);
        int nRows = 0;
        int nRowsDeTrended = 0;
        try {
            while (readH.doProcess) {
                HashMap h;
                readH.nextRecord();
                exp.inData = h = readH.outData;
                exp.process();
                if (exp.areAllEquals || exp.differents <= 2) continue;
                HashMap<Integer, double[]> d = exp.outDistances;
                HashMap<Integer, double[]> v = exp.outExperimentalVariogram;
                int j = 0;
                for (Map.Entry<Integer, double[]> tt : d.entrySet()) {
                    distance[j] = distance[j] + tt.getValue()[0];
                    variance[j] = variance[j] + v.get(tt.getKey())[0];
                    ++j;
                }
                if (this.doDetrended) {
                    stations.inData = h;
                    stations.execute();
                    int n1 = stations.hStationInitialSet.length - 1;
                    int[] idStations = Arrays.copyOfRange(stations.idStationInitialSet, 0, n1);
                    ResidualsEvaluator rEvaluator = new ResidualsEvaluator();
                    rEvaluator.doDetrended = true;
                    rEvaluator.hStations = Arrays.copyOfRange(stations.hStationInitialSet, 0, n1);
                    rEvaluator.zStations = Arrays.copyOfRange(stations.zStationInitialSet, 0, n1);
                    rEvaluator.process();
                    if (rEvaluator.isPValueOk) {
                        double[] hVal = rEvaluator.hResiduals;
                        ++nRowsDeTrended;
                        HashMap<Integer, double[]> hRes = new HashMap<Integer, double[]>();
                        for (int i = 0; i < n1; ++i) {
                            hRes.put(stations.idStationInitialSet[i], new double[]{hVal[i]});
                        }
                        exp.inData = hRes;
                        exp.process();
                        if (exp.differents > 2) {
                            d = exp.outDistances;
                            v = exp.outExperimentalVariogram;
                            j = 0;
                            for (Map.Entry<Integer, double[]> tt : d.entrySet()) {
                                distanceDeTrended[j] = distanceDeTrended[j] + tt.getValue()[0];
                                varianceDeTrended[j] = varianceDeTrended[j] + v.get(tt.getKey())[0];
                                ++j;
                            }
                        }
                    }
                }
                ++nRows;
            }
            readH.close();
            if (distance.length != variance.length) {
                throw new IllegalArgumentException();
            }
            for (int i = 0; i < distance.length; ++i) {
                distance[i] = distance[i] / (double)nRows;
                variance[i] = variance[i] / (double)nRows;
            }
            VariogramParamsEvaluator vEvaluator = new VariogramParamsEvaluator();
            vEvaluator.pSemivariogramType = this.pSemivariogramType;
            vEvaluator.x = distance;
            vEvaluator.y = variance;
            vEvaluator.proces();
            this.pNugGlogal = vEvaluator.nugget;
            this.pSGlogal = vEvaluator.sill;
            this.pAGlobal = vEvaluator.range;
            this.outGlobalVariogramType = vEvaluator.outSemivariogramType;
            this.pm.message("Global value for nugget: " + this.pNugGlogal + " sill:" + this.pSGlogal + " range: " + this.pAGlobal + "  semivariogram type:" + vEvaluator.outSemivariogramType);
            if (this.doDetrended) {
                if (distanceDeTrended.length != varianceDeTrended.length) {
                    throw new IllegalArgumentException();
                }
                if (nRowsDeTrended > 0) {
                    for (int i = 0; i < distanceDeTrended.length; ++i) {
                        distanceDeTrended[i] = distanceDeTrended[i] / (double)nRowsDeTrended;
                        varianceDeTrended[i] = varianceDeTrended[i] / (double)nRowsDeTrended;
                    }
                    vEvaluator = new VariogramParamsEvaluator();
                    vEvaluator.pSemivariogramType = this.pSemivariogramType;
                    vEvaluator.x = distanceDeTrended;
                    vEvaluator.y = varianceDeTrended;
                    vEvaluator.proces();
                    this.pNugGlogalDeTrended = vEvaluator.nugget;
                    this.pSGlogalDeTrended = vEvaluator.sill;
                    this.pAGlobalDeTrended = vEvaluator.range;
                    this.outGlobalDetrendedVariogramType = vEvaluator.outSemivariogramType;
                    this.pm.message("Global value with TREND for nugget: " + this.pNugGlogalDeTrended + " sill:" + this.pSGlogalDeTrended + " range: " + this.pAGlobalDeTrended + "  semivariogram type:" + vEvaluator.outSemivariogramType);
                } else {
                    this.pm.message("no trend has been found, so no parameters has been evauate");
                }
            }
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }
}

