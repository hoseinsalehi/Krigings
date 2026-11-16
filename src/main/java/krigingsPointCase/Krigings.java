/* This file is part of JGrasstools (http://www.jgrasstools.org)
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * JGrasstools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/* Source reconstructed from Kriging-2023-01-11-0.9.9.jar (CFR 0.152 decompilation, cleaned 2026; see UPSTREAM_DIFF.md). */
package krigingsPointCase;

import experimentalVariogram.ExperimentalVariogram;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Set;
import krigingsPointCase.GlobalParameterEvaluator;
import krigingsPointCase.ResidualsEvaluator;
import krigingsPointCase.StationsSelection;
import krigingsPointCase.Utility;
import krigingsPointCase.VariogramParamsEvaluator;
import linearSistemSolver.SimpleLinearSystemSolverFactory;
import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.feature.SchemaException;
import org.hortonmachine.gears.libs.exceptions.ModelsRuntimeException;
import org.hortonmachine.gears.libs.modules.HMConstants;
import org.hortonmachine.gears.libs.modules.HMModel;
import org.hortonmachine.gears.libs.monitor.IHMProgressMonitor;
import org.hortonmachine.gears.libs.monitor.LogProgressMonitor;
import org.hortonmachine.gears.utils.math.matrixes.ColumnVector;
import org.hortonmachine.hmachine.i18n.HortonMessageHandler;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import theoreticalVariogram.TheoreticalVariogram;

@Description(value="Ordinary kriging algorithm.")
@Documentation(value="Kriging.html")
@Author(name="Giuseppe Formetta, Daniele Andreis, Silvia Franceschi, Andrea Antonello, Marialaura Bancheri & Francesco Serafin")
@Keywords(value="Kriging, Hydrology")
@Label(value="")
@Name(value="kriging")
@Status
@License(value="General Public License Version 3 (GPLv3)")
public class Krigings
extends HMModel {
    @Description(value="The .shp of the measurement point, containing the position of the stations.")
    @In
    public SimpleFeatureCollection inStations = null;
    @Description(value="The field of the vector of stations, defining the id.")
    @In
    public String fStationsid = null;
    @Description(value="The field of the vector of stations, defining the elevation.")
    @In
    public String fStationsZ = null;
    @Description(value="The type of theoretical semivariogram: exponential, gaussian, spherical, pentasphericallinear, circular, bessel, periodic, hole, logaritmic, power, spline")
    @In
    public String pSemivariogramType = null;
    @Description(value="The HM with the measured data to be interpolated.")
    @In
    public HashMap<Integer, double[]> inData = null;
    @Description(value="The vector of the points in which the data have to be interpolated.")
    @In
    public SimpleFeatureCollection inInterpolate = null;
    @Description(value="The field of the interpolated vector points, defining the id.")
    @In
    public String fInterpolateid = null;
    @Description(value="The field of the interpolated vector points, defining the elevation.")
    @In
    public String fPointZ = null;
    @Description(value="The progress monitor.")
    @In
    public IHMProgressMonitor pm = new LogProgressMonitor();
    @Description(value="Include zeros in computations (default is true).")
    @In
    public boolean doIncludeZero = true;
    @Description(value="Deprecated spelling of doIncludeZero (GitHub master API); false here also excludes zeros.")
    @In
    public boolean doIncludezero = true;
    @Description(value="Nugget, sill and range for the current time step as an OMS time series (GitHub master API); overrides nugget, sill, range.")
    @In
    public HashMap<Integer, double[]> inParamTS;
    @Description(value="The range if the models runs with the gaussian variogram.")
    @In
    public double range;
    @Description(value="The sill if the models runs with the gaussian variogram.")
    @In
    public double sill;
    @Description(value="Is the nugget if the models runs with the gaussian variogram.")
    @In
    public double nugget;
    @Description(value="In the case of kriging with neighbor, maxdist is the maximum distance within the algorithm has to consider the stations")
    @In
    public double maxdist;
    @Description(value="In the case of kriging with neighbor, inNumCloserStations is the number of stations the algorithm has to consider")
    @In
    public int inNumCloserStations;
    @Description(value="Switch for detrended mode.")
    @In
    public boolean doDetrended;
    @Description(value="Degree of polynomial regression, default is 1")
    @In
    private static final int REGRESSION_ORDER = 1;
    @Description(value="Type of linear system solver")
    @In
    public String linearSystemSolverType = "default";
    @Description(value="Specified cutoff")
    @In
    public double cutoffInput;
    // ------------------------------------------------------------------ local PSO variogram option (2026, default off)
    @Description(value="New option (default false = 0.9.9 behaviour): per-target local spherical variogram on the nVario nearest stations fitted by PSO (numpy-compatible random stream), ordinary kriging on the nKrig nearest (Po engine, krigingsCore.local)")
    @In
    public boolean doLocalVariogram = false;
    @Description(value="local variogram: stations defining each local variogram (and local lapse)")
    @In
    public int nVario = 50;
    @Description(value="local variogram: kriging neighbours")
    @In
    public int nKrig = 7;
    @Description(value="local variogram: trend removal none | local (per-target lapse on the nVario stations) | global (one regression on all stations of the time step)")
    @In
    public String localTrend = "none";
    @Description(value="local variogram: estimate 0 when all nKrig neighbours are <= 0 (precipitation)")
    @In
    public boolean doZeroShortcut = false;
    @Description(value="local variogram: random stream numpy | java; seeds = BLAKE2b(seedTag|step, chunk offset)")
    @In
    public String rngMode = "numpy";
    @Description(value="local variogram: tag of the random-stream seeds (the time-step counter is appended)")
    @In
    public String seedTag = "krigings";
    @Description(value="local variogram: worker threads over chunks of targets (results independent of the number)")
    @In
    public int threads = 1;
    @Description(value="local variogram: kriging variance per interpolation point (other modes: not filled)")
    @Out
    public HashMap<Integer, double[]> outKrigingVariance = null;
    @Description(value="Number of bins to consider in the anlysis")
    @In
    public int cutoffDivide;
    @Description(value="Distances input path")
    @In
    public String inHValuesPath;
    @In
    public boolean doLogarithmic = false;
    @In
    public boolean boundedToZero = false;
    @Description(value="The hashmap withe the interpolated results")
    @Out
    public HashMap<Integer, double[]> outData = null;
    @Description(value="The hashmap withe the parameter for each time step: 0=>nugget, 1=>sill, 2=> range, 3 => is local (0 if used local parameter 1 if use global, 4 => is detrended (0 no trend, 1 with trend)")
    @Out
    public HashMap<Integer, double[]> outVariogramParams = null;
    private static final double TOLL = 1.0E-7;
    private int step = 0;
    private HortonMessageHandler msg = HortonMessageHandler.getInstance();
    int id;
    @In
    public double pAGlobal;
    @In
    public double pSGlobal;
    @In
    public double pNugGlobal;
    @In
    public double pAGlobalDeTrended;
    @In
    public double pSGlobalDeTrended;
    @In
    public double pNugGlobalDeTrended;
    @In
    public String fileNoValue = "-9999";
    @In
    public String tStart = null;
    @In
    public int tTimeStep = 60;
    private double aGlobal;
    private double sGlobal;
    private double nugGlobal;
    private double aGlobalDeTrended;
    private double sGlobalDeTrended;
    private double nugGlobalDeTrended;
    private boolean noGlobalParams = false;
    @In
    public String globalDetrendedVariogramType;
    @In
    public String globalVariogramType;
    private String actualSemivariogramType;
    private boolean noEvaluationGlobalParams;
    private boolean fixedParams = false;
    private double fixedNugget, fixedSill, fixedRange;

    @Execute
    public void executeKriging() throws Exception {
        if (this.doLocalVariogram) {
            executeLocalVariogram();
            return;
        }
        if (this.step == 0) {
            try {
                this.verifyInput();
                this.createDefaulParams();
            }
            catch (Exception e) {
                this.pm.errorMessage(e.toString());
            }
        }
        if (this.step == 0) {
            // 0.9.9 reset nugget/sill/range to 0 at every step, so fixed parameters given in the sim were never used
            // (although verifyInput accepts them). They are now honoured when all are given with pSemivariogramType.
            this.fixedParams = this.pSemivariogramType != null && !this.pSemivariogramType.isEmpty()
                    && !(this.nugget == 0.0 && this.sill == 0.0 && this.range == 0.0);
            this.fixedNugget = this.nugget; this.fixedSill = this.sill; this.fixedRange = this.range;
        }
        if (this.inParamTS != null) {
            this.fixedParams = true;
            this.fixedNugget = this.inParamTS.get(0)[0]; this.fixedSill = this.inParamTS.get(1)[0]; this.fixedRange = this.inParamTS.get(2)[0];
        }
        if (this.fixedParams) {
            this.nugget = this.fixedNugget; this.sill = this.fixedSill; this.range = this.fixedRange;
        } else {
            this.nugget = 0.0;
            this.sill = 0.0;
            this.range = 0.0;
        }
        LinkedHashMap<Integer, Coordinate> pointsToInterpolateId2Coordinates = null;
        pointsToInterpolateId2Coordinates = this.getCoordinate(0, this.inInterpolate, this.fInterpolateid);
        Set<Integer> pointsToInterpolateIdSet = pointsToInterpolateId2Coordinates.keySet();
        Iterator<Integer> idIterator = pointsToInterpolateIdSet.iterator();
        int j = 0;
        double[] result = new double[pointsToInterpolateId2Coordinates.size()];
        int[] idArray = new int[pointsToInterpolateId2Coordinates.size()];
        StationsSelection stations = new StationsSelection();
        stations.inStations = this.inStations;
        stations.doIncludezero = this.doIncludeZero && this.doIncludezero;
        stations.maxdist = this.maxdist;
        stations.fStationsid = this.fStationsid;
        stations.fStationsZ = this.fStationsZ;
        stations.doLogarithmic = this.doLogarithmic;
        ++this.step;
        stations.inData = this.inData;
        stations.execute();
        double[] xStations = stations.xStationInitialSet;
        double[] yStations = stations.yStationInitialSet;
        double[] zStations = stations.zStationInitialSet;
        double[] hStations = stations.hStationInitialSet;
        int n1 = xStations.length - 1;
        double[] hResiduals = new double[n1];
        boolean areAllEquals = stations.areAllEquals;
        double trendCoeff = 0.0;
        double trendIntercept = 0.0;
        if (!areAllEquals && n1 > 1) {
            ResidualsEvaluator residualsEvaluator = this.getResidualsEvaluator(Arrays.copyOfRange(zStations, 0, n1), Arrays.copyOfRange(hStations, 0, n1));
            hResiduals = residualsEvaluator.hResiduals;
            trendCoeff = residualsEvaluator.trend_coefficient;
            trendIntercept = residualsEvaluator.trend_intercept;
            this.setKrigingParams(stations, residualsEvaluator);
        } else {
            this.outVariogramParams = new HashMap();
            this.outVariogramParams.put(0, new double[]{-9999.0});
            this.outVariogramParams.put(1, new double[]{-9999.0});
            this.outVariogramParams.put(2, new double[]{-9999.0});
            this.outVariogramParams.put(3, new double[]{-9999.0});
            this.outVariogramParams.put(4, new double[]{-9999.0});
            this.outVariogramParams.put(4, new double[]{-9999.0});
        }
        while (idIterator.hasNext()) {
            double sum = 0.0;
            idArray[j] = this.id = idIterator.next().intValue();
            Coordinate coordinate = pointsToInterpolateId2Coordinates.get(this.id);
            if (this.maxdist > 0.0 || this.inNumCloserStations > 0) {
                if (this.inNumCloserStations > 0) {
                    stations.inNumCloserStations = this.inNumCloserStations;
                }
                if (this.maxdist > 0.0) {
                    stations.maxdist = this.maxdist;
                }
                stations.idx = coordinate.x;
                stations.idy = coordinate.y;
                stations.execute();
                xStations = stations.xStationInitialSet;
                yStations = stations.yStationInitialSet;
                zStations = stations.zStationInitialSet;
                hStations = stations.hStationInitialSet;
                n1 = xStations.length - 1;
            }
            xStations[n1] = coordinate.x;
            yStations[n1] = coordinate.y;
            if (n1 != 0) {
                if (!areAllEquals && n1 > 1) {
                    this.pm.beginTask(this.msg.message("kriging.working"), pointsToInterpolateId2Coordinates.size());
                    double h0 = 0.0;
                    double[][] covarianceMatrix = this.covMatrixCalculating(xStations, yStations, zStations, n1);
                    double[] knownTerm = this.knownTermsCalculation(xStations, yStations, zStations, n1);
                    ColumnVector solution = SimpleLinearSystemSolverFactory.solve(knownTerm, covarianceMatrix, this.linearSystemSolverType);
                    double[] moltiplicativeFactor = solution.copyValues1D();
                    for (int k = 0; k < n1; ++k) {
                        h0 += moltiplicativeFactor[k] * hResiduals[k];
                        sum += moltiplicativeFactor[k];
                    }
                    if (this.doDetrended) {
                        double trend = coordinate.z * trendCoeff + trendIntercept;
                        h0 += trend;
                    }
                    result[j] = h0;
                    ++j;
                    if (Math.abs(sum - 1.0) >= 1.0E-7) {
                        throw new ModelsRuntimeException("Error in the coffeicients calculation", this.getClass().getSimpleName());
                    }
                    this.pm.worked(1);
                } else if (n1 == 1 || areAllEquals) {
                    double tmp;
                    result[j] = tmp = hResiduals[0];
                    ++j;
                    n1 = 0;
                    this.pm.worked(1);
                }
                this.pm.done();
            } else {
                this.pm.errorMessage("No value for this time step");
                double[] value = this.inData.values().iterator().next();
                result[j] = value[0];
                ++j;
            }
            if (HMConstants.isNovalue(result[j - 1])) continue;
            if (this.doLogarithmic) {
                result[j - 1] = Utility.getInverseLog(result[j - 1]);
            }
            if (!this.boundedToZero || !(result[j - 1] < 0.0)) continue;
            result[j - 1] = 0.0;
        }
        this.storeResult(result, idArray);
    }

    private ResidualsEvaluator getResidualsEvaluator(double[] zStations, double[] hStations) {
        ResidualsEvaluator residualsEvaluator = new ResidualsEvaluator();
        residualsEvaluator.doDetrended = this.doDetrended;
        residualsEvaluator.hStations = hStations;
        residualsEvaluator.zStations = zStations;
        residualsEvaluator.regressionOrder = 1;
        residualsEvaluator.process();
        return residualsEvaluator;
    }

    private void setKrigingParams(StationsSelection stations, ResidualsEvaluator residualsEvaluator) throws Exception {
        int n1 = stations.hStationInitialSet.length - 1;
        double[] hResiduals = residualsEvaluator.hResiduals;
        if (this.inNumCloserStations != 0) {
            stations.inNumCloserStations = 0;
            stations.execute();
            n1 = stations.hStationInitialSet.length - 1;
            residualsEvaluator = this.getResidualsEvaluator(Arrays.copyOfRange(stations.zStationInitialSet, 0, n1), Arrays.copyOfRange(stations.hStationInitialSet, 0, n1));
            residualsEvaluator.process();
        }
        double isLocal = 0.0;
        double isTrend = 0.0;
        int[] idStations = Arrays.copyOfRange(stations.idStationInitialSet, 0, n1);
        hResiduals = residualsEvaluator.hResiduals;
        if (this.nugget == 0.0 && this.sill == 0.0 && this.range == 0.0) {
            VariogramParamsEvaluator variogramParamsEvaluator = new VariogramParamsEvaluator();
            variogramParamsEvaluator.expVar = this.getExperimentalVariogram(hResiduals, idStations);
            variogramParamsEvaluator.pSemivariogramType = this.pSemivariogramType;
            variogramParamsEvaluator.proces();
            if (!residualsEvaluator.isPValueOk) {
                if (variogramParamsEvaluator.nugget >= 0.0 && variogramParamsEvaluator.sill > 0.0 && variogramParamsEvaluator.range > 0.0 && variogramParamsEvaluator.isFitGood || this.noGlobalParams) {
                    this.nugget = variogramParamsEvaluator.nugget;
                    this.sill = variogramParamsEvaluator.sill;
                    this.range = variogramParamsEvaluator.range;
                    this.actualSemivariogramType = variogramParamsEvaluator.outSemivariogramType;
                } else {
                    this.nugget = this.nugGlobal;
                    this.sill = this.sGlobal;
                    this.range = this.aGlobal;
                    this.actualSemivariogramType = this.globalVariogramType;
                    isLocal = 1.0;
                }
            } else if (residualsEvaluator.isPValueOk) {
                isTrend = 1.0;
                if (variogramParamsEvaluator.nugget >= 0.0 && variogramParamsEvaluator.sill > 0.0 && variogramParamsEvaluator.range > 0.0 && variogramParamsEvaluator.isFitGood || this.noGlobalParams) {
                    this.nugget = variogramParamsEvaluator.nugget;
                    this.sill = variogramParamsEvaluator.sill;
                    this.range = variogramParamsEvaluator.range;
                    this.actualSemivariogramType = variogramParamsEvaluator.outSemivariogramType;
                } else {
                    this.nugget = this.nugGlobalDeTrended;
                    this.sill = this.sGlobalDeTrended;
                    this.range = this.aGlobalDeTrended;
                    this.actualSemivariogramType = this.globalDetrendedVariogramType;
                    isLocal = 1.0;
                }
            }
        } else {
            this.actualSemivariogramType = this.pSemivariogramType;
        }
        this.outVariogramParams = new HashMap();
        this.outVariogramParams.put(0, new double[]{this.nugget});
        this.outVariogramParams.put(1, new double[]{this.sill});
        this.outVariogramParams.put(2, new double[]{this.range});
        this.outVariogramParams.put(3, new double[]{isLocal});
        this.outVariogramParams.put(4, new double[]{isTrend});
        this.outVariogramParams.put(5, new double[]{Utility.getVariogramCode(this.actualSemivariogramType)});
    }

    private ExperimentalVariogram getExperimentalVariogram(double[] hresiduals, int[] idArray) {
        ExperimentalVariogram expVariogram = Utility.getExperimentalVariogram(this.fStationsid, this.inStations, this.doIncludeZero && this.doIncludezero, this.cutoffDivide, this.cutoffInput, 0);
        HashMap<Integer, double[]> tmpInData = new HashMap<Integer, double[]>();
        for (int i = 0; i < idArray.length; ++i) {
            tmpInData.put(idArray[i], new double[]{hresiduals[i]});
        }
        expVariogram.inData = tmpInData;
        return expVariogram;
    }

    private void verifyInput() {
        if (this.inData == null || this.inStations == null) {
            throw new NullPointerException(this.msg.message("kriging.stationProblem"));
        }
        if (this.doDetrended) {
            if (this.fPointZ == null || this.fStationsZ == null) {
                throw new NullPointerException("z field not found");
            }
            int ff = ((SimpleFeatureType)this.inStations.getSchema()).indexOf(this.fStationsZ);
            int ff2 = ((SimpleFeatureType)this.inInterpolate.getSchema()).indexOf(this.fPointZ);
            if (ff < 0 || ff2 < 0) {
                throw new NullPointerException("check if the z field name is correct");
            }
        }
        if (!(this.nugget == 0.0 && this.sill == 0.0 && this.range == 0.0 || this.pSemivariogramType != null && !this.pSemivariogramType.isEmpty())) {
            throw new NullPointerException("You provide incomplete fixed parameter check: nugget, sill, range or pSemivariogramType");
        }
    }

    private void createDefaulParams() throws Exception {
        boolean gloablZeroDeTrended;
        this.nugGlobal = this.pNugGlobal;
        this.aGlobal = this.pAGlobal;
        this.sGlobal = this.pSGlobal;
        this.nugGlobalDeTrended = this.pNugGlobalDeTrended;
        this.sGlobalDeTrended = this.pSGlobalDeTrended;
        this.aGlobalDeTrended = this.pAGlobalDeTrended;
        this.noEvaluationGlobalParams = this.inHValuesPath == null || this.inHValuesPath.isEmpty();
        boolean gloablZero = this.nugGlobal == 0.0 && this.aGlobal == 0.0 && this.sGlobal == 0.0;
        boolean bl = gloablZeroDeTrended = this.nugGlobalDeTrended == 0.0 && this.aGlobalDeTrended == 0.0 && this.sGlobalDeTrended == 0.0;
        if (this.noEvaluationGlobalParams && !gloablZero && this.globalVariogramType == null) {
            throw new NullPointerException("You provide incomplete fixed parameter for global: nugGlobal, sGlobal, rGlobal or globalVariogramType");
        }
        if (this.noEvaluationGlobalParams && !gloablZeroDeTrended && this.globalDetrendedVariogramType == null) {
            throw new NullPointerException("You provide incomplete fixed parameter for global de trended: nugGlobal, sGlobal, rGlobal or globalVariogramType");
        }
        if (this.doDetrended) {
            this.noGlobalParams = this.noEvaluationGlobalParams && gloablZero || this.globalVariogramType == null;
            this.noEvaluationGlobalParams = this.noEvaluationGlobalParams || !gloablZero && this.globalVariogramType != null;
        } else {
            this.noGlobalParams = this.noEvaluationGlobalParams && (gloablZero || this.globalVariogramType == null || gloablZeroDeTrended || this.globalDetrendedVariogramType == null);
            boolean bl2 = this.noEvaluationGlobalParams = this.noEvaluationGlobalParams || !gloablZero && this.globalVariogramType != null && !gloablZeroDeTrended && this.globalDetrendedVariogramType != null;
        }
        if (this.nugget == 0.0 && this.range == 0.0 && this.nugget == 0.0 && this.noGlobalParams) {
            this.pm.message("Warning: nugget,range and sill will be evaluate only at each timestep");
        }
        if (!this.noEvaluationGlobalParams) {
            GlobalParameterEvaluator gParam = new GlobalParameterEvaluator();
            gParam.inStations = this.inStations;
            gParam.doIncludeZero = this.doIncludeZero && this.doIncludezero;
            gParam.maxdist = this.maxdist;
            gParam.fStationsid = this.fStationsid;
            gParam.fStationsZ = this.fStationsZ;
            gParam.doLogarithmic = this.doLogarithmic;
            gParam.cutoffDivide = this.cutoffDivide;
            gParam.inHValuesPath = this.inHValuesPath;
            gParam.fStationsid = "ID";
            gParam.fileNoValue = this.fileNoValue;
            gParam.tStart = this.tStart;
            gParam.tTimeStep = this.tTimeStep;
            gParam.pSemivariogramType = this.pSemivariogramType;
            gParam.doDetrended = this.doDetrended;
            gParam.cutoffInput = this.cutoffInput;
            gParam.execute();
            this.nugGlobal = gParam.pNugGlogal;
            this.aGlobal = gParam.pAGlobal;
            this.sGlobal = gParam.pSGlogal;
            this.globalVariogramType = gParam.outGlobalVariogramType;
            this.nugGlobalDeTrended = gParam.pNugGlogalDeTrended;
            this.sGlobalDeTrended = gParam.pSGlogalDeTrended;
            this.aGlobalDeTrended = gParam.pAGlobalDeTrended;
            this.globalDetrendedVariogramType = gParam.outGlobalDetrendedVariogramType;
        }
    }

    public static double round(double value, int places) {
        if (places < 0) {
            throw new IllegalArgumentException();
        }
        long factor = (long)Math.pow(10.0, places);
        long tmp = Math.round(value *= (double)factor);
        return (double)tmp / (double)factor;
    }

    private LinkedHashMap<Integer, Coordinate> getCoordinate(int nStaz, SimpleFeatureCollection collection, String idField) throws Exception {
        LinkedHashMap<Integer, Coordinate> id2CoordinatesMcovarianceMatrix = new LinkedHashMap<Integer, Coordinate>();
        Coordinate coordinate = null;
        try (SimpleFeatureIterator iterator = collection.features();){
            while (iterator.hasNext()) {
                SimpleFeature feature = (SimpleFeature)iterator.next();
                int name = ((Number)feature.getAttribute(idField)).intValue();
                coordinate = ((Geometry)feature.getDefaultGeometry()).getCentroid().getCoordinate();
                double z = 0.0;
                if (this.fPointZ != null) {
                    try {
                        z = ((Number)feature.getAttribute(this.fPointZ)).doubleValue();
                    }
                    catch (NullPointerException e) {
                        this.pm.errorMessage(this.msg.message("kriging.noPointZ"));
                        throw new Exception(this.msg.message("kriging.noPointZ"));
                    }
                }
                coordinate.z = z;
                id2CoordinatesMcovarianceMatrix.put(name, coordinate);
            }
        }
        return id2CoordinatesMcovarianceMatrix;
    }

    private double[][] covMatrixCalculating(double[] x, double[] y, double[] z, int n) {
        return krigingsCore.KrigingSystem.covMatrix(x, y, z, n, this.actualSemivariogramType, this.nugget, this.range, this.sill);
    }

    private double[] knownTermsCalculation(double[] x, double[] y, double[] z, int n) {
        return krigingsCore.KrigingSystem.knownTerms(x, y, z, n, this.actualSemivariogramType, this.nugget, this.range, this.sill);
    }

    private double variogram(double nug, double range, double sill, double rx, double ry, double rz) {
        return krigingsCore.KrigingSystem.variogram(this.actualSemivariogramType, nug, range, sill, rx, ry, rz);
    }


    /**
     * Local-variogram mode (doLocalVariogram = true): for every interpolation point a spherical variogram is fitted by PSO on
     * the nVario nearest stations with a value at this time step and ordinary kriging uses the nKrig nearest (same kernels as
     * pokriging.PoKriging; METHOD_SPEC section 3).
     */
    private void executeLocalVariogram() throws Exception {
        if (this.inData == null || this.inStations == null) throw new NullPointerException(this.msg.message("kriging.stationProblem"));
        java.util.List<double[]> st = new java.util.ArrayList<>();
        try (SimpleFeatureIterator it = this.inStations.features()) {
            while (it.hasNext()) {
                SimpleFeature f = (SimpleFeature) it.next();
                int sid = ((Number) f.getAttribute(this.fStationsid)).intValue();
                double[] h = this.inData.get(sid);
                if (h == null || HMConstants.isNovalue(h[0])) continue;
                double v = this.doLogarithmic ? Utility.getLog(h[0]) : h[0];
                if (!(this.doIncludeZero && this.doIncludezero) && Math.abs(v) <= 0.0) continue;
                Coordinate co = ((Geometry) f.getDefaultGeometry()).getCentroid().getCoordinate();
                double z = 0.0;
                if (this.fStationsZ != null && f.getAttribute(this.fStationsZ) != null) z = ((Number) f.getAttribute(this.fStationsZ)).doubleValue();
                st.add(new double[]{co.x, co.y, z, v});
            }
        }
        LinkedHashMap<Integer, Coordinate> pts = this.getCoordinate(0, this.inInterpolate, this.fInterpolateid);
        int S = st.size(), n = pts.size();
        double[] sx = new double[S], sy = new double[S], sz = new double[S], y = new double[S];
        for (int i = 0; i < S; i++) { double[] r = st.get(i); sx[i] = r[0]; sy[i] = r[1]; sz[i] = r[2]; y[i] = r[3]; }
        float[] sd = new float[S * S];
        for (int i = 0; i < S; i++) for (int j = 0; j < S; j++) { double dx = sx[i] - sx[j], dy = sy[i] - sy[j]; sd[i * S + j] = (float) (Math.sqrt(dx * dx + dy * dy) / 1000.0); }
        int[] act = new int[S]; for (int i = 0; i < S; i++) act[i] = i;
        double[] tx = new double[n], ty = new double[n], tz = new double[n]; int[] self = new int[n]; long[] seed = new long[n]; int[] ids = new int[n];
        int k = 0;
        for (java.util.Map.Entry<Integer, Coordinate> e : pts.entrySet()) { ids[k] = e.getKey(); tx[k] = e.getValue().x; ty[k] = e.getValue().y; tz[k] = e.getValue().z; self[k] = -1; seed[k] = k; k++; }
        this.outData = new HashMap<Integer, double[]>(); this.outKrigingVariance = new HashMap<Integer, double[]>();
        String tag = this.seedTag + "|" + this.step; ++this.step;
        if (S < 3) { for (int i = 0; i < n; i++) { this.outData.put(ids[i], new double[]{HMConstants.doubleNovalue}); this.outKrigingVariance.put(ids[i], new double[]{HMConstants.doubleNovalue}); } return; }
        krigingsCore.local.KrigingSettings ks = new krigingsCore.local.KrigingSettings();
        ks.nVario = Math.min(this.nVario, S - 1); ks.nKrig = Math.min(this.nKrig, ks.nVario); ks.rngMode = this.rngMode;
        double[] yk = y; double gb0 = 0, gb1 = 0;
        if ("global".equals(this.localTrend)) {
            double zm = 0, ym = 0; for (int i = 0; i < S; i++) { zm += sz[i]; ym += y[i]; } zm /= S; ym /= S;
            double szz = 0, szy = 0; for (int i = 0; i < S; i++) { szz += (sz[i] - zm) * (sz[i] - zm); szy += (sz[i] - zm) * (y[i] - ym); }
            gb1 = szz > 1e-8 ? szy / szz : 0.0; gb0 = ym - gb1 * zm;
            yk = new double[S]; for (int i = 0; i < S; i++) yk[i] = y[i] - (gb0 + gb1 * sz[i]);
        }
        boolean lapse = "local".equals(this.localTrend);
        double[] gpar = krigingsCore.local.DayKriging.globalVariogram(ks, sd, S, act, lapse ? krigingsCore.local.DayKriging.detrend(yk, sz, act) : yk, tag);
        krigingsCore.local.DayKriging dk = new krigingsCore.local.DayKriging(ks, sx, sy, sz, sd, act, yk, this.doZeroShortcut,
                lapse ? krigingsCore.local.DayKriging.Lapse.LOCAL : krigingsCore.local.DayKriging.Lapse.NONE, Double.NaN, tag, gpar);
        krigingsCore.local.TargetSet ts = new krigingsCore.local.TargetSet(tx, ty, tz, self, seed);
        krigingsCore.local.TargetResult r = new krigingsCore.local.TargetResult(n);
        java.util.concurrent.ExecutorService ex = java.util.concurrent.Executors.newFixedThreadPool(Math.max(1, this.threads));
        try {
            java.util.List<java.util.concurrent.Future<?>> fs = new java.util.ArrayList<>();
            for (int c = 0; c < dk.nChunks(ts); c++) { final int cc = c; fs.add(ex.submit(() -> dk.chunk(ts, cc, r))); }
            for (java.util.concurrent.Future<?> f : fs) f.get();
        } finally { ex.shutdown(); }
        for (int i = 0; i < n; i++) {
            double v = r.est[i] + gb0 + gb1 * tz[i];
            if (!Double.isFinite(v)) v = HMConstants.doubleNovalue;
            else { if (this.doLogarithmic) v = Utility.getInverseLog(v); if (this.boundedToZero && v < 0.0) v = 0.0; }
            this.outData.put(ids[i], new double[]{v});
            this.outKrigingVariance.put(ids[i], new double[]{Double.isFinite(r.kvar[i]) ? r.kvar[i] : HMConstants.doubleNovalue});
        }
    }

    private void storeResult(double[] result, int[] id) throws SchemaException {
        this.outData = new HashMap();
        for (int i = 0; i < result.length; ++i) {
            this.outData.put(id[i], new double[]{result[i]});
        }
    }
}

