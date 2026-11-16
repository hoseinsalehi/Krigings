/*
 * GNU GPL v3 License (GEOframe Krigings, version 0.9.9 of 2023-01-11).
 */
/* Source reconstructed from Kriging-2023-01-11-0.9.9.jar (CFR 0.152 decompilation, cleaned 2026; see UPSTREAM_DIFF.md). */
package oldPointCase;

import flanagan.analysis.Regression;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Set;
import krigingsPointCase.StationsSelection;
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
    public boolean doIncludezero = true;
    @Description(value="The range if the models runs with the gaussian variogram.")
    @In
    public double range;
    @Description(value="The sill if the models runs with the gaussian variogram.")
    @In
    public double sill;
    @Description(value="Is the nugget if the models runs with the gaussian variogram.")
    @In
    public double nugget;
    @Description(value="Nugget, sill and range for the current time step. The parameters are provided as an OMS time series.")
    @In
    public HashMap<Integer, double[]> inParamTS;
    @Description(value="In the case of kriging with neighbor, maxdist is the maximum distance within the algorithm has to consider the stations")
    @In
    public double maxdist;
    @Description(value="In the case of kriging with neighbor, inNumCloserStations is the number of stations the algorithm has to consider")
    @In
    public int inNumCloserStations;
    @Description(value="Switch for detrended mode.")
    @In
    public boolean doDetrended;
    @Description(value="The double value of the trend")
    @In
    public double trend_intercept;
    @Description(value="The double value of the trend")
    @In
    public double trend_coefficient;
    @Description(value="Degree of polynomial regression, default is 1")
    @In
    public int regressionOrder = 1;
    @Description(value="Type of linear system solver")
    @In
    public String linearSystemSolverType = "default";
    @Description(value="The hashmap withe the interpolated results")
    @Out
    public HashMap<Integer, double[]> outData = null;
    private static final double TOLL = 1.0E-7;
    private HortonMessageHandler msg = HortonMessageHandler.getInstance();
    int id;

    @Execute
    public void executeKriging() throws Exception {
        this.verifyInput();
        LinkedHashMap<Integer, Coordinate> pointsToInterpolateId2Coordinates = null;
        pointsToInterpolateId2Coordinates = this.getCoordinate(0, this.inInterpolate, this.fInterpolateid);
        Set<Integer> pointsToInterpolateIdSet = pointsToInterpolateId2Coordinates.keySet();
        Iterator<Integer> idIterator = pointsToInterpolateIdSet.iterator();
        int j = 0;
        double[] result = new double[pointsToInterpolateId2Coordinates.size()];
        int[] idArray = new int[pointsToInterpolateId2Coordinates.size()];
        if (this.inParamTS != null) {
            this.nugget = this.inParamTS.get(0)[0];
            this.sill = this.inParamTS.get(1)[0];
            this.range = this.inParamTS.get(2)[0];
        }
        while (idIterator.hasNext()) {
            double sum = 0.0;
            idArray[j] = this.id = idIterator.next().intValue();
            Coordinate coordinate = pointsToInterpolateId2Coordinates.get(this.id);
            StationsSelection stations = new StationsSelection();
            stations.idx = coordinate.x;
            stations.idy = coordinate.y;
            stations.inStations = this.inStations;
            stations.inData = this.inData;
            stations.doIncludezero = this.doIncludezero;
            stations.maxdist = this.maxdist;
            stations.inNumCloserStations = this.inNumCloserStations;
            stations.fStationsid = this.fStationsid;
            stations.fStationsZ = this.fStationsZ;
            stations.execute();
            double[] xStations = stations.xStationInitialSet;
            double[] yStations = stations.yStationInitialSet;
            double[] zStations = stations.zStationInitialSet;
            double[] hStations = stations.hStationInitialSet;
            boolean areAllEquals = stations.areAllEquals;
            int n1 = xStations.length - 1;
            xStations[n1] = coordinate.x;
            yStations[n1] = coordinate.y;
            zStations[n1] = coordinate.z;
            double[] hresiduals = hStations;
            if (this.doDetrended) {
                Regression r = new Regression();
                r = new Regression(zStations, hStations);
                r.polynomial(this.regressionOrder);
                this.trend_intercept = r.getBestEstimates()[0];
                this.trend_coefficient = r.getBestEstimates()[1];
                hresiduals = r.getResiduals();
            }
            if (n1 != 0) {
                if (!areAllEquals && n1 > 1) {
                    this.pm.beginTask(this.msg.message("kriging.working"), pointsToInterpolateId2Coordinates.size());
                    double h0 = 0.0;
                    double[][] covarianceMatrix = this.covMatrixCalculating(xStations, yStations, zStations, n1);
                    double[] knownTerm = this.knownTermsCalculation(xStations, yStations, zStations, n1);
                    ColumnVector solution = SimpleLinearSystemSolverFactory.solve(knownTerm, covarianceMatrix, this.linearSystemSolverType);
                    double[] moltiplicativeFactor = solution.copyValues1D();
                    for (int k = 0; k < n1; ++k) {
                        h0 += moltiplicativeFactor[k] * hresiduals[k];
                        sum += moltiplicativeFactor[k];
                    }
                    double trend = this.doDetrended ? coordinate.z * this.trend_coefficient + this.trend_intercept : 0.0;
                    result[j] = h0 += trend;
                    ++j;
                    if (Math.abs(sum - 1.0) >= 1.0E-7) {
                        throw new ModelsRuntimeException("Error in the coffeicients calculation", this.getClass().getSimpleName());
                    }
                    this.pm.worked(1);
                } else if (n1 == 1 || areAllEquals) {
                    double tmp = hresiduals[0];
                    this.pm.message(this.msg.message("kriging.setequalsvalue"));
                    this.pm.beginTask(this.msg.message("kriging.working"), pointsToInterpolateId2Coordinates.size());
                    result[j] = tmp;
                    ++j;
                    n1 = 0;
                    this.pm.worked(1);
                }
                this.pm.done();
                continue;
            }
            this.pm.errorMessage("No value for this time step");
            double[] value = this.inData.values().iterator().next();
            result[j] = value[0];
            ++j;
        }
        this.storeResult(result, idArray);
    }

    private void verifyInput() {
        if (this.inData == null || this.inStations == null) {
            throw new NullPointerException(this.msg.message("kriging.stationProblem"));
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
        double[][] covarianceMatrix = new double[n + 1][n + 1];
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                double rx = x[i] - x[j];
                double ry = y[i] - y[j];
                double rz = z[i] - z[j];
                covarianceMatrix[j][i] = this.variogram(this.nugget, this.range, this.sill, rx, ry, rz);
                covarianceMatrix[i][j] = this.variogram(this.nugget, this.range, this.sill, rx, ry, rz);
            }
        }
        for (int i = 0; i < n; ++i) {
            covarianceMatrix[i][n] = 1.0;
            covarianceMatrix[n][i] = 1.0;
        }
        covarianceMatrix[n][n] = 0.0;
        return covarianceMatrix;
    }

    private double[] knownTermsCalculation(double[] x, double[] y, double[] z, int n) {
        double[] gamma = new double[n + 1];
        for (int i = 0; i < n; ++i) {
            double rx = x[i] - x[n];
            double ry = y[i] - y[n];
            double rz = z[i] - z[n];
            gamma[i] = this.variogram(this.nugget, this.range, this.sill, rx, ry, rz);
        }
        gamma[n] = 1.0;
        return gamma;
    }

    private double variogram(double nug, double range, double sill, double rx, double ry, double rz) {
        double vgmResult;
        double h2;
        if (HMConstants.isNovalue(rz)) {
            rz = 0.0;
        }
        if ((h2 = Math.sqrt(rx * rx + rz * rz + ry * ry)) != 0.0) {
            TheoreticalVariogram vgm = new TheoreticalVariogram();
            vgmResult = vgm.calculateVGM(this.pSemivariogramType, h2, sill, range, nug);
        } else {
            vgmResult = 0.0;
        }
        return vgmResult;
    }

    private void storeResult(double[] result, int[] id) throws SchemaException {
        this.outData = new HashMap();
        for (int i = 0; i < result.length; ++i) {
            this.outData.put(id[i], new double[]{result[i]});
        }
    }
}

