/*
 * GNU GPL v3 License
 *
 * Copyright 2015 Marialaura Bancheri
 *
 * This program is free software: you can redistribute it and/or modify
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
package experimentalVariogram;

import java.util.HashMap;
import krigingsPointCase.StationsSelection;
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
import org.geotools.feature.SchemaException;
import org.hortonmachine.gears.libs.modules.HMModel;
import org.hortonmachine.gears.libs.monitor.IHMProgressMonitor;
import org.hortonmachine.gears.libs.monitor.LogProgressMonitor;

@Description(value="Experimental semivariogram algorithm.")
@Documentation(value="Experimental semivariogram")
@Author(name="Giuseppe Formetta, Francesco Adami, Silvia Franceschi & Marialaura Bancheri")
@Keywords(value="Experimental semivariogram, Kriging, Hydrology")
@Label(value="HortonMachine/Statistics")
@Name(value="variogram")
@Status(value=40)
@License(value="General Public License Version 3 (GPLv3)")
public class ExperimentalVariogram
extends HMModel {
    @Description(value="The vector of the measurement point, containing the position of the stations.")
    @In
    public SimpleFeatureCollection inStations = null;
    @Description(value="The field of the vector of stations, defining the id.")
    @In
    public String fStationsid = null;
    @Description(value="The file with the measured data, to be interpolated.")
    @In
    public HashMap<Integer, double[]> inData = null;
    @Description(value="Include zeros in computations (default is true).")
    @In
    public boolean doIncludezero = true;
    @Description(value="Specified cutoff")
    @In
    public double Cutoffinput;
    @Description(value="In the case of kriging with neighbor, inNumCloserStations is the number of stations the algorithm has to consider")
    @In
    public int inNumCloserStations;
    @Description(value="Number of bins to consider in the anlysis")
    @In
    public int Cutoff_divide;
    @Description(value="The Experimental Distances.")
    @Out
    public HashMap<Integer, double[]> outDistances;
    @Description(value="The Experimental Variogram.")
    @Out
    public HashMap<Integer, double[]> outExperimentalVariogram;
    @Out
    public int differents;
    @Out
    public boolean areAllEquals;
    @Description(value="The progress monitor.")
    @In
    public IHMProgressMonitor pm = new LogProgressMonitor();

    @Execute
    public void process() throws Exception {
        StationsSelection stations = new StationsSelection();
        stations.inStations = this.inStations;
        stations.inData = this.inData;
        stations.doIncludezero = this.doIncludezero;
        stations.fStationsid = this.fStationsid;
        stations.execute();
        this.differents = stations.n1;
        this.areAllEquals = stations.areAllEquals;
        double[] xStations = stations.xStationInitialSet;
        double[] yStations = stations.yStationInitialSet;
        double[] hStations = stations.hStationInitialSet;
        int[] idStations = stations.idStationInitialSet;
        if (this.differents > 2) {
            double[][] outResult = this.processAlgorithm(xStations, yStations, hStations, this.Cutoffinput);
            this.storeResult(outResult);
        } else {
            System.out.println("Only 1 data >0 or All the data are equal. Variogram is not running");
        }
    }

    private double[][] processAlgorithm(double[] xStation, double[] yStation, double[] hStation, double Cutoffinput) {
        int i;
        double mean = 0.0;
        double maxDistance = 0.0;
        int iCount = xStation.length;
        double[][] distanceMatrix = new double[iCount][iCount];
        double x_max = xStation[0];
        double y_max = yStation[0];
        double x_min = xStation[0];
        double y_min = yStation[0];
        for (i = 1; i < iCount - 1; ++i) {
            x_min = Math.min(x_min, xStation[i]);
            y_min = Math.min(y_min, yStation[i]);
            x_max = Math.max(x_max, xStation[i]);
            y_max = Math.max(y_max, yStation[i]);
        }
        double diagonal = Math.sqrt((x_max - x_min) * (x_max - x_min) + (y_max - y_min) * (y_max - y_min));
        double Cutoff = Cutoffinput == 0.0 ? diagonal / 3.0 : Cutoffinput;
        for (i = 0; i < iCount - 1; ++i) {
            double x1 = xStation[i];
            double y1 = yStation[i];
            double value = hStation[i];
            mean += value;
            for (int j = 0; j < iCount - 1; ++j) {
                double x2 = xStation[j];
                double y2 = yStation[j];
                double dDifX = x2 - x1;
                double dDifY = y2 - y1;
                distanceMatrix[i][j] = Math.sqrt(dDifX * dDifX + dDifY * dDifY);
                maxDistance = Math.max(maxDistance, distanceMatrix[i][j]);
            }
        }
        return this.calculate(Cutoff, distanceMatrix, hStation, mean /= (double)iCount, maxDistance);
    }

    public double[][] calculate(double cutoff, double[][] distanceMatrix, double[] hStation, double mean, double maxDistance) {
        this.Cutoff_divide = this.Cutoff_divide == 0 ? 15 : this.Cutoff_divide;
        double binAmplitude = cutoff / (double)this.Cutoff_divide;
        int iClasses = (int)(maxDistance / binAmplitude + 2.0);
        double[] m_dSemivar = new double[iClasses];
        int[] iPointsInClass = new int[iClasses];
        double[] m_ddist = new double[iClasses];
        int contaNONzero = 0;
        for (int i = 0; i < distanceMatrix.length; ++i) {
            double value1 = hStation[i];
            for (int j = i + 1; j < distanceMatrix.length; ++j) {
                int iClass;
                if (!(distanceMatrix[i][j] > 0.0) || !(distanceMatrix[i][j] < cutoff)) continue;
                int n = iClass = (int)Math.floor(distanceMatrix[i][j] / binAmplitude);
                iPointsInClass[n] = iPointsInClass[n] + 1;
                double value2 = hStation[j];
                double dSemivar = Math.pow(value1 - value2, 2.0);
                int n2 = iClass;
                m_dSemivar[n2] = m_dSemivar[n2] + dSemivar;
                int n3 = iClass;
                m_ddist[n3] = m_ddist[n3] + distanceMatrix[i][j];
            }
        }
        double[][] result = new double[this.Cutoff_divide][2];
        for (int i = 0; i < this.Cutoff_divide; ++i) {
            contaNONzero = iPointsInClass[i] == 0 ? contaNONzero : contaNONzero + 1;
            m_dSemivar[i] = iPointsInClass[i] == 0 ? 0.0 : m_dSemivar[i] / (2.0 * (double)iPointsInClass[i]);
            m_ddist[i] = iPointsInClass[i] == 0 ? 0.0 : m_ddist[i] / (double)iPointsInClass[i];
            result[i][0] = m_ddist[i];
            result[i][1] = m_dSemivar[i];
        }
        return result;
    }

    private void storeResult(double[][] result) throws SchemaException {
        this.outDistances = new HashMap();
        this.outExperimentalVariogram = new HashMap();
        for (int i = 0; i < result.length; ++i) {
            this.outDistances.put(i, new double[]{result[i][0]});
            this.outExperimentalVariogram.put(i, new double[]{result[i][1]});
        }
    }
}

