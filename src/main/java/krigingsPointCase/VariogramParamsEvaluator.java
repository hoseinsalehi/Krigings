/*
 * GNU GPL v3 License (GEOframe Krigings, version 0.9.9 of 2023-01-11).
 */
/* Source reconstructed from Kriging-2023-01-11-0.9.9.jar (CFR 0.152 decompilation, cleaned 2026; see UPSTREAM_DIFF.md). */
package krigingsPointCase;

import curvefitter.VariogramFitter;
import curvefitter.VariogramFunction;
import experimentalVariogram.ExperimentalVariogram;
import java.util.ArrayList;
import java.util.Map;
import java.util.Set;
import krigingsPointCase.Utility;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import theoreticalVariogram.SimpleModelFactory;

public class VariogramParamsEvaluator {
    private static final double TRESHOLD = 0.5;
    @In
    public String pSemivariogramType;
    @In
    public ExperimentalVariogram expVar;
    @In
    public double[] x;
    @In
    public double[] y;
    @Out
    public double sill;
    @Out
    public double nugget;
    @Out
    public double range;
    @Out
    public double rmse;
    @Out
    public boolean isFitGood;
    @Out
    public double relError;
    @Out
    public String outSemivariogramType;

    @Execute
    public void proces() {
        try {
            ArrayList<WeightedObservedPoint> points = new ArrayList<WeightedObservedPoint>();
            if (this.expVar != null) {
                this.expVar.process();
                Set<Map.Entry<Integer, double[]>> entrySet = this.expVar.outDistances.entrySet();
                for (Map.Entry<Integer, double[]> entry : entrySet) {
                    Integer ID = entry.getKey();
                    double x = this.expVar.outDistances.get(ID)[0];
                    double y = this.expVar.outExperimentalVariogram.get(ID)[0];
                    if (y == -9999.0) continue;
                    WeightedObservedPoint point = new WeightedObservedPoint(1.0, x, y);
                    points.add(point);
                }
            } else if (this.x != null && this.y != null) {
                for (int i = 0; i < this.x.length; ++i) {
                    double distance = this.x[i];
                    double variance = this.y[i];
                    WeightedObservedPoint point = new WeightedObservedPoint(1.0, distance, variance);
                    points.add(point);
                }
            }
            String[] variogramType = this.pSemivariogramType != null ? new String[]{this.pSemivariogramType} : Utility.availableTheorethicalVariogra;
            this.performEvaluation(variogramType, points);
        }
        catch (Exception e) {
            this.isFitGood = false;
        }
    }

    private void performEvaluation(String[] variogramType, ArrayList<WeightedObservedPoint> points) {
        this.relError = Double.MAX_VALUE;
        for (int j = 0; j < variogramType.length; ++j) {
            try {
                VariogramFunction variogramFunction = new VariogramFunction(variogramType[j]);
                VariogramFitter fitter = new VariogramFitter(variogramFunction);
                ArrayList<WeightedObservedPoint> filteresPoints = variogramFunction.filterPoint(points);
                double[] coeffs = fitter.fit(filteresPoints);
                double distance = 0.0;
                int count = 0;
                for (int i = 0; i < points.size(); ++i) {
                    double actualValue = points.get(i).getY();
                    if (actualValue == 0.0) continue;
                    distance += (SimpleModelFactory.createModel(variogramType[j], points.get(i).getX(), coeffs[0], coeffs[1], coeffs[2]).computeSemivariance() - actualValue) / actualValue;
                    ++count;
                }
                double tmpError = Math.abs(distance / (double)points.size());
                System.out.println(variogramType[j] + " errore " + tmpError);
                if (j != 0 && !(tmpError < this.relError)) continue;
                this.relError = Math.abs(distance / (double)points.size());
                this.sill = coeffs[0];
                this.range = coeffs[1];
                this.nugget = coeffs[2];
                this.rmse = fitter.getRMS();
                this.isFitGood = this.relError < 0.5;
                this.outSemivariogramType = variogramType[j];
                continue;
            }
            catch (Exception e) {
                System.out.println(variogramType[j]);
                System.out.println(e.getMessage());
            }
        }
    }
}

