/*
 * GNU GPL v3 License (GEOframe Krigings, version 0.9.9 of 2023-01-11).
 */
/* Source reconstructed from Kriging-2023-01-11-0.9.9.jar (CFR 0.152 decompilation, cleaned 2026; see UPSTREAM_DIFF.md). */
package curvefitter;

import java.util.ArrayList;
import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import theoreticalVariogram.Model;
import theoreticalVariogram.SimpleModelFactory;

public class VariogramFunction
implements ParametricUnivariateFunction {
    String type;

    public VariogramFunction(String type) {
        this.type = type;
    }

    public double value(double x, double ... parameters) {
        Model variogram = SimpleModelFactory.createModel(this.type, x, parameters[0], parameters[1], parameters[2]);
        return variogram.computeSemivariance();
    }

    public double[] gradient(double x, double ... parameters) {
        Model variogram = SimpleModelFactory.createModel(this.type, x, parameters[0], parameters[1], parameters[2]);
        return variogram.computeGradient();
    }

    public ArrayList<WeightedObservedPoint> filterPoint(ArrayList<WeightedObservedPoint> points) {
        ArrayList<WeightedObservedPoint> newPoint = new ArrayList<WeightedObservedPoint>();
        if (this.type == "logarithmic" || this.type == "power") {
            for (WeightedObservedPoint point : points) {
                if (point.getX() == 0.0) continue;
                newPoint.add(point);
            }
            return newPoint;
        }
        return points;
    }
}

