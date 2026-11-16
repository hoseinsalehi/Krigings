/*
 * GNU GPL v3 License (GEOframe Krigings, version 0.9.9 of 2023-01-11).
 */
/* Source reconstructed from Kriging-2023-01-11-0.9.9.jar (CFR 0.152 decompilation, cleaned 2026; see UPSTREAM_DIFF.md). */
package curvefitter;

import curvefitter.KrigingParamValidator;
import java.util.Arrays;
import java.util.Collection;
import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.fitting.AbstractCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.fitting.leastsquares.EvaluationRmsChecker;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.ParameterValidator;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.ConvergenceChecker;

public class VariogramFitter
extends AbstractCurveFitter {
    private ParametricUnivariateFunction myFunction = null;
    private LeastSquaresProblem problem = null;
    private int nEvaluation = 0;
    private int nIteration = 0;
    private double rms = 0.0;
    private RealVector residuals;
    double tol = Math.pow(0.1, 9.0);

    private VariogramFitter() {
    }

    public VariogramFitter(ParametricUnivariateFunction myFun) {
        this.myFunction = myFun;
    }

    protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> points) {
        int len = points.size();
        double[] target = new double[len];
        double[] x = new double[len];
        int i = 0;
        for (WeightedObservedPoint point : points) {
            target[i] = point.getY();
            x[i] = point.getX();
            ++i;
        }
        double initialSill = 0.9 * Arrays.stream(target).max().getAsDouble();
        double initialNugget = 0.0;
        double[] semivarianceMinusSill = Arrays.stream(target).map(v -> Math.abs(v - initialSill)).toArray();
        int index = 0;
        double min = semivarianceMinusSill[index];
        for (int j = 1; j < semivarianceMinusSill.length; ++j) {
            if (!(semivarianceMinusSill[j] <= min)) continue;
            min = semivarianceMinusSill[j];
            index = j;
        }
        double initialRange = x[index];
        double[] initialGuess = new double[]{initialSill, initialRange, initialNugget};
        AbstractCurveFitter.TheoreticalValuesFunction model = new AbstractCurveFitter.TheoreticalValuesFunction(this.myFunction, points);
        this.problem = new LeastSquaresBuilder().maxEvaluations(1000).maxIterations(1000).start(initialGuess).target(target).model(model.getModelFunction(), model.getModelFunctionJacobian()).checker((ConvergenceChecker)new EvaluationRmsChecker(this.tol, this.tol)).parameterValidator((ParameterValidator)new KrigingParamValidator(new double[]{1.0E-6, 1.0E-6, 0.0})).lazyEvaluation(false).build();
        return this.problem;
    }

    public double[] fit(Collection<WeightedObservedPoint> points) {
        double[] optimalValues = null;
        LevenbergMarquardtOptimizer lmo = new LevenbergMarquardtOptimizer();
        LeastSquaresOptimizer.Optimum lsoo = lmo.withCostRelativeTolerance(1.0E-12).withParameterRelativeTolerance(1.0E-12).optimize(this.getProblem(points));
        optimalValues = lsoo.getPoint().toArray();
        this.nIteration = lsoo.getIterations();
        this.nEvaluation = lsoo.getEvaluations();
        this.rms = lsoo.getRMS();
        this.residuals = lsoo.getResiduals();
        return optimalValues;
    }

    public int getNumberOfIteration() {
        return this.nIteration;
    }

    public int getNumberOfEvaluation() {
        return this.nEvaluation;
    }

    public double getRMS() {
        return this.rms;
    }

    public RealVector getResiduals() {
        return this.residuals;
    }
}

