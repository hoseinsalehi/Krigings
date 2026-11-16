/*
 * GNU GPL v3 License (GEOframe Krigings, version 0.9.9 of 2023-01-11).
 */
/* Source reconstructed from Kriging-2023-01-11-0.9.9.jar (CFR 0.152 decompilation, cleaned 2026; see UPSTREAM_DIFF.md). */
package curvefitter;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.fitting.leastsquares.ParameterValidator;
import org.apache.commons.math3.linear.RealVector;

public class KrigingParamValidator
implements ParameterValidator {
    final double[] lowerBounds_;

    public KrigingParamValidator(double[] lowerBounds) {
        if (lowerBounds.length != 3) {
            throw new DimensionMismatchException(lowerBounds.length, 3);
        }
        this.lowerBounds_ = lowerBounds;
    }

    public RealVector validate(RealVector params) {
        if (params.getDimension() != 3) {
            throw new DimensionMismatchException(params.getDimension(), 3);
        }
        RealVector n = params.copy();
        for (int i = 0; i < 3; ++i) {
            if (!(n.getEntry(i) < this.lowerBounds_[i])) continue;
            n.setEntry(i, this.lowerBounds_[i]);
        }
        return n;
    }
}

