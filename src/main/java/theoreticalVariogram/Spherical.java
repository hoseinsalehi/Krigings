/*
 * GNU GPL v3 License
 *
 * Copyright 2016 Marialaura Bancheri
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
package theoreticalVariogram;

import theoreticalVariogram.Model;

public class Spherical
implements Model {
    double dist;
    double sill;
    double range;
    double nug;
    boolean isOk = false;

    public Spherical(double dist, double sill, double range, double nug) {
        this.dist = dist;
        this.sill = sill;
        this.range = range;
        this.nug = nug;
        this.isOk = nug >= 0.0 && sill >= 0.0 && range >= 0.0;
    }

    @Override
    public double computeSemivariance() {
        double result = Double.MAX_VALUE;
        if (this.isOk) {
            double hr = this.dist / this.range;
            if (this.dist < this.range) {
                result = this.nug + this.sill * hr * (1.5 - 0.5 * hr * hr);
            } else if (this.dist >= this.range) {
                result = this.sill + this.nug;
            }
        }
        return result;
    }

    @Override
    public double[] computeGradient() {
        double[] gradient = new double[]{Double.NaN, Double.NaN, Double.NaN};
        if (this.isOk) {
            double hr = this.dist / this.range;
            if (this.dist < this.range) {
                gradient = new double[]{hr * (1.5 - 0.5 * hr * hr), 1.5 * this.sill * hr / this.range * (1.0 - hr * hr), 1.0};
            } else if (this.dist >= this.range) {
                gradient = new double[]{1.0, 0.0, 1.0};
            }
        }
        return gradient;
    }

    /**
     * Spherical semivariance in the operation order of the Python Po engine (po_core.gamma_sph):
     * nug + c * (1.5 rr - 0.5 rr^3), rr = min(h / max(a, 1e-12), 1); gamma(h &lt; 1e-9) = 0. Algebraically equal to
     * {@link #computeSemivariance()} (sill = partial sill c); kept separate so that the local PSO kriging is bitwise
     * reproducible against the Python reference. Added 2026 (see UPSTREAM_DIFF.md).
     */
    public static double semivariancePoOrder(double h, double nug, double c, double a) {
        if (h < 1e-9) return 0.0;
        double r = h / Math.max(a, 1e-12);
        double rr = r < 1.0 ? r : (Double.isNaN(r) ? r : 1.0);
        return nug + c * (1.5 * rr - 0.5 * StrictMath.pow(rr, 3));
    }
}
