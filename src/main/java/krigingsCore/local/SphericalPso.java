package krigingsCore.local;

import utilities.random.*;

/**
 * Batched particle-swarm fit of a spherical variogram (nugget, partial sill, range) to empirical variograms,
 * identical to Python {@code po_numba.pso_spherical}: same bounds, weights, update rule, random-stream order
 * (U(B,P,3) random, N(B,P,3) normal(0,0.06), then per iteration r1(B,P,3), r2(B,P,3)) and numpy pairwise
 * summation of the 12-term objective.
 */
public final class SphericalPso {
    public static final class Settings {
        public int particles = 24, iterations = 35; public double inertia = 0.72, c1 = 1.45, c2 = 1.45;
        public boolean pairCountWeights = true;
    }

    private SphericalPso() { }

    /**
     * Fit B variograms with one random stream (rng). Output params[b*3 + {0,1,2}] = nugget, partial sill, range;
     * sse[b]; boundHit[b] = range >= 0.975 * upper range bound.
     */
    public static void fit(Variogram.Empirical[] ev, int B, Rng rng, Settings s,
                           double[] params, double[] sse, boolean[] boundHit) {
        int P = s.particles, L = ev[0].lag.length, D = 3;
        double[] lo = new double[B * D], hi = new double[B * D];
        double[] h = new double[B * L], g0 = new double[B * L], ww = new double[B * L];
        for (int b = 0; b < B; b++) {
            Variogram.Empirical e = ev[b];
            double maxg = Double.NaN, minh = Double.NaN, maxh = Double.NaN;
            for (int l = 0; l < L; l++) if (e.mask[l]) {
                if (Double.isNaN(maxg) || e.gam[l] > maxg) maxg = e.gam[l];
                if (Double.isNaN(minh) || e.lag[l] < minh) minh = e.lag[l];
                if (Double.isNaN(maxh) || e.lag[l] > maxh) maxh = e.lag[l];
            }
            maxg = NpMath.maximum(maxg, 1e-8); minh = NpMath.maximum(minh, 1e-4); maxh = NpMath.maximum(maxh, minh * 1.01);
            lo[b * D] = 0.0; lo[b * D + 1] = 1e-12; lo[b * D + 2] = NpMath.maximum(0.25 * minh, 1e-4);
            hi[b * D] = 1.5 * maxg; hi[b * D + 1] = 3 * maxg; hi[b * D + 2] = 2 * maxh;
            double tot = 0.0;
            if (s.pairCountWeights) { tot = NpMath.pairwiseSum(e.cnt, 0, L); }
            else { for (int l = 0; l < L; l++) if (e.mask[l]) tot += 1.0; }
            double den = Math.max(tot, 1);
            for (int l = 0; l < L; l++) {
                ww[b * L + l] = (s.pairCountWeights ? e.cnt[l] : (e.mask[l] ? 1.0 : 0.0)) / den;
                g0[b * L + l] = e.mask[l] ? e.gam[l] : 0.0;
                h[b * L + l] = e.mask[l] ? e.lag[l] : 1.0;
            }
        }
        int n = B * P * D;
        double[] pos = new double[n], vel = new double[n], pb = new double[n], pbs = new double[B * P], gb = new double[B * D], gbs = new double[B];
        rng.fill(pos);                                                       // U
        for (int i = 0; i < n; i++) vel[i] = rng.normal(0, .06);            // N
        double[] buf = new double[L];
        for (int b = 0; b < B; b++) {
            for (int p = 0; p < P; p++) for (int k = 0; k < D; k++) {
                int i = (b * P + p) * D + k; double span = hi[b * D + k] - lo[b * D + k];
                pos[i] = lo[b * D + k] + pos[i] * span; vel[i] = vel[i] * span; pb[i] = pos[i];
            }
            int gi = 0;
            for (int p = 0; p < P; p++) {
                int i = (b * P + p) * D;
                pbs[b * P + p] = obj(pos[i], pos[i + 1], pos[i + 2], h, g0, ww, b * L, L, buf);
                if (pbs[b * P + p] < pbs[b * P + gi]) gi = p;
            }
            for (int k = 0; k < D; k++) gb[b * D + k] = pos[(b * P + gi) * D + k];
            gbs[b] = pbs[b * P + gi];
        }
        double[] r1 = new double[n], r2 = new double[n];
        for (int t = 0; t < s.iterations; t++) {
            rng.fill(r1);
            rng.fill(r2);
            for (int b = 0; b < B; b++) {
                for (int p = 0; p < P; p++) for (int k = 0; k < D; k++) {
                    int i = (b * P + p) * D + k;
                    double v = s.inertia * vel[i];
                    v = v + (s.c1 * r1[i]) * (pb[i] - pos[i]);
                    v = v + (s.c2 * r2[i]) * (gb[b * D + k] - pos[i]);
                    vel[i] = v;
                    double x = pos[i] + v;
                    if (x < lo[b * D + k]) x = lo[b * D + k];
                    if (x > hi[b * D + k]) x = hi[b * D + k];
                    pos[i] = x;
                }
                for (int p = 0; p < P; p++) {
                    int i = (b * P + p) * D;
                    double sc = obj(pos[i], pos[i + 1], pos[i + 2], h, g0, ww, b * L, L, buf);
                    if (sc < pbs[b * P + p]) { pbs[b * P + p] = sc; pb[i] = pos[i]; pb[i + 1] = pos[i + 1]; pb[i + 2] = pos[i + 2]; }
                }
                int gi = 0;
                for (int p = 1; p < P; p++) if (pbs[b * P + p] < pbs[b * P + gi]) gi = p;
                if (pbs[b * P + gi] < gbs[b]) {
                    gbs[b] = pbs[b * P + gi];
                    for (int k = 0; k < D; k++) gb[b * D + k] = pb[(b * P + gi) * D + k];
                }
            }
        }
        for (int b = 0; b < B; b++) {
            for (int k = 0; k < D; k++) params[b * D + k] = gb[b * D + k];
            sse[b] = gbs[b];
            boundHit[b] = gb[b * D + 2] >= 0.975 * hi[b * D + 2];
        }
    }

    private static double obj(double x0, double x1, double x2, double[] h, double[] g0, double[] ww, int off, int L, double[] buf) {
        if (L == 12) return obj12(x0, x1, x2, h, g0, ww, off);
        double a = Math.max(x2, 1e-12);
        for (int l = 0; l < L; l++) {
            double r = h[off + l] / a; double rr = r < 1.0 ? r : 1.0;
            double v = rr * rr; v = v * rr; v = v * -0.5; v = v + 1.5 * rr;
            v = v * x1; v = v + x0; v = v - g0[off + l]; v = v * v; v = v * ww[off + l];
            buf[l] = v;
        }
        return NpMath.pairwiseSum(buf, 0, L);
    }

    private static double term(double x0, double x1, double a, double h, double g0, double ww) {
        double r = h / a; double rr = r < 1.0 ? r : 1.0;
        double v = rr * rr; v = v * rr; v = v * -0.5; v = v + 1.5 * rr;
        v = v * x1; v = v + x0; v = v - g0; v = v * v; return v * ww;
    }

    /** 12-bin objective, bitwise equal to obj() with numpy pairwise summation (n = 12: r0..r7, then + b8 + b9 + b10 + b11). */
    private static double obj12(double x0, double x1, double x2, double[] h, double[] g0, double[] ww, int o) {
        double a = Math.max(x2, 1e-12);
        double b0 = term(x0, x1, a, h[o], g0[o], ww[o]), b1 = term(x0, x1, a, h[o + 1], g0[o + 1], ww[o + 1]);
        double b2 = term(x0, x1, a, h[o + 2], g0[o + 2], ww[o + 2]), b3 = term(x0, x1, a, h[o + 3], g0[o + 3], ww[o + 3]);
        double b4 = term(x0, x1, a, h[o + 4], g0[o + 4], ww[o + 4]), b5 = term(x0, x1, a, h[o + 5], g0[o + 5], ww[o + 5]);
        double b6 = term(x0, x1, a, h[o + 6], g0[o + 6], ww[o + 6]), b7 = term(x0, x1, a, h[o + 7], g0[o + 7], ww[o + 7]);
        double s = ((b0 + b1) + (b2 + b3)) + ((b4 + b5) + (b6 + b7));
        s += term(x0, x1, a, h[o + 8], g0[o + 8], ww[o + 8]); s += term(x0, x1, a, h[o + 9], g0[o + 9], ww[o + 9]);
        s += term(x0, x1, a, h[o + 10], g0[o + 10], ww[o + 10]); s += term(x0, x1, a, h[o + 11], g0[o + 11], ww[o + 11]);
        return s;
    }
}
