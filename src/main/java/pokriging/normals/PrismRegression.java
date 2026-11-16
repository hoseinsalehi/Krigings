package pokriging.normals;

import krigingsCore.local.KdTree;
import krigingsCore.local.NpMath;

/**
 * PRISM-type local weighted regression of a monthly normal on elevation (Python s2_normals_grid.predict): K nearest
 * stations, bisquare distance weights, elevation-difference (and inactive facet/TPI) weights, slope bounded, target
 * elevation capped (precip, wetfreq) or above-network lapse substitution (temp).
 */
public final class PrismRegression {
    public enum Kind { PRECIP, WETFREQ, TEMP }

    public static final class Params {
        public int K; public double dz0, facetW = 1.0, tpi0 = 1e9, extrapDz = 1e9, extrapLapse = 0.0, slopeLo, slopeHi; public boolean useTpi;
    }

    /** Station set of one month. */
    public static final class Stations {
        final double[] x, y, z, tpi, val; final int[] facet; final KdTree tree;
        public Stations(double[] x, double[] y, double[] z, double[] tpi, int[] facet, double[] val) {
            this.x = x; this.y = y; this.z = z; this.tpi = tpi; this.facet = facet; this.val = val; tree = new KdTree(x, y);
        }
    }

    /**
     * Predict at M targets. self[t] = index (in the station set) of the target's own station or -1;
     * fixedSlope = per-target imposed slope or null. Returns {estimate[M], slope[M]}.
     */
    public static double[][] predict(Kind kind, Params p, double[] tx, double[] ty, double[] tz, double[] ttpi, int[] tfacet, int[] self,
                                     Stations st, double[] fixedSlope, int threads) throws InterruptedException {
        int M = tx.length, K = p.K, kq = Math.min(K + 1, st.x.length), P = kind == Kind.TEMP ? 3 : 2;
        int Kk = Math.min(K, kq);
        double tscale = p.useTpi ? 1.0 : 0.0;
        int[][] il = new int[M][]; double[][] w = new double[M][], xc1 = new double[M][], xc2 = new double[M][];
        double[][] A = new double[M][], b = new double[M][]; double[] zmax = new double[M], tzEff = new double[M];
        Thread[] th = new Thread[Math.max(1, threads)];
        final int nt = th.length;
        for (int ti = 0; ti < nt; ti++) {
            final int t0 = ti;
            th[ti] = new Thread(() -> {
                KdTree.Query q = new KdTree.Query(kq);
                for (int t = t0; t < M; t += nt) {
                    int got = st.tree.knn(tx[t], ty[t], kq, q);
                    int[] idx = new int[Kk]; double[] d = new double[Kk]; int c = 0;
                    for (int j = 0; j < got && c < Kk; j++) if (q.outIdx[j] != self[t]) { idx[c] = q.outIdx[j]; d[c] = q.outDist[j]; c++; }
                    for (int j = 0; j < got && c < Kk; j++) if (q.outIdx[j] == self[t]) { idx[c] = q.outIdx[j]; d[c] = q.outDist[j]; c++; }
                    double h = d[Kk - 1] * 1.01 + 1.0;
                    double[] ww = new double[Kk]; double zm = Double.NEGATIVE_INFINITY;
                    for (int j = 0; j < Kk; j++) {
                        double r = d[j] / h; double a = 1 - r * r; double wj = a * a;
                        double e = Math.abs(st.z[idx[j]] - tz[t]) / p.dz0; wj *= 1.0 / (1.0 + e * e);
                        boolean same = st.facet[idx[j]] == tfacet[t] || st.facet[idx[j]] == 8 || tfacet[t] == 8;
                        wj *= same ? 1.0 : p.facetW;
                        if (kind == Kind.TEMP) { double u = Math.abs(st.tpi[idx[j]] - ttpi[t]) / p.tpi0; wj *= 1.0 / (1.0 + u * u); }
                        ww[j] = wj; if (st.z[idx[j]] > zm) zm = st.z[idx[j]];
                    }
                    double ze = kind == Kind.TEMP ? tz[t] : NpMath.minimum(tz[t], zm + p.extrapDz);
                    double[] c1 = new double[Kk], c2 = new double[Kk];
                    double[] AA = new double[P * P], bb = new double[P];
                    for (int j = 0; j < Kk; j++) {
                        c1[j] = st.z[idx[j]] - ze;
                        c2[j] = st.tpi[idx[j]] * tscale - ttpi[t] * tscale;
                        double[] xr = P == 3 ? new double[]{1, c1[j], c2[j]} : new double[]{1, c1[j]};
                        for (int r = 0; r < P; r++) { bb[r] += xr[r] * ww[j] * st.val[idx[j]]; for (int s = 0; s < P; s++) AA[r * P + s] += xr[r] * ww[j] * xr[s]; }
                    }
                    il[t] = idx; w[t] = ww; xc1[t] = c1; xc2[t] = c2; A[t] = AA; b[t] = bb; zmax[t] = zm; tzEff[t] = ze;
                }
            });
            th[ti].start();
        }
        for (Thread x : th) x.join();
        double[] lam = new double[P];
        for (int j = 1; j < P; j++) { double s = 0; for (int t = 0; t < M; t++) s += A[t][j * P + j]; lam[j] = 1e-9 * Math.max(s / M, 1e-12); }
        double[] est = new double[M], slope = new double[M];
        for (int t = 0; t < M; t++) {
            double[] AA = A[t].clone(); for (int j = 0; j < P; j++) AA[j * P + j] += lam[j];
            double[] beta = solve(AA, b[t].clone(), P);
            double sl = fixedSlope == null ? NpMath.clip(beta[1], p.slopeLo, p.slopeHi) : fixedSlope[t];
            slope[t] = sl;
            int[] idx = il[t]; double[] ww = w[t];
            if (kind == Kind.TEMP) {
                double a00 = 0, a01 = 0, a11 = 0, r0 = 0, r1 = 0;
                for (int j = 0; j < idx.length; j++) {
                    double yr = st.val[idx[j]] - sl * xc1[t][j];
                    a00 += ww[j]; a01 += ww[j] * xc2[t][j]; a11 += xc2[t][j] * ww[j] * xc2[t][j]; r0 += ww[j] * yr; r1 += xc2[t][j] * ww[j] * yr;
                }
                double[] s2 = solve(new double[]{a00, a01, a01, a11 + lam[2]}, new double[]{r0, r1}, 2);
                double e = s2[0];
                if (p.extrapLapse < 0) { double up = Math.max(tz[t] - zmax[t], 0.0); e = e - sl * up + Math.min(sl, p.extrapLapse) * up; }
                est[t] = e;
            } else {
                double num = 0, den = 0;
                for (int j = 0; j < idx.length; j++) { num += (st.val[idx[j]] - sl * xc1[t][j]) * ww[j]; den += ww[j]; }
                est[t] = num / den;
            }
        }
        return new double[][]{est, slope};
    }

    static double[] solve(double[] A, double[] b, int n) {
        for (int c = 0; c < n; c++) {
            int p = c; for (int r = c + 1; r < n; r++) if (Math.abs(A[r * n + c]) > Math.abs(A[p * n + c])) p = r;
            if (A[p * n + c] == 0) { double[] x = new double[n]; java.util.Arrays.fill(x, Double.NaN); return x; }
            if (p != c) { for (int j = 0; j < n; j++) { double t = A[c * n + j]; A[c * n + j] = A[p * n + j]; A[p * n + j] = t; } double t = b[c]; b[c] = b[p]; b[p] = t; }
            for (int r = c + 1; r < n; r++) { double f = A[r * n + c] / A[c * n + c]; for (int j = c; j < n; j++) A[r * n + j] -= f * A[c * n + j]; b[r] -= f * b[c]; }
        }
        double[] x = new double[n];
        for (int r = n - 1; r >= 0; r--) { double s = b[r]; for (int j = r + 1; j < n; j++) s -= A[r * n + j] * x[j]; x[r] = s / A[r * n + r]; }
        return x;
    }
}
