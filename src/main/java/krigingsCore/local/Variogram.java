package krigingsCore.local;


/**
 * Empirical semivariogram with equal-width lag bins from 0 to the largest pair distance
 * (Python {@code po_numba._empirical} + {@code empirical_batch}, operation by operation).
 */
public final class Variogram {
    /** Per-target empirical variogram (nLags bins). */
    public static final class Empirical {
        public final double[] lag, gam, cnt; public final boolean[] mask; public boolean valid; public long npairs;
        public Empirical(int nLags) { lag = new double[nLags]; gam = new double[nLags]; cnt = new double[nLags]; mask = new boolean[nLags]; }
    }

    private Variogram() { }

    /**
     * Empirical variogram of values yv[0..n) located at stations nb[0..n), with station distances sdmat
     * (float32 km, row-major S x S).
     * @param minPairs minimum pairs per lag (MIN_PAIRS_PER_LAG)
     */
    public static void empirical(double[] yv, int[] nb, int n, float[] sdmat, int S, int nLags, int minPairs, Empirical out) {
        double md = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            double yi = yv[i]; int row = nb[i] * S;
            for (int j = i + 1; j < n; j++) {
                double d = sdmat[row + nb[j]]; double df = yi - yv[j]; double sv = 0.5 * (df * df);
                if (isFinite(d) && isFinite(sv) && d > 0) { if (d > md) md = d; }
            }
        }
        accumulate(yv, nb, n, sdmat, S, nLags, minPairs, md, out);
    }

    // shared second pass
    private static void accumulate(double[] yv, int[] nb, int n, float[] sdmat, int S, int nLags, int minPairs, double md, Empirical out) {
        if (!(isFinite(md) && md > 0)) md = 1.0;
        double stop = md * (1 + 1e-10), step = stop / nLags;
        double[] lag = out.lag, gam = out.gam, cnt = out.cnt;
        java.util.Arrays.fill(lag, 0); java.util.Arrays.fill(gam, 0); java.util.Arrays.fill(cnt, 0);
        double[] edges = new double[nLags + 1];
        for (int e = 1; e <= nLags; e++) edges[e] = (e == nLags) ? stop : e * step;
        long np = 0;
        for (int i = 0; i < n; i++) {
            double yi = yv[i]; int row = nb[i] * S;
            for (int j = i + 1; j < n; j++) {
                double d = sdmat[row + nb[j]]; double df = yi - yv[j]; double sv = 0.5 * (df * df);
                if (!(isFinite(d) && isFinite(sv) && d > 0)) continue;
                np++;
                int k = 0;
                for (int e = 1; e <= nLags; e++) if (d >= edges[e]) k++;
                if (k > nLags - 1) k = nLags - 1;
                cnt[k] += 1.0; lag[k] += d; gam[k] += sv;
            }
        }
        finish(out, np, nLags, minPairs);
    }

    private static void finish(Empirical out, long np, int nLags, int minPairs) {
        int nm = 0; double gmax = Double.NaN;
        for (int l = 0; l < nLags; l++) {
            boolean m = out.cnt[l] >= minPairs;
            out.mask[l] = m;
            if (m) { out.lag[l] = out.lag[l] / out.cnt[l]; out.gam[l] = out.gam[l] / out.cnt[l]; nm++;
                     if (Double.isNaN(gmax) || out.gam[l] > gmax) gmax = out.gam[l]; }
            else { out.lag[l] = Double.NaN; out.gam[l] = Double.NaN; out.cnt[l] = 0.0; }
        }
        double g = Double.isNaN(gmax) ? 0.0 : gmax;
        out.npairs = np;
        out.valid = np >= Math.max(20, minPairs * 3) && nm >= 3 && g > 1e-12;
    }

    /** Neutralise an invalid row before PSO exactly as Python (bins 0-2: lag 1,2,3, gamma 1, count 1). */
    public static void neutralise(Empirical e) {
        for (int l = 0; l < e.mask.length; l++) e.mask[l] = false;
        for (int l = 0; l < 3; l++) { e.mask[l] = true; e.lag[l] = l + 1.0; e.gam[l] = 1.0; e.cnt[l] = 1.0; }
        // counts of bins >= 3 are kept as they are (they still enter the pair-count weights, as in Python)
    }

    static boolean isFinite(double v) { return !Double.isNaN(v) && !Double.isInfinite(v); }

    // ---------------------------------------------------------------- spherical model
    /** 1.5 r - 0.5 r^3 with r = min(r, 1), operation order of po_core.spherical_core (no out argument). */
    public static double sphericalCore(double r) {
        double rr = r < 1.0 ? r : (Double.isNaN(r) ? r : 1.0);
        return 1.5 * rr - 0.5 * StrictMath.pow(rr, 3);
    }

    /** gamma(h) = nug + c * sph(h / max(a, 1e-12)); gamma(h < 1e-9) = 0. */
    public static double gamma(double h, double nug, double c, double a) {
        return theoreticalVariogram.Spherical.semivariancePoOrder(h, nug, c, a);
    }
}
