package krigingsCore.local;


/**
 * Floating-point reductions in the exact operation order of numpy, so that values feeding the PSO fit are
 * bitwise identical to the Python engine.
 */
public final class NpMath {
    private NpMath() { }

    /** numpy pairwise summation of a[off .. off+n) (contiguous reduction, float64, identity-initialised). */
    public static double pairwiseSum(double[] a, int off, int n) {
        if (n < 8) {
            double s = 0.0;               // numpy: res = 0.; res += a[i]  (keeps -0.0 semantics of numpy)
            for (int i = 0; i < n; i++) s += a[off + i];
            return s;
        } else if (n <= 128) {
            double r0 = a[off], r1 = a[off + 1], r2 = a[off + 2], r3 = a[off + 3];
            double r4 = a[off + 4], r5 = a[off + 5], r6 = a[off + 6], r7 = a[off + 7];
            int i = 8, lim = n - (n % 8);
            for (; i < lim; i += 8) {
                r0 += a[off + i]; r1 += a[off + i + 1]; r2 += a[off + i + 2]; r3 += a[off + i + 3];
                r4 += a[off + i + 4]; r5 += a[off + i + 5]; r6 += a[off + i + 6]; r7 += a[off + i + 7];
            }
            double s = ((r0 + r1) + (r2 + r3)) + ((r4 + r5) + (r6 + r7));
            for (; i < n; i++) s += a[off + i];
            return s;
        } else {
            int n2 = n / 2;
            n2 -= n2 % 8;
            return pairwiseSum(a, off, n2) + pairwiseSum(a, off + n2, n - n2);
        }
    }

    public static double pairwiseSum(double[] a) { return pairwiseSum(a, 0, a.length); }

    /** numpy pairwise summation in float32 arithmetic. */
    public static float pairwiseSumF(float[] a, int off, int n) {
        if (n < 8) {
            float s = 0.0f;
            for (int i = 0; i < n; i++) s += a[off + i];
            return s;
        } else if (n <= 128) {
            float r0 = a[off], r1 = a[off + 1], r2 = a[off + 2], r3 = a[off + 3];
            float r4 = a[off + 4], r5 = a[off + 5], r6 = a[off + 6], r7 = a[off + 7];
            int i = 8, lim = n - (n % 8);
            for (; i < lim; i += 8) {
                r0 += a[off + i]; r1 += a[off + i + 1]; r2 += a[off + i + 2]; r3 += a[off + i + 3];
                r4 += a[off + i + 4]; r5 += a[off + i + 5]; r6 += a[off + i + 6]; r7 += a[off + i + 7];
            }
            float s = ((r0 + r1) + (r2 + r3)) + ((r4 + r5) + (r6 + r7));
            for (; i < n; i++) s += a[off + i];
            return s;
        } else {
            int n2 = n / 2;
            n2 -= n2 % 8;
            return pairwiseSumF(a, off, n2) + pairwiseSumF(a, off + n2, n - n2);
        }
    }

    /** numpy.nanmax semantics (NaN if no finite... i.e. all NaN). */
    public static double nanMax(double[] a, int off, int n) {
        double m = Double.NaN;
        for (int i = 0; i < n; i++) { double v = a[off + i]; if (!Double.isNaN(v) && (Double.isNaN(m) || v > m)) m = v; }
        return m;
    }

    public static double nanMin(double[] a, int off, int n) {
        double m = Double.NaN;
        for (int i = 0; i < n; i++) { double v = a[off + i]; if (!Double.isNaN(v) && (Double.isNaN(m) || v < m)) m = v; }
        return m;
    }

    /** numpy.rint (round half to even). */
    public static double rint(double x) { return Math.rint(x); }

    /** numpy.maximum: NaN propagates. */
    public static double maximum(double a, double b) { return (Double.isNaN(a) || Double.isNaN(b)) ? Double.NaN : (a >= b ? a : b); }

    /** numpy.minimum: NaN propagates. */
    public static double minimum(double a, double b) { return (Double.isNaN(a) || Double.isNaN(b)) ? Double.NaN : (a <= b ? a : b); }

    /** numpy.clip(x, lo, hi) = minimum(maximum(x, lo), hi). */
    public static double clip(double x, double lo, double hi) { return minimum(maximum(x, lo), hi); }
}
