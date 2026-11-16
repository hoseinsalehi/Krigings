package krigingsCore.local;


/**
 * Ordinary kriging with k neighbours and a spherical variogram (Python {@code po_core.ok_solve}) plus the
 * weight diagnostics of {@code po_core.weight_diagnostics}. Not thread-safe (work arrays); one per thread.
 */
public final class OrdinaryKriging {
    private final int k;
    private final double[] A, rhs, sol; private final int[] piv;
    public final double[] w, g0;
    public double mu, var; public boolean ok;
    // diagnostics
    public int nNeg; public double neff, rpos12, rpos23, rpos13, rabs90, sumAbsW;

    public OrdinaryKriging(int k) {
        this.k = k; int m = k + 1;
        A = new double[m * m]; rhs = new double[m]; sol = new double[m]; piv = new int[m];
        w = new double[k]; g0 = new double[k];
    }

    /**
     * @param nx,ny coordinates of the k neighbours [m]; dk distances target-neighbour [km]
     */
    public void solve(double[] nx, double[] ny, double[] dk, double nug, double c, double a) {
        int m = k + 1;
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                double dx = nx[i] - nx[j], dy = ny[i] - ny[j];
                double dd = Math.sqrt(dx * dx + dy * dy) / 1000.0;
                A[i * m + j] = Variogram.gamma(dd, nug, c, a);
            }
            A[i * m + k] = 1.0; A[k * m + i] = 1.0;
            g0[i] = Variogram.gamma(dk[i], nug, c, a); rhs[i] = g0[i];
        }
        A[k * m + k] = 0.0; rhs[k] = 1.0;
        ok = luSolve(A, rhs, sol, piv, m);
        double s = 0.0;
        for (int i = 0; i < k; i++) { w[i] = ok ? sol[i] : Double.NaN; }
        mu = ok ? sol[k] : Double.NaN;
        for (int i = 0; i < k; i++) s += w[i] * g0[i];
        var = NpMath.maximum(s + mu, 0.0);
        if (ok) for (int i = 0; i <= k; i++) if (!Variogram.isFinite(sol[i])) { ok = false; break; }
    }

    /** Gaussian elimination with partial pivoting (LAPACK dgesv semantics: exact zero pivot -> singular). */
    static boolean luSolve(double[] A, double[] b, double[] x, int[] piv, int n) {
        for (int i = 0; i < n; i++) x[i] = b[i];
        for (int c = 0; c < n; c++) {
            int p = c; double best = Math.abs(A[c * n + c]);
            for (int r = c + 1; r < n; r++) { double v = Math.abs(A[r * n + c]); if (v > best) { best = v; p = r; } }
            if (best == 0.0) return false;
            if (p != c) {
                for (int j = 0; j < n; j++) { double t = A[c * n + j]; A[c * n + j] = A[p * n + j]; A[p * n + j] = t; }
                double t = x[c]; x[c] = x[p]; x[p] = t;
            }
            double inv = A[c * n + c];
            for (int r = c + 1; r < n; r++) {
                double f = A[r * n + c] / inv;
                if (f != 0.0) {
                    for (int j = c + 1; j < n; j++) A[r * n + j] -= f * A[c * n + j];
                    x[r] -= f * x[c];
                }
            }
        }
        for (int r = n - 1; r >= 0; r--) {
            double s = x[r];
            for (int j = r + 1; j < n; j++) s -= A[r * n + j] * x[j];
            x[r] = s / A[r * n + r];
        }
        return true;
    }

    /** Estimate sum(w_i y_i) (sequential, numpy sum of 7 terms). */
    public double estimate(double[] yk) {
        double s = 0.0;
        for (int i = 0; i < k; i++) s += w[i] * yk[i];
        return s;
    }

    /** Weight diagnostics; d = neighbour distances [km], ascending. */
    public void diagnostics(double[] d, double negTol) {
        double s = 0.0; nNeg = 0;
        for (int i = 0; i < k; i++) { s += Math.abs(w[i]); if (w[i] < negTol) nNeg++; }
        sumAbsW = s;
        double pp = 0.0, mx = Double.NaN;
        for (int i = 0; i < k; i++) { double p = Math.abs(w[i]) / s; pp += p * p; }
        neff = 1.0 / pp;
        for (int i = 0; i < k; i++) if (w[i] > 0 && (Double.isNaN(mx) || w[i] > mx)) mx = w[i];
        rpos12 = rmax(d, 0.5 * mx); rpos23 = rmax(d, (2.0 / 3.0) * mx); rpos13 = rmax(d, (1.0 / 3.0) * mx);
        double cum = 0.0; int i90 = 0;
        for (int i = 0; i < k; i++) { cum += Math.abs(w[i]) / s; if (cum >= 0.9) { i90 = i; break; } }
        rabs90 = d[i90];
    }

    private double rmax(double[] d, double thr) {
        double r = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < k; i++) if (w[i] >= thr && d[i] > r) r = d[i];
        return r;
    }
}
