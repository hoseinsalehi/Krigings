package krigingsCore.local;


/**
 * Row-wise OLS T = b0 + b1 E on finite pairs (Python {@code po_core.lapse_batch}), with numpy pairwise sums.
 * Not thread-safe (work arrays).
 */
public final class LocalLapse {
    public double b0, b1, r2, span; public int n;
    private double[] z, t, dz, dt, zz, zt, tt;

    public LocalLapse(int nmax) { alloc(nmax); }

    private void alloc(int m) { z = new double[m]; t = new double[m]; dz = new double[m]; dt = new double[m]; zz = new double[m]; zt = new double[m]; tt = new double[m]; }

    /** Fit on E[0..m), T[0..m) (NaN = missing). Invalid fit (n < 3 or Szz <= 1e-8) -> all NaN. */
    public void fit(double[] E, double[] T, int m) {
        if (z.length < m) alloc(m);
        int cnt = 0; double mx = Double.NEGATIVE_INFINITY, mn = Double.POSITIVE_INFINITY;
        for (int i = 0; i < m; i++) {
            boolean ok = Variogram.isFinite(E[i]) && Variogram.isFinite(T[i]);
            z[i] = ok ? E[i] : 0.0; t[i] = ok ? T[i] : 0.0;
            if (ok) { cnt++; if (E[i] > mx) mx = E[i]; if (E[i] < mn) mn = E[i]; }
        }
        double zm = NpMath.pairwiseSum(z, 0, m) / cnt, tm = NpMath.pairwiseSum(t, 0, m) / cnt;
        for (int i = 0; i < m; i++) {
            boolean ok = Variogram.isFinite(E[i]) && Variogram.isFinite(T[i]);
            dz[i] = ok ? E[i] - zm : 0.0; dt[i] = ok ? T[i] - tm : 0.0;
            zz[i] = dz[i] * dz[i]; zt[i] = dz[i] * dt[i]; tt[i] = dt[i] * dt[i];
        }
        double szz = NpMath.pairwiseSum(zz, 0, m), szt = NpMath.pairwiseSum(zt, 0, m), stt = NpMath.pairwiseSum(tt, 0, m);
        double bb1 = szt / szz, bb0 = tm - bb1 * zm;
        double rr2 = stt > 0 ? 1 - (stt - bb1 * szt) / stt : Double.NaN;
        n = cnt;
        boolean good = cnt >= 3 && szz > 1e-8;
        b0 = good ? bb0 : Double.NaN; b1 = good ? bb1 : Double.NaN; r2 = good ? rr2 : Double.NaN; span = good ? mx - mn : Double.NaN;
    }
}
