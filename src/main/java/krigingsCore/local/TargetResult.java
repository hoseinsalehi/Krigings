package krigingsCore.local;


import java.util.Arrays;

/** Per-target outputs of {@link DayKriging} (Python estimate_targets dict). */
public final class TargetResult {
    public final int n;
    public final double[] est, kvar, nug, psill, range, sse, b1, r2, span, nNeg, neff, rpos12, rpos23, rpos13, rabs90, sumAbsW, d1, d7, d50;
    /** components of est (est = resid + trend): kriged residual and lapse trend at the target (0 without lapse) */
    public final double[] resid, trend;
    public final int[] flags;

    public TargetResult(int n) {
        this.n = n;
        est = nan(n); kvar = nan(n); nug = nan(n); psill = nan(n); range = nan(n); sse = nan(n); b1 = nan(n); r2 = nan(n);
        span = nan(n); nNeg = nan(n); neff = nan(n); rpos12 = nan(n); rpos23 = nan(n); rpos13 = nan(n); rabs90 = nan(n);
        sumAbsW = nan(n); d1 = nan(n); d7 = nan(n); d50 = nan(n); flags = new int[n];
        resid = nan(n); trend = nan(n);
    }

    private static double[] nan(int n) { double[] a = new double[n]; Arrays.fill(a, Double.NaN); return a; }
}
