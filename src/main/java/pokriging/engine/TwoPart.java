package pokriging.engine;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** Occurrence calibration of the two-part method (Python s3_run._thresholds, _calibrate, _mix_std, enc_tp/dec_tp). */
public final class TwoPart {
    private TwoPart() { }

    public static final int PI_FILL = 65535, AMT_FILL = -32768;

    /** enc_tp("pi"): u2, scale 1e-4. Returned as int 0..65535. */
    public static int encPi(double pi) {
        double q = Math.rint(pi / 1e-4);
        if (Double.isNaN(q) || Double.isInfinite(q)) return PI_FILL;
        return (int) Math.max(0, Math.min(65534, q));
    }

    /** enc_tp("amount"/"samount"): i2, scale 0.1. */
    public static short encAmt(double a) {
        double q = Math.rint(a / 0.1);
        if (Double.isNaN(q) || Double.isInfinite(q)) return (short) AMT_FILL;
        return (short) Math.max(-32767, Math.min(32767, q));
    }

    /** dec_tp("pi") followed by .astype(float32). */
    public static float decPi(int q) { return q == PI_FILL ? Float.NaN : (float) (q * 1e-4); }

    public static double decAmt(short q) { return q == AMT_FILL ? Double.NaN : q * 0.1; }

    /** hurdle standard deviation in float64 (grid). */
    public static double mixStd(double pi, double amt, double samt) {
        double p = krigingsCore.local.NpMath.clip(pi, 0, 1);
        double v = p * (samt * samt) + (p * (1 - p)) * (amt * amt);
        return Math.sqrt(krigingsCore.local.NpMath.maximum(v, 0.0));
    }

    /** hurdle standard deviation in float32 arithmetic (LOO / points tables are float32 in Python). */
    public static float mixStdF(float pi, float amt, float samt) {
        float p = Float.isNaN(pi) ? Float.NaN : Math.min(Math.max(pi, 0f), 1f);
        float v = p * (samt * samt) + (p * (1f - p)) * (amt * amt);
        if (Float.isNaN(v)) return Float.NaN;
        return (float) Math.sqrt(Math.max(v, 0f));
    }

    /**
     * _thresholds for one target: pi values (NaN allowed) and wanted wet days k (double) -> threshold t such that
     * pi >= t marks about k days wet (t = +inf if k = 0; t <= 0 -> 1e-9).
     */
    public static double threshold(float[] pi, int n, double k) {
        float[] f = new float[n]; int m = 0;
        for (int i = 0; i < n; i++) if (!Float.isNaN(pi[i])) f[m++] = pi[i];
        long kk = (long) Math.rint(k);
        kk = Math.max(0, Math.min(m, kk));
        double t = Double.POSITIVE_INFINITY;
        if (kk > 0) { Arrays.sort(f, 0, m); t = f[m - (int) kk]; }
        return t <= 0 ? 1e-9 : t;
    }

    /**
     * _calibrate for a long table (LOO stations or points): rows with month (1-12), dateKey, group key, pi (float32),
     * p_target (float32), in table order. Returns wet flags and the mode ("monthly" | "pooled").
     */
    public static String calibrate(int[] month, int[] dateKey, int[] key, float[] pi, float[] pTarget, int calibMinDays, boolean[] wet) {
        int n = pi.length;
        java.util.Set<Integer>[] dm = new java.util.HashSet[12];
        for (int i = 0; i < 12; i++) dm[i] = new java.util.HashSet<>();
        for (int i = 0; i < n; i++) dm[month[i] - 1].add(dateKey[i]);
        int nmin = Integer.MAX_VALUE; for (int i = 0; i < 12; i++) nmin = Math.min(nmin, dm[i].size());
        boolean monthly = nmin >= calibMinDays;
        Map<Long, List<Integer>> groups = new HashMap<>();
        for (int i = 0; i < n; i++) {
            long g = monthly ? ((long) key[i] * 16 + month[i]) : key[i];
            groups.computeIfAbsent(g, z -> new ArrayList<>()).add(i);
        }
        for (List<Integer> rows : groups.values()) {
            // k = round(sum of p_target over rows with finite pi); pandas sums float32 with compensation, then rounds half-even
            double s = 0, comp = 0;
            List<Integer> fin = new ArrayList<>();
            for (int i : rows) if (!Float.isNaN(pi[i])) {
                fin.add(i);
                double y = pTarget[i] - comp, t = s + y; comp = (t - s) - y; s = t;
            }
            double k = Math.rint((float) s);
            fin.sort((a, b) -> { int c = pi[b] > pi[a] ? 1 : (pi[b] < pi[a] ? -1 : 0); return c != 0 ? c : Integer.compare(a, b); });   // rank 'first', descending
            for (int r = 0; r < fin.size(); r++) { int i = fin.get(r); wet[i] = (r + 1) <= k && pi[i] > 0; }
            for (int i : rows) if (Float.isNaN(pi[i])) wet[i] = false;
        }
        return monthly ? "monthly" : "pooled";
    }
}
