package krigingsCore.local;

import utilities.random.*;

/**
 * Kriging of one day at a set of targets: Python {@code po_core.estimate_targets} and {@code global_variogram}.
 * <p>
 * A {@code DayKriging} holds the day's active stations, the kriged quantity and the k-d tree; {@link #chunk}
 * processes targets [c0, c1) as one random-stream unit and may be called concurrently for different chunks.
 */
public final class DayKriging {
    public enum Lapse { NONE, LOCAL }

    final KrigingSettings ks;
    final double[] sx, sy, sz; final float[] sdmat; final int S;
    final int[] act; final double[] y; final KdTree tree;
    final boolean zeroShortcut; final Lapse lapse; final double lapseExtrapM; final String dayTag; final double[] globalParams;

    /**
     * @param act active station indices (ascending); y kriged quantity at every station (index = station);
     * @param zeroShortcut precipitation without lapse (direct, climatology): all nKrig neighbours <= zeroTol -> 0
     * @param lapseExtrapM NaN = no clip of the trend elevation (direct mode), else clip (anomaly mode)
     * @param globalParams the day's global variogram (nugget, psill, range) or null
     */
    public DayKriging(KrigingSettings ks, double[] sx, double[] sy, double[] sz, float[] sdmat, int[] act, double[] y,
                      boolean zeroShortcut, Lapse lapse, double lapseExtrapM, String dayTag, double[] globalParams) {
        this.ks = ks; this.sx = sx; this.sy = sy; this.sz = sz; this.sdmat = sdmat; this.S = sx.length;
        this.act = act; this.y = y; this.zeroShortcut = zeroShortcut; this.lapse = lapse; this.lapseExtrapM = lapseExtrapM;
        this.dayTag = dayTag; this.globalParams = globalParams;
        double[] ax = new double[act.length], ay = new double[act.length];
        for (int i = 0; i < act.length; i++) { ax[i] = sx[act[i]]; ay[i] = sy[act[i]]; }
        tree = new KdTree(ax, ay);
    }

    /** Number of chunks for a target set. */
    public int nChunks(TargetSet t) { return (t.n + ks.chunk - 1) / ks.chunk; }

    /** Process chunk ci of the target set into res (disjoint ranges -> thread-safe). */
    public void chunk(TargetSet t, int ci, TargetResult res) {
        int c0 = ci * ks.chunk, c1 = Math.min(t.n, c0 + ks.chunk), B = c1 - c0;
        int A = act.length, kq = Math.min(ks.nVario + 1, A), nKeep = Math.min(ks.nVario, kq), K = ks.nKrig;
        KdTree.Query q = new KdTree.Query(Math.max(kq, Math.min(ks.lapseFallbackN + 1, A)));
        int[][] nb = new int[B][nKeep]; double[][] dkm = new double[B][nKeep]; double[][] Y = new double[B][nKeep];
        double[] trend = new double[B];
        boolean[] zero = new boolean[B];
        LocalLapse ll = new LocalLapse(Math.max(nKeep, ks.lapseFallbackN + 1));
        double[] E = new double[Math.max(nKeep, ks.lapseFallbackN + 1)], T = new double[E.length];
        for (int b = 0; b < B; b++) {
            int ti = c0 + b, self = t.self[ti];
            int got = tree.knn(t.x[ti], t.y[ti], kq, q);
            // stable move of the target's own station to the end, keep nKeep
            int w = 0;
            for (int j = 0; j < got && w < nKeep; j++) { int s = act[q.outIdx[j]]; if (s != self) { nb[b][w] = s; dkm[b][w] = q.outDist[j] / 1000.0; w++; } }
            for (int j = 0; j < got && w < nKeep; j++) { int s = act[q.outIdx[j]]; if (s == self) { nb[b][w] = s; dkm[b][w] = q.outDist[j] / 1000.0; w++; } }
            res.d1[ti] = dkm[b][0]; res.d7[ti] = dkm[b][K - 1]; res.d50[ti] = dkm[b][nKeep - 1];
            for (int j = 0; j < nKeep; j++) Y[b][j] = y[nb[b][j]];
            if (lapse == Lapse.LOCAL) {
                for (int j = 0; j < nKeep; j++) { E[j] = sz[nb[b][j]]; T[j] = Y[b][j]; }
                ll.fit(E, T, nKeep);
                double b0 = ll.b0, b1 = ll.b1, r2 = ll.r2, span = ll.span;
                boolean fb = !Variogram.isFinite(b1) || span < ks.lapseMinSpanM || ll.n < ks.lapseMinStations;
                if (fb) {
                    int k2 = Math.min(ks.lapseFallbackN + 1, A);
                    int g2 = tree.knn(t.x[ti], t.y[ti], k2, q);
                    for (int j = 0; j < g2; j++) { int s = act[q.outIdx[j]]; E[j] = sz[s]; T[j] = (s == self) ? Double.NaN : y[s]; }
                    ll.fit(E, T, g2);
                    b0 = ll.b0; b1 = ll.b1; r2 = ll.r2; span = ll.span;
                }
                b0 = nanToNum(b0); b1 = nanToNum(b1);
                double zmin = Double.NaN, zmax = Double.NaN;
                for (int j = 0; j < nKeep; j++) {
                    double e = sz[nb[b][j]];
                    Y[b][j] = Y[b][j] - (b0 + b1 * e);
                    if (Variogram.isFinite(Y[b][j])) { if (Double.isNaN(zmin) || e < zmin) zmin = e; if (Double.isNaN(zmax) || e > zmax) zmax = e; }
                }
                double ze = t.z[ti];
                if (!Double.isNaN(lapseExtrapM)) ze = NpMath.clip(ze, zmin - lapseExtrapM, zmax + lapseExtrapM);
                trend[b] = b0 + b1 * ze;
                res.b1[ti] = b1; res.r2[ti] = r2; res.span[ti] = span;
                res.flags[ti] |= (fb ? KrigingSettings.FLAG_LAPSE_FALLBACK : 0) | (b1 > 0 ? KrigingSettings.FLAG_POSITIVE_LAPSE : 0);
            }
            if (zeroShortcut) {
                boolean z = true;
                for (int j = 0; j < K; j++) if (!(Y[b][j] <= ks.zeroTol)) { z = false; break; }
                zero[b] = z;
            }
        }
        // variograms + PSO for the targets that need one (chunk order)
        int nq = 0; for (int b = 0; b < B; b++) if (!zero[b]) nq++;
        double[] params = new double[B * 3]; double[] sse = new double[B]; boolean[] bh = new boolean[B], vfb = new boolean[B];
        for (int b = 0; b < B; b++) { params[b * 3] = 0.0; params[b * 3 + 1] = 1.0; params[b * 3 + 2] = 100.0; sse[b] = Double.NaN; }
        if (nq > 0) {
            Variogram.Empirical[] ev = new Variogram.Empirical[nq]; int[] qi = new int[nq]; boolean[] valid = new boolean[nq];
            int k = 0;
            for (int b = 0; b < B; b++) if (!zero[b]) {
                Variogram.Empirical e = new Variogram.Empirical(ks.nLags);
                Variogram.empirical(Y[b], nb[b], nKeep, sdmat, S, ks.nLags, ks.minPairsPerLag, e);
                valid[k] = e.valid;
                if (!e.valid) Variogram.neutralise(e);
                ev[k] = e; qi[k] = b; k++;
            }
            double[] p_ = new double[nq * 3], s_ = new double[nq]; boolean[] b_ = new boolean[nq];
            SphericalPso.fit(ev, nq, newRng(dayTag, t.seedIndex[c0]), ks.pso, p_, s_, b_);
            for (int j = 0; j < nq; j++) {
                int b = qi[j];
                if (globalParams != null && !valid[j]) { System.arraycopy(globalParams, 0, p_, j * 3, 3); vfb[b] = true; }
                System.arraycopy(p_, j * 3, params, b * 3, 3); sse[b] = s_[j]; bh[b] = b_[j];
            }
        }
        OrdinaryKriging ok = new OrdinaryKriging(K);
        double[] nx = new double[K], ny = new double[K], dk = new double[K], yk = new double[K];
        for (int b = 0; b < B; b++) {
            int ti = c0 + b;
            for (int j = 0; j < K; j++) { nx[j] = sx[nb[b][j]]; ny[j] = sy[nb[b][j]]; dk[j] = dkm[b][j]; yk[j] = Y[b][j]; }
            ok.solve(nx, ny, dk, params[b * 3], params[b * 3 + 1], params[b * 3 + 2]);
            double est = ok.estimate(yk), var = ok.var; boolean good = ok.ok;
            if (zero[b]) { est = 0.0; var = 0.0; good = true; }
            ok.diagnostics(dk, ks.negWeightTol);
            res.est[ti] = (good ? est : Double.NaN) + trend[b];
            res.resid[ti] = good ? est : Double.NaN; res.trend[ti] = trend[b];
            res.kvar[ti] = good ? var : Double.NaN;
            if (!zero[b]) {
                res.nug[ti] = params[b * 3]; res.psill[ti] = params[b * 3 + 1]; res.range[ti] = params[b * 3 + 2];
                res.nNeg[ti] = ok.nNeg; res.neff[ti] = ok.neff; res.rpos12[ti] = ok.rpos12; res.rpos23[ti] = ok.rpos23;
                res.rpos13[ti] = ok.rpos13; res.rabs90[ti] = ok.rabs90; res.sumAbsW[ti] = ok.sumAbsW;
            }
            res.sse[ti] = sse[b];
            res.flags[ti] |= (zero[b] ? KrigingSettings.FLAG_ZERO : 0)
                    | ((ok.sumAbsW > ks.oscSumAbsW && !zero[b]) ? KrigingSettings.FLAG_OSCILLATORY : 0)
                    | (bh[b] ? KrigingSettings.FLAG_RANGE_BOUND : 0) | (vfb[b] ? KrigingSettings.FLAG_GLOBAL_VARIOGRAM : 0);
        }
    }

    /** Process all chunks of a target set sequentially. */
    public TargetResult run(TargetSet t) {
        TargetResult r = new TargetResult(t.n);
        for (int c = 0; c < nChunks(t); c++) chunk(t, c, r);
        return r;
    }

    private Rng newRng(String tag, Object s0) {
        return Rng.create(ks.rngMode, Blake2b.seedFor(tag, s0));
    }

    static double nanToNum(double v) {
        if (Double.isNaN(v)) return 0.0;
        if (v == Double.POSITIVE_INFINITY) return Double.MAX_VALUE;
        if (v == Double.NEGATIVE_INFINITY) return -Double.MAX_VALUE;
        return v;
    }

    /**
     * Global variogram on all active stations (Python global_variogram): params or null if invalid.
     * @param yv value at every station (index = station)
     */
    public static double[] globalVariogram(KrigingSettings ks, float[] sdmat, int S, int[] act, double[] yv, String tag) {
        Variogram.Empirical e = new Variogram.Empirical(ks.nLags);
        double[] v = new double[act.length];
        for (int i = 0; i < act.length; i++) v[i] = yv[act[i]];
        Variogram.empirical(v, act, act.length, sdmat, S, ks.nLags, ks.minPairsPerLag, e);
        if (!e.valid) return null;
        double[] p = new double[3], s = new double[1]; boolean[] b = new boolean[1];
        SphericalPso.fit(new Variogram.Empirical[]{e}, 1, Rng.create(ks.rngMode, Blake2b.seedFor(tag, "global")), ks.pso, p, s, b);
        return p;
    }

    /**
     * y - polyval(polyfit(z[act], y[act], 1), z) at every station (detrending before the global variogram in the
     * lapse modes). Centred OLS instead of numpy's SVD least squares: identical up to ~1e-15 relative (the global
     * variogram is only a fallback for invalid local variograms).
     */
    public static double[] detrend(double[] yv, double[] z, int[] act) {
        int n = act.length; double zm = 0, ym = 0;
        for (int i : act) { zm += z[i]; ym += yv[i]; }
        zm /= n; ym /= n;
        double szz = 0, szy = 0;
        for (int i : act) { szz += (z[i] - zm) * (z[i] - zm); szy += (z[i] - zm) * (yv[i] - ym); }
        double p0 = szy / szz, p1 = ym - p0 * zm;
        double[] out = new double[yv.length];
        for (int i = 0; i < yv.length; i++) out[i] = yv[i] - (p0 * z[i] + p1);
        return out;
    }
}
