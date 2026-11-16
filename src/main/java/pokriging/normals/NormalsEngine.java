package pokriging.normals;

import krigingsCore.local.KdTree;
import krigingsCore.local.NpMath;
import pokriging.io.Grid;
import pokriging.io.NcOut;
import pokriging.io.StationData;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Monthly normals of the Po engine (Python s1_station_normals, s2_normals_grid, s1w_wetfreq_normals): station monthly
 * values with gap filling, accepted station normals, PRISM-type gridded normals and station LOO normals, wet-day
 * frequency normals (precipitation). Output files and columns = the Python ones, so either engine can feed PoKriging.
 */
public final class NormalsEngine {
    final NormalsConfig c; StationData sd; Grid g; int S, NT; boolean precip;
    public final Map<String, Object> info = new LinkedHashMap<>();
    static final double[] DIM_REP = {31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    public NormalsEngine(NormalsConfig c) { this.c = c; }

    public void run() throws Exception {
        long t0 = System.nanoTime();
        c.resolveDefaults(); precip = c.variable.equals("precip");
        sd = StationData.read(Paths.get(c.stationsFile), Paths.get(c.dataFile), c.missingValue, c.variable, c.tempMin, c.tempMax);
        if (c.targetGrid != null) g = c.targetGrid;
        else {
            if (c.demFile == null || c.demFile.isEmpty()) throw new IllegalArgumentException("demFile (or inTargets) is required");
            g = Grid.readAscii(Paths.get(c.demFile));
            if (c.gridCellsFile != null && !c.gridCellsFile.isEmpty()) g.readCellAttributes(Paths.get(c.gridCellsFile));
        }
        S = sd.S;
        Path out = Paths.get(c.outputDir); Files.createDirectories(out);
        // ---------------- s1: station monthly values, gap filling, normals
        LocalDate a = LocalDate.parse(c.normalsStart), b = LocalDate.parse(c.normalsEnd);
        int y0 = a.getYear(), m0 = a.getMonthValue();
        NT = (b.getYear() - y0) * 12 + b.getMonthValue() - m0 + 1;
        double[][] M = new double[S][NT];
        int[][] nval = new int[S][NT], nwet = new int[S][NT]; int[] ndays = new int[NT];
        float[][] sum = new float[S][NT];
        int[] first = new int[NT]; Arrays.fill(first, -1);
        for (int d = 0; d < sd.T; d++) {
            LocalDate dt = sd.dates[d]; if (dt.isBefore(a) || dt.isAfter(b)) continue;
            int t = (dt.getYear() - y0) * 12 + dt.getMonthValue() - m0; ndays[t]++; if (first[t] < 0) first[t] = d;
            for (int s = 0; s < S; s++) { float v = sd.values[d][s]; if (!Float.isNaN(v)) { nval[s][t]++; if (v >= c.occThreshold) nwet[s][t]++; } }
        }
        // numpy add.reduceat of float32 (Python monthly_matrix): first element + pairwise sum of the remaining ones, in float32
        float[] buf = new float[32];
        for (int t = 0; t < NT; t++) for (int s = 0; s < S; s++) {
            int n = ndays[t];
            for (int k = 0; k < n; k++) { float v = sd.values[first[t] + k][s]; buf[k] = Float.isNaN(v) ? 0f : v; }
            sum[s][t] = n == 0 ? 0f : buf[0] + NpMath.pairwiseSumF(buf, 1, n - 1);
        }
        for (int s = 0; s < S; s++) for (int t = 0; t < NT; t++)
            M[s][t] = nval[s][t] >= c.monthCompleteness * ndays[t] ? ((double) sum[s][t] / nval[s][t]) * (precip ? ndays[t] : 1.0) : Double.NaN;
        int[] cal = new int[NT]; for (int t = 0; t < NT; t++) cal[t] = (m0 - 1 + t) % 12;
        double[][] F = fillGaps(M, cal);
        double[][] normal = new double[S][12], obsOnly = new double[S][12]; int[][] nobs = new int[S][12], ntot = new int[S][12];
        for (int s = 0; s < S; s++) for (int m = 0; m < 12; m++) {
            double so = 0, sf = 0; int no = 0, nf = 0;
            double[] ao = new double[NT / 12 + 1], af = new double[NT / 12 + 1]; int k = 0;
            for (int t = 0; t < NT; t++) if (cal[t] == m) { ao[k] = Double.isNaN(M[s][t]) ? 0 : M[s][t]; af[k] = Double.isNaN(F[s][t]) ? 0 : F[s][t];
                if (!Double.isNaN(M[s][t])) no++; if (!Double.isNaN(F[s][t])) nf++; k++; }
            so = NpMath.pairwiseSum(ao, 0, k); sf = NpMath.pairwiseSum(af, 0, k);
            nobs[s][m] = no; ntot[s][m] = nf; obsOnly[s][m] = so / no; double fill = sf / nf;
            normal[s][m] = (no >= c.minObsYears && nf >= c.minTotalYears) ? fill : Double.NaN;
        }
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(out.resolve("station_normals_" + c.variable + ".csv")))) {
            w.println("station,month,normal,normal_obs_only,n_obs_years,n_total_years,ok" + (precip ? ",normal_daily_mm" : ""));
            for (int m = 0; m < 12; m++) for (int s = 0; s < S; s++) {
                boolean ok = !Double.isNaN(normal[s][m]);
                w.println(sd.ids[s] + "," + (m + 1) + "," + f(normal[s][m]) + "," + f(obsOnly[s][m]) + "," + nobs[s][m] + "," + ntot[s][m] + "," + (ok ? "True" : "False")
                        + (precip ? "," + f(normal[s][m] / DIM_REP[m]) : ""));
            }
        }
        int acc = 0; for (int s = 0; s < S; s++) for (int m = 0; m < 12; m++) if (!Double.isNaN(normal[s][m])) acc++;
        info.put("station_months_accepted", acc); info.put("station_months", S * 12);
        // ---------------- s2: gridded normals + station LOO normals
        PrismRegression.Params p = new PrismRegression.Params();
        p.K = c.regK; p.dz0 = c.regDz0; p.facetW = c.facetWeight; p.tpi0 = c.tpi0; p.useTpi = c.useTpi; p.extrapDz = c.regExtrapDz;
        p.extrapLapse = c.regExtrapLapse; p.slopeLo = c.slopeMin; p.slopeHi = c.slopeMax;
        PrismRegression.Kind kind = precip ? PrismRegression.Kind.PRECIP : PrismRegression.Kind.TEMP;
        double[][] ys = new double[S][12];
        for (int s = 0; s < S; s++) for (int m = 0; m < 12; m++) ys[s][m] = precip ? StrictMath.log(NpMath.maximum(normal[s][m], c.precipFloorMm)) : normal[s][m];
        double[][][] r = regression(kind, p, ys, c.slopeMode.equals("annual"));
        double[][] G = r[0], SL = r[1], LOO = r[2];
        for (int m = 0; m < 12; m++) for (int k = 0; k < g.n; k++) { if (precip) G[m][k] = StrictMath.exp(G[m][k]); SL[m][k] *= 1000; }
        for (int s = 0; s < S; s++) for (int m = 0; m < 12; m++) if (precip) LOO[s][m] = StrictMath.exp(LOO[s][m]);
        Map<String, double[][]> fields = new LinkedHashMap<>(); fields.put("normal", G); fields.put("elevation_slope", SL);
        Map<String, String[]> ul = new LinkedHashMap<>();
        ul.put("normal", new String[]{precip ? "mm month-1" : "degC", c.normalsStart.substring(0, 4) + "-" + c.normalsEnd.substring(0, 4) + " monthly normal of " + c.variable});
        ul.put("elevation_slope", new String[]{precip ? "ln(mm) km-1" : "degC km-1", "local regression slope on elevation"});
        Map<String, String> at = new LinkedHashMap<>(); at.put("method", "PRISM-type local weighted regression (Java PoNormals)");
        at.put("slope_mode", c.slopeMode); at.put("period", c.normalsStart + ".." + c.normalsEnd);
        String ncp = (c.nccopy == null || c.nccopy.isEmpty()) ? null : c.nccopy;
        NcOut.writeMonthly(out.resolve("normals_" + c.variable + ".nc"), g, fields, ul, at, c.crs, ncp);
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(out.resolve("station_normals_loo_" + c.variable + ".csv")))) {
            w.println("station,month,normal_station,normal_grid_loo");
            for (int s = 0; s < S; s++) for (int m = 0; m < 12; m++) w.println(sd.ids[s] + "," + (m + 1) + "," + f(normal[s][m]) + "," + f(LOO[s][m]));
        }
        // ---------------- s1w: wet-day frequency normals
        if (precip && c.wetFrequency) {
            double[][] pst = new double[S][12]; int[][] nm = new int[S][12];
            for (int s = 0; s < S; s++) for (int m = 0; m < 12; m++) {
                long sw = 0, sv = 0; int n = 0;
                for (int t = 0; t < NT; t++) if (cal[t] == m && nval[s][t] >= c.monthCompleteness * ndays[t]) { sw += nwet[s][t]; sv += nval[s][t]; n++; }
                nm[s][m] = n; pst[s][m] = n < c.minObsYears ? Double.NaN : sw / (double) sv;
            }
            double[][] yl = new double[S][12];
            for (int s = 0; s < S; s++) for (int m = 0; m < 12; m++) { double q = NpMath.clip(pst[s][m], c.wfPmin, c.wfPmax); yl[s][m] = StrictMath.log(q / (1 - q)); }
            PrismRegression.Params pw = new PrismRegression.Params();
            pw.K = c.wfK; pw.dz0 = c.wfDz0; pw.facetW = c.facetWeight; pw.extrapDz = c.wfExtrapDz; pw.slopeLo = c.wfSlopeMin; pw.slopeHi = c.wfSlopeMax;
            double[][][] rw = regression(PrismRegression.Kind.WETFREQ, pw, yl, true);
            double[][] GW = rw[0], SW = rw[1], LW = rw[2];
            for (int m = 0; m < 12; m++) for (int k = 0; k < g.n; k++) { GW[m][k] = 1.0 / (1.0 + StrictMath.exp(-GW[m][k])); SW[m][k] *= 1000; }
            Map<String, double[][]> fw = new LinkedHashMap<>(); fw.put("wetfreq", GW); fw.put("elevation_slope", SW);
            Map<String, String[]> uw = new LinkedHashMap<>(); uw.put("wetfreq", new String[]{"1", "wet-day frequency (days >= " + c.occThreshold + " mm / days)"});
            uw.put("elevation_slope", new String[]{"logit km-1", "local slope of logit(frequency) on elevation"});
            NcOut.writeMonthly(out.resolve("wetfreq_precip.nc"), g, fw, uw, at, c.crs, ncp);
            try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(out.resolve("station_wetfreq_loo_precip.csv")))) {
                w.println("station,month,wetfreq_station,wetfreq_grid_loo,n_months");
                for (int s = 0; s < S; s++) for (int m = 0; m < 12; m++) w.println(sd.ids[s] + "," + (m + 1) + "," + f(pst[s][m]) + "," + f(1.0 / (1.0 + StrictMath.exp(-LW[s][m]))) + "," + nm[s][m]);
            }
        }
        info.put("seconds", (System.nanoTime() - t0) / 1e9);
    }

    /**
     * Grid estimates [12][ncell], grid slopes [12][ncell] and station LOO estimates [S][12] for y[S][12] (NaN = no normal).
     * annual: slope = mean of the 12 bounded monthly slopes (grid and, leave-one-out, stations), then imposed.
     */
    double[][][] regression(PrismRegression.Kind kind, PrismRegression.Params p, double[][] y, boolean annual) throws InterruptedException {
        double[] tz = g.z, ttpi = g.tpi; int[] tfac = g.facet; int[] gself = new int[g.n]; Arrays.fill(gself, -1);
        double[] gs = new double[g.n], ss = new double[S];
        PrismRegression.Stations[] sets = new PrismRegression.Stations[12]; int[][] pos = new int[12][];
        for (int m = 0; m < 12; m++) {
            List<Integer> idx = new ArrayList<>(); for (int s = 0; s < S; s++) if (!Double.isNaN(y[s][m])) idx.add(s);
            int n = idx.size(); double[] x = new double[n], yy = new double[n], z = new double[n], tp = new double[n], v = new double[n]; int[] fc = new int[n];
            pos[m] = new int[S]; Arrays.fill(pos[m], -1);
            for (int i = 0; i < n; i++) { int s = idx.get(i); x[i] = sd.x[s]; yy[i] = sd.y[s]; z[i] = sd.z[s]; tp[i] = sd.tpi[s]; fc[i] = sd.facet[s]; v[i] = y[s][m]; pos[m][s] = i; }
            sets[m] = new PrismRegression.Stations(x, yy, z, tp, fc, v);
            if (annual) {
                double[] gsl = PrismRegression.predict(kind, p, g.x, g.y, tz, ttpi, tfac, gself, sets[m], null, c.threads)[1];
                double[] ssl = PrismRegression.predict(kind, p, sd.x, sd.y, sd.z, sd.tpi, sd.facet, pos[m], sets[m], null, c.threads)[1];
                for (int k = 0; k < g.n; k++) gs[k] += gsl[k] / 12; for (int s = 0; s < S; s++) ss[s] += ssl[s] / 12;
            }
        }
        double[][] G = new double[12][], SL = new double[12][], L = new double[S][12];
        for (int m = 0; m < 12; m++) {
            double[] le = PrismRegression.predict(kind, p, sd.x, sd.y, sd.z, sd.tpi, sd.facet, pos[m], sets[m], annual ? ss : null, c.threads)[0];
            for (int s = 0; s < S; s++) L[s][m] = le[s];
            double[][] ge = PrismRegression.predict(kind, p, g.x, g.y, tz, ttpi, tfac, gself, sets[m], annual ? gs : null, c.threads);
            G[m] = ge[0]; SL[m] = ge[1];
        }
        return new double[][][]{G, SL, L};
    }

    /** s1 gap filling (Python s1_station_normals.fill_gaps). */
    double[][] fillGaps(double[][] M, int[] cal) {
        double[][] mu = new double[S][12], sdv = new double[S][12];
        for (int s = 0; s < S; s++) for (int m = 0; m < 12; m++) {
            double[] a = new double[NT / 12 + 1]; int k = 0, n = 0;
            for (int t = 0; t < NT; t++) if (cal[t] == m) { a[k++] = Double.isNaN(M[s][t]) ? 0 : M[s][t]; if (!Double.isNaN(M[s][t])) n++; }
            double avg = NpMath.pairwiseSum(a, 0, k) / n; mu[s][m] = avg;
            double[] q = new double[k]; int j = 0;
            for (int t = 0; t < NT; t++) if (cal[t] == m) { double v = Double.isNaN(M[s][t]) ? 0 : M[s][t] - avg; q[j] = v * v; j++; }
            sdv[s][m] = Math.sqrt(NpMath.pairwiseSum(q, 0, k) / n);
        }
        double[][] Z = new double[S][NT];
        for (int s = 0; s < S; s++) for (int t = 0; t < NT; t++) { double sdd = sdv[s][cal[t]] > 0 ? sdv[s][cal[t]] : Double.NaN; Z[s][t] = (M[s][t] - mu[s][cal[t]]) / sdd; }
        KdTree tree = new KdTree(sd.x, sd.y); int kq = Math.min(c.fillMaxCandidates + 1, S); KdTree.Query qq = new KdTree.Query(kq);
        double[][] F = new double[S][]; for (int s = 0; s < S; s++) F[s] = M[s].clone();
        double[] a = new double[NT], bb = new double[NT], cv = new double[NT], va = new double[NT], vb = new double[NT];
        for (int s = 0; s < S; s++) {
            boolean anyMiss = false; for (int t = 0; t < NT; t++) if (Double.isNaN(M[s][t])) { anyMiss = true; break; }
            if (!anyMiss) continue;
            int got = tree.knn(sd.x[s], sd.y[s], kq, qq);
            List<Integer> cand = new ArrayList<>(); List<Double> rr = new ArrayList<>();
            for (int j = 0; j < got; j++) {
                double d = qq.outDist[j]; if (!(d > 0 && d <= c.fillRadiusKm * 1000)) continue;
                int o = qq.outIdx[j]; int nov = 0;
                for (int t = 0; t < NT; t++) { boolean both = !Double.isNaN(Z[s][t]) && !Double.isNaN(Z[o][t]); a[t] = both ? Z[s][t] : 0; bb[t] = both ? Z[o][t] : 0; if (both) nov++; }
                double am = NpMath.pairwiseSum(a, 0, NT) / nov, bm = NpMath.pairwiseSum(bb, 0, NT) / nov;
                for (int t = 0; t < NT; t++) { boolean both = !Double.isNaN(Z[s][t]) && !Double.isNaN(Z[o][t]);
                    cv[t] = both ? (Z[s][t] - am) * (Z[o][t] - bm) : 0; va[t] = both ? (Z[s][t] - am) * (Z[s][t] - am) : 0; vb[t] = both ? (Z[o][t] - bm) * (Z[o][t] - bm) : 0; }
                double r = NpMath.pairwiseSum(cv, 0, NT) / Math.sqrt(NpMath.pairwiseSum(va, 0, NT) * NpMath.pairwiseSum(vb, 0, NT));
                if (nov >= c.fillMinOverlapMonths && r >= c.fillMinR) { cand.add(o); rr.add(r); }
            }
            if (cand.isEmpty()) continue;
            Integer[] ord = new Integer[cand.size()]; for (int i = 0; i < ord.length; i++) ord[i] = i;
            Arrays.sort(ord, (x, y) -> Double.compare(rr.get(y), rr.get(x)));
            int nc = ord.length; int[] cd = new int[nc]; double[] rv = new double[nc];
            for (int i = 0; i < nc; i++) { cd[i] = cand.get(ord[i]); rv[i] = rr.get(ord[i]); }
            double[][] rel = new double[nc][12];
            for (int i = 0; i < nc; i++) for (int m = 0; m < 12; m++) {
                double[] xs = new double[NT / 12 + 1], xn = new double[NT / 12 + 1]; int k = 0, nb = 0;
                for (int t = 0; t < NT; t++) if (cal[t] == m) {
                    boolean both = !Double.isNaN(M[s][t]) && !Double.isNaN(M[cd[i]][t]);
                    if (precip) { xs[k] = both ? M[s][t] : 0; xn[k] = both ? M[cd[i]][t] : 0; } else xs[k] = both ? M[s][t] - M[cd[i]][t] : 0;
                    if (both) nb++; k++;
                }
                double v = precip ? NpMath.pairwiseSum(xs, 0, k) / NpMath.pairwiseSum(xn, 0, k) : NpMath.pairwiseSum(xs, 0, k) / nb;
                rel[i][m] = nb < 5 ? Double.NaN : v;
            }
            for (int t = 0; t < NT; t++) {
                if (!Double.isNaN(M[s][t])) continue;
                double num = 0, den = 0; int used = 0;
                for (int i = 0; i < nc; i++) {
                    double e = precip ? M[cd[i]][t] * rel[i][cal[t]] : M[cd[i]][t] + rel[i][cal[t]];
                    boolean fin = !Double.isNaN(e) && !Double.isInfinite(e);
                    if (!fin) continue;
                    used++; if (used > c.fillTopK) break;
                    double w = rv[i] * rv[i]; num += e * w; den += w;
                }
                double val = num / den;
                if (!Double.isNaN(val) && !Double.isInfinite(val)) F[s][t] = precip ? Math.max(val, 0.0) : val;
            }
        }
        return F;
    }

    static String f(double v) { return Double.isNaN(v) ? "" : Double.toString(v); }
}
