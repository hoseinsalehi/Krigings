package pokriging.uncertainty;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;

import pokriging.io.Csv;

/**
 * Calibration of the kriging standard deviation with the leave-one-out errors (port of PoGrid_v3_code
 * s9_uncertainty_calibration.py:main; method: docs/UNCERTAINTY_CALIBRATION.md). OMS-independent.
 * <p>
 * Reads loo_&lt;var&gt;_&lt;Y&gt;.csv[.gz] (date,station,obs,pred,kriging_std,...) for the years of the period, and writes
 * uncertainty_&lt;var&gt;_&lt;start&gt;_&lt;end&gt;.csv (the s9 table = input of PoKriging uncertaintyTable) and
 * uncertainty_summary_&lt;var&gt;_&lt;start&gt;_&lt;end&gt;.csv (coverage of several nominal levels, station-half cross-validation).
 * Arithmetic follows pandas/numpy on float32 columns: e, |e|/sigma and 1.96 sigma in float32, quantiles (numpy linear) and
 * k in float64.
 */
public final class UncertaintyCalibration {

    /** Settings (= @In keys of pokriging.PoUncertainty). */
    public static final class Config {
        public String variable;                         // precip | temp
        public String looDir;                           // folder with loo_<var>_<Y>.csv[.gz]
        public String periodStart = "1991-01-01", periodEnd = "2020-12-31";
        public String stationsFile;                     // id,...,elevation (temperature classes; optional for precip)
        public String outputDir;                        // default looDir
        public String nominalLevels = "0.5,0.8,0.9,0.95";
        public int minCount = 50;                       // k needs more than minCount values (s9: 50)
    }

    public static final String[] SEASONS = {"DJF", "MAM", "JJA", "SON"};
    static final int[] SEASON_OF_MONTH = {0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 0};   // Jan..Dec -> DJF,MAM,JJA,SON
    public static final String[] PCLS = {"0", "0-1", "1-5", "5-20", "20-50", ">50", "nan"};
    public static final String[] TCLS = {"<500", "500-1000", "1000-1500", "1500-2000", ">=2000", "nan"};

    final Config c; final boolean precip; final String[] clsNames;
    // NaN-free LOO rows (dropna), in file order
    int n; byte[] season, cls, half; float[] e, sd;
    String first, last; int nFiles, nStations, nDropped, nNoElev;
    public Path outTable, outSummary;
    public final List<Map<String, Object>> table = new ArrayList<>();

    public UncertaintyCalibration(Config c) {
        this.c = c;
        if (!"precip".equals(c.variable) && !"temp".equals(c.variable)) throw new IllegalArgumentException("variable must be precip or temp");
        precip = c.variable.equals("precip"); clsNames = precip ? PCLS : TCLS;
    }

    /** class index of a LOO row (pd.cut, right-closed bins; outside -> "nan"). */
    static int precipClass(double p) {
        if (!(p > -1) || !(p <= 1e9)) return 6;
        return p <= 0 ? 0 : p <= 1 ? 1 : p <= 5 ? 2 : p <= 20 ? 3 : p <= 50 ? 4 : 5;
    }

    static int tempClass(double z) {
        if (!(z > -1e3) || !(z <= 9e3)) return 5;
        return z <= 500 ? 0 : z <= 1000 ? 1 : z <= 1500 ? 2 : z <= 2000 ? 3 : 4;
    }

    // ------------------------------------------------------------------ input
    void read() throws IOException {
        Path dir = Paths.get(c.looDir);
        Map<String, Double> elev = new HashMap<>();
        if (c.stationsFile != null && !c.stationsFile.isEmpty()) {
            Csv s = Csv.read(Paths.get(c.stationsFile));
            for (int i = 0; i < s.rows.size(); i++) elev.put(s.get(i, "id"), s.num(i, "elevation"));
        } else if (!precip) throw new IllegalArgumentException("temp: stationsFile (elevation classes) is required");
        int y0 = Integer.parseInt(c.periodStart.substring(0, 4)), y1 = Integer.parseInt(c.periodEnd.substring(0, 4));
        int cap = 1 << 20; season = new byte[cap]; cls = new byte[cap]; half = new byte[cap]; e = new float[cap]; sd = new float[cap];
        Map<String, Integer> stIndex = new HashMap<>(); Map<String, Integer> stClass = new HashMap<>();
        for (int y = y0; y <= y1; y++) {
            Path f = dir.resolve("loo_" + c.variable + "_" + y + ".csv");
            if (!Files.exists(f)) f = dir.resolve("loo_" + c.variable + "_" + y + ".csv.gz");
            if (!Files.exists(f)) continue;
            nFiles++;
            try (BufferedReader r = Csv.open(f)) {
                String[] h = Csv.split(r.readLine().replace("\uFEFF", ""));
                int iD = -1, iS = -1, iO = -1, iP = -1, iK = -1;
                for (int i = 0; i < h.length; i++) switch (h[i].trim()) {
                    case "date": iD = i; break; case "station": iS = i; break; case "obs": iO = i; break;
                    case "pred": iP = i; break; case "kriging_std": iK = i; break; default:
                }
                if (iD < 0 || iS < 0 || iO < 0 || iP < 0 || iK < 0) throw new IOException(f + ": needs columns date,station,obs,pred,kriging_std");
                int nc = Math.max(Math.max(Math.max(iD, iS), Math.max(iO, iP)), iK) + 1;
                String[] fld = new String[nc]; String l;
                while ((l = r.readLine()) != null) {
                    if (l.isEmpty()) continue;
                    int k = 0, s = 0;
                    for (int i = 0; i < l.length() && k < nc; i++) if (l.charAt(i) == ',') { fld[k++] = l.substring(s, i); s = i + 1; }
                    if (k < nc) fld[k++] = l.substring(s);
                    String d = fld[iD];
                    if (d.compareTo(c.periodStart) < 0 || d.compareTo(c.periodEnd) > 0) continue;
                    float o = num(fld[iO]), p = num(fld[iP]), q = num(fld[iK]);
                    if (Float.isNaN(o) || Float.isNaN(p) || Float.isNaN(q) || d.isEmpty() || fld[iS].isEmpty()) { nDropped++; continue; }
                    String st = fld[iS].trim();
                    Integer si = stIndex.get(st);
                    if (si == null) {                      // pd.factorize: index of first appearance
                        si = stIndex.size(); stIndex.put(st, si);
                        Double z = elev.get(st); if (!precip && z == null) nNoElev++;
                        stClass.put(st, precip ? -1 : tempClass(z == null ? Double.NaN : z));
                    }
                    if (n == e.length) grow();
                    int m = Integer.parseInt(d.substring(5, 7));
                    season[n] = (byte) SEASON_OF_MONTH[m - 1];
                    cls[n] = (byte) (precip ? precipClass(p) : stClass.get(st));
                    half[n] = (byte) (si % 2);
                    e[n] = o - p; sd[n] = q;               // float32 as pandas
                    if (first == null || d.compareTo(first) < 0) first = d;
                    if (last == null || d.compareTo(last) > 0) last = d;
                    n++;
                }
            }
        }
        nStations = stIndex.size();
        if (n == 0) throw new IOException("no LOO rows of " + c.variable + " in " + dir + " for " + c.periodStart + " .. " + c.periodEnd);
    }

    static float num(String s) {
        s = s.trim();
        if (s.isEmpty() || s.equalsIgnoreCase("nan")) return Float.NaN;
        return Float.parseFloat(s);
    }

    void grow() {
        int m = e.length * 2;
        season = Arrays.copyOf(season, m); cls = Arrays.copyOf(cls, m); half = Arrays.copyOf(half, m); e = Arrays.copyOf(e, m); sd = Arrays.copyOf(sd, m);
    }

    // ------------------------------------------------------------------ statistics
    /** numpy.quantile(method='linear') of a sorted float32 sample (numpy _lerp evaluated in float64). */
    public static double quantile(float[] x, int len, double p) {
        if (len == 0) return Double.NaN;
        double v = (len - 1) * p; double fl = Math.floor(v);
        int j = (int) fl; if (j < 0) j = 0; if (j > len - 1) j = len - 1;
        int jn = Math.min(j + 1, len - 1);
        double g = v - j; float a = x[j], b = x[jn]; float d = b - a;
        return g >= 0.5 ? b - (double) d * (1 - g) : a + (double) d * g;
    }

    static float[] sorted(float[] v, int len) { float[] s = Arrays.copyOf(v, len); Arrays.sort(s); return s; }

    static double zOf(double p) {
        if (Math.abs(p - 0.95) < 1e-12) return 1.96;          // s9 convention
        return new org.apache.commons.math3.distribution.NormalDistribution(0, 1).inverseCumulativeProbability((1 + p) / 2);
    }

    /** One group: rows idx[0..m). Returns the s9 row; adds the summary rows to sum. */
    Map<String, Object> group(String sName, String cName, int[] idx, int m, double[] levels, List<Map<String, Object>> sum) {
        float[] zAll = new float[m], z0 = new float[m], z1 = new float[m], eAll = new float[m], e0 = new float[m], e1 = new float[m];
        int np = 0, n0 = 0, n1 = 0, ne0 = 0, ne1 = 0, nzero = 0;
        for (int t = 0; t < m; t++) {
            int i = idx[t]; float ae = Math.abs(e[i]); eAll[t] = e[i];
            if (half[i] == 0) e0[ne0++] = e[i]; else e1[ne1++] = e[i];
            if (sd[i] > 0) { float z = ae / sd[i]; zAll[np++] = z; if (half[i] == 0) z0[n0++] = z; else z1[n1++] = z; } else nzero++;
        }
        float[] sz = sorted(zAll, np), sz0 = sorted(z0, n0), sz1 = sorted(z1, n1), se = sorted(eAll, m), se0 = sorted(e0, ne0), se1 = sorted(e1, ne1);
        Map<String, Object> row = new java.util.LinkedHashMap<>();
        row.put("variable", c.variable); row.put("season", sName); row.put("cls", cName); row.put("n", m);
        row.put("share_zero_std", nzero / (double) m);
        for (double lev : levels) {
            double zp = zOf(lev); float zpf = (float) zp;
            long cr = 0; for (int t = 0; t < m; t++) { int i = idx[t]; if (sd[i] > 0 && Math.abs(e[i]) <= zpf * sd[i]) cr++; }
            double covRaw = np > 0 ? cr / (double) np : Double.NaN;
            double k = np > c.minCount ? quantile(sz, np, lev) / zp : Double.NaN;
            double[] kh = new double[2]; Arrays.fill(kh, Double.NaN);
            List<Double> cov = new ArrayList<>(), covE = new ArrayList<>();
            for (int h = 0; h < 2; h++) {
                int ntr = h == 0 ? n0 : n1, nte = h == 0 ? n1 : n0;
                if (ntr > c.minCount && nte > c.minCount) {
                    kh[h] = quantile(h == 0 ? sz0 : sz1, ntr, lev) / zp;
                    float bnd = (float) (zp * kh[h]); long cc = 0;
                    for (int t = 0; t < m; t++) { int i = idx[t]; if (sd[i] > 0 && half[i] != h && Math.abs(e[i]) <= bnd * sd[i]) cc++; }
                    cov.add(cc / (double) nte);
                }
                int etr = h == 0 ? ne0 : ne1, ete = h == 0 ? ne1 : ne0;
                if (etr > c.minCount && ete > c.minCount) {
                    float[] st = h == 0 ? se0 : se1;
                    double lo = quantile(st, etr, (1 - lev) / 2), hi = quantile(st, etr, (1 + lev) / 2); long cc = 0;
                    for (int t = 0; t < m; t++) { int i = idx[t]; if (half[i] != h && e[i] >= lo && e[i] <= hi) cc++; }
                    covE.add(cc / (double) ete);
                }
            }
            double covCal = cov.isEmpty() ? Double.NaN : cov.size() == 1 ? cov.get(0) : (cov.get(0) + cov.get(1)) / 2;
            double covEmp = covE.isEmpty() ? Double.NaN : covE.size() == 1 ? covE.get(0) : (covE.get(0) + covE.get(1)) / 2;
            Map<String, Object> s = new java.util.LinkedHashMap<>();
            s.put("variable", c.variable); s.put("season", sName); s.put("cls", cName); s.put("nominal", lev); s.put("z", zp);
            s.put("n", m); s.put("n_pos", np); s.put("coverage_raw", covRaw); s.put("scale_k", k); s.put("k_half0", kh[0]); s.put("k_half1", kh[1]);
            s.put("coverage_cal_crossval", covCal); s.put("coverage_empirical_crossval", covEmp);
            sum.add(s);
            if (Math.abs(lev - 0.95) < 1e-12) { row.put("coverage95_raw", covRaw); row.put("scale_k", k); row.put("coverage95_cal_crossval", covCal); }
        }
        if (!row.containsKey("scale_k")) throw new IllegalStateException("nominalLevels must contain 0.95");
        double[] qs = {0.05, 0.25, 0.5, 0.75, 0.95}; String[] qn = {"e_q05", "e_q25", "e_q50", "e_q75", "e_q95"};
        for (int j = 0; j < qs.length; j++) row.put(qn[j], quantile(se, m, (qs[j] * 100.0) / 100.0));   // pandas: percentile(100 p)
        // s9 column order
        Map<String, Object> o = new java.util.LinkedHashMap<>();
        for (String key : new String[]{"variable", "season", "cls", "n", "share_zero_std", "coverage95_raw", "scale_k", "coverage95_cal_crossval",
                "e_q05", "e_q25", "e_q50", "e_q75", "e_q95"}) o.put(key, row.get(key));
        return o;
    }

    // ------------------------------------------------------------------ run
    public void run() throws IOException {
        long t0 = System.currentTimeMillis();
        read();
        double[] levels = Arrays.stream(c.nominalLevels.split(",")).map(String::trim).filter(s -> !s.isEmpty()).mapToDouble(Double::parseDouble).toArray();
        // groups (season, class) that occur, sorted by the strings as pandas groupby, then (all, all)
        int nc = clsNames.length; int[] cnt = new int[4 * nc];
        for (int i = 0; i < n; i++) cnt[season[i] * nc + cls[i]]++;
        TreeMap<String, Integer> order = new TreeMap<>((a, b) -> {
            String[] x = a.split("\u0000"), y = b.split("\u0000"); int r = x[0].compareTo(y[0]); return r != 0 ? r : x[1].compareTo(y[1]); });
        for (int g = 0; g < cnt.length; g++) if (cnt[g] > 0) order.put(SEASONS[g / nc] + "\u0000" + clsNames[g % nc], g);
        int[][] members = new int[cnt.length][]; int[] fill = new int[cnt.length];
        for (int g = 0; g < cnt.length; g++) members[g] = new int[cnt[g]];
        for (int i = 0; i < n; i++) { int g = season[i] * nc + cls[i]; members[g][fill[g]++] = i; }
        List<Map<String, Object>> sum = new ArrayList<>();
        for (Map.Entry<String, Integer> en : order.entrySet()) {
            int g = en.getValue();
            table.add(group(SEASONS[g / nc], clsNames[g % nc], members[g], cnt[g], levels, sum));
        }
        int[] all = new int[n]; for (int i = 0; i < n; i++) all[i] = i;
        table.add(group("all", "all", all, n, levels, sum));
        Path od = Paths.get(c.outputDir == null || c.outputDir.isEmpty() ? c.looDir : c.outputDir); Files.createDirectories(od);
        outTable = od.resolve("uncertainty_" + c.variable + "_" + first + "_" + last + ".csv");
        outSummary = od.resolve("uncertainty_summary_" + c.variable + "_" + first + "_" + last + ".csv");
        write(outTable, table); write(outSummary, sum);
        System.out.printf(Locale.ROOT, "UncertaintyCalibration %s: %d LOO files, %d rows (%d with missing values dropped), %d stations%s, %s .. %s, %.1f s%n",
                c.variable, nFiles, n, nDropped, nStations, nNoElev > 0 ? " (" + nNoElev + " without elevation -> class nan)" : "", first, last,
                (System.currentTimeMillis() - t0) / 1000.0);
        System.out.println("  season cls        n          cov95_raw  k        cov95_cal_cv");
        for (Map<String, Object> r : table)
            System.out.printf(Locale.ROOT, "  %-6s %-10s %-10d %-10.4f %-8.4f %.4f%n", r.get("season"), r.get("cls"), r.get("n"),
                    (Double) r.get("coverage95_raw"), (Double) r.get("scale_k"), (Double) r.get("coverage95_cal_crossval"));
        System.out.println("  -> " + outTable + "\n  -> " + outSummary);
    }

    static void write(Path p, List<Map<String, Object>> rows) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(p, StandardCharsets.UTF_8))) {
            w.println(String.join(",", rows.get(0).keySet()));
            for (Map<String, Object> r : rows) {
                StringBuilder sb = new StringBuilder();
                for (Object v : r.values()) {
                    if (sb.length() > 0) sb.append(',');
                    if (v instanceof Double && ((Double) v).isNaN()) continue;   // pandas writes NaN as an empty field
                    sb.append(v);
                }
                w.println(sb);
            }
        }
    }

    /** java pokriging.uncertainty.UncertaintyCalibration key=value ... (keys of Config). */
    public static void main(String[] args) throws Exception {
        Config c = new Config();
        for (String a : args) {
            int k = a.indexOf('='); Field f = Config.class.getField(a.substring(0, k)); String v = a.substring(k + 1);
            if (f.getType() == int.class) f.setInt(c, Integer.parseInt(v)); else f.set(c, v);
        }
        new UncertaintyCalibration(c).run();
    }
}
