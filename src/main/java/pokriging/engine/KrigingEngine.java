package pokriging.engine;

import krigingsCore.local.DayKriging;
import krigingsCore.local.KdTree;
import krigingsCore.local.KrigingSettings;
import krigingsCore.local.NpMath;
import krigingsCore.local.TargetResult;
import krigingsCore.local.TargetSet;
import pokriging.io.Csv;
import pokriging.io.Grid;
import pokriging.io.NcOut;
import pokriging.io.NormalsIO;
import pokriging.io.StationData;

import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Daily interpolation and leave-one-out of the Po engine (Python s3_run.process_day / twopart_day /
 * finalize_twopart), OMS-independent. Multithreaded over (day, target chunk) units; every chunk has its own random
 * stream, so results do not depend on the number of threads.
 */
public final class KrigingEngine {
    private static final Logger LOG = Logger.getLogger("pokriging");
    private static final LocalDate EPOCH = LocalDate.of(1991, 1, 1);

    final KrigingConfig c; final KrigingSettings ks = new KrigingSettings();
    StationData sd; float[] sdmat; boolean[] usable; int S;
    boolean precip, twopart, clim, lapseLocal, lapseGlobal, doGrid, doPoints, doLoo;
    String modeTag;
    Grid grid; int[] cells; TargetSet tsGrid;
    String[] pid; TargetSet tsPts; int[] pCell;
    double[][] NS, NL, PSF, PLF; boolean[][] nsFallback; double[][] NG, PGF;
    float[][] thrTable;                       // [cell][month] for calibration = table
    double[][] ukTable;                       // uncertainty k [season][class]
    List<Integer> days = new ArrayList<>();
    Path out;
    public final Map<String, Object> info = new LinkedHashMap<>();

    public KrigingEngine(KrigingConfig c) { this.c = c; }

    // ======================================================================================== set-up
    void init() throws IOException {
        long t0 = System.currentTimeMillis();
        if (c.variable == null || !(c.variable.equals("precip") || c.variable.equals("temp"))) throw new IllegalArgumentException("variable must be precip|temp");
        precip = c.variable.equals("precip");
        if (!c.variogramModel.equals("spherical")) throw new IllegalArgumentException("only variogramModel = spherical");
        if (precip && !Arrays.asList("direct", "climatology", "twopart").contains(c.precipMethod)) throw new IllegalArgumentException("precipMethod: direct|climatology|twopart");
        if (!precip && !Arrays.asList("direct", "anomaly").contains(c.tempMethod)) throw new IllegalArgumentException("tempMethod: direct|anomaly");
        if (!Arrays.asList("local", "none", "global").contains(c.temperatureTrend)) throw new IllegalArgumentException("temperatureTrend: local|none|global");
        twopart = precip && c.precipMethod.equals("twopart");
        clim = precip ? !c.precipMethod.equals("direct") : c.tempMethod.equals("anomaly");
        modeTag = clim ? "anomaly" : "direct";
        // Python: temperature direct always uses the local lapse; anomaly uses temperatureTrend
        lapseLocal = !precip && (c.tempMethod.equals("direct") ? !c.temperatureTrend.equals("global") && !c.temperatureTrend.equals("none") || c.temperatureTrend.equals("local") : c.temperatureTrend.equals("local"));
        if (!precip && c.tempMethod.equals("direct") && c.temperatureTrend.equals("none")) lapseLocal = false;
        lapseGlobal = !precip && c.temperatureTrend.equals("global");
        doLoo = !c.mode.equals("interpolation"); doGrid = !c.mode.equals("loo") && c.target.equals("grid"); doPoints = !c.mode.equals("loo") && c.target.equals("points");
        ks.nVario = c.nVario; ks.nKrig = c.nKrig; ks.nLags = c.nLags; ks.minPairsPerLag = c.minPairsPerLag; ks.chunk = c.chunk;
        ks.lapseMinSpanM = c.lapseMinSpanM; ks.lapseMinStations = c.lapseMinStations; ks.lapseFallbackN = c.lapseFallbackN;
        ks.anomLapseMaxExtrapM = c.anomLapseMaxExtrapM; ks.zeroTol = c.zeroTol; ks.globalVariogramFallback = c.globalVariogramFallback;
        ks.rngMode = c.rngMode; ks.pso.particles = c.psoParticles; ks.pso.iterations = c.psoIterations; ks.pso.inertia = c.psoInertia;
        ks.pso.c1 = c.psoC1; ks.pso.c2 = c.psoC2; ks.pso.pairCountWeights = c.fitWeighting.equals("pair-count");
        if (c.minActiveStations < 0) c.minActiveStations = c.nVario + 1;
        out = Paths.get(c.outputDir); Files.createDirectories(out);

        sd = StationData.read(Paths.get(c.stationsFile), Paths.get(c.dataFile), c.missingValue, c.variable, c.tempMin, c.tempMax);
        S = sd.S; sdmat = sd.distanceMatrixKm(); usable = new boolean[S]; Arrays.fill(usable, true);
        if (c.excludeStationsFile != null && !c.excludeStationsFile.isEmpty()) {
            Set<String> ex = new HashSet<>(); for (String l : Files.readAllLines(Paths.get(c.excludeStationsFile))) if (!l.trim().isEmpty()) ex.add(l.trim());
            for (int s = 0; s < S; s++) if (ex.contains(sd.ids[s])) usable[s] = false;
        }
        boolean needGrid = doGrid || doPoints || clim;
        if (needGrid) {
            if (c.targetGrid != null) grid = c.targetGrid;
            else {
                if (c.demFile == null || c.demFile.isEmpty()) throw new IllegalArgumentException("demFile (or inTargets) is required");
                grid = Grid.readAscii(Paths.get(c.demFile));
                if (c.gridCellsFile != null && !c.gridCellsFile.isEmpty()) grid.readCellAttributes(Paths.get(c.gridCellsFile));
            }
        }
        nsFallback = new boolean[S][12];
        if (clim) {
            Path nd = Paths.get(c.normalsDir);
            double[][][] p = NormalsIO.stationPair(nd.resolve("station_normals_loo_" + c.variable + ".csv"), sd.ids, "normal_station", "normal_grid_loo");
            NS = p[0]; NL = p[1];
            for (int s = 0; s < S; s++) for (int m = 0; m < 12; m++) if (!Double.isFinite(NS[s][m])) { NS[s][m] = NL[s][m]; nsFallback[s][m] = true; }
            NG = NormalsIO.gridField(nd.resolve("normals_" + c.variable + ".nc"), "normal", grid);
            if (twopart) {
                double[][][] q = NormalsIO.stationPair(nd.resolve("station_wetfreq_loo_precip.csv"), sd.ids, "wetfreq_station", "wetfreq_grid_loo");
                PSF = q[0]; PLF = q[1];
                for (int s = 0; s < S; s++) for (int m = 0; m < 12; m++) if (!Double.isFinite(PSF[s][m])) { PSF[s][m] = PLF[s][m]; nsFallback[s][m] = true; }
                PGF = NormalsIO.gridField(nd.resolve("wetfreq_precip.nc"), "wetfreq", grid);
            }
        }
        if (doGrid) {
            if (c.targetCellsFile != null && !c.targetCellsFile.isEmpty()) {
                Csv tc = Csv.read(Paths.get(c.targetCellsFile)); cells = new int[tc.rows.size()];
                for (int i = 0; i < cells.length; i++) cells[i] = (int) tc.num(i, "cell");
            } else { cells = new int[grid.n]; for (int i = 0; i < grid.n; i++) cells[i] = i; }
            int n = cells.length; double[] x = new double[n], y = new double[n], z = new double[n]; int[] self = new int[n]; long[] seed = new long[n];
            for (int i = 0; i < n; i++) { int k = cells[i]; x[i] = grid.x[k]; y[i] = grid.y[k]; z[i] = grid.z[k]; self[i] = -1; seed[i] = k; }
            tsGrid = new TargetSet(x, y, z, self, seed);
        }
        if (doPoints) {
            Csv pc = Csv.read(Paths.get(c.pointsFile)); int n = pc.rows.size();
            pid = new String[n]; pCell = new int[n];
            double[] x = new double[n], y = new double[n], z = new double[n]; int[] self = new int[n]; long[] seed = new long[n];
            KdTree gt = new KdTree(grid.x, grid.y); KdTree.Query q = new KdTree.Query(1);
            for (int i = 0; i < n; i++) {
                pid[i] = pc.has("id") ? pc.get(i, "id") : pc.get(i, "name");
                x[i] = pc.num(i, "x"); y[i] = pc.num(i, "y");
                gt.knn(x[i], y[i], 1, q); pCell[i] = q.outIdx[0];
                double e = pc.has("elevation") ? pc.num(i, "elevation") : Double.NaN;
                if (!Double.isFinite(e)) { int k = grid.nearestCellRounded(x[i], y[i]); e = grid.z[k >= 0 ? k : pCell[i]]; }
                z[i] = e; self[i] = -1; seed[i] = i;
            }
            tsPts = new TargetSet(x, y, z, self, seed);
        }
        if (twopart && c.calibration.equals("table")) {
            Csv t = Csv.read(Paths.get(c.occurrenceThresholds)); thrTable = new float[grid.n][12];
            for (float[] r : thrTable) Arrays.fill(r, Float.NaN);
            int cc = t.col("cell"), cm = t.col("month"), ct = t.col("threshold");
            for (String[] r : t.rows) thrTable[Integer.parseInt(r[cc].trim())][Integer.parseInt(r[cm].trim()) - 1] = (float) Csv.num(r[ct]);
        }
        if (c.uncertaintyTable != null && !c.uncertaintyTable.isEmpty()) readUncertainty(Paths.get(c.uncertaintyTable));
        setupOutputs();
        // days
        LocalDate a = LocalDate.parse(c.periodStart), b = LocalDate.parse(c.periodEnd);
        Set<LocalDate> want = null;
        if (c.datesFile != null && !c.datesFile.isEmpty()) {
            Csv dc = Csv.read(Paths.get(c.datesFile)); want = new HashSet<>();
            String col = dc.has("date") ? "date" : "Date";
            for (int i = 0; i < dc.rows.size(); i++) want.add(LocalDate.parse(dc.get(i, col).substring(0, 10)));
        }
        for (int d = 0; d < sd.T; d++) {
            LocalDate dt = sd.dates[d];
            if (dt.isBefore(a) || dt.isAfter(b)) continue;
            if (want != null && !want.contains(dt)) continue;
            days.add(d);
        }
        info.put("setup_s", (System.currentTimeMillis() - t0) / 1000.0);
        LOG.info(String.format("PoKriging: variable=%s method=%s lapse=%s days=%d stations=%d grid targets=%d points=%d loo=%b threads=%d",
                c.variable, precip ? c.precipMethod : c.tempMethod, lapseLocal ? "local" : lapseGlobal ? "global" : "none", days.size(), S,
                tsGrid == null ? 0 : tsGrid.n, tsPts == null ? 0 : tsPts.n, doLoo, c.threads));
    }

    // ======================================================================================== one day
    final class Day {
        int di, m; LocalDate date; double nd; String tag; double[] vals; int[] act, actW; boolean skipped; double r0; long t0; double seconds;
        DayKriging k1, k2; double gb0, gb1, gr2 = Double.NaN;                         // k1 occurrence / main, k2 amount
        TargetSet tsLoo; TargetResult rG, rP, rL, rG2, rP2, rL2; List<Callable<Void>> tasks = new ArrayList<>();
    }

    private Day prepare(int di) {
        Day d = new Day(); d.t0 = System.nanoTime(); d.di = di; d.date = sd.dates[di]; d.m = d.date.getMonthValue() - 1; d.nd = d.date.lengthOfMonth();
        double[] vals = new double[S]; for (int s = 0; s < S; s++) vals[s] = sd.values[di][s]; d.vals = vals;
        if (twopart) prepareTwoPart(d); else prepareSingle(d);
        if (d.skipped) return d;
        if (doLoo) {
            int n = d.act.length; double[] x = new double[n], y = new double[n], z = new double[n]; long[] seed = new long[n];
            for (int i = 0; i < n; i++) { int s = d.act[i]; x[i] = sd.x[s]; y[i] = sd.y[s]; z[i] = sd.z[s]; seed[i] = i; }
            d.tsLoo = new TargetSet(x, y, z, d.act.clone(), seed);
        }
        if (tsGrid != null) { d.rG = new TargetResult(tsGrid.n); addTasks(d, d.k1, tsGrid, d.rG); if (d.k2 != null) { d.rG2 = new TargetResult(tsGrid.n); addTasks(d, d.k2, tsGrid, d.rG2); } }
        if (tsPts != null) { d.rP = new TargetResult(tsPts.n); addTasks(d, d.k1, tsPts, d.rP); if (d.k2 != null) { d.rP2 = new TargetResult(tsPts.n); addTasks(d, d.k2, tsPts, d.rP2); } }
        if (d.tsLoo != null) { d.rL = new TargetResult(d.tsLoo.n); addTasks(d, d.k1, d.tsLoo, d.rL); if (d.k2 != null) { d.rL2 = new TargetResult(d.tsLoo.n); addTasks(d, d.k2, d.tsLoo, d.rL2); } }
        return d;
    }

    private void addTasks(Day d, DayKriging k, TargetSet ts, TargetResult r) {
        for (int ci = 0; ci < k.nChunks(ts); ci++) { final int cc = ci; d.tasks.add(() -> { k.chunk(ts, cc, r); return null; }); }
    }

    private int[] active(double[] y, boolean requireFinite) {
        int n = 0; int[] a = new int[S];
        for (int s = 0; s < S; s++) if (usable[s] && Double.isFinite(y[s])) a[n++] = s;
        return Arrays.copyOf(a, n);
    }

    private void prepareSingle(Day d) {
        double[] y = new double[S];
        if (!clim) System.arraycopy(d.vals, 0, y, 0, S);
        else if (precip) for (int s = 0; s < S; s++) { double cd = NpMath.maximum(NS[s][d.m] / d.nd, c.minDailyNormalMm); y[s] = NpMath.minimum(d.vals[s] / cd, c.maxPrecipRatio); }
        else for (int s = 0; s < S; s++) y[s] = d.vals[s] - NS[s][d.m];
        d.tag = c.variable + "|" + modeTag + "|" + d.date;
        d.act = active(y, true);
        if (d.act.length < c.minActiveStations) { d.skipped = true; return; }
        double[] yk = y;
        if (lapseGlobal) {                                   // one daily OLS on all active stations (not in Python, §A-16)
            double zm = 0, ym = 0; for (int s : d.act) { zm += sd.z[s]; ym += y[s]; } zm /= d.act.length; ym /= d.act.length;
            double szz = 0, szy = 0; for (int s : d.act) { szz += (sd.z[s] - zm) * (sd.z[s] - zm); szy += (sd.z[s] - zm) * (y[s] - ym); }
            d.gb1 = szz > 1e-8 ? szy / szz : 0.0; d.gb0 = ym - d.gb1 * zm;
            double syy = 0; for (int s : d.act) syy += (y[s] - ym) * (y[s] - ym);
            d.gr2 = szz > 1e-8 && syy > 0 ? szy * szy / (szz * syy) : Double.NaN;   // diagnostics only
            yk = new double[S]; for (int s = 0; s < S; s++) yk[s] = y[s] - (d.gb0 + d.gb1 * sd.z[s]);
        }
        double[] gpar = null;
        if (c.globalVariogramFallback) gpar = DayKriging.globalVariogram(ks, sdmat, S, d.act, lapseLocal ? DayKriging.detrend(yk, sd.z, d.act) : yk, d.tag);
        d.k1 = new DayKriging(ks, sd.x, sd.y, sd.z, sdmat, d.act, yk, precip, lapseLocal ? DayKriging.Lapse.LOCAL : DayKriging.Lapse.NONE,
                (clim && lapseLocal) ? c.anomLapseMaxExtrapM : Double.NaN, d.tag, gpar);
    }

    private void prepareTwoPart(Day d) {
        d.tag = c.variable + "|twopart|" + d.date;
        double[] yOcc = new double[S], yAmt = new double[S];
        boolean[] fin = new boolean[S], wet = new boolean[S];
        for (int s = 0; s < S; s++) {
            double v = d.vals[s]; fin[s] = Double.isFinite(v); wet[s] = fin[s] && v >= c.occThreshold;
            double ps = PSF[s][d.m], ist = NS[s][d.m] / (NpMath.maximum(ps, c.pFloor) * d.nd);
            yOcc[s] = fin[s] ? (wet[s] ? 1.0 : 0.0) - ps : Double.NaN;
            yAmt[s] = wet[s] ? NpMath.minimum(v / NpMath.maximum(ist, c.minDailyNormalMm), c.maxPrecipRatio) : Double.NaN;
        }
        int n = 0; int[] a = new int[S]; for (int s = 0; s < S; s++) if (fin[s] && usable[s]) a[n++] = s; d.act = Arrays.copyOf(a, n);
        if (d.act.length < c.minActiveStations) { d.skipped = true; return; }
        int nw = 0; int[] w = new int[n]; for (int s : d.act) if (wet[s] && Double.isFinite(yAmt[s])) w[nw++] = s; d.actW = Arrays.copyOf(w, nw);
        double[] gOcc = c.globalVariogramFallback ? DayKriging.globalVariogram(ks, sdmat, S, d.act, yOcc, d.tag + "|occ") : null;
        double[] gAmt = (c.globalVariogramFallback && nw >= 10) ? DayKriging.globalVariogram(ks, sdmat, S, d.actW, yAmt, d.tag + "|amt") : null;
        d.k1 = new DayKriging(ks, sd.x, sd.y, sd.z, sdmat, d.act, yOcc, false, DayKriging.Lapse.NONE, Double.NaN, d.tag + "|occ", gOcc);
        if (nw >= c.nKrig + 1) d.k2 = new DayKriging(ks, sd.x, sd.y, sd.z, sdmat, d.actW, yAmt, false, DayKriging.Lapse.NONE, Double.NaN, d.tag + "|amt", gAmt);
        else if (nw > 0) { double[] v = new double[nw]; for (int i = 0; i < nw; i++) v[i] = yAmt[d.actW[i]]; Arrays.sort(v);
            d.r0 = nw % 2 == 1 ? v[nw / 2] : (v[nw / 2 - 1] + v[nw / 2]) / 2.0; }
        else d.r0 = 0.0;
    }

    // ======================================================================================== run
    public void run() throws Exception {
        init();
        ExecutorService ex = Executors.newFixedThreadPool(Math.max(1, c.threads));
        ExecutorService wr = Executors.newSingleThreadExecutor();      // outputs written in date order while the next block runs
        long t0 = System.nanoTime();
        try {
            Outputs o = new Outputs();
            int block = Math.max(4 * c.threads, 8);
            List<Future<?>> pending = new ArrayList<>();
            for (int b0 = 0; b0 < days.size(); b0 += block) {
                List<Integer> blk = days.subList(b0, Math.min(days.size(), b0 + block));
                List<Future<Day>> fd = new ArrayList<>();
                for (int di : blk) fd.add(ex.submit(() -> prepare(di)));
                List<Day> ds = new ArrayList<>(); for (Future<Day> f : fd) ds.add(f.get());
                List<Callable<Void>> all = new ArrayList<>(); for (Day d : ds) all.addAll(d.tasks);
                for (Future<Void> f : ex.invokeAll(all)) f.get();
                long tb = System.nanoTime();
                for (Day d : ds) d.seconds = (tb - d.t0) / 1e9;
                for (Future<?> f : pending) f.get();          // at most one block of outputs in flight (bounded memory)
                pending.clear();
                pending.add(wr.submit(() -> { for (Day d : ds) o.day(d); return null; }));
            }
            for (Future<?> f : pending) f.get();
            o.finish();
        } finally { ex.shutdown(); wr.shutdown(); }
        info.put("wall_s", (System.nanoTime() - t0) / 1e9);
        writeRunConfig();
    }


    // ======================================================================================== output selection
    /** selected value-file (selV) and diagnostics-file (selG) NetCDF variables, in canonical order */
    final Set<String> selV = new java.util.LinkedHashSet<>(), selG = new java.util.LinkedHashSet<>();
    final Map<String, NcOut.Var> encV = new LinkedHashMap<>();
    final List<String> diagAllowed = new ArrayList<>();
    boolean outDefault, needOccAnom, needAmtRatio;

    String finalName() { return precip ? "precipitation" : "temperature"; }

    static NcOut.Var f4(String n, String u, String l) { return new NcOut.Var(n, "f4", 1, 0, u, l); }

    void setupOutputs() {
        String vname = finalName(), vunits = precip ? "mm day-1" : "degC";
        String ty = c.quantise ? "i2" : "f4"; double sc = c.quantise ? (precip ? 0.1 : 0.01) : 1;
        encV.put(vname, new NcOut.Var(vname, ty, sc, -32768, vunits, "daily " + vname));
        encV.put("kriging_std", new NcOut.Var("kriging_std", ty, sc, -32768, vunits, twopart ? "standard deviation of the two-part estimate (occurrence x amount)" : "ordinary-kriging standard deviation (back-transformed)"));
        if (ukTable != null) encV.put(c.variable + "_uncertainty", new NcOut.Var(c.variable + "_uncertainty", ty, sc, -32768, vunits, "calibrated uncertainty = kriging_std x k(season, class)"));
        if (precip && !clim) encV.put("estimate_raw", f4("estimate_raw", "mm day-1", "ordinary-kriging estimate before negative values are set to 0 (precipitation = max(estimate_raw, 0))"));
        if (precip && clim && !twopart) {
            encV.put("normal_daily", f4("normal_daily", "mm day-1", "daily climatology used as divisor = max(monthly normal / days in month, minDailyNormalMm)"));
            encV.put("ratio", f4("ratio", "1", "kriged daily precipitation ratio before back-transform (precipitation = max(ratio, 0) x normal_daily)"));
        }
        if (twopart) {
            encV.put("wet_probability", f4("wet_probability", "1", "occurrence index pi = clip(wetfreq_normal + occurrence_anomaly, 0, 1) before calibration, rounded to 1e-4 exactly as compared with occurrence_threshold (float32)"));
            encV.put("occurrence_anomaly", f4("occurrence_anomaly", "1", "kriged occurrence anomaly (station wet indicator minus station wet-day frequency normal)"));
            encV.put("occurrence_threshold", f4("occurrence_threshold", "1", "calibrated occurrence threshold of the cell (per month if monthly calibration); wet_flag = wet_probability >= threshold"));
            encV.put("wet_flag", new NcOut.Var("wet_flag", "u1", 1, 255, "1", "wet day after the occurrence calibration (1 = wet, 0 = dry)"));
            encV.put("amount", new NcOut.Var("amount", ty, sc, -32768, "mm day-1", "conditional wet-day amount = max(amount_ratio, 0) x intensity_normal, 0.1 mm resolution, before the wet flag (precipitation = wet_flag x amount)"));
            encV.put("amount_std", new NcOut.Var("amount_std", ty, sc, -32768, "mm day-1", "kriging standard deviation of the wet-day amount = sqrt(kriging variance of amount_ratio) x intensity_normal"));
            encV.put("amount_ratio", f4("amount_ratio", "1", "kriged wet-day intensity ratio (station amount / station intensity normal); median of the wet-station ratios on days with fewer than nKrig+1 wet stations"));
            encV.put("wetfreq_normal", f4("wetfreq_normal", "1", "wet-day frequency normal p(x) of the month (wet = precipitation >= occThreshold)"));
            encV.put("intensity_normal", f4("intensity_normal", "mm day-1", "wet-day intensity normal i(x) = M / (max(p, pFloor) x days in month) of the month"));
        }
        if (!precip) {
            if (clim) {
                encV.put("normal_daily", f4("normal_daily", "degC", "monthly temperature normal of the month, used unchanged for every day (temperature = normal_daily + anomaly)"));
                encV.put("anomaly", f4("anomaly", "degC", "kriged daily anomaly from the monthly normal, lapse term included (anomaly = residual + lapse_trend)"));
            }
            encV.put("residual", f4("residual", "degC", "kriged residual after removing the lapse trend (equals the kriged " + (clim ? "anomaly" : "temperature") + " when there is no lapse)"));
            if (lapseLocal || lapseGlobal) encV.put("lapse_trend", f4("lapse_trend", "degC", lapseLocal
                    ? "local lapse trend at the cell b0 + b1 x z" + (clim ? " (z clipped to the neighbour elevation range +- anomLapseMaxExtrapM)" : "")
                    : "global daily lapse trend at the cell b0 + b1 x z (one OLS on all active stations)"));
        }
        for (String k : ENC.keySet()) {
            if (k.startsWith("lapse")) { if (!precip && (lapseLocal || lapseGlobal)) diagAllowed.add(k); continue; }
            if (k.endsWith("_amount")) { if (twopart) diagAllowed.add(k); continue; }
            if (k.equals("wet_probability")) continue;                 // value-file variable when listed explicitly
            diagAllowed.add(k);
        }
        String ov = c.outputVariables == null ? "default" : c.outputVariables.trim();
        outDefault = ov.isEmpty() || ov.equalsIgnoreCase("default");
        if (outDefault) {                                               // behaviour before outputVariables existed
            selV.add(vname); selV.add("kriging_std"); if (ukTable != null) selV.add(c.variable + "_uncertainty");
            if (c.diagnostics) for (String k : ENC.keySet()) {
                if (precip && k.startsWith("lapse")) continue;
                if (!twopart && (k.endsWith("_amount") || k.equals("wet_probability"))) continue;
                selG.add(k);
            }
        } else if (ov.equalsIgnoreCase("all")) { selV.addAll(encV.keySet()); selG.addAll(diagAllowed); }
        else {
            Set<String> want = new HashSet<>();
            for (String tok : ov.split("[,;\\s]+")) {
                String t = tok.trim(); if (t.isEmpty()) continue;
                if (t.equals(c.variable)) t = vname;
                if (!encV.containsKey(t) && !diagAllowed.contains(t))
                    throw new IllegalArgumentException("outputVariables: '" + tok.trim() + "' is not available for " + methodKey() + ". Allowed: default | all | "
                            + allowedList() + (ukTable == null && t.equals(c.variable + "_uncertainty") ? " (" + t + " needs uncertaintyTable)" : ""));
                want.add(t);
            }
            for (String k : encV.keySet()) if (want.contains(k)) selV.add(k);
            for (String k : diagAllowed) if (want.contains(k)) selG.add(k);
        }
        needOccAnom = twopart && selV.contains("occurrence_anomaly"); needAmtRatio = twopart && selV.contains("amount_ratio");
        if (!outDefault && !c.outputFormat.equals("netcdf")) LOG.warning("outputVariables applies to NetCDF grid output only; outputFormat = " + c.outputFormat + " writes the usual CSV columns");
        LOG.info("PoKriging outputs: " + c.variable + "_<Y>.nc " + selV + (selG.isEmpty() ? "" : ", diagnostics_" + c.variable + "_<Y>.nc " + selG));
    }

    String methodKey() {
        return "variable=" + c.variable + (precip ? " precipMethod=" + c.precipMethod : " tempMethod=" + c.tempMethod + " temperatureTrend=" + c.temperatureTrend);
    }

    String allowedList() {
        List<String> l = new ArrayList<>();
        for (String k : encV.keySet()) l.add(k.equals(finalName()) ? c.variable : k);
        l.addAll(diagAllowed);
        return String.join(",", l);
    }

    static Map<String, double[]> select(Map<String, double[]> f, Set<String> keys, int n) {
        Map<String, double[]> o = new LinkedHashMap<>(); for (String k : keys) { double[] a = f.get(k); o.put(k, a != null ? a : nan(n)); } return o;
    }

    /** intermediate fields of the single-field methods at the targets of r (climT = monthly normal, trendT = global trend). */
    Map<String, double[]> components(TargetResult r, double[] climT, double[] trendT, double nd) {
        int n = r.n; Map<String, double[]> f = new LinkedHashMap<>();
        if (precip && !clim) f.put("estimate_raw", r.est.clone());
        else if (precip) {
            double[] cd = new double[n]; for (int i = 0; i < n; i++) cd[i] = NpMath.maximum(climT[i] / nd, c.minDailyNormalMm);
            f.put("normal_daily", cd); f.put("ratio", r.est.clone());
        } else {
            double[] an = new double[n], lt = new double[n];
            for (int i = 0; i < n; i++) { double g = trendT != null ? trendT[i] : 0.0; an[i] = r.est[i] + g; lt[i] = r.trend[i] + g; }
            if (clim) { f.put("normal_daily", climT.clone()); f.put("anomaly", an); }
            f.put("residual", r.resid.clone()); f.put("lapse_trend", lt);
        }
        return f;
    }

    // ======================================================================================== outputs
    static final String[] DIAG = {"D1", "D7", "D50", "range", "nugget_fraction", "log10_sill", "lapse", "lapse_r2", "n_negative_weights", "n_eff",
            "R_pos_half", "R_pos_2_3", "R_pos_1_3", "R_abs90", "qc_flags"};
    static final Map<String, NcOut.Var> ENC = new LinkedHashMap<>();
    static {
        ENC.put("D1", new NcOut.Var("D1", "u2", 0.1, 65535, "km", "distance to nearest active station"));
        ENC.put("D7", new NcOut.Var("D7", "u2", 0.1, 65535, "km", "distance to 7th nearest active station (kriging neighbourhood radius)"));
        ENC.put("D50", new NcOut.Var("D50", "u2", 0.5, 65535, "km", "distance to 50th nearest active station (variogram neighbourhood radius)"));
        ENC.put("range", new NcOut.Var("range", "u2", 0.5, 65535, "km", "spherical variogram range"));
        ENC.put("nugget_fraction", new NcOut.Var("nugget_fraction", "u1", 0.01, 255, "1", "nugget / (nugget + partial sill)"));
        ENC.put("log10_sill", new NcOut.Var("log10_sill", "i2", 0.01, -32768, "log10(units^2)", "log10 of total sill"));
        ENC.put("lapse", new NcOut.Var("lapse", "i2", 0.01, -32768, "degC km-1", "local daily temperature-elevation slope"));
        ENC.put("lapse_r2", new NcOut.Var("lapse_r2", "u1", 0.01, 255, "1", "R2 of the local daily lapse regression"));
        ENC.put("n_negative_weights", new NcOut.Var("n_negative_weights", "u1", 1, 255, "1", "number of negative kriging weights"));
        ENC.put("n_eff", new NcOut.Var("n_eff", "u1", 0.05, 255, "1", "effective number of stations (inverse Herfindahl of |w|)"));
        ENC.put("R_pos_half", new NcOut.Var("R_pos_half", "u2", 0.5, 65535, "km", "farthest station with weight >= 1/2 of the largest positive weight"));
        ENC.put("R_pos_2_3", new NcOut.Var("R_pos_2_3", "u2", 0.5, 65535, "km", "farthest station with weight >= 2/3 of the largest positive weight"));
        ENC.put("R_pos_1_3", new NcOut.Var("R_pos_1_3", "u2", 0.5, 65535, "km", "farthest station with weight >= 1/3 of the largest positive weight"));
        ENC.put("R_abs90", new NcOut.Var("R_abs90", "u2", 0.5, 65535, "km", "distance holding 90% of the total absolute weight"));
        ENC.put("qc_flags", new NcOut.Var("qc_flags", "u1", 1, 0, "1", "bit flags: 1=zero_shortcut, 2=lapse_fallback, 4=positive_lapse, 8=oscillatory_weights, 16=range_bound_hit, 32=global_variogram_fallback, 64=elevation_filled, 128=station_normal_fallback"));
        ENC.put("range_amount", new NcOut.Var("range_amount", "u2", 0.5, 65535, "km", "spherical variogram range, wet-day intensity ratios"));
        ENC.put("nugget_fraction_amount", new NcOut.Var("nugget_fraction_amount", "u1", 0.01, 255, "1", "nugget / (nugget + partial sill), wet-day intensity ratios"));
        ENC.put("log10_sill_amount", new NcOut.Var("log10_sill_amount", "i2", 0.01, -32768, "log10(1)", "log10 of total sill, wet-day intensity ratios"));
        ENC.put("R_pos_half_amount", new NcOut.Var("R_pos_half_amount", "u2", 0.5, 65535, "km", "farthest wet station with weight >= 1/2 of the largest positive weight (amount kriging)"));
        ENC.put("wet_probability", new NcOut.Var("wet_probability", "u1", 1.0 / 250, 255, "1", "occurrence index pi = p(x) + kriged occurrence anomaly (before calibration)"));
    }

    /** diagnostic fields of a target result as in s3_run (python names). */
    Map<String, double[]> diagFields(TargetResult r, TargetResult r2, int[] cellIdx) {
        Map<String, double[]> f = new LinkedHashMap<>(); int n = r.n;
        double[] nf = new double[n], ls = new double[n], lp = new double[n], qc = new double[n];
        for (int i = 0; i < n; i++) {
            double tot = r.nug[i] + r.psill[i]; double t = tot > 0 ? tot : Double.NaN;
            nf[i] = r.nug[i] / t; ls[i] = Math.log10(t); lp[i] = r.b1[i] * 1000;
            qc[i] = r.flags[i] | ((cellIdx != null && grid.elevFilled[cellIdx[i]]) ? KrigingSettings.FLAG_ELEV_FILLED : 0);
        }
        f.put("D1", r.d1); f.put("D7", r.d7); f.put("D50", r.d50); f.put("range", r.range); f.put("nugget_fraction", nf); f.put("log10_sill", ls);
        if (!precip) { f.put("lapse", lp); f.put("lapse_r2", r.r2); }
        f.put("n_negative_weights", r.nNeg); f.put("n_eff", r.neff); f.put("R_pos_half", r.rpos12); f.put("R_pos_2_3", r.rpos23);
        f.put("R_pos_1_3", r.rpos13); f.put("R_abs90", r.rabs90); f.put("qc_flags", qc);
        return f;
    }

    void amountDiag(Map<String, double[]> f, TargetResult r2, int n) {
        double[] ra = nan(n), nfa = nan(n), lsa = nan(n), rpa = nan(n);
        if (r2 != null) for (int i = 0; i < n; i++) {
            double tot = r2.nug[i] + r2.psill[i]; double t = tot > 0 ? tot : Double.NaN;
            ra[i] = r2.range[i]; nfa[i] = r2.nug[i] / t; lsa[i] = Math.log10(t); rpa[i] = r2.rpos12[i];
        }
        f.put("range_amount", ra); f.put("nugget_fraction_amount", nfa); f.put("log10_sill_amount", lsa); f.put("R_pos_half_amount", rpa);
    }

    static double[] nan(int n) { double[] a = new double[n]; Arrays.fill(a, Double.NaN); return a; }

    /** back-transform of the single-field methods: returns {value, kstd}. */
    double[][] back(TargetResult r, double[] climT, double[] trendT, double nd) {
        int n = r.n; double[] v = new double[n], s = new double[n];
        for (int i = 0; i < n; i++) {
            double est = r.est[i] + (trendT != null ? trendT[i] : 0.0), ksd = Math.sqrt(r.kvar[i]);
            if (!clim) { v[i] = precip ? NpMath.maximum(est, 0.0) : est; s[i] = ksd; }
            else if (precip) { double cd = NpMath.maximum(climT[i] / nd, c.minDailyNormalMm); v[i] = NpMath.maximum(est, 0.0) * cd; s[i] = ksd * cd; }
            else { v[i] = est + climT[i]; s[i] = ksd; }
        }
        return new double[][]{v, s};
    }

    double[] trendAt(Day d, TargetSet ts) {
        if (!lapseGlobal) return null;
        double[] t = new double[ts.n]; for (int i = 0; i < ts.n; i++) t[i] = d.gb0 + d.gb1 * ts.z[i]; return t;
    }

    static final String[] SEAS = {"DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF"};
    static final String[] PCLS = {"0", "0-1", "1-5", "5-20", "20-50", ">50"}, TCLS = {"<500", "500-1000", "1000-1500", "1500-2000", ">=2000"};

    void readUncertainty(Path p) throws IOException {
        Csv t = Csv.read(p); ukTable = new double[5][6];
        for (double[] r : ukTable) Arrays.fill(r, Double.NaN);
        List<String> seas = Arrays.asList("DJF", "MAM", "JJA", "SON", "all");
        String[] cls = precip ? PCLS : TCLS;
        for (int i = 0; i < t.rows.size(); i++) {
            int si = seas.indexOf(t.get(i, "season")); String cl = t.get(i, "cls");
            int ci = cl.equals("all") ? 5 : Arrays.asList(cls).indexOf(cl);
            if (si >= 0 && ci >= 0) ukTable[si][ci] = t.num(i, "scale_k");
        }
    }

    double uncertaintyK(int month0, double est, double elev) {
        int si = Arrays.asList("DJF", "MAM", "JJA", "SON").indexOf(SEAS[month0]); int ci;
        if (precip) { if (!(est > -1)) return Double.NaN; ci = est <= 0 ? 0 : est <= 1 ? 1 : est <= 5 ? 2 : est <= 20 ? 3 : est <= 50 ? 4 : 5; }
        else { ci = elev <= 500 ? 0 : elev <= 1000 ? 1 : elev <= 1500 ? 2 : elev <= 2000 ? 3 : 4; }
        double k = ukTable[si][ci];
        return Double.isNaN(k) ? ukTable[4][5] : k;
    }

    /** Output writers and two-part pass-1 store. */
    final class Outputs {
        final String vname = precip ? "precipitation" : "temperature", vunits = precip ? "mm day-1" : "degC";
        final double vscale = precip ? 0.1 : 0.01;
        final Map<Integer, List<Integer>> yearDays = new TreeMap<>();
        NcOut ncVal, ncDiag; PrintWriter gridCsv; int curYear = -1;
        PrintWriter looCsv, ptsCsv, sumCsv;
        // two-part pass 1
        List<int[]> p1Pi = new ArrayList<>(); List<short[]> p1Amt = new ArrayList<>(), p1Samt = new ArrayList<>(); List<Integer> p1Day = new ArrayList<>();
        List<float[]> p1Occ = new ArrayList<>(), p1Rat = new ArrayList<>();   // optional outputs occurrence_anomaly, amount_ratio
        List<double[][]> p1Raw = new ArrayList<>();                // csv test output: unquantised pi, amount, samount, p_target + diag rows
        List<Map<String, double[]>> p1Diag = new ArrayList<>();
        // two-part LOO / points tables (float32 as Python)
        final List<Object[]> looRows = new ArrayList<>(), ptRows = new ArrayList<>();
        final Map<Integer, Map<String, Object>> summary = new TreeMap<>();

        Outputs() { for (int di : days) yearDays.computeIfAbsent(sd.dates[di].getYear(), z -> new ArrayList<>()).add(di); }

        int timeIndex(int di) { return yearDays.get(sd.dates[di].getYear()).indexOf(di); }

        void openYear(int y) throws IOException {
            closeYear();
            curYear = y;
            List<Integer> dl = yearDays.get(y); int[] t = new int[dl.size()];
            for (int i = 0; i < t.length; i++) t[i] = (int) (sd.dates[dl.get(i)].toEpochDay() - EPOCH.toEpochDay());
            Map<String, String> at = new LinkedHashMap<>();
            at.put("source", "Po River District 1 km daily " + vname + " (GEOframe PoKriging, Java port of the v2.0/v5 engine)");
            at.put("method", methodText()); at.put("crs", c.crs); at.put("created", java.time.LocalDateTime.now().withNano(0).toString());
            String nc = (c.nccopy == null || c.nccopy.isEmpty()) ? null : c.nccopy;
            if (doGrid && c.outputFormat.equals("netcdf")) {
                if (!outDefault) at.put("output_variables", String.join(",", selV) + (selG.isEmpty() ? "" : " | diagnostics: " + String.join(",", selG)));
                if ((!twopart || thrTable != null) && !selV.isEmpty()) {
                    List<NcOut.Var> vl = new ArrayList<>(); for (String k : selV) vl.add(encV.get(k));
                    ncVal = new NcOut(out.resolve(c.variable + "_" + y + ".nc"), grid.xc, grid.yc, t, vl, at, c.crs, nc);
                }
                if (!selG.isEmpty()) {
                    List<NcOut.Var> vl = new ArrayList<>(); for (String k : selG) vl.add(ENC.get(k));
                    ncDiag = new NcOut(out.resolve("diagnostics_" + c.variable + "_" + y + ".nc"), grid.xc, grid.yc, t, vl, at, c.crs, nc);
                }
            }
            if (doGrid && c.outputFormat.equals("csv")) gridCsv = writer(out.resolve("grid_" + c.variable + "_" + y + ".csv.gz"));
            if (doLoo && !twopart) { looCsv = writer(out.resolve("loo_" + c.variable + "_" + y + ".csv"));
                looCsv.println("date,station,obs,pred,kriging_std,D1,D7,qc_flags" + (ukTable != null ? ",uncertainty" : "")); }
            if (doPoints && !twopart) { ptsCsv = writer(out.resolve("points_" + c.variable + "_" + y + ".csv")); ptsCsv.println("id,date,value,kriging_std,D1,D7,qc_flags"); }
        }

        PrintWriter writer(Path p) throws IOException {
            java.io.OutputStream os = Files.newOutputStream(p);
            if (p.toString().endsWith(".gz")) os = new GZIPOutputStream(os, 1 << 16);
            return new PrintWriter(new BufferedWriter(new OutputStreamWriter(os, StandardCharsets.UTF_8), 1 << 16));
        }

        void closeYear() throws IOException {
            if (ncVal != null) { ncVal.close(); ncVal = null; }
            if (ncDiag != null) { ncDiag.close(); ncDiag = null; }
            if (gridCsv != null) { gridCsv.close(); gridCsv = null; }
            if (looCsv != null) { looCsv.close(); looCsv = null; }
            if (ptsCsv != null) { ptsCsv.close(); ptsCsv = null; }
            if (twopart && thrTable == null && curYear >= 0 && !p1Day.isEmpty()) spillPass1(curYear);
            if (curYear >= 0) writeSummary(curYear);
        }

        void day(Day d) throws IOException {
            int y = d.date.getYear();
            if (y != curYear) openYear(y);
            int ti = timeIndex(d.di);
            Map<String, Object> s = new LinkedHashMap<>(); s.put("date", d.date.toString()); s.put("n_active", d.act == null ? 0 : d.act.length);
            summary.put(d.di, s);
            if (twopart) dayTwoPart(d, ti, s); else daySingle(d, ti, s);
            s.put("seconds", d.seconds);
        }

        void writeCells(NcOut nc, Map<String, double[]> f, int ti) throws IOException {
            int[] iy = new int[cells.length], ix = new int[cells.length];
            for (int i = 0; i < cells.length; i++) { iy[i] = grid.iy[cells[i]]; ix[i] = grid.ix[cells[i]]; }
            for (Map.Entry<String, double[]> e : f.entrySet()) nc.writeSlice(e.getKey(), ti, e.getValue(), iy, ix);
        }

        void daySingle(Day d, int ti, Map<String, Object> s) throws IOException {
            if (doGrid) {
                double[][] vk; Map<String, double[]> dg; Map<String, double[]> comp = new LinkedHashMap<>();
                if (d.skipped) { vk = new double[][]{nan(cells.length), nan(cells.length)}; dg = null; }
                else {
                    double[] ct = null;
                    if (clim) { ct = new double[cells.length]; for (int i = 0; i < cells.length; i++) ct[i] = NG[d.m][cells[i]]; }
                    vk = back(d.rG, ct, trendAt(d, tsGrid), d.nd); dg = diagFields(d.rG, null, cells);
                    gridStats(d, s, vk[0]);
                    if (!outDefault && ncVal != null) comp = components(d.rG, ct, trendAt(d, tsGrid), d.nd);
                }
                double[] unc = null;
                if (ukTable != null) { unc = new double[cells.length]; for (int i = 0; i < cells.length; i++) unc[i] = vk[1][i] * uncertaintyK(d.m, vk[0][i], tsGrid.z[i]); }
                if (ncVal != null) {
                    Map<String, double[]> f = new LinkedHashMap<>(); f.put(vname, vk[0]); f.put("kriging_std", vk[1]); if (unc != null) f.put(c.variable + "_uncertainty", unc);
                    f.putAll(comp);
                    writeCells(ncVal, select(f, selV, cells.length), ti);
                }
                if (ncDiag != null) { if (dg == null) { dg = new LinkedHashMap<>(); for (String k : DIAG) if (!(precip && k.startsWith("lapse"))) dg.put(k, nan(cells.length)); }
                    Map<String, double[]> g = dg;
                    if (!outDefault && lapseGlobal && !d.skipped) {           // global lapse: the day's single regression in the lapse diagnostics
                        g = new LinkedHashMap<>(dg); double[] lp = new double[cells.length], lr = new double[cells.length];
                        Arrays.fill(lp, d.gb1 * 1000); Arrays.fill(lr, d.gr2); g.put("lapse", lp); g.put("lapse_r2", lr);
                    }
                    writeCells(ncDiag, select(g, selG, cells.length), ti); }
                if (gridCsv != null && !d.skipped) gridCsvRows(d, vk[0], vk[1], dg, null);
            }
            if (d.skipped) { s.put("skipped", true); return; }
            if (doPoints) {
                double[] ct = null;
                if (clim) { ct = new double[tsPts.n]; for (int i = 0; i < tsPts.n; i++) ct[i] = NG[d.m][pCell[i]]; }
                double[][] vk = back(d.rP, ct, trendAt(d, tsPts), d.nd);
                for (int i = 0; i < tsPts.n; i++) ptsCsv.println(pid[i] + "," + d.date + "," + fmt(vk[0][i]) + "," + fmt(vk[1][i]) + "," + fmt(d.rP.d1[i]) + "," + fmt(d.rP.d7[i]) + "," + d.rP.flags[i]);
            }
            if (doLoo) {
                int n = d.tsLoo.n; double[] ct = null;
                if (clim) { ct = new double[n]; for (int i = 0; i < n; i++) ct[i] = NL[d.act[i]][d.m]; }
                double[][] vk = back(d.rL, ct, trendAt(d, d.tsLoo), d.nd);
                double[] obs = new double[n];
                for (int i = 0; i < n; i++) {
                    int st = d.act[i]; obs[i] = d.vals[st];
                    int fl = d.rL.flags[i] | (nsFallback[st][d.m] ? KrigingSettings.FLAG_NORMAL_FALLBACK : 0);
                    looCsv.print(d.date + "," + sd.ids[st] + "," + (float) obs[i] + "," + (float) vk[0][i] + "," + (float) vk[1][i] + "," + (float) d.rL.d1[i] + "," + (float) d.rL.d7[i] + "," + fl);
                    if (ukTable != null) looCsv.print("," + (float) (vk[1][i] * uncertaintyK(d.m, vk[0][i], sd.z[st])));
                    looCsv.println();
                }
                looMetrics(s, vk[0], obs);
            }
        }

        void gridStats(Day d, Map<String, Object> s, double[] val) {
            double sum = 0; int n = 0; for (double v : val) if (Double.isFinite(v)) { sum += v; n++; }
            s.put("field_mean", n > 0 ? sum / n : Double.NaN);
        }

        void looMetrics(Map<String, Object> s, double[] pred, double[] obs) {
            int n = 0, nw = 0; double se = 0, sae = 0, sse = 0, swae = 0, swe = 0, swo = 0;
            List<Integer> ok = new ArrayList<>();
            for (int i = 0; i < pred.length; i++) {
                double e = pred[i] - obs[i]; if (!Double.isFinite(e)) continue;
                ok.add(i); n++; se += e; sae += Math.abs(e); sse += e * e;
                if (precip && obs[i] >= c.wetThresholdStats) { nw++; swae += Math.abs(e); swe += e; swo += obs[i]; }
            }
            s.put("loo_n", n); s.put("loo_bias", se / n); s.put("loo_mae", sae / n); s.put("loo_rmse", Math.sqrt(sse / n));
            double r = Double.NaN;
            if (n > 3) { double mp = 0, mo = 0; for (int i : ok) { mp += pred[i]; mo += obs[i]; } mp /= n; mo /= n;
                double spo = 0, spp = 0, soo = 0; for (int i : ok) { spo += (pred[i] - mp) * (obs[i] - mo); spp += (pred[i] - mp) * (pred[i] - mp); soo += (obs[i] - mo) * (obs[i] - mo); }
                if (soo > 0) r = spo / Math.sqrt(spp * soo); }
            s.put("loo_r", r);
            if (precip) { s.put("loo_wet_n", nw); s.put("loo_wet_mae", nw > 0 ? swae / nw : Double.NaN); s.put("loo_wet_pbias", nw > 0 && swo > 0 ? 100 * swe / swo : Double.NaN); }
        }

        void gridCsvRows(Day d, double[] val, double[] kstd, Map<String, double[]> dg, Map<String, double[]> extra) {
            if (gridCsvHeader) {
                StringBuilder h = new StringBuilder("date,cell,value,kstd"); for (String k : dg.keySet()) h.append(',').append(k);
                if (extra != null) for (String k : extra.keySet()) h.append(',').append(k);
                gridCsv.println(h); gridCsvHeader = false;
            }
            for (int i = 0; i < cells.length; i++) {
                StringBuilder r = new StringBuilder(); r.append(d.date).append(',').append(cells[i]).append(',').append(fmt(val[i])).append(',').append(fmt(kstd[i]));
                for (double[] a : dg.values()) r.append(',').append(fmt(a[i]));
                if (extra != null) for (double[] a : extra.values()) r.append(',').append(fmt(a[i]));
                gridCsv.println(r);
            }
        }
        boolean gridCsvHeader = true;

        // ---------------------------------------------------------------- two-part
        void dayTwoPart(Day d, int ti, Map<String, Object> s) throws IOException {
            if (!d.skipped) { s.put("n_wet", d.actW.length); s.put("wet_share", d.actW.length / (double) d.act.length); }
            if (doGrid) {
                int n = cells.length; double[] pi = nan(n), amt = nan(n), samt = nan(n), pt = nan(n); Map<String, double[]> dg;
                double[] occA = nan(n), ratA = nan(n), itA = nan(n);
                if (!d.skipped) {
                    for (int i = 0; i < n; i++) {
                        int k = cells[i]; double p = PGF[d.m][k], it = NG[d.m][k] / (NpMath.maximum(p, c.pFloor) * d.nd); pt[i] = p;
                        occA[i] = d.rG.est[i]; ratA[i] = d.rG2 != null ? d.rG2.est[i] : d.r0; itA[i] = it;
                        pi[i] = NpMath.clip(p + d.rG.est[i], 0.0, 1.0);
                        if (d.rG2 != null) { amt[i] = NpMath.maximum(d.rG2.est[i], 0.0) * it; samt[i] = Math.sqrt(d.rG2.kvar[i]) * it; }
                        else { amt[i] = d.r0 * it; samt[i] = 0.0 * it; }
                    }
                    dg = diagFields(d.rG, null, cells); amountDiag(dg, d.rG2, n); dg.put("wet_probability", pi.clone());
                } else { dg = new LinkedHashMap<>(); for (String k : ENC.keySet()) if (!k.startsWith("lapse")) dg.put(k, nan(n)); }
                int[] qp = new int[n]; short[] qa = new short[n], qs = new short[n];
                for (int i = 0; i < n; i++) { qp[i] = TwoPart.encPi(pi[i]); qa[i] = TwoPart.encAmt(amt[i]); qs[i] = TwoPart.encAmt(samt[i]); }
                if (ncDiag != null) writeCells(ncDiag, select(dg, selG, n), ti);
                if (thrTable != null) {                        // thresholds known: final values now
                    double[] val = new double[n], ksd = new double[n], thr = new double[n], wet = new double[n];
                    for (int i = 0; i < n; i++) {
                        float pq = TwoPart.decPi(qp[i]); float t = thrTable[cells[i]][d.m]; thr[i] = t;
                        boolean w = pq >= t; wet[i] = w ? 1 : 0;
                        double a = TwoPart.decAmt(qa[i]), sa = TwoPart.decAmt(qs[i]);
                        val[i] = Float.isNaN(pq) ? Double.NaN : (w ? a : 0.0); ksd[i] = TwoPart.mixStd(pq, a, sa);
                    }
                    if (ncVal != null) { Map<String, double[]> f = new LinkedHashMap<>(); f.put(vname, val); f.put("kriging_std", ksd);
                        if (ukTable != null) f.put(c.variable + "_uncertainty", uncertainty(val, ksd, d.m));
                        if (!outDefault) tpComponents(f, qp, qa, qs, thr, d.skipped ? null : occA, d.skipped ? null : ratA, d.m, d.nd);
                        writeCells(ncVal, select(f, selV, n), ti); }
                    if (gridCsv != null && !d.skipped) gridCsvRows(d, val, ksd, dg, extraTp(pi, amt, samt, qp, qa, qs, pt, thr, wet));
                } else {
                    p1Day.add(d.di); p1Pi.add(qp); p1Amt.add(qa); p1Samt.add(qs);
                    if (needOccAnom) p1Occ.add(toFloat(occA)); if (needAmtRatio) p1Rat.add(toFloat(ratA));
                    if (gridCsv != null) { p1Raw.add(new double[][]{pi, amt, samt, pt}); p1Diag.add(d.skipped ? null : dg); }
                }
            }
            if (d.skipped) { s.put("skipped", true); return; }
            if (doPoints) addRows(d, d.rP, d.rP2, tsPts.n, ptRows, true);
            if (doLoo) addRows(d, d.rL, d.rL2, d.tsLoo.n, looRows, false);
        }

        /** two-part intermediate fields of one day (pass-1 values qp/qa/qs, thresholds thr per cell, optional occ / ratio). */
        void tpComponents(Map<String, double[]> f, int[] qp, short[] qa, short[] qs, double[] thr, double[] occ, double[] rat, int m, double nd) {
            int n = qp.length; double[] wp = new double[n], wf = new double[n], a = new double[n], sa = new double[n], th = new double[n], pn = new double[n], inn = new double[n];
            for (int i = 0; i < n; i++) {
                float pq = TwoPart.decPi(qp[i]); boolean fin = !Float.isNaN(pq); int k = cells[i];
                double p = PGF[m][k];
                wp[i] = pq; wf[i] = fin ? (pq >= thr[i] ? 1 : 0) : Double.NaN; a[i] = TwoPart.decAmt(qa[i]); sa[i] = TwoPart.decAmt(qs[i]);
                th[i] = thr[i]; pn[i] = fin ? p : Double.NaN; inn[i] = fin ? NG[m][k] / (NpMath.maximum(p, c.pFloor) * nd) : Double.NaN;
            }
            f.put("wet_probability", wp); if (occ != null) f.put("occurrence_anomaly", occ); f.put("occurrence_threshold", th); f.put("wet_flag", wf);
            f.put("amount", a); f.put("amount_std", sa); if (rat != null) f.put("amount_ratio", rat); f.put("wetfreq_normal", pn); f.put("intensity_normal", inn);
        }

        float[] toFloat(double[] a) { float[] f = new float[a.length]; for (int i = 0; i < a.length; i++) f[i] = (float) a[i]; return f; }
        double[] toDouble(float[] a) { double[] f = new double[a.length]; for (int i = 0; i < a.length; i++) f[i] = a[i]; return f; }

        double[] uncertainty(double[] val, double[] ksd, int m) {
            double[] u = new double[val.length]; for (int i = 0; i < val.length; i++) u[i] = ksd[i] * uncertaintyK(m, val[i], 0); return u;
        }

        Map<String, double[]> extraTp(double[] pi, double[] amt, double[] samt, int[] qp, short[] qa, short[] qs, double[] pt, double[] thr, double[] wet) {
            int n = pi.length; double[] piq = new double[n], aq = new double[n], sq = new double[n];
            for (int i = 0; i < n; i++) { piq[i] = TwoPart.decPi(qp[i]); aq[i] = TwoPart.decAmt(qa[i]); sq[i] = TwoPart.decAmt(qs[i]); }
            Map<String, double[]> e = new LinkedHashMap<>();
            e.put("pi", pi); e.put("amount", amt); e.put("samount", samt); e.put("pi_q", piq); e.put("amount_q", aq); e.put("samount_q", sq);
            e.put("p_target", pt); e.put("threshold", thr); e.put("wet", wet);
            return e;
        }

        /** LOO / points rows: key, month, date, obs, pi, amount, samount, p_target, D1, D7, flags (float32 as Python). */
        void addRows(Day d, TargetResult r1, TargetResult r2, int n, List<Object[]> rows, boolean points) {
            for (int i = 0; i < n; i++) {
                double p, it; int st = points ? -1 : d.act[i];
                if (points) { p = PGF[d.m][pCell[i]]; it = NG[d.m][pCell[i]] / (NpMath.maximum(p, c.pFloor) * d.nd); }
                else { p = PLF[st][d.m]; it = NL[st][d.m] / (NpMath.maximum(p, c.pFloor) * d.nd); }
                double pi = NpMath.clip(p + r1.est[i], 0.0, 1.0), amt, samt;
                if (r2 != null) { amt = NpMath.maximum(r2.est[i], 0.0) * it; samt = Math.sqrt(r2.kvar[i]) * it; } else { amt = d.r0 * it; samt = 0.0 * it; }
                int fl = r1.flags[i] | (!points && nsFallback[st][d.m] ? KrigingSettings.FLAG_NORMAL_FALLBACK : 0);
                rows.add(new Object[]{points ? i : st, d.date.getMonthValue(), d.di, (float) (points ? Double.NaN : d.vals[st]), (float) pi, (float) amt, (float) samt,
                        (float) p, (float) r1.d1[i], (float) r1.d7[i], fl});
            }
        }

        // pass-1 spill (gzip raw shorts) to keep memory bounded on long runs
        void spillPass1(int year) throws IOException {
            Path dir = out.resolve("_twopart"); Files.createDirectories(dir);
            try (DataOutputStream o = new DataOutputStream(new java.io.BufferedOutputStream(new GZIPOutputStream(Files.newOutputStream(dir.resolve("pass1_" + year + ".bin.gz")), 1 << 16)))) {
                o.writeInt(p1Day.size()); o.writeInt(cells.length);
                for (int k = 0; k < p1Day.size(); k++) {
                    o.writeInt(p1Day.get(k));
                    for (int v : p1Pi.get(k)) o.writeShort(v); for (short v : p1Amt.get(k)) o.writeShort(v); for (short v : p1Samt.get(k)) o.writeShort(v);
                }
            }
            if (needOccAnom || needAmtRatio)
                try (DataOutputStream o = new DataOutputStream(new java.io.BufferedOutputStream(new GZIPOutputStream(Files.newOutputStream(dir.resolve("pass1x_" + year + ".bin.gz")), 1 << 16)))) {
                    for (int k = 0; k < p1Day.size(); k++) {
                        if (needOccAnom) for (float v : p1Occ.get(k)) o.writeFloat(v);
                        if (needAmtRatio) for (float v : p1Rat.get(k)) o.writeFloat(v);
                    }
                }
            p1Day.clear(); p1Pi.clear(); p1Amt.clear(); p1Samt.clear(); p1Occ.clear(); p1Rat.clear();
        }

        /** optional pass-1 extras of a year: per day {occurrence_anomaly or null, amount_ratio or null}. */
        List<double[][]> readPass1x(int year, int nd, int n) throws IOException {
            List<double[][]> r = new ArrayList<>();
            Path f = out.resolve("_twopart").resolve("pass1x_" + year + ".bin.gz");
            try (DataInputStream in = new DataInputStream(new java.io.BufferedInputStream(new GZIPInputStream(Files.newInputStream(f), 1 << 16)))) {
                for (int k = 0; k < nd; k++) {
                    double[] oc = null, ra = null;
                    if (needOccAnom) { oc = new double[n]; for (int i = 0; i < n; i++) oc[i] = in.readFloat(); }
                    if (needAmtRatio) { ra = new double[n]; for (int i = 0; i < n; i++) ra[i] = in.readFloat(); }
                    r.add(new double[][]{oc, ra});
                }
            }
            return r;
        }

        List<Object[]> readPass1(int year, boolean piOnly) throws IOException {
            List<Object[]> r = new ArrayList<>();
            Path f = out.resolve("_twopart").resolve("pass1_" + year + ".bin.gz");
            try (DataInputStream in = new DataInputStream(new java.io.BufferedInputStream(new GZIPInputStream(Files.newInputStream(f), 1 << 16)))) {
                int nd = in.readInt(), n = in.readInt();
                for (int k = 0; k < nd; k++) {
                    int di = in.readInt(); int[] p = new int[n]; short[] a = new short[n], s = new short[n];
                    for (int i = 0; i < n; i++) p[i] = in.readShort() & 0xffff;
                    for (int i = 0; i < n; i++) a[i] = in.readShort(); for (int i = 0; i < n; i++) s[i] = in.readShort();
                    r.add(new Object[]{di, p, piOnly ? null : a, piOnly ? null : s});
                }
            }
            return r;
        }

        void finish() throws IOException {
            closeYear(); curYear = -1;
            if (twopart && doGrid && thrTable == null) finalizeGrid();
            if (twopart && doLoo) finalizeTable(looRows, false);
            if (twopart && doPoints) finalizeTable(ptRows, true);
            for (int y : yearDays.keySet()) writeSummary(y);
            if (twopart && !c.writeTwopartPass1) {
                Path dir = out.resolve("_twopart");
                if (Files.isDirectory(dir)) { try (java.util.stream.Stream<Path> st = Files.list(dir)) { for (Path p : (Iterable<Path>) st::iterator) Files.delete(p); } Files.delete(dir); }
            }
        }

        /**
         * Write the occurrence-calibration thresholds of this run as occurrence_thresholds_<start>_<end>.csv
         * (cell,month,threshold,wet_days_target,p_normal,calibration). A long monthly-calibrated run (e.g. 1991-2020)
         * produces the table that shorter runs read back with calibration = table. Pooled runs repeat their single
         * threshold for the 12 months (calibration = pooled); wet_days_target is then the pooled target of the run.
         */
        void writeThresholdTable(double[][] thr, boolean monthly, int[] mon, int[] perMonth, int D) throws IOException {
            java.util.List<Integer> fd = yearDays.get(yearDays.keySet().stream().min(Integer::compare).get());
            java.util.List<Integer> ld = yearDays.get(yearDays.keySet().stream().max(Integer::compare).get());
            String a = sd.dates[fd.get(0)].toString(), b = sd.dates[ld.get(ld.size() - 1)].toString();
            Path f = out.resolve("occurrence_thresholds_" + a + "_" + b + ".csv");
            try (PrintWriter w = writer(f)) {
                w.println("cell,month,threshold,wet_days_target,p_normal,calibration");
                for (int i = 0; i < cells.length; i++) {
                    int cell = cells[i]; double kPooled = 0;
                    if (!monthly) for (int k = 0; k < D; k++) kPooled += PGF[mon[k]][cell];
                    for (int m = 0; m < 12; m++) {
                        double t = monthly ? thr[i][m] : thr[i][0];
                        long target = Math.round(monthly ? PGF[m][cell] * perMonth[m] : kPooled);
                        w.println(cell + "," + (m + 1) + "," + (Double.isInfinite(t) ? "inf" : Double.toString((float) t)) + "," + target + ","
                                + (float) PGF[m][cell] + "," + (monthly ? "monthly" : "pooled"));
                    }
                }
            }
            info.put("threshold_table", f.getFileName().toString());
        }

        void finalizeGrid() throws IOException {
            int n = cells.length; List<Integer> dlist = new ArrayList<>(); List<int[]> pis = new ArrayList<>();
            for (int y : yearDays.keySet()) for (Object[] r : readPass1(y, true)) { dlist.add((Integer) r[0]); pis.add((int[]) r[1]); }
            int D = dlist.size(); int[] mon = new int[D]; int[] perMonth = new int[12];
            for (int k = 0; k < D; k++) { mon[k] = sd.dates[dlist.get(k)].getMonthValue() - 1; perMonth[mon[k]]++; }
            int mn = Integer.MAX_VALUE; for (int v : perMonth) mn = Math.min(mn, v);
            boolean monthly = c.calibration.equals("monthly") || (c.calibration.equals("auto") && mn >= c.calibMinDays);
            info.put("grid_calibration", monthly ? "monthly" : "pooled");
            double[][] thr = new double[n][monthly ? 12 : 1];
            float[] col = new float[D];
            for (int i = 0; i < n; i++) {
                int cell = cells[i];
                if (monthly) {
                    for (int m = 0; m < 12; m++) {
                        int q = 0, fin = 0;
                        for (int k = 0; k < D; k++) if (mon[k] == m) { col[q] = TwoPart.decPi(pis.get(k)[i]); if (!Float.isNaN(col[q])) fin++; q++; }
                        thr[i][m] = TwoPart.threshold(col, q, PGF[m][cell] * fin);
                    }
                } else {
                    double k_ = 0;
                    for (int k = 0; k < D; k++) { col[k] = TwoPart.decPi(pis.get(k)[i]); if (!Float.isNaN(col[k])) k_ += PGF[mon[k]][cell]; }
                    thr[i][0] = TwoPart.threshold(col, D, k_);
                }
            }
            pis = null;
            writeThresholdTable(thr, monthly, mon, perMonth, D);
            String ncp = (c.nccopy == null || c.nccopy.isEmpty()) ? null : c.nccopy;
            long wetN = 0, totN = 0; int pos = 0;
            for (int y : yearDays.keySet()) {
                List<Object[]> rows = readPass1(y, false);
                NcOut nc = null; PrintWriter gc = null; boolean ncMode = c.outputFormat.equals("netcdf");
                List<double[][]> ext = (ncMode && (needOccAnom || needAmtRatio)) ? readPass1x(y, rows.size(), n) : null;
                if (ncMode && !selV.isEmpty()) {
                    List<Integer> dl = yearDays.get(y); int[] t = new int[dl.size()];
                    for (int i = 0; i < t.length; i++) t[i] = (int) (sd.dates[dl.get(i)].toEpochDay() - EPOCH.toEpochDay());
                    String ty = c.quantise ? "i2" : "f4"; double sc = c.quantise ? 0.1 : 1;
                    List<NcOut.Var> vl = new ArrayList<>(); for (String k : selV) vl.add(encV.get(k));
                    Map<String, String> at = new LinkedHashMap<>(); at.put("method", methodText()); at.put("crs", c.crs); at.put("occurrence_calibration", monthly ? "monthly" : "pooled");
                    if (!outDefault) at.put("output_variables", String.join(",", selV));
                    nc = new NcOut(out.resolve(c.variable + "_" + y + ".nc"), grid.xc, grid.yc, t, vl, at, c.crs, ncp);
                } else if (!ncMode) gc = writer(out.resolve("grid_" + c.variable + "_" + y + ".csv.gz"));
                boolean header = true; int rk = -1;
                for (Object[] r : rows) { rk++;
                    int di = (Integer) r[0]; int[] qp = (int[]) r[1]; short[] qa = (short[]) r[2], qs = (short[]) r[3]; int m = sd.dates[di].getMonthValue() - 1;
                    double[] val = new double[n], ksd = new double[n], th = new double[n], wet = new double[n];
                    for (int i = 0; i < n; i++) {
                        float pq = TwoPart.decPi(qp[i]); double t = thr[i][monthly ? m : 0]; th[i] = t;
                        boolean w = pq >= t; wet[i] = w ? 1 : 0; if (!Float.isNaN(pq)) { totN++; if (w) wetN++; }
                        double a = TwoPart.decAmt(qa[i]), sa = TwoPart.decAmt(qs[i]);
                        val[i] = Float.isNaN(pq) ? Double.NaN : (w ? a : 0.0); ksd[i] = TwoPart.mixStd(pq, a, sa);
                    }
                    if (ncMode) {
                        if (nc == null) continue;
                        int[] iy = new int[n], ix = new int[n]; for (int i = 0; i < n; i++) { iy[i] = grid.iy[cells[i]]; ix[i] = grid.ix[cells[i]]; }
                        int ti = yearDays.get(y).indexOf(di);
                        Map<String, double[]> f = new LinkedHashMap<>(); f.put(vname, val); f.put("kriging_std", ksd);
                        if (ukTable != null) f.put(c.variable + "_uncertainty", uncertainty(val, ksd, m));
                        if (!outDefault) tpComponents(f, qp, qa, qs, th, ext != null ? ext.get(rk)[0] : null, ext != null ? ext.get(rk)[1] : null, m, sd.dates[di].lengthOfMonth());
                        for (Map.Entry<String, double[]> e : select(f, selV, n).entrySet()) nc.writeSlice(e.getKey(), ti, e.getValue(), iy, ix);
                    } else {
                        int k = p1Day.indexOf(di);   // never: spilled
                        double[][] raw = p1Raw.get(pos); Map<String, double[]> dg = p1Diag.get(pos); pos++;
                        if (dg == null) continue;
                        Map<String, double[]> ex = extraTp(raw[0], raw[1], raw[2], qp, qa, qs, raw[3], th, wet);
                        if (header) { StringBuilder h = new StringBuilder("date,cell,value,kstd"); for (String kk : dg.keySet()) h.append(',').append(kk);
                            for (String kk : ex.keySet()) h.append(',').append(kk); gc.println(h); header = false; }
                        for (int i = 0; i < n; i++) {
                            StringBuilder sb = new StringBuilder(); sb.append(sd.dates[di]).append(',').append(cells[i]).append(',').append(fmt(val[i])).append(',').append(fmt(ksd[i]));
                            for (double[] a : dg.values()) sb.append(',').append(fmt(a[i])); for (double[] a : ex.values()) sb.append(',').append(fmt(a[i]));
                            gc.println(sb);
                        }
                    }
                }
                if (nc != null) nc.close(); if (gc != null) gc.close();
            }
            info.put("grid_wet_fraction", totN > 0 ? wetN / (double) totN : Double.NaN);
        }

        void finalizeTable(List<Object[]> rows, boolean points) throws IOException {
            int n = rows.size(); int[] month = new int[n], dk = new int[n], key = new int[n]; float[] pi = new float[n], pt = new float[n];
            for (int i = 0; i < n; i++) { Object[] r = rows.get(i); key[i] = (Integer) r[0]; month[i] = (Integer) r[1]; dk[i] = (Integer) r[2]; pi[i] = (Float) r[4]; pt[i] = (Float) r[7]; }
            boolean[] wet = new boolean[n];
            String mode = TwoPart.calibrate(month, dk, key, pi, pt, c.calibMinDays, wet);
            info.put(points ? "points_calibration" : "loo_calibration", mode);
            Map<Integer, PrintWriter> w = new TreeMap<>(); Map<Integer, List<double[]>> perDay = new TreeMap<>();
            long wp = 0, wo = 0; double mae = 0; int nn = 0;
            for (int i = 0; i < n; i++) {
                Object[] r = rows.get(i); int di = dk[i]; int y = sd.dates[di].getYear();
                PrintWriter pw = w.get(y);
                if (pw == null) {
                    pw = writer(out.resolve((points ? "points_" : "loo_") + c.variable + "_" + y + ".csv")); w.put(y, pw);
                    pw.println(points ? "id,date,value,kriging_std,D1,D7,qc_flags,pi,amount" : "date,station,obs,pred,kriging_std,D1,D7,qc_flags,pi,amount,samount,p_target,wet_pred");
                }
                float amt = (Float) r[5], samt = (Float) r[6];
                float pred = wet[i] ? amt : 0f; float ksd = TwoPart.mixStdF(pi[i], amt, samt);
                if (points) pw.println(pid[key[i]] + "," + sd.dates[di] + "," + pred + "," + ksd + "," + r[8] + "," + r[9] + "," + r[10] + "," + pi[i] + "," + amt);
                else {
                    float obs = (Float) r[3];
                    pw.println(sd.dates[di] + "," + sd.ids[key[i]] + "," + obs + "," + pred + "," + ksd + "," + r[8] + "," + r[9] + "," + r[10] + "," + pi[i] + "," + amt + "," + samt + "," + pt[i] + "," + (wet[i] ? 1 : 0));
                    perDay.computeIfAbsent(di, z -> new ArrayList<>()).add(new double[]{pred, obs});
                    if (pred >= c.occThreshold) wp++; if (obs >= c.occThreshold) wo++; mae += Math.abs(pred - obs); nn++;
                }
            }
            for (PrintWriter pw : w.values()) pw.close();
            if (!points) {
                for (Map.Entry<Integer, List<double[]>> e : perDay.entrySet()) {
                    List<double[]> l = e.getValue(); double[] p = new double[l.size()], o = new double[l.size()];
                    for (int i = 0; i < l.size(); i++) { p[i] = l.get(i)[0]; o[i] = l.get(i)[1]; }
                    looMetrics(summary.get(e.getKey()), p, o);
                }
                info.put("loo_wet_fraction_pred", wp / (double) nn); info.put("loo_wet_fraction_obs", wo / (double) nn); info.put("loo_mae", mae / nn);
            }
        }

        void writeSummary(int y) throws IOException {
            List<Map<String, Object>> rows = new ArrayList<>();
            for (int di : yearDays.get(y)) if (summary.containsKey(di)) rows.add(summary.get(di));
            if (rows.isEmpty()) return;
            Set<String> keys = new java.util.LinkedHashSet<>(); for (Map<String, Object> r : rows) keys.addAll(r.keySet());
            try (PrintWriter pw = writer(out.resolve("daily_summary_" + c.variable + "_" + y + ".csv"))) {
                pw.println(String.join(",", keys));
                for (Map<String, Object> r : rows) {
                    StringBuilder sb = new StringBuilder(); boolean first = true;
                    for (String k : keys) { if (!first) sb.append(','); first = false; Object v = r.get(k); if (v != null) sb.append(v instanceof Double ? fmt((Double) v) : v.toString()); }
                    pw.println(sb);
                }
            }
        }
    }

    static String fmt(double v) {
        if (Double.isNaN(v)) return "";
        if (v == Double.POSITIVE_INFINITY) return "inf"; if (v == Double.NEGATIVE_INFINITY) return "-inf";
        if (v == Math.rint(v) && Math.abs(v) < 1e15) return Long.toString((long) v) + (v == 0 && 1 / v < 0 ? "" : "");
        return Double.toString(v);
    }

    String methodText() {
        String base = "local spherical variogram on the " + c.nVario + " nearest active stations, PSO " + c.psoParticles + "x" + c.psoIterations
                + ", ordinary kriging from the " + c.nKrig + " nearest; ";
        if (twopart) return "two-part method: occurrence = kriged occurrence anomalies (wet >= " + c.occThreshold + " mm) added to the wet-day frequency normal and "
                + "calibrated to it; amount = kriged wet-day intensity ratios times the intensity normal M/(p days); " + base;
        return base + (lapseLocal ? "per-cell daily lapse on the same stations; " : lapseGlobal ? "one daily lapse regression on all active stations; " : "")
                + (clim ? "anomaly method on monthly normals" : "direct interpolation");
    }

    void writeRunConfig() throws IOException {
        StringBuilder sb = new StringBuilder("{\n");
        for (java.lang.reflect.Field f : KrigingConfig.class.getFields()) {
            try { Object v = f.get(c); sb.append("  \"").append(f.getName()).append("\": ").append(v instanceof String ? "\"" + v.toString().replace("\\", "/") + "\"" : String.valueOf(v)).append(",\n"); }
            catch (IllegalAccessException e) { /* skip */ }
        }
        for (Map.Entry<String, Object> e : info.entrySet()) sb.append("  \"").append(e.getKey()).append("\": ").append(e.getValue() instanceof String ? "\"" + e.getValue() + "\"" : e.getValue()).append(",\n");
        sb.append("  \"method\": \"").append(methodText()).append("\",\n  \"version\": \"pokriging-java 1.0\"\n}\n");
        Files.write(out.resolve("run_config.json"), sb.toString().getBytes(StandardCharsets.UTF_8));
    }
}
