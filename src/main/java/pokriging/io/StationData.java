package pokriging.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.time.LocalDate;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Stations (stations.csv: id,x,y,elevation[,tpi20,facet,...]) and their daily values (wide data.csv: date + one
 * column per station id). Engine station order = data.csv column order (METHOD_SPEC §0).
 */
public final class StationData {
    public String[] ids; public double[] x, y, z, tpi; public int[] facet;
    public LocalDate[] dates; public float[][] values;   // [day][station], NaN = missing
    public int S, T;

    /**
     * @param missing missing code (in addition to the empty field); values numerically close to it (numpy isclose)
     *                are missing as well
     * @param variable precip (negative -> missing) or temp (outside [tMin, tMax] -> missing)
     */
    public static StationData read(Path stationsCsv, Path dataCsv, String missing, String variable, double tMin, double tMax) throws IOException {
        Csv st = Csv.read(stationsCsv);
        Map<String, Integer> row = new HashMap<>();
        for (int i = 0; i < st.rows.size(); i++) row.put(st.get(i, "id"), i);
        StationData d = new StationData();
        double miss = Double.NaN;
        try { if (missing != null && !missing.trim().isEmpty()) miss = Double.parseDouble(missing.trim()); } catch (NumberFormatException e) { /* string code */ }
        try (BufferedReader r = Csv.open(dataCsv)) {
            String[] h = Csv.split(Csv.nextLine(r));
            if (!h[0].trim().equalsIgnoreCase("date")) throw new IOException("first column of " + dataCsv + " must be 'date'");
            int S = h.length - 1; d.S = S; d.ids = new String[S];
            d.x = new double[S]; d.y = new double[S]; d.z = new double[S]; d.tpi = new double[S]; d.facet = new int[S];
            for (int j = 0; j < S; j++) {
                String id = h[j + 1].trim(); d.ids[j] = id;
                Integer i = row.get(id);
                if (i == null) throw new IOException("data column " + id + " has no row in " + stationsCsv);
                d.x[j] = st.num(i, "x"); d.y[j] = st.num(i, "y"); d.z[j] = st.num(i, "elevation");
                if (!Double.isFinite(d.z[j])) throw new IOException("station " + id + ": elevation missing");
                d.tpi[j] = st.has("tpi20") ? st.num(i, "tpi20") : 0.0;
                d.facet[j] = st.has("facet") ? (int) st.num(i, "facet") : 8;
            }
            List<LocalDate> dates = new ArrayList<>(); List<float[]> vals = new ArrayList<>();
            String l;
            boolean prec = variable.equals("precip");
            while ((l = Csv.nextLine(r)) != null) {
                String[] f = Csv.split(l);
                LocalDate dt = LocalDate.parse(f[0].trim());
                if (!dates.isEmpty() && !dt.equals(dates.get(dates.size() - 1).plusDays(1)))
                    throw new IOException("dates of " + dataCsv + " are not strictly daily at " + dt);
                float[] v = new float[S];
                for (int j = 0; j < S; j++) {
                    String s = j + 1 < f.length ? f[j + 1].trim() : "";
                    float x;
                    if (s.isEmpty() || s.equals(missing) || s.equalsIgnoreCase("nan")) x = Float.NaN;
                    else {
                        double dv = Double.parseDouble(s);              // pandas: parse float64 then cast to float32
                        x = (float) dv;
                        if (!Double.isNaN(miss) && Math.abs(x - miss) <= 1e-8 + 1e-5 * Math.abs(miss)) x = Float.NaN;
                    }
                    if (prec) { if (x < 0) x = Float.NaN; } else if (x < tMin || x > tMax) x = Float.NaN;
                    v[j] = x;
                }
                dates.add(dt); vals.add(v);
            }
            d.T = dates.size(); d.dates = dates.toArray(new LocalDate[0]); d.values = vals.toArray(new float[0][]);
        }
        return d;
    }

    /** float32 station-station distance matrix in km (po_common.station_distance_km), row-major. */
    public float[] distanceMatrixKm() {
        float[] m = new float[S * S];
        for (int i = 0; i < S; i++) for (int j = 0; j < S; j++) {
            double dx = x[i] - x[j], dy = y[i] - y[j];
            m[i * S + j] = (float) (Math.sqrt(dx * dx + dy * dy) / 1000.0);
        }
        return m;
    }

    public int dayIndex(LocalDate d) { long k = d.toEpochDay() - dates[0].toEpochDay(); return (k < 0 || k >= T) ? -1 : (int) k; }
}
