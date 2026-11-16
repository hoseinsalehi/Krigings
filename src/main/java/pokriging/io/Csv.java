package pokriging.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/** Small CSV reader (comma, no quoting, header line, '#' comment lines). Gzip if the name ends with .gz. */
public final class Csv {
    public final String[] header; public final List<String[]> rows = new ArrayList<>();
    private final Map<String, Integer> col = new HashMap<>();

    private Csv(String[] header) { this.header = header; for (int i = 0; i < header.length; i++) col.put(header[i].trim(), i); }

    public static BufferedReader open(Path p) throws IOException {
        if (p.toString().endsWith(".gz"))
            return new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(p), 1 << 16), StandardCharsets.UTF_8), 1 << 16);
        return Files.newBufferedReader(p, StandardCharsets.UTF_8);
    }

    public static Csv read(Path p) throws IOException {
        try (BufferedReader r = open(p)) {
            String h = nextLine(r);
            if (h == null) throw new IOException("empty CSV " + p);
            Csv c = new Csv(split(h));
            String l;
            while ((l = nextLine(r)) != null) c.rows.add(split(l));
            return c;
        }
    }

    static String nextLine(BufferedReader r) throws IOException {
        String l;
        while ((l = r.readLine()) != null) { if (l.isEmpty() || l.startsWith("#")) continue; if (l.charAt(0) == '\uFEFF') l = l.substring(1); return l; }
        return null;
    }

    /** split on commas keeping empty trailing fields. */
    public static String[] split(String l) {
        int n = 1; for (int i = 0; i < l.length(); i++) if (l.charAt(i) == ',') n++;
        String[] out = new String[n]; int k = 0, s = 0;
        for (int i = 0; i < l.length(); i++) if (l.charAt(i) == ',') { out[k++] = l.substring(s, i); s = i + 1; }
        out[k] = l.substring(s);
        return out;
    }

    public boolean has(String name) { return col.containsKey(name); }

    public int col(String name) {
        Integer i = col.get(name);
        if (i == null) throw new IllegalArgumentException("column '" + name + "' missing (have " + String.join(",", header) + ")");
        return i;
    }

    public String get(int row, String name) { String[] r = rows.get(row); int c = col(name); return c < r.length ? r[c].trim() : ""; }

    /** parse a double; empty / nan / NaN -> NaN. */
    public static double num(String s) {
        s = s.trim();
        if (s.isEmpty() || s.equalsIgnoreCase("nan")) return Double.NaN;
        if (s.equalsIgnoreCase("inf")) return Double.POSITIVE_INFINITY;
        if (s.equalsIgnoreCase("-inf")) return Double.NEGATIVE_INFINITY;
        if (s.equalsIgnoreCase("true")) return 1; if (s.equalsIgnoreCase("false")) return 0;
        return Double.parseDouble(s);
    }

    public double num(int row, String name) { return num(get(row, name)); }
}
