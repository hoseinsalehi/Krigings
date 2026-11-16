package pokriging.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Target grid: ESRI ASCII DEM (cell centres, NODATA outside the domain), valid cells in row-major order from the
 * northern row (= Python engine cell index), optional per-cell elevation-filled flag from grid_cells.csv.
 */
public final class Grid {
    public int nx, ny; public double[] xc, yc;                 // cell-centre coordinates (x ascending, y descending)
    public int n; public int[] iy, ix; public double[] x, y, z, tpi; public int[] facet; public boolean[] elevFilled;
    public int[][] cellIndex;                                  // [iy][ix] -> cell or -1

    public static Grid readAscii(Path asc) throws IOException {
        Grid g = new Grid();
        try (BufferedReader r = Files.newBufferedReader(asc)) {
            double xll = 0, yll = 0, cs = 0, nod = -9999; boolean centre = false;
            for (int k = 0; k < 6; k++) {
                String[] t = r.readLine().trim().split("\\s+"); String key = t[0].toLowerCase(); double v = Double.parseDouble(t[1]);
                switch (key) {
                    case "ncols": g.nx = (int) v; break; case "nrows": g.ny = (int) v; break;
                    case "xllcorner": xll = v; break; case "yllcorner": yll = v; break;
                    case "xllcenter": xll = v; centre = true; break; case "yllcenter": yll = v; centre = true; break;
                    case "cellsize": cs = v; break; case "nodata_value": nod = v; break;
                    default: throw new IOException("unexpected ASCII grid header " + key);
                }
            }
            double x0 = centre ? xll : xll + cs / 2, ytop = (centre ? yll : yll + cs / 2) + (g.ny - 1) * cs;
            g.xc = new double[g.nx]; g.yc = new double[g.ny];
            for (int i = 0; i < g.nx; i++) g.xc[i] = x0 + i * cs;
            for (int j = 0; j < g.ny; j++) g.yc[j] = ytop - j * cs;
            List<double[]> cells = new ArrayList<>();
            g.cellIndex = new int[g.ny][g.nx];
            int row = 0; String l; int c = 0;
            while (row < g.ny && (l = r.readLine()) != null) {
                String[] t = l.trim().split("\\s+");
                if (t.length == 1 && t[0].isEmpty()) continue;
                for (int i = 0; i < g.nx; i++) {
                    double v = Double.parseDouble(t[i]);
                    if (v == nod || Double.isNaN(v)) { g.cellIndex[row][i] = -1; continue; }
                    g.cellIndex[row][i] = c++; cells.add(new double[]{row, i, v});
                }
                row++;
            }
            g.n = cells.size(); g.iy = new int[g.n]; g.ix = new int[g.n]; g.x = new double[g.n]; g.y = new double[g.n]; g.z = new double[g.n];
            g.tpi = new double[g.n]; g.facet = new int[g.n]; g.elevFilled = new boolean[g.n];
            for (int k = 0; k < g.n; k++) {
                double[] e = cells.get(k); g.iy[k] = (int) e[0]; g.ix[k] = (int) e[1]; g.x[k] = g.xc[g.ix[k]]; g.y[k] = g.yc[g.iy[k]]; g.z[k] = e[2]; g.facet[k] = 8;
            }
        }
        return g;
    }

    /** Read elev_filled (and optional tpi20/facet) per cell from grid_cells.csv (id = cell index), checking x/y. */
    public void readCellAttributes(Path gridCellsCsv) throws IOException {
        Csv c = Csv.read(gridCellsCsv);
        if (c.rows.size() != n) throw new IOException(gridCellsCsv + ": " + c.rows.size() + " cells, DEM has " + n);
        for (int r = 0; r < n; r++) {
            int id = (int) c.num(r, "id");
            if (Math.abs(c.num(r, "x") - x[id]) > 1e-6 || Math.abs(c.num(r, "y") - y[id]) > 1e-6)
                throw new IOException("grid_cells.csv cell " + id + " does not match the DEM cell centre");
            if (c.has("elev_filled")) elevFilled[id] = c.num(r, "elev_filled") != 0;
            if (c.has("tpi20")) tpi[id] = c.num(r, "tpi20");
            if (c.has("facet")) facet[id] = (int) c.num(r, "facet");
        }
    }

    /** nearest grid cell of a point (rounding to the cell, clipped to the grid; -1 outside the valid domain). */
    public int nearestCellRounded(double px, double py) {
        double dx = xc.length > 1 ? xc[1] - xc[0] : 1, dy = yc.length > 1 ? yc[1] - yc[0] : -1;
        int i = (int) Math.max(0, Math.min(nx - 1, Math.rint((px - xc[0]) / dx)));
        int j = (int) Math.max(0, Math.min(ny - 1, Math.rint((py - yc[0]) / dy)));
        return cellIndex[j][i];
    }

    /**
     * Build the target grid from cell centres (any order) with their elevation, e.g. the cells of a basin or of the whole
     * domain read from a point or polygon shapefile. The cells must lie on a regular square lattice; the cell size is the
     * smallest spacing found. Cells absent from the list are outside the domain (NODATA). Valid cells are numbered row-major
     * from the northern row, as with readAscii, so a grid built here and the same grid written with writeAscii and read back
     * are identical.
     */
    public static Grid fromCells(double[] px, double[] py, double[] pz) throws IOException {
        int m = px.length; if (m == 0) throw new IOException("no target cells");
        double cs = Math.min(minSpacing(px), minSpacing(py));
        if (!(cs > 0) || !Double.isFinite(cs)) cs = 1000;      // a single cell or a single row/column: fall back to 1 km
        double xmin = Double.POSITIVE_INFINITY, ymax = Double.NEGATIVE_INFINITY, xmax = Double.NEGATIVE_INFINITY, ymin = Double.POSITIVE_INFINITY;
        for (int k = 0; k < m; k++) { xmin = Math.min(xmin, px[k]); xmax = Math.max(xmax, px[k]); ymin = Math.min(ymin, py[k]); ymax = Math.max(ymax, py[k]); }
        Grid g = new Grid();
        g.nx = (int) Math.rint((xmax - xmin) / cs) + 1; g.ny = (int) Math.rint((ymax - ymin) / cs) + 1;
        g.xc = new double[g.nx]; g.yc = new double[g.ny];
        for (int i = 0; i < g.nx; i++) g.xc[i] = xmin + i * cs;
        for (int j = 0; j < g.ny; j++) g.yc[j] = ymax - j * cs;
        double[][] zz = new double[g.ny][g.nx]; for (double[] r : zz) java.util.Arrays.fill(r, Double.NaN);
        for (int k = 0; k < m; k++) {
            double fi = (px[k] - xmin) / cs, fj = (ymax - py[k]) / cs; int i = (int) Math.rint(fi), j = (int) Math.rint(fj);
            if (Math.abs(fi - i) > 1e-6 || Math.abs(fj - j) > 1e-6)
                throw new IOException(String.format("cell (%.3f, %.3f) is not on the %.3f m lattice of the other cells", px[k], py[k], cs));
            if (!Double.isNaN(zz[j][i])) throw new IOException(String.format("two cells at (%.3f, %.3f)", px[k], py[k]));
            if (!Double.isFinite(pz[k])) throw new IOException(String.format("cell (%.3f, %.3f) has no elevation", px[k], py[k]));
            zz[j][i] = pz[k];
        }
        g.cellIndex = new int[g.ny][g.nx]; g.n = m;
        g.iy = new int[m]; g.ix = new int[m]; g.x = new double[m]; g.y = new double[m]; g.z = new double[m];
        g.tpi = new double[m]; g.facet = new int[m]; g.elevFilled = new boolean[m];
        int c = 0;
        for (int j = 0; j < g.ny; j++) for (int i = 0; i < g.nx; i++) {
            if (Double.isNaN(zz[j][i])) { g.cellIndex[j][i] = -1; continue; }
            g.cellIndex[j][i] = c; g.iy[c] = j; g.ix[c] = i; g.x[c] = g.xc[i]; g.y[c] = g.yc[j]; g.z[c] = zz[j][i]; g.facet[c] = 8; c++;
        }
        return g;
    }

    private static double minSpacing(double[] v) {
        double[] s = v.clone(); java.util.Arrays.sort(s); double best = Double.POSITIVE_INFINITY;
        for (int k = 1; k < s.length; k++) { double d = s[k] - s[k - 1]; if (d > 1e-6 && d < best) best = d; }
        return best;
    }

    /**
     * Target grid from an OMS/GeoTools feature collection (e.g. the output of OmsShapefileFeatureReader): point features
     * (cell centres) or polygon features (cells; the centroid is used). zField holds the cell elevation.
     */
    public static Grid fromFeatures(org.geotools.data.simple.SimpleFeatureCollection fc, String zField, String filledField) throws IOException {
        List<double[]> p = new ArrayList<>(); boolean useFlag = filledField != null && !filledField.isEmpty();
        try (org.geotools.data.simple.SimpleFeatureIterator it = fc.features()) {
            while (it.hasNext()) {
                org.opengis.feature.simple.SimpleFeature f = it.next();
                Object a = f.getAttribute(zField);
                if (a == null) throw new IOException("feature " + f.getID() + ": field " + zField + " missing or empty");
                org.locationtech.jts.geom.Point ct = ((org.locationtech.jts.geom.Geometry) f.getDefaultGeometry()).getCentroid();
                Object fl = useFlag ? f.getAttribute(filledField) : null;
                p.add(new double[]{ct.getX(), ct.getY(), ((Number) a).doubleValue(), fl == null ? 0 : ((Number) fl).doubleValue()});
            }
        }
        double[] x = new double[p.size()], y = new double[p.size()], z = new double[p.size()];
        for (int k = 0; k < x.length; k++) { x[k] = p.get(k)[0]; y[k] = p.get(k)[1]; z[k] = p.get(k)[2]; }
        Grid g = fromCells(x, y, z);
        for (double[] e : p) if (e[3] != 0) g.elevFilled[g.nearestCellRounded(e[0], e[1])] = true;
        return g;
    }

    /** Write the grid as ESRI ASCII DEM (cell corners, NODATA = -9999), readable by readAscii and by GIS software. */
    public void writeAscii(Path asc) throws IOException {
        double cs = nx > 1 ? xc[1] - xc[0] : (ny > 1 ? yc[0] - yc[1] : 1000);
        try (java.io.PrintWriter w = new java.io.PrintWriter(Files.newBufferedWriter(asc))) {
            w.printf(java.util.Locale.ROOT, "ncols %d%nnrows %d%nxllcorner %s%nyllcorner %s%ncellsize %s%nNODATA_value -9999%n", nx, ny,
                    fmt(xc[0] - cs / 2), fmt(yc[ny - 1] - cs / 2), fmt(cs));
            StringBuilder b = new StringBuilder();
            for (int j = 0; j < ny; j++) {
                b.setLength(0);
                for (int i = 0; i < nx; i++) { if (i > 0) b.append(' '); int k = cellIndex[j][i]; b.append(k < 0 ? "-9999" : fmt(z[k])); }
                w.println(b);
            }
        }
    }

    private static String fmt(double v) { return v == Math.rint(v) && Math.abs(v) < 1e15 ? Long.toString((long) v) : Double.toString(v); }
}
