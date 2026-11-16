package pokriging.io;

import ucar.ma2.Array;
import ucar.nc2.NetcdfFile;
import ucar.nc2.NetcdfFiles;
import ucar.nc2.Variable;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/** Readers for the normals written by PoNormals or by the Python engine (same file names and columns). */
public final class NormalsIO {
    private NormalsIO() { }

    /** Station table (station, month, colA, colB) -> {A[S][12], B[S][12]} in engine station order (NaN where absent). */
    public static double[][][] stationPair(Path csv, String[] ids, String colA, String colB) throws IOException {
        Csv c = Csv.read(csv);
        Map<String, Integer> pos = new HashMap<>();
        for (int i = 0; i < ids.length; i++) pos.put(ids[i], i);
        double[][] a = nan(ids.length), b = nan(ids.length);
        int cs = c.col("station"), cm = c.col("month"), ca = c.col(colA), cb = c.col(colB);
        for (String[] r : c.rows) {
            Integer i = pos.get(r[cs].trim());
            if (i == null) continue;
            int m = Integer.parseInt(r[cm].trim()) - 1;
            a[i][m] = Csv.num(r[ca]); b[i][m] = Csv.num(r[cb]);
        }
        return new double[][][]{a, b};
    }

    private static double[][] nan(int n) { double[][] x = new double[n][12]; for (double[] r : x) Arrays.fill(r, Double.NaN); return x; }

    /**
     * Gridded monthly field (month, y, x) of a NetCDF file at the valid cells: [12][ncell] as float values widened
     * to double (the Python engine stores float32 normals).
     */
    public static double[][] gridField(Path nc, String var, Grid g) throws IOException {
        try (NetcdfFile f = NetcdfFiles.open(nc.toString())) {
            Variable v = f.findVariable(var);
            if (v == null) throw new IOException(nc + ": variable " + var + " missing");
            int[] sh = v.getShape();
            if (sh.length != 3 || sh[0] != 12 || sh[1] != g.ny || sh[2] != g.nx)
                throw new IOException(nc + ": " + var + " shape " + Arrays.toString(sh) + " does not match grid 12x" + g.ny + "x" + g.nx);
            Array a = v.read();
            double[][] out = new double[12][g.n];
            for (int m = 0; m < 12; m++) for (int k = 0; k < g.n; k++)
                out[m][k] = a.getFloat((m * g.ny + g.iy[k]) * g.nx + g.ix[k]);
            return out;
        }
    }
}
