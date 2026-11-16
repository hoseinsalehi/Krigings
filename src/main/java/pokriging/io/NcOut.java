package pokriging.io;

import ucar.ma2.Array;
import ucar.ma2.ArrayByte;
import ucar.ma2.ArrayDouble;
import ucar.ma2.ArrayFloat;
import ucar.ma2.ArrayInt;
import ucar.ma2.ArrayShort;
import ucar.ma2.DataType;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Attribute;
import ucar.nc2.write.NetcdfFileFormat;
import ucar.nc2.write.NetcdfFormatWriter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

/**
 * CF grid file (time, y, x) written as NetCDF-3 (64-bit offset) by netcdf-java, one time slice at a time; on close it
 * is optionally converted to NetCDF-4 with zlib level 6 and (1, ny, nx) chunks by the netCDF C utility nccopy
 * (netcdf-java cannot write NetCDF-4 without the C library). Packed integer variables follow the Python encoding:
 * value = stored * scale_factor, unsigned types use the _Unsigned = "true" convention.
 */
public final class NcOut implements AutoCloseable {
    public static final class Var {
        final String name, type, units, longName; final double scale; final long fill;
        public Var(String name, String type, double scale, long fill, String units, String longName) {
            this.name = name; this.type = type; this.scale = scale; this.fill = fill; this.units = units; this.longName = longName;
        }
        DataType dt() { switch (type) { case "i2": case "u2": return DataType.SHORT; case "u1": case "i1": return DataType.BYTE;
            case "i4": return DataType.INT; case "f4": return DataType.FLOAT; default: throw new IllegalArgumentException(type); } }
        boolean unsigned() { return type.startsWith("u"); }
        /** numpy-style encoding: rint(x/scale) clipped to [min+1 (signed), max-1 (unsigned)], NaN -> fill. */
        long encode(double x) {
            double q = Math.rint(x / scale);
            if (Double.isNaN(q) || Double.isInfinite(q)) return fill;
            double lo, hi;
            switch (type) {
                case "i2": lo = -32767; hi = 32767; break; case "u2": lo = 0; hi = 65534; break;
                case "u1": lo = 0; hi = 254; break; case "i1": lo = -127; hi = 127; break;
                default: lo = -2147483647; hi = 2147483647;
            }
            return (long) Math.max(lo, Math.min(hi, q));
        }
    }

    private final Path finalPath, partPath; private final NetcdfFormatWriter w; private final int ny, nx;
    private final Map<String, Var> vars = new LinkedHashMap<>(); private final String nccopy;

    /**
     * @param days days since 1991-01-01 of the time steps; nccopy path or null (keep NetCDF-3)
     */
    public NcOut(Path file, double[] xc, double[] yc, int[] days, List<Var> vl, Map<String, String> attrs, String crs, String nccopy) throws IOException {
        this.finalPath = file; this.partPath = file.resolveSibling(file.getFileName() + ".part"); this.ny = yc.length; this.nx = xc.length;
        this.nccopy = nccopy;
        Files.createDirectories(file.toAbsolutePath().getParent());
        NetcdfFormatWriter.Builder b = NetcdfFormatWriter.builder().setNewFile(true).setFormat(NetcdfFileFormat.NETCDF3_64BIT_OFFSET)
                .setLocation(partPath.toString()).setFill(false);
        b.addDimension("time", days.length); b.addDimension("y", ny); b.addDimension("x", nx);
        b.addVariable("time", DataType.INT, "time").addAttribute(new Attribute("units", "days since 1991-01-01"))
                .addAttribute(new Attribute("calendar", "standard")).addAttribute(new Attribute("standard_name", "time"));
        b.addVariable("y", DataType.DOUBLE, "y").addAttribute(new Attribute("units", "m")).addAttribute(new Attribute("standard_name", "projection_y_coordinate"));
        b.addVariable("x", DataType.DOUBLE, "x").addAttribute(new Attribute("units", "m")).addAttribute(new Attribute("standard_name", "projection_x_coordinate"));
        b.addVariable("crs", DataType.INT, "").addAttribute(new Attribute("spatial_ref", crs)).addAttribute(new Attribute("epsg_code", crs));
        for (Var v : vl) {
            vars.put(v.name, v);
            ucar.nc2.Variable.Builder<?> vb = b.addVariable(v.name, v.dt(), "time y x");
            if (v.unsigned()) vb.addAttribute(new Attribute("_Unsigned", "true"));
            vb.addAttribute(fillAttr(v));
            if (v.scale != 1) vb.addAttribute(new Attribute("scale_factor", v.scale));
            vb.addAttribute(new Attribute("units", v.units)).addAttribute(new Attribute("long_name", v.longName))
              .addAttribute(new Attribute("grid_mapping", "crs"));
        }
        for (Map.Entry<String, String> e : attrs.entrySet()) b.addAttribute(new Attribute(e.getKey(), e.getValue()));
        b.addAttribute(new Attribute("Conventions", "CF-1.8"));
        w = b.build();
        try {
            ArrayInt.D1 t = new ArrayInt.D1(days.length, false); for (int i = 0; i < days.length; i++) t.set(i, days[i]); w.write("time", t);
            ArrayDouble.D1 ya = new ArrayDouble.D1(ny); for (int i = 0; i < ny; i++) ya.set(i, yc[i]); w.write("y", ya);
            ArrayDouble.D1 xa = new ArrayDouble.D1(nx); for (int i = 0; i < nx; i++) xa.set(i, xc[i]); w.write("x", xa);
        } catch (InvalidRangeException e) { throw new IOException(e); }
    }

    private static Attribute fillAttr(Var v) {
        switch (v.dt()) {
            case SHORT: return new Attribute("_FillValue", (short) v.fill);
            case BYTE: return new Attribute("_FillValue", (byte) v.fill);
            case INT: return new Attribute("_FillValue", (int) v.fill);
            default: return new Attribute("_FillValue", Float.NaN);
        }
    }

    /** Write one time slice of variable name from cell values (cells at iy/ix; other pixels = fill). */
    public synchronized void writeSlice(String name, int t, double[] cellValues, int[] iy, int[] ix) throws IOException {
        Var v = vars.get(name);
        Array a;
        switch (v.dt()) {
            case SHORT: { short[] d = new short[ny * nx]; java.util.Arrays.fill(d, (short) v.fill);
                for (int k = 0; k < iy.length; k++) d[iy[k] * nx + ix[k]] = (short) v.encode(cellValues[k]);
                a = Array.factory(DataType.SHORT, new int[]{1, ny, nx}, d); break; }
            case BYTE: { byte[] d = new byte[ny * nx]; java.util.Arrays.fill(d, (byte) v.fill);
                for (int k = 0; k < iy.length; k++) d[iy[k] * nx + ix[k]] = (byte) v.encode(cellValues[k]);
                a = Array.factory(DataType.BYTE, new int[]{1, ny, nx}, d); break; }
            case INT: { int[] d = new int[ny * nx]; java.util.Arrays.fill(d, (int) v.fill);
                for (int k = 0; k < iy.length; k++) d[iy[k] * nx + ix[k]] = (int) v.encode(cellValues[k]);
                a = Array.factory(DataType.INT, new int[]{1, ny, nx}, d); break; }
            default: { float[] d = new float[ny * nx]; java.util.Arrays.fill(d, Float.NaN);
                for (int k = 0; k < iy.length; k++) d[iy[k] * nx + ix[k]] = (float) cellValues[k];
                a = Array.factory(DataType.FLOAT, new int[]{1, ny, nx}, d); }
        }
        try { w.write(name, new int[]{t, 0, 0}, a); } catch (InvalidRangeException e) { throw new IOException(e); }
    }

    @Override public void close() throws IOException {
        w.close();
        if (nccopy != null && new File(nccopy).canExecute()) {
            Path tmp = finalPath.resolveSibling(finalPath.getFileName() + ".nc4.part");
            ProcessBuilder pb = new ProcessBuilder(nccopy, "-k", "nc4", "-d", "6", "-s", "-c", "time/1,y/" + ny + ",x/" + nx,
                    partPath.toString(), tmp.toString()).redirectErrorStream(true);
            try {
                Process p = pb.start(); String out = new String(p.getInputStream().readAllBytes());
                if (!p.waitFor(600, TimeUnit.SECONDS) || p.exitValue() != 0) throw new IOException("nccopy failed: " + out);
                Files.move(tmp, finalPath, java.nio.file.StandardCopyOption.REPLACE_EXISTING);
                Files.delete(partPath);
                return;
            } catch (InterruptedException e) { throw new IOException(e); }
        }
        Files.move(partPath, finalPath, java.nio.file.StandardCopyOption.REPLACE_EXISTING);
    }

    /** Monthly normals file (month, y, x) as float32, used by PoNormals. */
    public static void writeMonthly(Path file, Grid g, Map<String, double[][]> fields, Map<String, String[]> unitsLong,
                                    Map<String, String> attrs, String crs, String nccopy) throws IOException {
        Path part = file.resolveSibling(file.getFileName() + ".part");
        Files.createDirectories(file.toAbsolutePath().getParent());
        NetcdfFormatWriter.Builder b = NetcdfFormatWriter.builder().setNewFile(true).setFormat(NetcdfFileFormat.NETCDF3_64BIT_OFFSET).setLocation(part.toString());
        b.addDimension("month", 12); b.addDimension("y", g.ny); b.addDimension("x", g.nx);
        b.addVariable("month", DataType.INT, "month");
        b.addVariable("y", DataType.DOUBLE, "y").addAttribute(new Attribute("units", "m"));
        b.addVariable("x", DataType.DOUBLE, "x").addAttribute(new Attribute("units", "m"));
        for (String k : fields.keySet()) b.addVariable(k, DataType.FLOAT, "month y x").addAttribute(new Attribute("_FillValue", Float.NaN))
                .addAttribute(new Attribute("units", unitsLong.get(k)[0])).addAttribute(new Attribute("long_name", unitsLong.get(k)[1]));
        for (Map.Entry<String, String> e : attrs.entrySet()) b.addAttribute(new Attribute(e.getKey(), e.getValue()));
        b.addAttribute(new Attribute("crs", crs));
        try (NetcdfFormatWriter w = b.build()) {
            ArrayInt.D1 m = new ArrayInt.D1(12, false); for (int i = 0; i < 12; i++) m.set(i, i + 1); w.write("month", m);
            ArrayDouble.D1 ya = new ArrayDouble.D1(g.ny); for (int i = 0; i < g.ny; i++) ya.set(i, g.yc[i]); w.write("y", ya);
            ArrayDouble.D1 xa = new ArrayDouble.D1(g.nx); for (int i = 0; i < g.nx; i++) xa.set(i, g.xc[i]); w.write("x", xa);
            for (Map.Entry<String, double[][]> e : fields.entrySet()) {
                float[] d = new float[12 * g.ny * g.nx]; java.util.Arrays.fill(d, Float.NaN);
                for (int mm = 0; mm < 12; mm++) for (int k = 0; k < g.n; k++) d[(mm * g.ny + g.iy[k]) * g.nx + g.ix[k]] = (float) e.getValue()[mm][k];
                w.write(e.getKey(), Array.factory(DataType.FLOAT, new int[]{12, g.ny, g.nx}, d));
            }
        } catch (InvalidRangeException e) { throw new IOException(e); }
        if (nccopy != null && new File(nccopy).canExecute()) {
            Path tmp = file.resolveSibling(file.getFileName() + ".nc4.part");
            try {
                Process p = new ProcessBuilder(nccopy, "-k", "nc4", "-d", "4", part.toString(), tmp.toString()).redirectErrorStream(true).start();
                String out = new String(p.getInputStream().readAllBytes());
                if (p.waitFor() != 0) throw new IOException("nccopy failed: " + out);
                Files.move(tmp, file, java.nio.file.StandardCopyOption.REPLACE_EXISTING); Files.delete(part); return;
            } catch (InterruptedException e) { throw new IOException(e); }
        }
        Files.move(part, file, java.nio.file.StandardCopyOption.REPLACE_EXISTING);
    }
}
