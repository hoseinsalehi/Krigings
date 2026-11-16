package pokrigingTests;

import utilities.random.Blake2b;
import utilities.random.NumpyRandom;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

/** Compares NumpyRandom and Blake2b.seedFor with numpy draws (tests/data, from tools/make_rng_testdata.py). */
public class RngTest {
    public static void main(String[] a) throws Exception {
        String dir = a.length > 0 ? a[0] : "tests/data";
        List<String> lines = Files.readAllLines(Paths.get(dir, "rng_numpy.csv"));
        long cur = -1; String kind = ""; NumpyRandom g = null; int n = 0, bad = 0;
        for (String l : lines.subList(1, lines.size())) {
            String[] f = l.split(",");
            long seed = Long.parseLong(f[0]); int i = Integer.parseInt(f[2]);
            if (seed != cur || !f[1].equals(kind)) { cur = seed; kind = f[1]; g = new NumpyRandom(seed); }
            double want = Double.parseDouble(f[3].replace("0x", "0x").contains("p") ? hexToJava(f[3]) : f[3]);
            double got;
            if (kind.equals("random")) got = g.random();
            else if (kind.equals("normal")) got = g.standardNormal();
            else got = (i < 30 || i >= 60) ? g.random() : g.normal(0, .06);
            n++;
            if (Double.doubleToLongBits(got) != Double.doubleToLongBits(want)) {
                if (bad < 5) System.out.println("MISMATCH seed=" + seed + " " + kind + " i=" + i + " got=" + got + " want=" + want);
                bad++;
            }
        }
        System.out.println("RNG draws compared: " + n + ", bitwise mismatches: " + bad);
        List<String> sl = Files.readAllLines(Paths.get(dir, "seed_for.csv"));
        int sb = 0;
        for (String l : sl.subList(1, sl.size())) {
            int k = l.lastIndexOf(',');
            long want = Long.parseLong(l.substring(k + 1)), got = Blake2b.seedFor(l.substring(0, k));
            if (got != want) { sb++; System.out.println("SEED MISMATCH " + l + " got " + got); }
        }
        System.out.println("seed_for compared: " + (sl.size() - 1) + ", mismatches: " + sb);
        if (bad > 0 || sb > 0) System.exit(1);
    }

    /** Python float.hex() -> Java hex float literal ("0x1.8p+3" is valid for Double.parseDouble). */
    private static String hexToJava(String s) {
        if (s.equals("0x0.0p+0")) return "0.0";
        return s;
    }
}
