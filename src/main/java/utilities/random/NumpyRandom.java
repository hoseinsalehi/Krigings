package utilities.random;

/**
 * Bit-exact port of {@code numpy.random.default_rng(seed)} for the two draws used by the engine:
 * {@code random()} (53-bit doubles) and {@code normal(loc, scale)} (ziggurat {@code random_standard_normal}).
 * <ul>
 * <li>SeedSequence(entropy = non-negative integer, pool_size 4) -> generate_state(4, uint64)
 *     (numpy/random/bit_generator.pyx);</li>
 * <li>PCG64 = pcg_setseq_128_xsl_rr_64 seeded with pcg64_set_seed(state[0..1], state[2..3])
 *     (numpy/random/src/pcg64);</li>
 * <li>standard normal: numpy's 256-layer ziggurat with the numpy tables ({@link ZigguratTables}).</li>
 * </ul>
 * numpy is distributed under the BSD-3-Clause licence (Copyright (c) 2005-2025 NumPy Developers); the
 * algorithms and constants are reproduced here for numerical equivalence with the Python engine.
 * Not thread-safe: one instance per stream.
 */
public final class NumpyRandom {
    private static final int INIT_A = 0x43b0d7e5, MULT_A = 0x931e8875, INIT_B = 0x8b51f9dd, MULT_B = 0x58f38ded;
    private static final int MIX_MULT_L = 0xca01f9dd, MIX_MULT_R = 0x4973f715;
    private static final long MUL_HI = 2549297995355413924L, MUL_LO = 4865540595714422341L;

    private long stHi, stLo, incHi, incLo;

    /** Generator equivalent to numpy.random.default_rng(seed), seed >= 0. */
    public NumpyRandom(long seed) {
        if (seed < 0) throw new IllegalArgumentException("seed must be >= 0");
        int[] entropy;
        if (seed == 0) entropy = new int[]{0};
        else if ((seed >>> 32) == 0) entropy = new int[]{(int) seed};
        else entropy = new int[]{(int) seed, (int) (seed >>> 32)};
        long[] st = seedSequenceState64(entropy, 4);
        setSeed(st[0], st[1], st[2], st[3]);
    }

    // ---------------------------------------------------------------- SeedSequence
    private static int hashmix(int value, int[] hc) {
        value ^= hc[0];
        hc[0] *= MULT_A;
        value *= hc[0];
        value ^= value >>> 16;
        return value;
    }

    private static int mix(int x, int y) {
        int r = MIX_MULT_L * x - MIX_MULT_R * y;
        r ^= r >>> 16;
        return r;
    }

    /** SeedSequence(entropy).generate_state(nWords64, uint64), pool size 4, no spawn key. */
    static long[] seedSequenceState64(int[] entropy, int nWords64) {
        int pool = 4;
        int[] mixer = new int[pool];
        int[] hc = {INIT_A};
        for (int i = 0; i < pool; i++) mixer[i] = hashmix(i < entropy.length ? entropy[i] : 0, hc);
        for (int s = 0; s < pool; s++)
            for (int d = 0; d < pool; d++)
                if (s != d) mixer[d] = mix(mixer[d], hashmix(mixer[s], hc));
        for (int s = pool; s < entropy.length; s++)
            for (int d = 0; d < pool; d++) mixer[d] = mix(mixer[d], hashmix(entropy[s], hc));
        int n32 = 2 * nWords64;
        int[] state = new int[n32];
        int h = INIT_B;
        for (int i = 0; i < n32; i++) {
            int v = mixer[i % pool];
            v ^= h;
            h *= MULT_B;
            v *= h;
            v ^= v >>> 16;
            state[i] = v;
        }
        long[] out = new long[nWords64];
        for (int i = 0; i < nWords64; i++)
            out[i] = (state[2 * i] & 0xffffffffL) | ((long) state[2 * i + 1] << 32);   // little-endian view
        return out;
    }

    // ---------------------------------------------------------------- PCG64 (128-bit LCG, XSL-RR output)
    private void setSeed(long seedHi, long seedLo, long seqHi, long seqLo) {
        incHi = (seqHi << 1) | (seqLo >>> 63);
        incLo = (seqLo << 1) | 1L;
        stHi = 0; stLo = 0;
        step();
        long lo = stLo + seedLo;
        long carry = Long.compareUnsigned(lo, stLo) < 0 ? 1 : 0;
        stHi = stHi + seedHi + carry; stLo = lo;
        step();
    }

    private void step() {
        // state = state * MUL + inc (mod 2^128)
        long lo = stLo * MUL_LO;
        long hi = Math.multiplyHigh(stLo, MUL_LO) + ((stLo >> 63) & MUL_LO) + ((MUL_LO >> 63) & stLo);  // unsigned high
        hi += stHi * MUL_LO + stLo * MUL_HI;
        long nlo = lo + incLo;
        long carry = Long.compareUnsigned(nlo, lo) < 0 ? 1 : 0;
        stHi = hi + incHi + carry; stLo = nlo;
    }

    /** next_uint64 of numpy's PCG64. */
    public long nextLong() {
        step();
        long x = stHi ^ stLo;
        int rot = (int) (stHi >>> 58);
        return Long.rotateRight(x, rot);
    }

    /** Generator.random(): (next_uint64 >> 11) * 2^-53. */
    public double random() {
        return (nextLong() >>> 11) * (1.0 / 9007199254740992.0);
    }

    /** Generator.standard_normal() (ziggurat, numpy distributions.c random_standard_normal). */
    public double standardNormal() {
        for (;;) {
            long r = nextLong();
            int idx = (int) (r & 0xff);
            r >>>= 8;
            int sign = (int) (r & 1);
            long rabs = (r >>> 1) & 0x000fffffffffffffL;
            double x = rabs * ZigguratTables.WI[idx];
            if (sign == 1) x = -x;
            if (rabs < ZigguratTables.KI[idx]) return x;
            if (idx == 0) {
                for (;;) {
                    double xx = -ZigguratTables.ZIGGURAT_NOR_INV_R * StrictMath.log1p(-random());
                    double yy = -StrictMath.log1p(-random());
                    if (yy + yy > xx * xx)
                        return ((rabs >>> 8) & 1) == 1 ? -(ZigguratTables.ZIGGURAT_NOR_R + xx) : ZigguratTables.ZIGGURAT_NOR_R + xx;
                }
            } else {
                if (((ZigguratTables.FI[idx - 1] - ZigguratTables.FI[idx]) * random() + ZigguratTables.FI[idx])
                        < StrictMath.exp(-0.5 * x * x))
                    return x;
            }
        }
    }

    /** Generator.normal(loc, scale) = loc + scale * standard_normal. */
    public double normal(double loc, double scale) {
        return loc + scale * standardNormal();
    }

    /** Fill {@code out} with random() draws (C order). */
    public void random(double[] out) {
        for (int i = 0; i < out.length; i++) out[i] = random();
    }
}
