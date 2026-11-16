package utilities.random;

/** Random stream used by the PSO: numpy-exact ({@link NumpyRandom}) or Java SplittableRandom ("java" mode). */
public interface Rng {
    double random();
    double normal(double loc, double scale);

    /** Fill out[] with random() draws (C order). */
    default void fill(double[] out) { for (int i = 0; i < out.length; i++) out[i] = random(); }

    static Rng create(String mode, long seed) {
        if ("java".equals(mode)) {
            java.util.SplittableRandom r = new java.util.SplittableRandom(seed);
            final java.util.Random g = new java.util.Random(seed ^ 0x5DEECE66DL);
            return new Rng() {
                public double random() { return r.nextDouble(); }
                public double normal(double loc, double scale) { return loc + scale * g.nextGaussian(); }
            };
        }
        final NumpyRandom n = new NumpyRandom(seed);
        return new Rng() {
            public double random() { return n.random(); }
            public double normal(double loc, double scale) { return n.normal(loc, scale); }
            public void fill(double[] out) { n.random(out); }
        };
    }
}
