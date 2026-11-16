package krigingsCore.local;


/** Interpolation targets of one call: coordinates [m], elevation [m], own station index (LOO) or -1, seed index. */
public final class TargetSet {
    public final double[] x, y, z; public final int[] self; public final long[] seedIndex; public final int n;

    public TargetSet(double[] x, double[] y, double[] z, int[] self, long[] seedIndex) {
        this.x = x; this.y = y; this.z = z; this.self = self; this.seedIndex = seedIndex; this.n = x.length;
    }
}
