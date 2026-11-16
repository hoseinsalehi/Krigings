package krigingsCore.local;


import java.util.Arrays;

/**
 * Static 2-D k-d tree for k-nearest-neighbour queries (Euclidean). Results are sorted by (distance, point index),
 * a deterministic tie rule (scipy cKDTree has no ties on the engine station sets). Distances are
 * sqrt(dx*dx + dy*dy) as scipy computes them. Thread-safe for queries after construction.
 */
public final class KdTree {
    private final double[] x, y;
    private final int[] idx;          // permutation of 0..n-1 in tree order
    private final int n;
    private static final int LEAF = 8;

    public KdTree(double[] x, double[] y) {
        this.x = x; this.y = y; n = x.length;
        idx = new int[n];
        for (int i = 0; i < n; i++) idx[i] = i;
        build(0, n, 0);
    }

    public int size() { return n; }

    private void build(int lo, int hi, int depth) {
        if (hi - lo <= LEAF) return;
        int mid = (lo + hi) >>> 1;
        select(lo, hi - 1, mid, depth & 1);
        build(lo, mid, depth + 1);
        build(mid + 1, hi, depth + 1);
    }

    private double coord(int p, int axis) { return axis == 0 ? x[p] : y[p]; }

    private void select(int lo, int hi, int k, int axis) {
        while (hi > lo) {
            double pivot = coord(idx[(lo + hi) >>> 1], axis);
            int i = lo, j = hi;
            while (i <= j) {
                while (coord(idx[i], axis) < pivot) i++;
                while (coord(idx[j], axis) > pivot) j--;
                if (i <= j) { int t = idx[i]; idx[i] = idx[j]; idx[j] = t; i++; j--; }
            }
            if (k <= j) hi = j; else if (k >= i) lo = i; else return;
        }
    }

    /** Per-thread query buffers. */
    public static final class Query {
        double[] hd; int[] hi; int cnt, k;
        public final int[] outIdx; public final double[] outDist;
        public Query(int kmax) { hd = new double[kmax]; hi = new int[kmax]; outIdx = new int[kmax]; outDist = new double[kmax]; }
    }

    /**
     * k nearest points of (qx, qy); writes point indices and distances (ascending) into q.outIdx / q.outDist and
     * returns the number found (min(k, n)).
     */
    public int knn(double qx, double qy, int k, Query q) {
        k = Math.min(k, n);
        if (q.hd.length < k) throw new IllegalArgumentException("query buffer too small");
        q.cnt = 0; q.k = k;
        search(0, n, 0, qx, qy, q);
        // heap -> sorted by (squared distance, index)
        int c = q.cnt;
        double[] od = q.outDist; int[] oi = q.outIdx;
        for (int i = 0; i < c; i++) {                 // insertion sort (k <= ~100)
            double d = q.hd[i]; int p = q.hi[i]; int j = i - 1;
            while (j >= 0 && worse(od[j], oi[j], d, p)) { od[j + 1] = od[j]; oi[j + 1] = oi[j]; j--; }
            od[j + 1] = d; oi[j + 1] = p;
        }
        for (int i = 0; i < c; i++) od[i] = Math.sqrt(od[i]);
        return c;
    }

    // max-heap on (d2, index) of the current k best
    private static boolean worse(double da, int ia, double db, int ib) { return da > db || (da == db && ia > ib); }

    private void offer(Query q, double d2, int p) {
        if (q.cnt < q.k) {
            int i = q.cnt++;
            q.hd[i] = d2; q.hi[i] = p;
            while (i > 0) { int par = (i - 1) >>> 1; if (worse(q.hd[i], q.hi[i], q.hd[par], q.hi[par])) { swap(q, i, par); i = par; } else break; }
        } else if (worse(q.hd[0], q.hi[0], d2, p)) {
            q.hd[0] = d2; q.hi[0] = p;
            int i = 0;
            for (;;) {
                int l = 2 * i + 1, r = l + 1, m = i;
                if (l < q.cnt && worse(q.hd[l], q.hi[l], q.hd[m], q.hi[m])) m = l;
                if (r < q.cnt && worse(q.hd[r], q.hi[r], q.hd[m], q.hi[m])) m = r;
                if (m == i) break;
                swap(q, i, m); i = m;
            }
        }
    }

    private static void swap(Query q, int a, int b) {
        double t = q.hd[a]; q.hd[a] = q.hd[b]; q.hd[b] = t; int u = q.hi[a]; q.hi[a] = q.hi[b]; q.hi[b] = u;
    }

    private void search(int lo, int hi, int depth, double qx, double qy, Query q) {
        if (hi - lo <= LEAF) {
            for (int i = lo; i < hi; i++) { int p = idx[i]; double dx = qx - x[p], dy = qy - y[p]; offer(q, dx * dx + dy * dy, p); }
            return;
        }
        int mid = (lo + hi) >>> 1, axis = depth & 1, p = idx[mid];
        double diff = (axis == 0 ? qx - x[p] : qy - y[p]);
        { double dx = qx - x[p], dy = qy - y[p]; offer(q, dx * dx + dy * dy, p); }
        int nlo, nhi, flo, fhi;
        if (diff < 0) { nlo = lo; nhi = mid; flo = mid + 1; fhi = hi; } else { nlo = mid + 1; nhi = hi; flo = lo; fhi = mid; }
        search(nlo, nhi, depth + 1, qx, qy, q);
        if (q.cnt < q.k || diff * diff <= q.hd[0]) search(flo, fhi, depth + 1, qx, qy, q);
    }
}
