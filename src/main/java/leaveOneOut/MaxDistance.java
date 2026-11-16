/*
 * GNU GPL v3 License. Original code: Marialaura Bancheri, Giuseppe Formetta et al. (GEOframe Krigings).
 * Refactored 2026: thin package wrapper over the shared krigingsCore implementation (see UPSTREAM_DIFF.md).
 */
package leaveOneOut;

/** The Class MaxDistance (implementation in krigingsCore.MaxDistance). */
public class MaxDistance extends krigingsCore.MaxDistance implements Model {

	public MaxDistance(double[] distanceVector, double maxdist) {
		super(distanceVector, maxdist);
	}
}
