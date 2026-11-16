/*
 * GNU GPL v3 License. Original code: Marialaura Bancheri, Giuseppe Formetta et al. (GEOframe Krigings).
 * Refactored 2026: thin package wrapper over the shared krigingsCore implementation (see UPSTREAM_DIFF.md).
 */
package krigingsRasterCase;

/** The Class MaxDistance (implementation in krigingsCore.MaxDistance). */
public class MaxDistance extends krigingsCore.MaxDistance implements Model {

	public MaxDistance(double[] distanceVector, double maxdist) {
		// upstream krigingsRasterCase.MaxDistance never assigned maxdist (stays 0); kept for compatibility
		super(distanceVector, 0.0);
	}
}
