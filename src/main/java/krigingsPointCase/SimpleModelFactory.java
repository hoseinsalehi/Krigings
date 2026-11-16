/*
 * GNU GPL v3 License. Original code: Marialaura Bancheri, Giuseppe Formetta et al. (GEOframe Krigings).
 * Refactored 2026: thin package wrapper over the shared krigingsCore implementation (see UPSTREAM_DIFF.md).
 */
package krigingsPointCase;

/** A simple design factory for creating Model objects (same rule as krigingsCore.SimpleModelFactory). */
public class SimpleModelFactory {

	public static Model createModel(double[] distanceVector, int inNumCloserStations, double maxdist) {
		Model model = null;
		if (maxdist > 0) {
			model = new MaxDistance(distanceVector, maxdist);
		} else if (inNumCloserStations > 0) {
			model = new NumberOfStations(inNumCloserStations);
		}
		return model;
	}
}
