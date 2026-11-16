/*
 * GNU GPL v3 License. Original code: Marialaura Bancheri, Giuseppe Formetta et al. (GEOframe Krigings).
 * Refactored 2026: thin package wrapper over the shared krigingsCore implementation (see UPSTREAM_DIFF.md).
 */
package krigingsPointCase;

/** The Class StationsSelection (implementation in krigingsCore.StationsSelection). */
public class StationsSelection extends krigingsCore.StationsSelection {

	@Override
	protected krigingsCore.Model createModel(double[] distanceVector, int inNumCloserStations, double maxdist) {
		return SimpleModelFactory.createModel(distanceVector, inNumCloserStations, maxdist);
	}
}
