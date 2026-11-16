/*
 * GNU GPL v3 License. Original code: Marialaura Bancheri, Giuseppe Formetta et al. (GEOframe Krigings).
 * Refactored 2026: thin package wrapper over the shared krigingsCore implementation (see UPSTREAM_DIFF.md).
 */
package krigingsPointCase;

/** The Class NumberOfStations (implementation in krigingsCore.NumberOfStations). */
public class NumberOfStations extends krigingsCore.NumberOfStations implements Model {

	public NumberOfStations(int inNumCloserStations) {
		super(inNumCloserStations);
	}
}
