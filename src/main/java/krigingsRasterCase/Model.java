/*
 * GNU GPL v3 License. Original code: Marialaura Bancheri, Giuseppe Formetta et al. (GEOframe Krigings).
 * Refactored 2026: thin package wrapper over the shared krigingsCore implementation (see UPSTREAM_DIFF.md).
 */
package krigingsRasterCase;

/** The Interface Model (kept for API compatibility; see krigingsCore.Model). */
public interface Model extends krigingsCore.Model {

	/**
	 * Number of stations.
	 *
	 * @return the int number of the stations
	 */
	abstract public int numberOfStations();
}
