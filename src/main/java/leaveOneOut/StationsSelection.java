/*
 * GNU GPL v3 License. Original code: Marialaura Bancheri, Giuseppe Formetta et al. (GEOframe Krigings).
 * Refactored 2026: thin package wrapper over the shared krigingsCore implementation (see UPSTREAM_DIFF.md).
 */
package leaveOneOut;

/** Station selection without the station idOut (leave-one-out); implementation in krigingsCore.StationsSelection. */
public class StationsSelection extends krigingsCore.StationsSelection {

	/** Id of the station left out. */
	public int idOut;

	boolean stationOut;

	@Override
	protected boolean include(int id) {
		return id != idOut;
	}

	@Override
	protected krigingsCore.Model createModel(double[] distanceVector, int inNumCloserStations, double maxdist) {
		return SimpleModelFactory.createModel(distanceVector, inNumCloserStations, maxdist);
	}
}
