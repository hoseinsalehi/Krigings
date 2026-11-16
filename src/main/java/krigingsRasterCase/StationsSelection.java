/*
 * GNU GPL v3 License. Original code: Marialaura Bancheri, Giuseppe Formetta et al. (GEOframe Krigings).
 * Refactored 2026: thin package wrapper over the shared krigingsCore implementation (see UPSTREAM_DIFF.md).
 */
package krigingsRasterCase;

/**
 * Raster-case station selection (upstream variant, not used by krigingsRasterCase.Krigings, which uses
 * krigingsPointCase.StationsSelection). Legacy behaviour kept through krigingsCore.StationsSelection.rasterLegacy.
 */
public class StationsSelection extends krigingsCore.StationsSelection {

	public StationsSelection() {
		rasterLegacy = true;
	}

	@Override
	protected krigingsCore.Model createModel(double[] distanceVector, int inNumCloserStations, double maxdist) {
		return SimpleModelFactory.createModel(distanceVector, inNumCloserStations, maxdist);
	}
}
