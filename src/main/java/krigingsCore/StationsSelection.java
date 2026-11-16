package krigingsCore;

import static org.hortonmachine.gears.libs.modules.HMConstants.isNovalue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.hortonmachine.gears.libs.modules.ModelsEngine;
import org.hortonmachine.gears.libs.monitor.IHMProgressMonitor;
import org.hortonmachine.gears.libs.monitor.LogProgressMonitor;
import org.hortonmachine.gears.utils.sorting.QuickSortAlgorithm;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.opengis.feature.simple.SimpleFeature;

/**
 * Shared station selection of the three kriging packages (code of krigingsPointCase.StationsSelection, unchanged
 * in behaviour). Package-specific variants are subclasses: leaveOneOut excludes one station id
 * ({@link #include(int)}), krigingsRasterCase keeps its legacy neighbour arrays ({@link #rasterLegacy}).
 */
public class StationsSelection {

	/** The .shp of the measurement point, containing the position of the stations. */
	public SimpleFeatureCollection inStations = null;

	/** The HM with the measured data to be interpolated. */
	public HashMap<Integer, double[]> inData = null;

	/** Include zeros in computations (default is true). */
	public boolean doIncludezero = true;

	/** In the case of kriging with neighbor, maxdist is the maximum distance within the algorithm has to consider the stations */
	public double maxdist;

	/** In the case of kriging with neighbor, inNumCloserStations is the number of stations the algorithm has to consider */
	public int inNumCloserStations;

	/** The field of the vector of stations, defining the id.*/
	public String fStationsid = null;

	/** The pm. */
	public IHMProgressMonitor pm = new LogProgressMonitor();

	public double[] xStationInitialSet;
	public double[] yStationInitialSet;
	public double[] zStationInitialSet;
	public double[] hStationInitialSet;
	public int[] idStationInitialSet;

	/** The flag is true if all the stations are equals */
	public boolean areAllEquals = true;

	/** The number of different stations. */
	public int n1 = 0;

	public double idx;

	public double idy;

	public String fStationsZ = null;

	/** Interpolate log values (0.9.9 point case): h -> Utility.getLog(h) when the station is read. */
	public boolean doLogarithmic = false;

	/** The model selection for the choice of the stations. */
	protected Model modelSelection;

	/**
	 * krigingsRasterCase legacy behaviour: distances measured from xStationInitialSet[n1] (not idx/idy), no cap of
	 * inNumCloserStations at the number of stations, neighbour arrays of length dim (no empty target slot).
	 */
	protected boolean rasterLegacy = false;

	/** Hook for subclasses: false = the station is skipped (leave-one-out). */
	protected boolean include(int id) {
		return true;
	}

	public void execute() throws Exception {

		List<Double> xStationList = new ArrayList<Double>();
		List<Double> yStationList = new ArrayList<Double>();
		List<Double> zStationList = new ArrayList<Double>();
		List<Double> hStationList = new ArrayList<Double>();
		List<Integer> idStationList = new ArrayList<Integer>();

		// Store the station coordinates and measured data in the array. Skip data for non existing stations and novalues.
		FeatureIterator<SimpleFeature> stationsIter = inStations.features();
		try {
			while (stationsIter.hasNext()) {
				SimpleFeature feature = stationsIter.next();
				int id = ((Number) feature.getAttribute(fStationsid)).intValue();

				double z = 0;
				if (fStationsZ != null) {
					try {
						z = ((Number) feature.getAttribute(fStationsZ)).doubleValue();
					} catch (NullPointerException e) {
					}
				}

				Coordinate coordinate = ((Geometry) feature.getDefaultGeometry()).getCentroid().getCoordinate();
				double[] h = inData.get(id);
				if (h == null || isNovalue(h[0])) {
					continue;
				}
				if (doLogarithmic) {
					h[0] = krigingsPointCase.Utility.getLog(h[0]);
				}
				if (!include(id)) {
					continue;
				}

				if (doIncludezero) {
					if (Math.abs(h[0]) >= 0.0) { // TOLL
						xStationList.add(coordinate.x);
						yStationList.add(coordinate.y);
						zStationList.add(z);
						hStationList.add(h[0]);
						idStationList.add(id);
						n1 = n1 + 1;
					}
				} else {
					if (Math.abs(h[0]) > 0.0) { // TOLL
						xStationList.add(coordinate.x);
						yStationList.add(coordinate.y);
						zStationList.add(z);
						hStationList.add(h[0]);
						idStationList.add(id);
						n1 = n1 + 1;
					}
				}
			}
		} finally {
			stationsIter.close();
		}

		int nStaz = xStationList.size();

		/*
		 * xStationInitialSet has the dimensions of the coordinates of the measurements points plus 1 (the station where
		 * it is going to interpolate)
		 */
		xStationInitialSet = new double[nStaz + 1];
		yStationInitialSet = new double[nStaz + 1];
		zStationInitialSet = new double[nStaz + 1];
		hStationInitialSet = new double[nStaz + 1];
		idStationInitialSet = new int[nStaz + 1];

		if (nStaz != 0) {
			xStationInitialSet[0] = xStationList.get(0);
			yStationInitialSet[0] = yStationList.get(0);
			zStationInitialSet[0] = zStationList.get(0);
			hStationInitialSet[0] = hStationList.get(0);
			idStationInitialSet[0] = idStationList.get(0);
			double previousValue = hStationInitialSet[0];

			for (int i = 0; i < nStaz; i++) {
				double xTmp = xStationList.get(i);
				double yTmp = yStationList.get(i);
				double zTmp = zStationList.get(i);
				double hTmp = hStationList.get(i);
				int idTmp = idStationList.get(i);

				boolean doubleStation = ModelsEngine.verifyDoubleStation(xStationInitialSet, yStationInitialSet, zStationInitialSet,
						hStationInitialSet, xTmp, yTmp, zTmp, hTmp, i, false);
				if (!doubleStation) {
					xStationInitialSet[i] = xTmp;
					yStationInitialSet[i] = yTmp;
					zStationInitialSet[i] = zTmp;
					hStationInitialSet[i] = hTmp;
					idStationInitialSet[i] = idTmp;
					if (areAllEquals && hStationInitialSet[i] != previousValue) {
						areAllEquals = false;
					}
					previousValue = hStationInitialSet[i];
				}
			}
		}

		// kriging with neighbours: distances from the target, sorted
		if (inNumCloserStations > 0 || maxdist > 0) {

			double x2, y2;
			double dDifX, dDifY;
			double distanceVector[] = new double[xStationInitialSet.length];
			double pos[] = new double[xStationInitialSet.length];

			for (int jj = 0; jj < xStationInitialSet.length; jj++) {
				x2 = xStationInitialSet[jj];
				y2 = yStationInitialSet[jj];
				if (rasterLegacy) {
					dDifX = xStationInitialSet[n1] - x2;
					dDifY = yStationInitialSet[n1] - y2;
				} else {
					dDifX = idx - x2;
					dDifY = idy - y2;
				}
				distanceVector[jj] = Math.sqrt(dDifX * dDifX + dDifY * dDifY);
				pos[jj] = jj;
			}

			QuickSortAlgorithm t = new QuickSortAlgorithm(pm);
			t.sort(distanceVector, pos);

			if (!rasterLegacy) inNumCloserStations = (inNumCloserStations > nStaz) ? nStaz : inNumCloserStations;

			modelSelection = createModel(distanceVector, inNumCloserStations, maxdist);
			int dim = modelSelection.numberOfStations();
			int len = rasterLegacy ? dim : dim + 1;

			double[] xStationWithNeighbour = new double[len];
			double[] yStationWithNeighbour = new double[len];
			int[] idStationWithNeighbour = new int[len];
			double[] zStationWithNeighbour = new double[len];
			double[] hWithNeighbour = new double[len];

			for (int i = 0; i < dim; i++) {
				xStationWithNeighbour[i] = xStationInitialSet[(int) pos[i]];
				yStationWithNeighbour[i] = yStationInitialSet[(int) pos[i]];
				zStationWithNeighbour[i] = zStationInitialSet[(int) pos[i]];
				idStationWithNeighbour[i] = idStationInitialSet[(int) pos[i]];
				hWithNeighbour[i] = hStationInitialSet[(int) pos[i]];
			}

			xStationInitialSet = xStationWithNeighbour;
			yStationInitialSet = yStationWithNeighbour;
			zStationInitialSet = zStationWithNeighbour;
			idStationInitialSet = idStationWithNeighbour;
			hStationInitialSet = hWithNeighbour;
		}
	}

	/** Station-selection rule; subclasses return their package's Model implementation. */
	protected Model createModel(double[] distanceVector, int inNumCloserStations, double maxdist) {
		return SimpleModelFactory.createModel(distanceVector, inNumCloserStations, maxdist);
	}
}
