package krigingsCore;

/** Factory of the station-selection rule: maxdist > 0 -> MaxDistance, else inNumCloserStations > 0 -> NumberOfStations. */
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
