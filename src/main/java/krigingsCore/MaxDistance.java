package krigingsCore;

/** Stations within maxdist of the target (at least the nearest 3 when only one is found), as in the original packages. */
public class MaxDistance implements Model {

	/** The maxdist. */
	protected double maxdist;

	/** The distance vector (sorted ascending). */
	protected double[] distanceVector;

	public MaxDistance(double[] distanceVector, double maxdist) {
		this.distanceVector = distanceVector;
		this.maxdist = maxdist;
	}

	@Override
	public int numberOfStations() {
		// posDist is the number of the stations within the distance
		int posDist = distanceVector.length;
		for (int k = 0; k < distanceVector.length; k++) {
			if (distanceVector[k] > maxdist) {
				posDist = k;
				break;
			}
		}
		// in case there are no stations within the distance, the algorithm considers at least the nearest 3
		posDist = (posDist == 1) ? posDist += 3 : posDist;
		return posDist;
	}
}
