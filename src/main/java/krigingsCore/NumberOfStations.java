package krigingsCore;

/** Fixed number of closest stations (returns inNumCloserStations + 1, as in the original packages). */
public class NumberOfStations implements Model {

	/** The in num closer stations. */
	protected int inNumCloserStations;

	public NumberOfStations(int inNumCloserStations) {
		this.inNumCloserStations = inNumCloserStations;
	}

	@Override
	public int numberOfStations() {
		return inNumCloserStations + 1;
	}
}
