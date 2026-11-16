package krigingsCore;

/**
 * Rule that decides how many of the distance-sorted stations enter a kriging system (shared by krigingsPointCase,
 * krigingsRasterCase and leaveOneOut; each package keeps its own {@code Model} interface extending this one).
 */
public interface Model {
	public int numberOfStations();
}
