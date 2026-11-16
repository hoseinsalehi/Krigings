package krigingsCore;

import static org.hortonmachine.gears.libs.modules.HMConstants.isNovalue;

import theoreticalVariogram.TheoreticalVariogram;

/**
 * Ordinary-kriging system of the original Krigings components (covariance matrix, known terms and variogram of a
 * 3-D distance), shared by krigingsPointCase, krigingsRasterCase and leaveOneOut. Code and operation order are the
 * original ones.
 */
public final class KrigingSystem {

	private KrigingSystem() {
	}

	/** Matrix of semivariances among the n stations, bordered by ones (Lagrange row/column). */
	public static double[][] covMatrix(double[] x, double[] y, double[] z, int n, String type, double nugget, double range, double sill) {
		double[][] covarianceMatrix = new double[n + 1][n + 1];
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				double rx = x[i] - x[j];
				double ry = y[i] - y[j];
				double rz = z[i] - z[j];
				covarianceMatrix[j][i] = variogram(type, nugget, range, sill, rx, ry, rz);
				covarianceMatrix[i][j] = variogram(type, nugget, range, sill, rx, ry, rz);
			}
		}
		for (int i = 0; i < n; i++) {
			covarianceMatrix[i][n] = 1.0;
			covarianceMatrix[n][i] = 1.0;
		}
		covarianceMatrix[n][n] = 0;
		return covarianceMatrix;
	}

	/** Semivariances between the n stations and the target stored at index n, plus 1. */
	public static double[] knownTerms(double[] x, double[] y, double[] z, int n, String type, double nugget, double range, double sill) {
		double[] gamma = new double[n + 1];
		for (int i = 0; i < n; i++) {
			double rx = x[i] - x[n];
			double ry = y[i] - y[n];
			double rz = z[i] - z[n];
			gamma[i] = variogram(type, nugget, range, sill, rx, ry, rz);
		}
		gamma[n] = 1.0;
		return gamma;
	}

	/** Semivariance at the 3-D separation (rx, ry, rz); novalue rz counts as 0; gamma(0) = 0. */
	public static double variogram(String type, double nug, double range, double sill, double rx, double ry, double rz) {
		if (isNovalue(rz)) {
			rz = 0;
		}
		double h2 = Math.sqrt(rx * rx + rz * rz + ry * ry);
		double vgmResult;
		if (h2 != 0) {
			TheoreticalVariogram vgm = new TheoreticalVariogram();
			vgmResult = vgm.calculateVGM(type, h2, sill, range, nug);
		} else {
			vgmResult = 0;
		}
		return vgmResult;
	}
}
