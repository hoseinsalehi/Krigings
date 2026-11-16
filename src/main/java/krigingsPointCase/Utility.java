/*
 * GNU GPL v3 License (GEOframe Krigings, version 0.9.9 of 2023-01-11).
 */
/* Source reconstructed from Kriging-2023-01-11-0.9.9.jar (CFR 0.152 decompilation, cleaned 2026; see UPSTREAM_DIFF.md). */
package krigingsPointCase;

import experimentalVariogram.ExperimentalVariogram;
import org.geotools.data.simple.SimpleFeatureCollection;

public class Utility {
    public static String[] availableTheorethicalVariogra = new String[]{"exponential", "linear", "power", "spherical"};

    public static final int getVariogramCode(String name) {
        for (int i = 0; i < availableTheorethicalVariogra.length; ++i) {
            if (availableTheorethicalVariogra[i] != name) continue;
            return i;
        }
        return -9999;
    }

    public static final double getLog(double h) {
        return Math.log(h + 1.0);
    }

    public static final double getInverseLog(double h) {
        return Math.exp(h) - 1.0;
    }

    public static final ExperimentalVariogram getExperimentalVariogram(String fStationsid, SimpleFeatureCollection inStations, boolean doIncludeZero, int cutoffDivide, double cutoffInput, int numCloseStation) {
        ExperimentalVariogram expVariogram = new ExperimentalVariogram();
        expVariogram.fStationsid = fStationsid;
        expVariogram.inStations = inStations;
        expVariogram.doIncludezero = doIncludeZero;
        expVariogram.inNumCloserStations = numCloseStation;
        if (cutoffDivide > 0) {
            expVariogram.Cutoff_divide = cutoffDivide;
        }
        if (cutoffInput > 0.0) {
            expVariogram.Cutoffinput = cutoffInput;
        }
        return expVariogram;
    }
}

