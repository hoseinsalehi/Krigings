package pokriging;

import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;
import pokriging.uncertainty.UncertaintyCalibration;

/**
 * OMS3 component: calibration of the kriging standard deviation with the leave-one-out errors of a PoKriging run (port of
 * the Python step s9_uncertainty_calibration.py; method in docs/UNCERTAINTY_CALIBRATION.md). Reads the LOO tables
 * loo_&lt;variable&gt;_&lt;Y&gt;.csv of the period and writes uncertainty_&lt;variable&gt;_&lt;start&gt;_&lt;end&gt;.csv (variable,season,cls,
 * n,share_zero_std,coverage95_raw,scale_k,coverage95_cal_crossval,e_q05..e_q95; the input of PoKriging uncertaintyTable)
 * and uncertainty_summary_&lt;variable&gt;_&lt;start&gt;_&lt;end&gt;.csv (coverage of the nominal levels, station-half cross-validation).
 * Chain: PoKriging.outLooDir -&gt; PoUncertainty.looDir.
 */
@Description("LOO calibration of the kriging standard deviation: k(season, class) table and coverage summary")
@Author(name = "GEOframe Krigings contributors; Po extension 2026")
@Keywords("Kriging, uncertainty, leave-one-out, calibration, coverage")
@Name("PoUncertainty")
@Status(Status.EXPERIMENTAL)
@License("General Public License Version 3 (GPLv3)")
public class PoUncertainty {

    @Description("false = do nothing (outputs null)") @In public boolean enabled = true;
    @Description("precip | temp") @In public String variable;
    @Description("folder with loo_<variable>_<Y>.csv[.gz] (from PoKriging.outLooDir or set directly)") @In public String looDir;
    @Description("first LOO day used (YYYY-MM-DD)") @In public String periodStart = "1991-01-01";
    @Description("last LOO day used (YYYY-MM-DD)") @In public String periodEnd = "2020-12-31";
    @Description("stations.csv with id and elevation (required for temp: elevation classes)") @In public String stationsFile;
    @Description("output folder (empty = looDir)") @In public String outputDir = "";
    @Description("nominal central-interval levels of the summary (must contain 0.95)") @In public String nominalLevels = "0.5,0.8,0.9,0.95";
    @Description("a factor k needs more than this many values (s9: 50)") @In public int minCount = 50;

    @Description("uncertainty_<variable>_<start>_<end>.csv (give it to PoKriging.uncertaintyTable of later runs)") @Out public String outUncertaintyTable;
    @Description("uncertainty_summary_<variable>_<start>_<end>.csv") @Out public String outUncertaintySummary;

    @Execute
    public void execute() throws Exception {
        if (!enabled) { System.out.println("PoUncertainty: disabled"); return; }
        if (looDir == null || looDir.isEmpty()) throw new IllegalArgumentException("PoUncertainty: looDir is empty (PoKriging mode must be loo or both)");
        UncertaintyCalibration.Config c = new UncertaintyCalibration.Config();
        c.variable = variable; c.looDir = looDir; c.periodStart = periodStart; c.periodEnd = periodEnd; c.stationsFile = stationsFile;
        c.outputDir = outputDir; c.nominalLevels = nominalLevels; c.minCount = minCount;
        UncertaintyCalibration u = new UncertaintyCalibration(c);
        u.run();
        outUncertaintyTable = u.outTable.toString(); outUncertaintySummary = u.outSummary.toString();
    }
}
