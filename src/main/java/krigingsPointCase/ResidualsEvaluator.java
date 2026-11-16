/*
 * GNU GPL v3 License (GEOframe Krigings, version 0.9.9 of 2023-01-11).
 */
/* Source reconstructed from Kriging-2023-01-11-0.9.9.jar (CFR 0.152 decompilation, cleaned 2026; see UPSTREAM_DIFF.md). */
package krigingsPointCase;

import flanagan.analysis.Regression;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;

public class ResidualsEvaluator {
    @Description(value="Switch for detrended mode.")
    @In
    public boolean doDetrended;
    @Description(value="Switch for detrended mode.")
    @In
    public int regressionOrder = 1;
    @In
    public double[] zStations;
    @In
    public double[] hStations;
    @Out
    public double trend_intercept;
    @Out
    public double trend_coefficient;
    @In
    @Out
    public double[] hResiduals;
    @Out
    public boolean isPValueOk = false;

    @Execute
    public void process() {
        this.hResiduals = this.hStations;
        this.trend_intercept = 0.0;
        this.trend_coefficient = 0.0;
        this.isPValueOk = false;
        try {
            if (this.doDetrended && this.zStations.length > 2) {
                Regression r = new Regression();
                r = new Regression(this.zStations, this.hStations);
                r.polynomial(this.regressionOrder);
                if (r.getPvalues()[1] < 0.05) {
                    this.isPValueOk = true;
                    this.trend_intercept = r.getBestEstimates()[0];
                    this.trend_coefficient = r.getBestEstimates()[1];
                    this.hResiduals = r.getResiduals();
                } else {
                    this.trend_intercept = 0.0;
                    this.trend_coefficient = 0.0;
                    this.hResiduals = this.hStations;
                }
            }
        }
        catch (Exception e) {
            this.isPValueOk = false;
            this.hResiduals = this.hStations;
            this.trend_intercept = 0.0;
            this.trend_coefficient = 0.0;
        }
    }
}

