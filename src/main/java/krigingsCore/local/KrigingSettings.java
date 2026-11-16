package krigingsCore.local;


/** Method parameters of the daily kriging; defaults = locked Po values (po_config.py). */
public final class KrigingSettings {
    public int nVario = 50, nKrig = 7, nLags = 12, minPairsPerLag = 5, chunk = 2048;
    public double lapseMinSpanM = 300, anomLapseMaxExtrapM = 300; public int lapseMinStations = 10, lapseFallbackN = 100;
    public double zeroTol = 0.0, oscSumAbsW = 2.0, negWeightTol = -1e-10;
    public boolean globalVariogramFallback = true;
    /** numpy: PCG64/SeedSequence stream seeded per chunk with BLAKE2b(dayTag + seedIndex of the chunk's first target). */
    public String rngMode = "numpy";
    public final SphericalPso.Settings pso = new SphericalPso.Settings();

    public static final int FLAG_ZERO = 1, FLAG_LAPSE_FALLBACK = 2, FLAG_POSITIVE_LAPSE = 4, FLAG_OSCILLATORY = 8,
            FLAG_RANGE_BOUND = 16, FLAG_GLOBAL_VARIOGRAM = 32, FLAG_ELEV_FILLED = 64, FLAG_NORMAL_FALLBACK = 128;
}
