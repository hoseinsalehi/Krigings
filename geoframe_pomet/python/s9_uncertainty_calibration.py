#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 9 - calibration of the published uncertainty with the leave-one-out errors of the production run.

The kriging standard deviation is a model quantity; its 95 % interval does not cover 95 % of the LOO errors
(precipitation especially). For each class (season x elevation band for temperature; season x predicted
intensity for precipitation) this writes
  coverage95_raw     share of LOO errors inside +-1.96 kriging_std
  scale_k            factor so that +-1.96 * k * kriging_std covers 95 % (k = q95(|e| / std) / 1.96)
  coverage95_cal     coverage after scaling, evaluated on the other half of the stations (split by station,
                     k fitted on even, tested on odd station index and vice versa) - an honest check
  e_q05 ... e_q95    empirical quantiles of obs - estimate (a distribution-free alternative)
Usage:  python s9_uncertainty_calibration.py --variable precip|temp [--run <var>_anomaly]
Output: <OUT_ROOT>/summary/uncertainty_calibration_<var>.csv
"""
import argparse, glob, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import po_config as C, po_common as P

SEAS = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--variable", choices=["precip", "temp"], required=True); ap.add_argument("--run", default=None)
    a = ap.parse_args(); v = a.variable; run = a.run or f"{v}_anomaly"
    fs = sorted(glob.glob(str(C.OUT_ROOT / run / f"loo_{v}_*.parquet")))
    if not fs: sys.exit("no LOO files")
    L = pd.concat([pd.read_parquet(f, columns=["date", "station", "obs", "pred", "kriging_std"]) for f in fs], ignore_index=True).dropna()
    L["station"] = L.station.astype(str); L["season"] = pd.to_datetime(L.date).dt.month.map(SEAS)
    st = P.load_station_data(v); el = pd.Series(st["elev"], index=st["ids"]); L["elev"] = L.station.map(el)
    if v == "precip":
        L["cls"] = pd.cut(L.pred, [-1, 0.0, 1, 5, 20, 50, 1e9], labels=["0", "0-1", "1-5", "5-20", "20-50", ">50"]).astype(str)
    else:
        L["cls"] = pd.cut(L.elev, [-1e3, 500, 1000, 1500, 2000, 9e3], labels=["<500", "500-1000", "1000-1500", "1500-2000", ">=2000"]).astype(str)
    L["e"] = L.obs - L.pred; L["half"] = pd.factorize(L.station)[0] % 2
    rows = []
    for (s, c), g in list(L.groupby(["season", "cls"])) + [(("all", "all"), L)]:
        pos = g.kriging_std > 0; z = (g.e.abs() / g.kriging_std.where(pos)).dropna()
        k = float(np.quantile(z, 0.95) / 1.96) if len(z) > 50 else np.nan
        cov_cal = []
        for h in (0, 1):
            tr = g[pos & (g.half == h)]; te = g[pos & (g.half != h)]
            if len(tr) > 50 and len(te) > 50:
                kh = np.quantile(tr.e.abs() / tr.kriging_std, 0.95) / 1.96; cov_cal.append(np.mean(te.e.abs() <= 1.96 * kh * te.kriging_std))
        q = g.e.quantile([.05, .25, .5, .75, .95]).to_numpy()
        rows.append(dict(variable=v, season=s, cls=c, n=len(g), share_zero_std=float((~pos).mean()),
                         coverage95_raw=float(np.mean(g.e[pos].abs() <= 1.96 * g.kriging_std[pos])) if pos.any() else np.nan, scale_k=k,
                         coverage95_cal_crossval=float(np.mean(cov_cal)) if cov_cal else np.nan,
                         e_q05=q[0], e_q25=q[1], e_q50=q[2], e_q75=q[3], e_q95=q[4]))
    T = pd.DataFrame(rows); od = C.OUT_ROOT / "summary"; od.mkdir(exist_ok=True); T.to_csv(od / f"uncertainty_calibration_{v}.csv", index=False)
    print(T.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
