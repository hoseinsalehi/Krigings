#!/usr/bin/env python3
"""Java (Krigings working copy, pokriging) vs Python reference/: grid and LOO differences per method.
Usage: compare_equivalence.py <method> <java_output_dir> [--json out.json]"""
import sys, json, numpy as np, pandas as pd
from pathlib import Path
B = Path(__file__).resolve().parents[1]
m, jd = sys.argv[1], Path(sys.argv[2])
ref = B / "reference" / "daily" / m
res = {"method": m}
COLS = ["value", "kstd", "D1", "D7", "D50", "range", "nugget_fraction", "log10_sill", "lapse", "lapse_r2", "n_negative_weights",
        "n_eff", "R_pos_half", "R_pos_2_3", "R_pos_1_3", "R_abs90", "qc_flags", "pi", "amount", "samount", "wet", "threshold",
        "range_amount", "R_pos_half_amount", "wet_probability"]
def stats(a, b, name):
    a = np.asarray(a, float); b = np.asarray(b, float)
    both_nan = np.isnan(a) & np.isnan(b); fin = np.isfinite(a) & np.isfinite(b)
    inf_eq = np.isinf(a) & np.isinf(b) & (np.sign(a) == np.sign(b))
    mism_nan = ~(both_nan | fin | inf_eq)
    d = np.abs(a[fin] - b[fin]); rel = d / np.maximum(np.abs(b[fin]), 1e-12)
    return {"n": int(len(a)), "nan_mismatch": int(mism_nan.sum()), "max_abs": float(d.max()) if len(d) else 0.0,
            "p99_abs": float(np.percentile(d, 99)) if len(d) else 0.0, "frac_rel_le_1e-6": float(np.mean(rel <= 1e-6)) if len(d) else 1.0,
            "frac_exact_10sig": float(np.mean(rel <= 5e-10)) if len(d) else 1.0}
g = pd.read_csv(ref / "grid.csv.gz"); gj = pd.concat([pd.read_csv(f) for f in sorted(jd.glob("grid_*.csv.gz"))], ignore_index=True)
mg = g.merge(gj, on=["date", "cell"], suffixes=("_py", "_j"), how="outer", indicator=True)
res["grid_rows"] = {"python": len(g), "java": len(gj), "matched": int((mg._merge == "both").sum())}
mg = mg[mg._merge == "both"]
res["grid"] = {c: stats(mg[c + "_j"], mg[c + "_py"], c) for c in COLS if c + "_py" in mg and c + "_j" in mg}
if (ref / "loo.csv.gz").exists() and list(jd.glob("loo_*.csv")):
    l = pd.read_csv(ref / "loo.csv.gz", dtype={"station": str}); lj = pd.concat([pd.read_csv(f, dtype={"station": str}) for f in sorted(jd.glob("loo_*.csv"))])
    ml = l.merge(lj, on=["date", "station"], suffixes=("_py", "_j"))
    res["loo_rows"] = {"python": len(l), "java": len(lj), "matched": len(ml)}
    res["loo"] = {c: stats(ml[c + "_j"], ml[c + "_py"], c) for c in ["obs", "pred", "kriging_std", "D1", "D7", "qc_flags", "pi", "amount", "wet_pred"] if c + "_py" in ml}
    def metrics(o, p, precip):
        ok = np.isfinite(o) & np.isfinite(p); o, p = o[ok], p[ok]
        r = np.corrcoef(o, p)[0, 1]; kge = 1 - np.sqrt((r - 1) ** 2 + (p.std() / o.std() - 1) ** 2 + (p.mean() / o.mean() - 1) ** 2)
        out = {"n": int(ok.sum()), "KGE": float(kge), "MAE": float(np.mean(np.abs(p - o))), "bias": float(np.mean(p - o))}
        if precip:
            out["wet_freq_obs_pct"] = float(100 * np.mean(o >= 1)); out["wet_freq_pred_pct"] = float(100 * np.mean(p >= 1))
            out["intensity_obs"] = float(o[o >= 1].mean()); out["intensity_pred"] = float(p[p >= 1].mean())
        return out
    pr = m.startswith("precip")
    res["loo_metrics_python"] = metrics(ml.obs_py.to_numpy(float), ml.pred_py.to_numpy(float), pr)
    res["loo_metrics_java"] = metrics(ml.obs_j.to_numpy(float), ml.pred_j.to_numpy(float), pr)
txt = json.dumps(res, indent=1)
if "--json" in sys.argv: Path(sys.argv[sys.argv.index("--json") + 1]).write_text(txt)
print(txt)
