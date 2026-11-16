"""Station-wise LOO metrics of the final precipitation product (two-part, v2.0) from <V5>/precip_anomaly/loo_precip_*.parquet.
Same metric definitions as s14_paper_tables.py (KGE Gupta 2009; wet day >= 0.1 mm). Writes ../data/loo_station[_year]_precip.csv"""
import sys, glob, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path.home()/"Downloads/PoGrid_v3_code"))
import importlib.util
spec = importlib.util.spec_from_file_location("s14", str(Path.home()/"Downloads/PoGrid_v3_code/s14_paper_tables.py")); s14 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s14)
V5 = Path.home()/"Downloads/PoGrid_v5"; OD = Path(__file__).resolve().parent.parent/"data"
fs = sorted(glob.glob(str(V5/"precip_anomaly"/"loo_precip_*.parquet")))
L = pd.concat([pd.read_parquet(f, columns=["date","station","obs","pred"]) for f in fs]).dropna(); L["station"] = L.station.astype(str)
meta = pd.read_csv(Path.home()/"Downloads/PoGrid_v4/paper_tables/stations_precip.csv", dtype={"id": str})
sm = s14.station_metrics(L, "precip").merge(meta[["id","x","y","lon","lat","elevation","group"]], left_on="station", right_on="id").drop(columns="id")
sm.to_csv(OD/"loo_station_precip.csv", index=False)
L["year"] = pd.to_datetime(L.date).dt.year
sy = s14.station_metrics(L, "precip", by=("station","year"), min_n=200); sy.to_csv(OD/"loo_station_year_precip.csv", index=False)
print("LOO stations", L.station.nunique(), "station-days", len(L))
print(sm.groupby("days")[["kge","cc","mae","rmse","bias"]].mean().round(3)); print(sm.groupby("days").size())
