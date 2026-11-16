"""
s14_paper_tables.py -- station-level tables behind the manuscript figures of v2.0 (no plotting).

Writes to OUT_ROOT/paper_tables/:
  stations_<var>.csv                station metadata used in the run: id, x, y (EPSG:3035), lon, lat, elevation, group
                                    (inside district / buffer / MeteoSwiss), has_normal, valid days per year (wide)
  loo_station_<var>.csv             station-wise LOO metrics (KGE Gupta 2009, CC, MAE, RMSE, bias; precip all + wet days)
  loo_station_year_<var>.csv        same per station and calendar year (>= 200 valid days)
  dayshift_station_<var>.csv        station-wise LOO metrics without / with day-label harmonisation (if the
                                    dayshift_off_/dayshift_on_ LOO runs exist) + shift class per station
Station observations themselves are not written.
"""
import os
import glob, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd, geopandas as gpd
import po_config as C

OD = C.OUT_ROOT / "paper_tables"; OD.mkdir(parents=True, exist_ok=True)


def metrics(o, p, v):
    if len(o) < 2 or o.std() == 0 or p.std() == 0: return None
    r = np.corrcoef(o, p)[0, 1]; a = p.std() / o.std()
    if v == "precip":
        b = p.mean() / o.mean() if o.mean() > 0 else np.nan; kge = 1 - np.sqrt((r - 1) ** 2 + (a - 1) ** 2 + (b - 1) ** 2); bias = 100 * (b - 1)
    else:
        kge = 1 - np.sqrt((r - 1) ** 2 + (a - 1) ** 2 + ((p.mean() - o.mean()) / o.std()) ** 2); bias = p.mean() - o.mean()
    return dict(n=len(o), kge=kge, cc=r, alpha=a, mae=np.abs(p - o).mean(), rmse=np.sqrt(((p - o) ** 2).mean()), bias=bias)


def station_metrics(d, v, by=("station",), min_n=365):
    rows = []
    for k, g in d.groupby(list(by)):
        o = g.obs.to_numpy(float); p = g.pred.to_numpy(float)
        if len(o) < min_n: continue
        m = metrics(o, p, v)
        if m is None: continue
        row = dict(zip(by, k if isinstance(k, tuple) else (k,))); row.update(days="all", **m); rows.append(row)
        if v == "precip":
            w = o >= 0.1
            if w.sum() >= min_n // 4:
                mw = metrics(o[w], p[w], v)
                if mw: row2 = dict(zip(by, k if isinstance(k, tuple) else (k,))); row2.update(days="wet", **mw); rows.append(row2)
    return pd.DataFrame(rows)


def load_loo(run, v):
    fs = sorted(glob.glob(str(C.OUT_ROOT / run / f"loo_{v}_*.parquet")))
    if not fs: return None
    d = pd.concat([pd.read_parquet(f, columns=["date", "station", "obs", "pred"]) for f in fs]).dropna()
    d["station"] = d.station.astype(str); return d


def main():
    import po_common as P
    B = gpd.read_file(os.environ.get("PO_INPUT_ROOT", "data") + "/Po_District.shp").to_crs(3035).union_all()
    for v in ("precip", "temp"):
        st = gpd.read_file(C.DATA[v]["stations"]).to_crs(3035); st["ID"] = st.ID.astype(str)
        d = P.load_station_data(v); ids = list(d["ids"]); V = pd.DataFrame(d["values"], index=d["dates"], columns=ids)
        st = st.set_index("ID").loc[ids]; ll = st.to_crs(4326)
        inside = st.within(B).to_numpy()
        grp = np.where(pd.Index(ids).str.startswith("CH_"), "MeteoSwiss", np.where(inside, "inside district", "buffer (20 km)"))
        nd = pd.read_csv(C.OUT_ROOT / "normals" / f"station_normals_{v}.csv", dtype={"station": str}); hasn = nd.groupby("station").ok.all()
        meta = pd.DataFrame(dict(id=ids, x=st.geometry.x.to_numpy(), y=st.geometry.y.to_numpy(), lon=ll.geometry.x.to_numpy(), lat=ll.geometry.y.to_numpy(),
                                 elevation=pd.to_numeric(st[C.STATION_ELEV_FIELD], errors="coerce").to_numpy(), group=grp,
                                 has_normal=pd.Series(ids).map(hasn).fillna(False).to_numpy()))
        avail = V.notna().groupby(V.index.year).sum().T; avail.columns = [f"days_{c}" for c in avail.columns]
        meta = meta.join(avail.reset_index(drop=True)); meta.to_csv(OD / f"stations_{v}.csv", index=False)
        L = load_loo(f"{v}_anomaly", v)
        sm = station_metrics(L, v).merge(meta[["id", "x", "y", "lon", "lat", "elevation", "group"]], left_on="station", right_on="id").drop(columns="id")
        sm.to_csv(OD / f"loo_station_{v}.csv", index=False)
        L["year"] = pd.to_datetime(L.date).dt.year
        sy = station_metrics(L, v, by=("station", "year"), min_n=200); sy.to_csv(OD / f"loo_station_year_{v}.csv", index=False)
        on, off = load_loo(f"dayshift_on_{v}", v), load_loo(f"dayshift_off_{v}", v)
        if on is not None and off is not None:
            sh = pd.read_csv(C.OUT_ROOT / "harmonisation" / f"day_shift_{v}.csv", dtype={"station": str})
            cls = sh.groupby("station").shift_days.agg(lambda s: "unchanged" if (s == 0).all() else ("-1 day" if (s == -1).all() else ("+1 day" if (s == 1).all() else "changed over time")))
            a = station_metrics(off, v); b = station_metrics(on, v)
            m = a.merge(b, on=["station", "days"], suffixes=("_original", "_harmonised"))
            m["shift_class"] = m.station.map(cls).fillna("unchanged")
            m = m.merge(meta[["id", "x", "y", "lon", "lat", "elevation", "group"]], left_on="station", right_on="id").drop(columns="id")
            m.to_csv(OD / f"dayshift_station_{v}.csv", index=False)
            print(v, "dayshift stations", len(m), m.groupby("shift_class").size().to_dict())
        print(v, "stations", len(meta), "LOO stations", sm.station.nunique(), "station-years", len(sy))
    print("->", OD)


if __name__ == "__main__":
    main()
