"""
s00_add_meteoswiss.py -- add the official MeteoSwiss stations (Open Data, ch.meteoschweiz.ogd-smn automatic network and
ch.meteoschweiz.ogd-nime manual precipitation network) located inside the Po River District or within 20 km of its
boundary (the same rule the existing network follows) and not already present in the network.

Steps
  1. read daily series 1991-2020 from ~/Downloads/MeteoSwiss_stations
       precip SMN : rka150d0 (00-00 UTC); days missing it are filled with rre150d0 (06-06 UTC, start-day label)
       precip NIME: rre150d0 (06-06 UTC)  -- the day-alignment step s0 tests/harmonises the label per station-year
       temp   SMN : tre200d0 (24-h mean); tre200dx / tre200dn used for consistency QC only
  2. pre-QC (MeteoSwiss data are already quality controlled; these are gross-error checks only)
       precip : < 0 or > 500 mm/day; >= 5 consecutive identical non-zero values
       temp   : outside -45..+40 C; Tmean outside [Tmin-0.5, Tmax+0.5]; >= 5 consecutive identical values
       neighbour-based checks are left to s7_station_qc.py, which treats new and old stations identically
  3. duplicate test against every network station within 15 km (raw files, lags -1/0/+1 to allow for the
     original day labelling): identical-value fraction (precip on days with max >= 0.5 mm, |diff| <= 0.1 mm;
     temp |diff| <= 0.05 C). >= 50 % identical -> same station already in the network -> not added.
     Co-location: a network gauge within 1 km whose 1991-2020 record is >= 50 % as long -> not added (same site,
     typically the same gauge stored with another day window, ~20 % identical values). If the network record of that
     co-located gauge is < 50 % as long and >= 15 % of values are identical, the full MeteoSwiss series replaces it.
  4. coverage: >= 365 valid days in 1991-2020
  5. write the extended network to OUT (originals untouched):
       Final_precip_plusCH.csv, STATIONS_plusCH_precip.shp, Final_temp_plusCH.csv, STATIONS_plusCH_temp.shp,
       meteoswiss_added_stations.csv (one row per candidate with the decision and all diagnostics),
       meteoswiss_preqc_flags.csv
IDs of added stations: CH_<abbr> (e.g. CH_GSB). Run: python s00_add_meteoswiss.py
"""
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd, geopandas as gpd
import po_config as C

MS = Path(os.environ.get("PO_INPUT_ROOT", "data") + "/MeteoSwiss_stations")
OUT = Path(os.environ.get("PO_INPUT_ROOT", "data") + "/PoGrid_v4_input")
BOUNDARY = Path(os.environ.get("PO_INPUT_ROOT", "data") + "/Po_District.shp")
BUFFER_M, DUP_RADIUS_M, DUP_FRAC, MIN_DAYS, COLOC_KM = 20000.0, 15000.0, 0.50, 365, 1.0
DATES = pd.date_range("1991-01-01", "2020-12-31", freq="D")
RAW = {"precip": dict(obs=C.MG / "Precip_New/Final_precip.csv", stations=C.MG / "Final_Data/STATIONS_FINAL_CLEANED_V3.shp", fmt="%d/%m/%y"),
       "temp": dict(obs=C.MG / "Final_Data/Final_Merged_Temp_Timeseries_Cleaned_V2_Filtered.csv",
                    stations=C.MG / "Final_Data/Temp_new_all_Filtered.shp", fmt="%d/%m/%y %H:%M")}


def runs_identical(s, n=5, nonzero=False):
    v = s.to_numpy(float); same = np.r_[False, np.isclose(v[1:], v[:-1])] & np.isfinite(v)
    if nonzero: same &= np.r_[False, v[1:] > 0]
    grp = np.cumsum(~same); ln = pd.Series(1, index=s.index).groupby(grp).transform("sum").to_numpy()
    return pd.Series((ln >= n) & np.isfinite(v), index=s.index)      # run length counted from the first value of the run


def read_ms(net, abbr):
    fs = [MS / f"ogd-{net}_{abbr.lower()}_d_{k}.csv" for k in ("historical", "recent")]
    d = pd.concat([pd.read_csv(f, sep=";", na_values=["", "-"]) for f in fs if f.exists()])
    d["date"] = pd.to_datetime(d["reference_timestamp"], format="%d.%m.%Y %H:%M").dt.floor("D")
    return d.drop_duplicates("date").set_index("date").sort_index().reindex(DATES)


def candidates():
    B = gpd.read_file(BOUNDARY).to_crs(C.METRIC_CRS).union_all(); rows = []
    for net in ("smn", "nime"):
        m = pd.read_csv(MS / f"ogd-{net}_meta_stations.csv", sep=";", encoding="latin1")
        g = gpd.GeoDataFrame(m, geometry=gpd.points_from_xy(m.station_coordinates_wgs84_lon, m.station_coordinates_wgs84_lat), crs=4326).to_crs(C.METRIC_CRS)
        g["dist_boundary_km"] = g.distance(B) / 1e3; g["inside"] = g.within(B); g["net"] = net
        rows.append(g[g.distance(B) <= BUFFER_M])
    g = pd.concat(rows, ignore_index=True)
    return gpd.GeoDataFrame(g, geometry="geometry", crs=C.METRIC_CRS)


def build_series(c):
    P, T, flags = {}, {}, []
    for _, r in c.iterrows():
        a = r.station_abbr; d = read_ms(r.net, a)
        if r.net == "smn":
            if "rka150d0" in d and "rre150d0" in d: p = d["rka150d0"].fillna(d["rre150d0"]); src = "rka150d0+rre150d0"
            elif "rre150d0" in d: p = d["rre150d0"]; src = "rre150d0"
            else: p = None
            if "tre200d0" in d:
                t = d["tre200d0"].astype(float).copy(); tx, tn = d.get("tre200dx"), d.get("tre200dn")
                bad = {"range": (t < -45) | (t > 40), "stuck5": runs_identical(t)}
                if tx is not None and tn is not None: bad["mean_outside_min_max"] = (t < tn - 0.5) | (t > tx + 0.5)
                for k, m in bad.items():
                    m = m.fillna(False).astype(bool)
                    for dt in t.index[m]: flags.append(dict(variable="temp", station=f"CH_{a}", date=dt, test=k, value=t[dt]))
                    t[m] = np.nan
                if t.notna().sum() > 0: T[a] = (t, r)
        else:
            if a in P: continue                               # SMN record of the same station preferred
            p = d.get("rre150d0"); src = "rre150d0"
        if p is not None:
            p = p.astype(float).copy()
            bad = {"range": (p < 0) | (p > 500), "stuck5_nonzero": runs_identical(p, nonzero=True)}
            for k, m in bad.items():
                m = m.fillna(False).astype(bool)
                for dt in p.index[m]: flags.append(dict(variable="precip", station=f"CH_{a}", date=dt, test=k, value=p[dt]))
                p[m] = np.nan
            if p.notna().sum() > 0: P[a] = (p, r, src)
    return P, T, pd.DataFrame(flags)


def load_net(v):
    g = gpd.read_file(RAW[v]["stations"]); gm = g.to_crs(C.METRIC_CRS)
    df = pd.read_csv(RAW[v]["obs"], low_memory=False); dates = pd.to_datetime(df["Date"], format=RAW[v]["fmt"]).dt.floor("D")
    df.index = dates; df = df.drop(columns="Date").apply(pd.to_numeric, errors="coerce").replace(C.MISSING_VALUE, np.nan)
    return g, gm, df.reindex(DATES)


def dup_test(v, s, xy, gm, net):
    d = np.hypot(gm.geometry.x - xy[0], gm.geometry.y - xy[1]); near = gm[d <= DUP_RADIUS_M].assign(d=d[d <= DUP_RADIUS_M]).sort_values("d")
    best = dict(dup_id=None, dup_frac=0.0, dup_lag=None, dup_dist_km=np.nan, dup_r=np.nan, nearest_id=None, nearest_km=np.nan)
    if len(near): best.update(nearest_id=str(near.iloc[0][C.STATION_ID_FIELD]), nearest_km=near.iloc[0]["d"] / 1e3)
    for _, n in near.iterrows():
        sid = str(n[C.STATION_ID_FIELD])
        if sid not in net.columns: continue
        o = net[sid]
        for L in (-1, 0, 1):
            oo = o.shift(L); ok = s.notna() & oo.notna()
            if v == "precip": ok &= np.maximum(s, oo) >= 0.5; tol = 0.1
            else: tol = 0.05
            if ok.sum() < 60: continue
            fr = float((np.abs(s[ok] - oo[ok]) <= tol).mean())
            if fr > best["dup_frac"]:
                best.update(dup_id=sid, dup_frac=fr, dup_lag=L, dup_dist_km=n["d"] / 1e3, dup_r=float(np.corrcoef(s[ok], oo[ok])[0, 1]))
    return best


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    c = candidates(); print(f"candidates within {BUFFER_M/1e3:.0f} km: {len(c)} ({(c.net=='smn').sum()} SMN, {(c.net=='nime').sum()} NIME)")
    P, T, flags = build_series(c); flags.to_csv(OUT / "meteoswiss_preqc_flags.csv", index=False)
    report = []
    for v, S in (("precip", P), ("temp", T)):
        g, gm, net = load_net(v); add = {}; drop = []
        for a, item in S.items():
            s, r = item[0], item[1]; src = item[2] if v == "precip" else "tre200d0"
            xy = (r.geometry.x, r.geometry.y); n_ok = int(s.notna().sum())
            dup = dup_test(v, s, xy, gm, net)
            nn = dup["nearest_id"]; nn_days = int(net[nn].notna().sum()) if nn in net.columns else 0
            dup["nearest_valid_days"] = nn_days
            if dup["dup_frac"] >= DUP_FRAC: dec = "duplicate_of_network"
            elif dup["nearest_km"] < COLOC_KM and nn_days >= 0.5 * n_ok: dec = "co_located_with_network"
            elif dup["nearest_km"] < COLOC_KM and dup["dup_frac"] >= 0.15 and n_ok >= MIN_DAYS:
                dec = "replaces_short_network_record"; add[f"CH_{a}"] = (s, r); drop.append(nn)
            elif n_ok < MIN_DAYS: dec = "too_short"
            else: dec = "added"; add[f"CH_{a}"] = (s, r)
            report.append(dict(variable=v, id=f"CH_{a}", abbr=a, name=r.station_name, network=r.net.upper(), source_field=src,
                               elev_m=r.station_height_masl, lon=r.station_coordinates_wgs84_lon, lat=r.station_coordinates_wgs84_lat,
                               inside_district=bool(r.inside), dist_boundary_km=round(r.dist_boundary_km, 2), valid_days_1991_2020=n_ok,
                               first=s.first_valid_index(), last=s.last_valid_index(),
                               preqc_removed=int(((flags.station == f"CH_{a}") & (flags.variable == v)).sum()) if len(flags) else 0,
                               decision=dec, **dup))
        # ---- write extended files
        df = pd.read_csv(RAW[v]["obs"], low_memory=False); dates = pd.to_datetime(df["Date"], format=RAW[v]["fmt"]).dt.floor("D")
        di = pd.DatetimeIndex(dates); print(f"{v}: network file {di[0].date()} .. {di[-1].date()} ({len(di)} days)")
        df = df.drop(columns=[c for c in df.columns if str(c) in set(drop)])
        newcols = pd.DataFrame({k: s.reindex(di).round(2).fillna(C.MISSING_VALUE).to_numpy() for k, (s, _) in add.items()})
        pd.concat([df, newcols], axis=1).to_csv(OUT / f"Final_{v}_plusCH.csv", index=False)
        gs = g.copy(); gs[C.STATION_ID_FIELD] = gs[C.STATION_ID_FIELD].astype(str); gs = gs[~gs[C.STATION_ID_FIELD].isin(drop)]
        newg = gpd.GeoDataFrame({C.STATION_ID_FIELD: list(add), C.STATION_ELEV_FIELD: [float(r.station_height_masl) for _, r in add.values()]},
                                geometry=gpd.GeoSeries([r.geometry for _, r in add.values()], crs=C.METRIC_CRS).to_crs(g.crs), crs=g.crs)
        for col in gs.columns:
            if col not in newg.columns and col != "geometry": newg[col] = np.nan
        if "x" in gs.columns and "y" in gs.columns: newg["x"] = newg.geometry.x; newg["y"] = newg.geometry.y
        out = pd.concat([gs, newg[gs.columns]], ignore_index=True)
        gpd.GeoDataFrame(out, geometry="geometry", crs=g.crs).to_file(OUT / f"STATIONS_plusCH_{v}.shp")
        print(f"{v}: network {len(g)} stations -> {len(out)} (+{len(add)} added, -{len(drop)} short duplicates replaced: {drop})")
    R = pd.DataFrame(report); R.to_csv(OUT / "meteoswiss_added_stations.csv", index=False)
    pd.set_option("display.width", 250)
    print(R.groupby(["variable", "decision"]).size().to_string())
    print(R[["variable", "id", "name", "elev_m", "inside_district", "valid_days_1991_2020", "decision", "dup_id", "dup_frac", "dup_r", "nearest_km", "nearest_valid_days"]]
          .sort_values(["variable", "decision", "id"]).round(2).to_string(index=False))
    if len(flags): print("pre-QC removed:\n", flags.groupby(["variable", "test"]).size().to_string())
    print("->", OUT)


if __name__ == "__main__":
    main()
