"""Fields and station statistics for Figs. 8 and 9 (run before fig08_09_comparison.py; env with geopandas, xarray, pyproj).
Precipitation P/F/I: same procedure as Table R12 of the reply letters (933 in-domain stations with >= 3650 days; product value
of the 1-km cell containing the station; references on their native grid, copied to the 1-km cells they contain, gaps <= 3 km
filled by nearest neighbour). Storm Alex: station totals and product values at the station cell from
figures_response_v5/extra/storm_alex_station_totals.csv (Table R13; 847 in-District stations with complete 2-5 Oct records).
Temperature: events chosen on the district-mean v2.0 series (clean boundary); stations inside the clean boundary with
>= 90 % of days 1991-2020 (mean annual), all 92 days (JJA 2003) or a value on the day; bias = product - station."""
import os, sys, glob, warnings; warnings.filterwarnings("ignore")
os.environ["PO_STATION_SET"] = "v4"; os.environ["PO_OUT_ROOT"] = os.environ.get("PO_INPUT_ROOT", "data") + "/PoGrid_v5"
sys.path.insert(0, os.environ.get("PO_INPUT_ROOT", "data") + "/PoGrid_v3_code")
from pathlib import Path
import numpy as np, pandas as pd, geopandas as gpd, xarray as xr, shapely
from scipy.spatial import cKDTree
from pyproj import Transformer
import po_common as P, s8_reference_comparison as S8

H = Path.home() / "Downloads"; V4 = H / "PoGrid_v4"; V5 = H / "PoGrid_v5"; O = V4 / "reviewer2_2026-10-20/ref_FI/out"
OUTD = V5 / "manuscript_revision/figures_v2"; (OUTD / "data").mkdir(exist_ok=True)
Gp = np.load(V4 / "grid_predictors.npz"); EXT = [Gp["x"].min()/1e3-0.5, Gp["x"].max()/1e3+0.5, Gp["y"].min()/1e3-0.5, Gp["y"].max()/1e3+0.5]
a = np.full(Gp["dem2d"].shape, np.nan); a[Gp["iy"], Gp["ix"]] = 1; MASK = np.isfinite(a)
grid = P.load_grid()
poly = gpd.read_file(H / "Po_District_clean.shp").to_crs(3035).union_all()
inC = shapely.contains_xy(poly, grid["xy"][:, 0], grid["xy"][:, 1])
CMASK = np.zeros(MASK.shape, bool); CMASK[grid["iy"], grid["ix"]] = inC


def bboxof(cfg):
    tr = Transformer.from_crs("EPSG:4326", cfg["crs"], always_xy=True)
    bx, by = tr.transform([6.0, 6.0, 14.5, 14.5], [43.5, 47.2, 43.5, 47.2]); return (min(bx), min(by), max(bx), max(by))


def to1km(A, arr):
    Ac = A.tocoo(); g = np.full(A.shape[1], np.nan); g[Ac.col] = arr[Ac.row]
    ok = np.isfinite(g); dd, j = cKDTree(grid["xy"][ok]).query(grid["xy"][~ok]); w = dd <= 3000
    idx = np.flatnonzero(~ok)[w]; g[idx] = g[ok][j[w]]
    m = np.full(MASK.shape, np.nan); m[grid["iy"], grid["ix"]] = g; return m


# ---------------- precipitation climatologies (Table R12 procedure)
R = dict(np.load(V4 / "loo_three_methods/freq_intensity_maps_cache.npz"))      # tp_* = v2.0 (two-part) P, F, I
PR = ["IH-GAR", "EMO", "ARCIS31", "RAPD", "E-OBS"]
for p in PR:
    cfg = S8.PRODUCTS["precip"][p]; tt, ra, x, y = S8.read_year(cfg, 2005, bboxof(cfg)); A, iyc, ixc, meta, info = S8.build_mapping(cfg, grid, x, y)
    c = pd.read_parquet(O / f"cells_{p}.parquet"); c = c[c.group == "all"].reset_index(drop=True); good = c.n > 3000
    for v, arr in dict(P=np.where(good, c.tot/c.n*365.25, np.nan), F=np.where(good, c.w1/c.n*100, np.nan), I=np.where(good, c.s1/c.w1, np.nan)).items():
        R[f"{p}_{v}"] = to1km(A, arr)
L = pd.concat([pd.read_parquet(f, columns=["station", "obs"]) for f in sorted(glob.glob(str(V5 / "precip_anomaly/loo_precip_*.parquet")))]).dropna()
L["station"] = L.station.astype(str)
So = L.groupby("station").obs.agg(n="size", P=lambda v: 365.25*v.mean(), F=lambda v: 100*(v >= 1).mean(), I=lambda v: v[v >= 1].mean())
st = pd.read_csv(V4 / "paper_tables/stations_precip.csv", dtype={"id": str}).set_index("id"); So = So[So.n >= 3650].join(st[["x", "y", "elevation"]], how="inner")
xs = np.sort(np.unique(Gp["xy"][:, 0])); ys = np.unique(Gp["xy"][:, 1])[::-1]


def cellidx(x, y):
    ci = np.abs(xs[None, :] - x[:, None]).argmin(1); cj = np.abs(ys[None, :] - y[:, None]).argmin(1)
    return ci, cj, MASK[cj, ci] & (np.abs(xs[ci] - x) <= 500) & (np.abs(ys[cj] - y) <= 500)


ci, cj, inside = cellidx(So.x.values, So.y.values); Sin = So[inside]; ci, cj = ci[inside], cj[inside]
KEYS = ["tp"] + PR; stats = {}
for k in KEYS:
    for v in "PFI":
        g = R[f"{k}_{v}"][cj, ci]; ok = np.isfinite(g); rel = 100*(g[ok]/Sin[v].values[ok]-1)
        stats[(k, v)] = dict(bias=np.median(rel), mae=np.median(np.abs(rel)), n=int(ok.sum()))

# ---------------- Storm Alex fields
F = {}; A0, A1 = "2020-10-02", "2020-10-05"
with xr.open_dataset(V5 / "precip_anomaly/precip_2020.nc") as ds:
    F[("tp", "Alex")] = ds["precipitation"].sel(time=slice(A0, A1)).sum("time", skipna=False).values
for p in ["IH-GAR", "ARCIS31", "EMO", "E-OBS"]:
    cfg = S8.PRODUCTS["precip"][p]; tt, ra, x, y = S8.read_year(cfg, 2020, bboxof(cfg)); A, iyc, ixc, meta, info = S8.build_mapping(cfg, grid, x, y)
    sel = (tt >= A0) & (tt <= A1); assert sel.sum() == 4; F[(p, "Alex")] = to1km(A, ra[sel][:, iyc, ixc].sum(0))
sa = pd.read_csv(V5 / "figures_response_v5/extra/storm_alex_station_totals.csv", dtype={"id": str}); SA = sa[sa["v2.0"].notna()].reset_index(drop=True)

# ---------------- temperature: events from the district-mean v2.0 series
gy, gx = grid["iy"][inC], grid["ix"][inC]; ser = []
for yr in range(1991, 2021):
    with xr.open_dataset(V5 / f"temp_anomaly/temp_{yr}.nc") as ds:
        ser.append(pd.Series(np.nanmean(ds["temperature"].values[:, gy, gx], 1), index=pd.DatetimeIndex(ds.time.values).floor("D")))
TS = pd.concat(ser); jja = TS[TS.index.month.isin([6, 7, 8])]; djf = TS[TS.index.month.isin([12, 1, 2])]
yJ = int(jja.groupby(jja.index.year).mean().idxmax()); dC = djf.idxmin().strftime("%Y-%m-%d"); dH = jja.idxmax().strftime("%Y-%m-%d")
print("warmest JJA", yJ, "coldest DJF day", dC, "hottest JJA day", dH)
EV = {"Tjja": (f"{yJ}-06-01", f"{yJ}-08-31"), "Tcold": (dC, dC), "Thot": (dH, dH)}
with xr.open_dataset(V5 / "summary/cell_climatology_temp_anomaly.nc") as ds: F[("v2", "Tann")] = ds.mean_temperature.values.astype(float)
for key, (a_, b_) in EV.items():
    with xr.open_dataset(V5 / f"temp_anomaly/temp_{a_[:4]}.nc") as ds: F[("v2", key)] = ds.temperature.sel(time=slice(a_, b_)).mean("time").values
for p in ["ARCIS_T", "EMO", "E-OBS"]:
    cfg = S8.PRODUCTS["temp"][p]; bb = bboxof(cfg); A = None
    for key, (a_, b_) in EV.items():
        tt, ra, x, y = S8.read_year(cfg, int(a_[:4]), bb)
        if A is None: A, iyc, ixc, meta, info = S8.build_mapping(cfg, grid, x, y)
        sel = (tt >= a_) & (tt <= b_); Y = ra[sel][:, iyc, ixc]
        F[(p, key)] = to1km(A, np.where(np.isfinite(Y).all(0), np.nanmean(Y, 0), np.nan))
    c = pd.read_parquet(V5 / f"reference_comparison/cells_temp_{p}.parquet")
    mm = meta[["ref_iy", "ref_ix"]].merge(c[["ref_iy", "ref_ix", "mean_ref"]], how="left", on=["ref_iy", "ref_ix"])
    F[(p, "Tann")] = to1km(A, mm.mean_ref.values)
for k in F: F[k] = np.where(CMASK, F[k], np.nan)

stT = P.load_station_data("temp"); sxy = stT["xy"]; ins = shapely.contains_xy(poly, sxy[:, 0], sxy[:, 1])
Vt = stT["values"][:, ins]; dts = stT["dates"]; frac = np.isfinite(Vt).mean(0); tx, ty = sxy[ins, 0], sxy[ins, 1]
tci, tcj, tok = cellidx(tx, ty); j = (dts >= EV["Tjja"][0]) & (dts <= EV["Tjja"][1])
STV = dict(Tann=np.where((frac >= 0.9) & tok, np.nanmean(Vt, 0), np.nan),
           Tjja=np.where(np.isfinite(Vt[j]).all(0) & tok, np.nanmean(Vt[j], 0), np.nan),
           Tcold=np.where(tok, Vt[dts == dC][0], np.nan), Thot=np.where(tok, Vt[dts == dH][0], np.nan))


def sample(fld, ci, cj):           # cell value; for cells just outside the clean mask, nearest masked cell within 1.5 cells
    v = fld[cj, ci].copy(); bad = ~np.isfinite(v)
    if bad.any():
        fy, fx = np.nonzero(np.isfinite(fld)); d, k = cKDTree(np.c_[fx, fy]).query(np.c_[ci[bad], cj[bad]])
        v[np.flatnonzero(bad)[d <= 1.5]] = fld[fy[k[d <= 1.5]], fx[k[d <= 1.5]]]
    return v


rows = []
for key in ["Tann", "Tjja", "Tcold", "Thot"]:
    s = STV[key]; m = np.isfinite(s); rows.append(dict(product="stations", column=key, n_stations=int(m.sum())))
    for p in ["v2", "ARCIS_T", "EMO", "E-OBS"]:
        g = sample(F[(p, key)], tci[m], tcj[m]); ok = np.isfinite(g); d = g[ok] - s[m][ok]
        rows.append(dict(product=p, column=key, n_stations=int(ok.sum()), median_bias=np.median(d), median_abs_error=np.median(np.abs(d)),
                         district_mean=np.nanmean(F[(p, key)])))
T9 = pd.DataFrame(rows)
rows = []
for v in "PFI":
    rows.append(dict(product="stations", column=v, n_stations=len(Sin), station_median=np.median(Sin[v])))
    for p in KEYS:
        s = stats[(p, v)]; rows.append(dict(product=p, column=v, n_stations=s["n"], median_bias_pct=s["bias"], median_abs_error_pct=s["mae"],
                                            district_mean=np.nanmean(np.where(CMASK, R[f"{p}_{v}"], np.nan))))
rows.append(dict(product="stations", column="Alex", n_stations=len(SA), station_mean=SA.obs.mean()))
for p, col in [("tp", "v2.0"), ("IH-GAR", "IH-GAR"), ("ARCIS31", "ARCIS31"), ("EMO", "EMO"), ("E-OBS", "E-OBS")]:
    mk = SA[col].notna(); o, g = SA.obs[mk], SA[col][mk]; rel = 100*(g/o-1)[o > 0]
    rows.append(dict(product=p, column="Alex", n_stations=int(mk.sum()), total_bias_pct=100*(g.mean()/o.mean()-1), mae_mm=(g-o).abs().mean(),
                     median_bias_pct=np.median(rel), median_abs_error_pct=np.median(np.abs(rel)), station_mean=o.mean(), district_mean=np.nanmean(F[(p, "Alex")])))
rows.append(dict(product="RAPD", column="Alex", n_stations=0))
T8 = pd.DataFrame(rows)

# ---------------- write
Fs = {f"{p}_{v}": np.where(CMASK, R[f"{p}_{v}"], np.nan).astype(np.float32) for p in KEYS for v in "PFI"}
Fs.update({f"{p}_{k}": a.astype(np.float32) for (p, k), a in F.items()}); Fs["CMASK"] = CMASK; Fs["EXT"] = np.array(EXT)
np.savez_compressed(OUTD / "data/fig08_09_fields.npz", **Fs)
Sin.reset_index().rename(columns={"index": "station"}).to_csv(OUTD / "data/fig08_station_climatology.csv", index=False)
SA[["id", "x", "y", "elev", "obs"]].to_csv(OUTD / "data/fig08_storm_alex_stations.csv", index=False)
pd.DataFrame(dict(id=np.array(stT["ids"])[ins], x=tx, y=ty, elev=stT["elev"][ins], frac_days=frac, **STV)).to_csv(OUTD / "data/fig09_station_values.csv", index=False)
ren = {"tp": "v2.0", "IH-GAR": "Manara et al. (2026)", "RAPD": "APGD", "ARCIS31": "ARCIS", "EMO": "EMO-1arcmin", "v2": "v2.0", "ARCIS_T": "ARCIS"}
cn8 = {"P": "mean_annual_precipitation_mm", "F": "wet_day_frequency_pct", "I": "wet_day_intensity_mm_per_wet_day", "Alex": "storm_alex_total_2_5_Oct_2020_mm"}
T8.assign(product=T8["product"].replace(ren), column=T8.column.map(cn8)).round(3).to_csv(OUTD / "stats_fig08.csv", index=False)
cn9 = {"Tann": "mean_annual_T_1991_2020", "Tjja": f"mean_T_JJA_{yJ}", "Tcold": f"daily_mean_T_{dC}", "Thot": f"daily_mean_T_{dH}"}
T9.assign(product=T9["product"].replace(ren), column=T9.column.map(cn9)).rename(columns={"median_bias": "median_bias_degC",
          "median_abs_error": "median_abs_error_degC", "district_mean": "district_mean_degC"}).round(3).to_csv(OUTD / "stats_fig09.csv", index=False)
print(T8.round(2).to_string()); print(T9.round(2).to_string())
