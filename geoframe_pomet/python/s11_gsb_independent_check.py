"""
s11_gsb_independent_check.py -- point check against MeteoSwiss Col du Grand St-Bernard (GSB, 2472 m, 45.8691N 7.1707E).

GSB is NOT in the precipitation network (nearest gauge 3056000, 2084 m, 4.4 km), so the precipitation check is
independent. For temperature the Italian-side station 3065 (2435 m, 1.3 km) is in the network, so that check is
only semi-independent (different instrument/owner, same site).

MeteoSwiss Open Data (ch.meteoschweiz.ogd-smn), daily:
  rka150d0  precipitation 00-00 UTC   (same convention as our harmonised calendar day)
  rre150d0  precipitation 06-06 UTC, labelled with the start day (used only for the lag check)
  tre200d0  24-h mean T 2 m;  tre200dx / tre200dn  daily max / min
Also gives the station estimate of (Tx+Tn)/2 - T24 used to interpret the ARCIS (Tmax+Tmin)/2 bias.

Extracts the nearest 1 km cell (and the 3x3 block mean) from PoGrid v3, the old released product and the
reference datasets (via s8 readers), 1991-2020. Elevation differences cell-station are reported, temperature is
also shown lapse-adjusted to the station height with -6.5 K/km (illustrative only).

Outputs (OUT_ROOT/gsb_check/):  gsb_daily.parquet, gsb_metrics.csv, gsb_monthly_offset_txtn.csv, gsb_lags.csv
"""
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd, xarray as xr
import po_config as C
import s8_reference_comparison as S8

GSB = dict(lat=45.869092, lon=7.170683, z=2472.0)
DL = Path(os.environ.get("PO_INPUT_ROOT", "data") + ""); MS = DL / "MeteoSwiss_GSB"
Y0, Y1 = 1991, 2020
OD = C.OUT_ROOT / "gsb_check"; OD.mkdir(parents=True, exist_ok=True)


def load_station():
    d = pd.concat([pd.read_csv(MS / f, sep=";", na_values=["", "-"]) for f in ("ogd-smn_gsb_d_historical.csv", "ogd-smn_gsb_d_recent.csv")])
    d["date"] = pd.to_datetime(d["reference_timestamp"], format="%d.%m.%Y %H:%M").dt.floor("D")
    d = d.drop_duplicates("date").set_index("date").sort_index()
    return d.loc[f"{Y0}-01-01":f"{Y1}-12-31", ["rka150d0", "rre150d0", "tre200d0", "tre200dx", "tre200dn"]]


def grid_cell():
    from pyproj import Transformer
    g = np.load(C.OUT_ROOT / "grid_predictors.npz")
    X, Y = Transformer.from_crs(4326, 3035, always_xy=True).transform(GSB["lon"], GSB["lat"])
    d = np.hypot(g["xy"][:, 0] - X, g["xy"][:, 1] - Y); k = int(np.argmin(d))
    iy, ix = int(g["iy"][k]), int(g["ix"][k])
    blk = [(j, i) for j in range(iy - 1, iy + 2) for i in range(ix - 1, ix + 2) if np.isfinite(g["dem2d"][j, i])]
    return dict(iy=iy, ix=ix, dist_km=d[k] / 1e3, elev=float(g["elev"][k]), blk=blk,
                blk_elev=float(np.nanmean([g["dem2d"][j, i] for j, i in blk])), X=X, Y=Y)


def read_own(pattern, var, years, cell):
    out = []
    for y in years:
        f = Path(pattern.format(y=y))
        if not f.exists(): continue
        with xr.open_dataset(f) as ds:
            a = ds[var].isel(y=slice(cell["iy"] - 1, cell["iy"] + 2), x=slice(cell["ix"] - 1, cell["ix"] + 2)).values.astype(float)
            t = pd.DatetimeIndex(ds["time"].values).floor("D")
        c = a[:, 1, 1]; m = np.array([np.nanmean(np.where(np.isfinite(a[i]) & (a[i] > -900), a[i], np.nan)) for i in range(len(t))])
        out.append(pd.DataFrame({"cell": c, "blk3": m}, index=t))
    return pd.concat(out) if out else None


def read_ref(cfg, years):
    """nearest reference cell to GSB in the reference grid's own CRS (no aggregation)."""
    from pyproj import Transformer
    X, Y = Transformer.from_crs(4326, cfg["crs"], always_xy=True).transform(GSB["lon"], GSB["lat"])
    pad = 0.3 if cfg["crs"] == "EPSG:4326" else 30000.0
    out = []
    for y in years:
        r = S8.read_year(cfg, y, (X - pad, Y - pad, X + pad, Y + pad))
        if r is None: continue
        t, a, x, yy = r; i = int(np.argmin(np.abs(x - X))); j = int(np.argmin(np.abs(yy - Y)))
        out.append(pd.Series(a[:, j, i].astype(float), index=t))
    return pd.concat(out) if out else None


def pmetrics(o, m):
    ok = o.notna() & m.notna(); o, m = o[ok], m[ok]
    yrs = ok.groupby(ok.index.year).sum(); full = yrs[yrs >= 330].index
    oy = o[o.index.year.isin(full)].groupby(o.index.year[o.index.year.isin(full)]).sum()
    my = m[m.index.year.isin(full)].groupby(m.index.year[m.index.year.isin(full)]).sum()
    return dict(n_days=int(ok.sum()), n_full_years=len(full), obs_annual_mm=oy.mean(), mod_annual_mm=my.mean(),
                bias_pct=100 * (my.mean() / oy.mean() - 1), r_daily=np.corrcoef(o, m)[0, 1],
                mae=float(np.abs(o - m).mean()), wetfreq_obs=float((o >= 1).mean()), wetfreq_mod=float((m >= 1).mean()),
                p99_obs=float(o.quantile(.99)), p99_mod=float(m.quantile(.99)),
                rx1_obs=float(o.groupby(o.index.year).max().mean()), rx1_mod=float(m.groupby(m.index.year).max().mean()))


def tmetrics(o, m):
    ok = o.notna() & m.notna(); o, m = o[ok], m[ok]; e = m - o
    return dict(n_days=int(ok.sum()), obs_mean=o.mean(), mod_mean=m.mean(), bias=e.mean(), mae=float(e.abs().mean()),
                rmse=float(np.sqrt((e ** 2).mean())), r_daily=np.corrcoef(o, m)[0, 1],
                r_dtd=np.corrcoef(o.diff().dropna(), m.diff().reindex(o.diff().dropna().index).fillna(0))[0, 1],
                bias_DJF=e[e.index.month.isin([12, 1, 2])].mean(), bias_JJA=e[e.index.month.isin([6, 7, 8])].mean())


def lags(o, m, name, var):
    rows = []
    for L in (-1, 0, 1):
        mm = m.shift(L, freq="D").reindex(o.index); ok = o.notna() & mm.notna()
        rows.append(dict(variable=var, product=name, lag_days=L, r=np.corrcoef(o[ok], mm[ok])[0, 1], n=int(ok.sum())))
    return rows


def main():
    years = range(Y0, Y1 + 1); st = load_station(); cell = grid_cell()
    print(f"GSB cell iy={cell['iy']} ix={cell['ix']} dist={cell['dist_km']:.2f} km  cell z={cell['elev']:.0f} m  3x3 z={cell['blk_elev']:.0f} m  station z={GSB['z']:.0f} m")
    P = {C.OUT_ROOT.name: read_own(str(C.OUT_ROOT / "precip_anomaly" / "precip_{y}.nc"), "precipitation", years, cell),
         "GEOframe_old": read_own(str(DL / "precipitation_data" / "precipitation_{y}.nc"), "precipitation", years, cell)}
    T = {C.OUT_ROOT.name: read_own(str(C.OUT_ROOT / "temp_anomaly" / "temp_{y}.nc"), "temperature", years, cell),
         "GEOframe_old": read_own(str(DL / "temperature_data" / "temperature_{y}.nc"), "temperature", years, cell)}
    for nm, cfg in S8.PRODUCTS["precip"].items():
        s = read_ref(cfg, years); print(nm, None if s is None else len(s))
        if s is not None: P[nm] = pd.DataFrame({"cell": s})
    for nm, cfg in S8.PRODUCTS["temp"].items():
        s = read_ref(cfg, years); print(nm, None if s is None else len(s))
        if s is not None: T[nm] = pd.DataFrame({"cell": s})

    rows, lagrows, daily = [], [], {"obs_rka150d0": st["rka150d0"], "obs_rre150d0": st["rre150d0"], "obs_tre200d0": st["tre200d0"],
                                    "obs_tre200dx": st["tre200dx"], "obs_tre200dn": st["tre200dn"]}
    for nm, df in P.items():
        for col in df.columns:
            m = df[col].reindex(st.index); daily[f"P_{nm}_{col}"] = m
            rows.append(dict(variable="precip", product=nm, sample=col, **pmetrics(st["rka150d0"], m)))
        lagrows += lags(st["rka150d0"], df["cell"], nm, "precip_rka150d0(00-00)")
        lagrows += lags(st["rre150d0"], df["cell"], nm, "precip_rre150d0(06-06,start-day label)")
    dz = cell["elev"] - GSB["z"]
    for nm, df in T.items():
        for col in df.columns:
            m = df[col].reindex(st.index); daily[f"T_{nm}_{col}"] = m
            rows.append(dict(variable="temp", product=nm, sample=col, **tmetrics(st["tre200d0"], m)))
            rows.append(dict(variable="temp_vs_TxTn2", product=nm, sample=col, **tmetrics((st["tre200dx"] + st["tre200dn"]) / 2, m)))
        lagrows += lags(st["tre200d0"], df["cell"], nm, "temp_tre200d0")
    R = pd.DataFrame(rows); R.insert(3, "cell_minus_station_elev_m", np.where(R["product"].isin([C.OUT_ROOT.name, "GEOframe_old"]), dz, np.nan))
    R.to_csv(OD / "gsb_metrics.csv", index=False); pd.DataFrame(lagrows).to_csv(OD / "gsb_lags.csv", index=False)
    pd.DataFrame(daily).to_parquet(OD / "gsb_daily.parquet")
    off = ((st["tre200dx"] + st["tre200dn"]) / 2 - st["tre200d0"]).dropna()
    O = off.groupby(off.index.month).agg(["mean", "std", "count"]); O.loc["annual"] = [off.mean(), off.std(), off.count()]
    O.to_csv(OD / "gsb_monthly_offset_txtn.csv")
    pd.set_option("display.width", 250)
    print(R.round(3).to_string(index=False)); print(pd.DataFrame(lagrows).pivot_table(index=["variable", "product"], columns="lag_days", values="r").round(3).to_string())
    print("(Tx+Tn)/2 - T24 at GSB:\n", O.round(2).to_string()); print("->", OD)


if __name__ == "__main__":
    main()
