"""
s13_zenodo_package.py -- build the Zenodo v2.0 package of the Po River District 1-km daily dataset.

  <OUT>/daily/precipitation_YYYY.nc, temperature_YYYY.nc   (same file and variable names as v1.0)
        precipitation (mm day-1) / temperature (degC) and <var>_uncertainty: calibrated 1-sigma
        = ordinary-kriging standard deviation x k(season, class), k from LOO (s9_uncertainty_calibration.py;
        precipitation classes by estimated daily amount, temperature classes by cell elevation), so that
        +-1.96 sigma covers 95 % of the leave-one-out errors. Packed int16 (0.1 mm / 0.01 degC), zlib.
  <OUT>/normals/normals_precipitation_1991-2020.nc, normals_temperature_1991-2020.nc   (monthly + annual)
  <OUT>/effective_resolution/effective_resolution_{precipitation,temperature}.nc
        per-cell median over days of R_pos_2_3, R_pos_half, R_pos_1_3 (farthest station whose kriging weight is
        >= 2/3, 1/2, 1/3 of the largest weight), D1 and D7 (distance to nearest / 7th station), per decade and 1991-2020
  <OUT>/tables/*.csv     aggregate validation tables (no station data)
  <OUT>/README.md, CHANGELOG.md, checksums.md5
Station observations and station lists are NOT included (data-provider licences).

Usage:  PO_OUT_ROOT=~/Downloads/PoGrid_v4 python s13_zenodo_package.py [--out ~/Downloads/PoGrid_v2.0_zenodo] [--years 1991-2020]
"""
import argparse, hashlib, shutil, sys, time
from datetime import datetime, timezone
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd, xarray as xr
import po_config as C

VERSION = "2.0"
CITED_DOI = "10.5281/zenodo.19207256"   # DOI cited for v1.0 in the submitted manuscript; override with --doi
TITLE = "A 30-year 1-km daily precipitation and air temperature dataset for the Po River District (Italy), 1991-2020"
PERIODS = {"1991-2000": (1991, 2000), "2001-2010": (2001, 2010), "2011-2020": (2011, 2020), "1991-2020": (1991, 2020)}
RAD = ["R_pos_2_3", "R_pos_half", "R_pos_1_3", "D1", "D7"]
SEASON = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}
NAMES = {"precip": ("precipitation", "mm day-1", "lwe_thickness_of_precipitation_amount", 0.1, "daily total precipitation, 00-24 UTC"),
         "temp": ("temperature", "degC", "air_temperature", 0.01, "daily mean 2-m air temperature, 00-24 UTC")}


def crs_var():
    from pyproj import CRS
    c = CRS.from_epsg(3035); a = c.to_cf(); a["spatial_ref"] = c.to_wkt(); a["crs_wkt"] = c.to_wkt(); a["epsg_code"] = "EPSG:3035"
    return xr.DataArray(np.int32(0), attrs=a)


def global_attrs(var_long, extra=None):
    a = dict(title=TITLE, product_version=VERSION, Conventions="CF-1.8", institution="University of Trento; Po River Basin District Authority (ADBPO)",
             source="daily station observations (regional ARPA networks, EEAR-Clim, MeteoSwiss Open Data), quality controlled and harmonised "
                    "to the 00-24 UTC calendar day; climatology-aided (anomaly) ordinary kriging on 1991-2020 PRISM-type monthly normals; "
                    "local spherical variograms (50 nearest stations), 7 kriging neighbours, local daily lapse rate",
             variable=var_long, crs="EPSG:3035 (ETRS89-LAEA), 1 km", time_coverage="1991-01-01/2020-12-31",
             doi_v1_record=CITED_DOI, license="CC BY 4.0", references="Salehi et al., Earth Syst. Sci. Data (in review)",
             date_created=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
             comment="Effective daily resolution is coarser than the 1-km grid (see effective_resolution/ and README).")
    a.update(extra or {}); return a


def k_lookup(v):
    u = pd.read_csv(C.OUT_ROOT / "summary" / f"uncertainty_calibration_{v}.csv")
    u = u[u.season != "all"]; return {(r.season, str(r.cls)): float(r.scale_k) for r in u.itertuples()}


def calibrated_std(v, std, est, months, elev2d, K):
    out = np.full_like(std, np.nan, dtype=np.float32)
    for t, m in enumerate(months):
        s = SEASON[int(m)]; e = est[t]; sd = std[t]; k = np.full(e.shape, np.nan, np.float32)
        if v == "precip":
            edges = [(-np.inf, 0.0, "0"), (0.0, 1, "0-1"), (1, 5, "1-5"), (5, 20, "5-20"), (20, 50, "20-50"), (50, np.inf, ">50")]
            for lo, hi, c in edges: k[(e > lo) & (e <= hi)] = K[(s, c)]
        else:
            edges = [(-np.inf, 500, "<500"), (500, 1000, "500-1000"), (1000, 1500, "1000-1500"), (1500, 2000, "1500-2000"), (2000, np.inf, ">=2000")]
            for lo, hi, c in edges: k[(elev2d > lo) & (elev2d <= hi)] = K[(s, c)]
        out[t] = sd * k
    return out


def enc(scale, name):
    return {name: dict(dtype="int16", scale_factor=scale, add_offset=0.0, _FillValue=np.int16(-32768), zlib=True, complevel=5, shuffle=True)}


def write_daily(v, years, od, grid):
    name, units, std_name, scale, long = NAMES[v]; K = k_lookup(v)
    src = C.OUT_ROOT / f"{v}_anomaly"; elev2d = grid["dem2d"]
    lat2d = np.full(elev2d.shape, np.nan, np.float32); lon2d = lat2d.copy()
    lat2d[grid["iy"], grid["ix"]] = grid["lat"]; lon2d[grid["iy"], grid["ix"]] = grid["lon"]
    for y in years:
        t0 = time.time(); f = src / f"{v}_{y}.nc"
        with xr.open_dataset(f) as ds:
            est = ds[NAMES[v][0]].values.astype(np.float32); std = ds["kriging_std"].values.astype(np.float32)
            t = ds["time"].values; x = ds["x"].values; yy = ds["y"].values
        cal = calibrated_std(v, std, est, pd.DatetimeIndex(t).month, elev2d, K)
        o = xr.Dataset({name: (("time", "y", "x"), est, dict(units=units, standard_name=std_name, long_name=long, grid_mapping="crs", coordinates="lat lon")),
                        f"{name}_uncertainty": (("time", "y", "x"), cal, dict(units=units, long_name=f"calibrated 1-sigma uncertainty of {name}",
                                                comment="ordinary-kriging standard deviation scaled by a leave-one-out factor k(season, class); "
                                                        "+-1.96 sigma covers 95 % of leave-one-out errors (cross-validated)", grid_mapping="crs", coordinates="lat lon")),
                        "crs": crs_var()},
                       coords=dict(time=("time", t, dict(standard_name="time", long_name="calendar day, 00-24 UTC")),
                                   y=("y", yy, dict(standard_name="projection_y_coordinate", units="m", axis="Y")),
                                   x=("x", x, dict(standard_name="projection_x_coordinate", units="m", axis="X")),
                                   lat=(("y", "x"), lat2d, dict(standard_name="latitude", units="degrees_north")),
                                   lon=(("y", "x"), lon2d, dict(standard_name="longitude", units="degrees_east"))),
                       attrs=global_attrs(long, dict(year=int(y))))
        e = {**enc(scale, name), **enc(scale, f"{name}_uncertainty"), "time": dict(units="days since 1991-01-01", calendar="standard", dtype="int32"),
             "lat": dict(dtype="float32", zlib=True), "lon": dict(dtype="float32", zlib=True)}
        out = od / "daily" / f"{name}_{y}.nc"; o.to_netcdf(out, encoding=e)
        print(f"  {out.name}  ({out.stat().st_size / 1e6:.0f} MB, {time.time() - t0:.0f} s)", flush=True)


def write_normals(od, grid):
    for v in ("precip", "temp"):
        name, _, std_name, _, _ = NAMES[v]
        n = xr.open_dataset(C.OUT_ROOT / "normals" / f"normals_{v}.nc")
        mon = n["normal"].astype("float32"); ann = mon.sum("month") if v == "precip" else mon.mean("month")
        u_m, u_a = ("mm month-1", "mm yr-1") if v == "precip" else ("degC", "degC")
        o = xr.Dataset({f"{name}_monthly": (("month", "y", "x"), mon.values, dict(units=u_m, long_name=f"1991-2020 monthly normal of {name}", grid_mapping="crs")),
                        f"{name}_annual": (("y", "x"), ann.values, dict(units=u_a, long_name=f"1991-2020 annual normal of {name}", grid_mapping="crs")),
                        "crs": crs_var()},
                       coords=dict(month=("month", np.arange(1, 13)), y=("y", n["y"].values, dict(standard_name="projection_y_coordinate", units="m")),
                                   x=("x", n["x"].values, dict(standard_name="projection_x_coordinate", units="m"))),
                       attrs=global_attrs(f"1991-2020 normals of {name}",
                                          dict(method="PRISM-type local weighted regression on station normals (ln P or T vs elevation, distance/elevation/"
                                                      "facet weights); no vertical extrapolation of precipitation above the highest local station; "
                                                      "temperature above the highest local station with a LOO-tuned lapse rate",
                                               parameters=n.attrs.get("parameters", ""))))
        out = od / "normals" / f"normals_{name}_1991-2020.nc"
        o.to_netcdf(out, encoding={k: dict(zlib=True, complevel=5, dtype="float32") for k in (f"{name}_monthly", f"{name}_annual")}); print("  ", out.name)


def write_effres(od, grid):
    for v in ("precip", "temp"):
        name = NAMES[v][0]; src = C.OUT_ROOT / f"{v}_anomaly"; data = {}
        import warnings
        for p, (a, b) in PERIODS.items():
            for r in RAD:                                  # one variable at a time: peak memory ~2 x (days x cells x 4 B)
                st_ = np.empty((0,) + grid["dem2d"].shape, np.float32)
                parts = []
                for y in range(a, b + 1):
                    with xr.open_dataset(src / f"diagnostics_{v}_{y}.nc") as d: parts.append(d[r].values.astype(np.float32))
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN cells outside the domain
                    data[(p, r)] = np.nanmedian(np.concatenate(parts, 0), axis=0).astype(np.float32)
                del parts
            print(f"   {name} {p}", flush=True)
        with xr.open_dataset(src / f"diagnostics_{v}_2020.nc") as d:
            x, y = d["x"].values, d["y"].values; lab = {r: d[r].attrs.get("long_name", r) for r in RAD}
        per = list(PERIODS)
        o = xr.Dataset({r: (("period", "y", "x"), np.stack([data[(p, r)] for p in per]), dict(units="km", long_name=f"median over days: {lab[r]}", grid_mapping="crs"))
                        for r in RAD} | {"crs": crs_var()},
                       coords=dict(period=("period", per), y=("y", y, dict(standard_name="projection_y_coordinate", units="m")),
                                   x=("x", x, dict(standard_name="projection_x_coordinate", units="m"))),
                       attrs=global_attrs(f"effective resolution of daily {name}",
                                          dict(comment="Distances describe how far the information behind each 1-km daily value comes from; "
                                                       "they are a measure of effective resolution, not of error.")))
        out = od / "effective_resolution" / f"effective_resolution_{name}.nc"
        o.to_netcdf(out, encoding={r: dict(zlib=True, complevel=5, dtype="float32") for r in RAD}); print("  ", out.name)


def copy_tables(od):
    s = C.OUT_ROOT; t = od / "tables"
    files = {"summary/loo_metrics_precip_anomaly.csv": "loo_metrics_precipitation.csv", "summary/loo_metrics_temp_anomaly.csv": "loo_metrics_temperature.csv",
             "summary/uncertainty_calibration_precip.csv": "uncertainty_k_precipitation.csv", "summary/uncertainty_calibration_temp.csv": "uncertainty_k_temperature.csv",
             "summary/support_radii_by_region_period_precip_anomaly.csv": "effective_resolution_by_region_period_precipitation.csv",
             "summary/support_radii_by_region_period_temp_anomaly.csv": "effective_resolution_by_region_period_temperature.csv",
             "summary/loo_stationwise_mean_like_paper.csv": "loo_stationwise_means.csv",
             "reference_comparison/summary_precip.csv": "reference_comparison_precipitation.csv",
             "reference_comparison/summary_temp.csv": "reference_comparison_temperature.csv",
             "isac_comparison/summary.csv": "normals_vs_isac_cnr_1961-1990.csv", "isac_comparison/gradients.csv": "normals_elevation_gradients_vs_isac_cnr.csv"}
    for a, b in files.items():
        if (s / a).exists(): shutil.copy(s / a, t / b)
        else: print("   missing (skipped):", a)


def readme(od):
    (od / "README.md").write_text(f"""# {TITLE} - version {VERSION}

Zenodo record of version 1.0 (as cited in the submitted manuscript): https://doi.org/{CITED_DOI}. Licence: CC BY 4.0.

## Contents
| Folder | Content |
|---|---|
| `daily/` | `precipitation_YYYY.nc` (mm day-1) and `temperature_YYYY.nc` (degC), 1991-2020, daily, 1 km, EPSG:3035; each with `<variable>_uncertainty` (calibrated 1-sigma) |
| `normals/` | 1991-2020 monthly and annual normals (the climatological layer of the interpolation) |
| `effective_resolution/` | per-cell median distances describing the effective resolution, per decade and 1991-2020 |
| `tables/` | aggregate leave-one-out, uncertainty-calibration, effective-resolution and inter-comparison tables |

Day definition: calendar day 00-24 UTC (station records were harmonised to this convention).
Grid: ETRS89-LAEA (EPSG:3035), 1 km, 494 x 330 cells, 88,074 valid cells; 2-D `lat`/`lon` are included.
Values are packed as int16 (precipitation 0.1 mm, temperature 0.01 degC); xarray/netCDF4 unpack them automatically.

## Uncertainty
`<variable>_uncertainty` is the ordinary-kriging standard deviation multiplied by a factor k estimated from
leave-one-out errors per season and class (precipitation: estimated daily amount; temperature: elevation), so that
the interval value +- 1.96 x uncertainty contains 95 % of leave-one-out errors (cross-validated on independent
station halves). Factors are listed in `tables/uncertainty_k_*.csv`. For precipitation the error distribution is
skewed; the interval is a calibrated spread, not a symmetric Gaussian error model.

## Limitations
* Effective daily resolution is about 10-20 km (median distance of the stations carrying 1/2 of the kriging
  weight: ~11 km for precipitation and ~13 km for temperature in 2011-2020, larger in the 1990s); 1-km detail
  comes mainly from the elevation-dependent normals.
* Kriging smooths daily extremes; heavy precipitation days are underestimated.
* Precipitation gauges are not corrected for wind-induced undercatch.
* Above the highest local station, precipitation normals are not extrapolated upward; isolated pass-top or summit
  gauges with exceptional totals (e.g. Gran San Bernardo) are underestimated when withheld.
* Station observations are not redistributed (data-provider licences).

## Citation
Salehi, H., Andreis, D., Baig, S., Roati, G., Brian, M., Tornatore, F., Formetta, G. and Rigon, R.: A 30-year 1-km
daily precipitation and air temperature dataset for the Po River District (Italy), Earth Syst. Sci. Data (in review).
""")
    (od / "CHANGELOG.md").write_text("""# Changelog

## 2.0
* Station day labelling harmonised to the 00-24 UTC calendar day (many records were labelled with the end day).
* Iterative leave-one-out-based quality control of station data.
* Network extended with official MeteoSwiss SwissMetNet and manual precipitation stations within 20 km of the district.
* New interpolation: climatology-aided (anomaly) ordinary kriging on 1991-2020 PRISM-type monthly normals, local
  variograms and local daily lapse rates (v1.0: domain-wide kriging with one daily lapse rate).
* New: calibrated uncertainty layer, monthly normals, effective-resolution maps, aggregate validation tables.
* File and variable names unchanged (`precipitation_YYYY.nc`, `temperature_YYYY.nc`).

## 1.0
* Initial release (submitted manuscript).
""")


def checksums(od):
    lines = []
    for f in sorted(p for p in od.rglob("*") if p.is_file() and p.name != "checksums.md5"):
        h = hashlib.md5()
        with open(f, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 20), b""): h.update(b)
        lines.append(f"{h.hexdigest()}  {f.relative_to(od)}")
    (od / "checksums.md5").write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default=str(Path.home() / "Downloads" / f"PoGrid_v{VERSION}_zenodo"))
    ap.add_argument("--years", default="1991-2020"); ap.add_argument("--skip-effres", action="store_true")
    ap.add_argument("--doi", default=None, help="DOI to write into metadata/README (default: the v1.0 DOI cited in the manuscript)")
    a = ap.parse_args(); od = Path(a.out)
    global CITED_DOI
    if a.doi: CITED_DOI = a.doi
    y0, y1 = (int(s) for s in a.years.split("-")); years = range(y0, y1 + 1)
    for sub in ("daily", "normals", "effective_resolution", "tables"): (od / sub).mkdir(parents=True, exist_ok=True)
    grid = np.load(C.OUT_ROOT / "grid_predictors.npz")
    print("daily files"); [write_daily(v, years, od, grid) for v in ("precip", "temp")]
    print("normals"); write_normals(od, grid)
    if not a.skip_effres: print("effective resolution (reads all daily diagnostics, ~10-20 min)"); write_effres(od, grid)
    copy_tables(od); readme(od); print("checksums"); checksums(od)
    tot = sum(p.stat().st_size for p in od.rglob("*") if p.is_file()) / 1e9; print(f"-> {od}  ({tot:.2f} GB)")


if __name__ == "__main__":
    main()
