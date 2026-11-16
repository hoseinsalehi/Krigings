#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 3 - daily 1 km fields (direct or anomaly method) with built-in leave-one-out.

Per day (one worker per day):
  1. active stations; value to interpolate at every station
       direct : precipitation as observed; temperature as observed (per-cell lapse removed below)
       anomaly: precipitation ratio P / Pclim_day; temperature difference T - Tclim
     (station climatology = step-1 normal; where a station has no accepted normal, the step-2
      grid normal predicted WITHOUT that station is used and the station is flagged)
  2. for every target (grid cells, optional points, and every active station for LOO, which
     never sees its own value or its own normal):
       50 nearest active stations -> [temperature: local daily lapse, 100-station fallback below
       300 m elevation span] -> empirical variogram -> PSO spherical (24 x 35) -> ordinary kriging
       from the 7 nearest -> estimate, kriging variance, weight diagnostics
     precipitation: estimate 0 when all 7 neighbours are dry; estimates bounded at 0
  3. back-transform: anomaly x / + target climatology (grid normal; LOO normal for stations)
Per year one set of files is written and nothing else is kept on disk:
  <run>/<var>_<Y>.nc          value + kriging_std              (int16, zlib, 1-day chunks)
  <run>/diagnostics_<var>_<Y>.nc  D1/D7/D50, variogram range/nugget fraction/log10 sill,
                                   lapse slope/R2, negative weights, N_eff, R_pos1/2, R_abs90, qc flags
  <run>/loo_<var>_<Y>.parquet  every active station-day predicted without itself
  <run>/daily_summary_<var>_<Y>.csv
  <run>/points_<var>_<Y>.csv  (if --points)
Years already complete are skipped (resume).

Examples
  python s3_run.py --variable precip --mode anomaly --workers 12
  python s3_run.py --variable temp --mode direct --loo-only --dates-csv selected_days.csv --tag loo_direct
"""
import argparse, json, sys, time, hashlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import multiprocessing as mp
import numpy as np, pandas as pd
import po_config as C, po_common as P, po_core as K

G = {}
QC_BITS = {"zero_shortcut": 1, "lapse_fallback": 2, "positive_lapse": 4, "oscillatory_weights": 8,
           "range_bound_hit": 16, "global_variogram_fallback": 32, "elevation_filled": 64, "station_normal_fallback": 128}


# ----------------------------------------------------------------------------- encoding
ENC = {  # name: (dtype, scale, fill, units, long_name)
    "D1": ("u2", 0.1, 65535, "km", "distance to nearest active station"),
    "D7": ("u2", 0.1, 65535, "km", "distance to 7th nearest active station (kriging neighbourhood radius)"),
    "D50": ("u2", 0.5, 65535, "km", "distance to 50th nearest active station (variogram neighbourhood radius)"),
    "range": ("u2", 0.5, 65535, "km", "spherical variogram range"),
    "nugget_fraction": ("u1", 1 / 100, 255, "1", "nugget / (nugget + partial sill)"),
    "log10_sill": ("i2", 1e-2, -32768, "log10(units^2)", "log10 of total sill"),
    "lapse": ("i2", 0.01, -32768, "degC km-1", "local daily temperature-elevation slope"),
    "lapse_r2": ("u1", 1 / 100, 255, "1", "R2 of the local daily lapse regression"),
    "n_negative_weights": ("u1", 1, 255, "1", "number of negative kriging weights"),
    "n_eff": ("u1", 1 / 20, 255, "1", "effective number of stations (inverse Herfindahl of |w|)"),
    "R_pos_half": ("u2", 0.5, 65535, "km", "farthest station with weight >= 1/2 of the largest positive weight"),
    "R_pos_2_3": ("u2", 0.5, 65535, "km", "farthest station with weight >= 2/3 of the largest positive weight"),
    "R_pos_1_3": ("u2", 0.5, 65535, "km", "farthest station with weight >= 1/3 of the largest positive weight"),
    "R_abs90": ("u2", 0.5, 65535, "km", "distance holding 90% of the total absolute weight"),
    "qc_flags": ("u1", 1, 0, "1", "bit flags: " + ", ".join(f"{v}={k}" for k, v in QC_BITS.items())),
    # two-part precipitation: variogram and support of the AMOUNT kriging (wet stations only); the unsuffixed
    # variogram/support fields above then refer to the OCCURRENCE kriging (all active stations)
    "range_amount": ("u2", 0.5, 65535, "km", "spherical variogram range, wet-day intensity ratios"),
    "nugget_fraction_amount": ("u1", 1 / 100, 255, "1", "nugget / (nugget + partial sill), wet-day intensity ratios"),
    "log10_sill_amount": ("i2", 1e-2, -32768, "log10(1)", "log10 of total sill, wet-day intensity ratios"),
    "R_pos_half_amount": ("u2", 0.5, 65535, "km", "farthest wet station with weight >= 1/2 of the largest positive weight (amount kriging)"),
    "wet_probability": ("u1", 1 / 250, 255, "1", "occurrence index pi = p(x) + kriged occurrence anomaly (before calibration)"),
}
TP = {  # two-part pass-1 fields kept until the occurrence calibration (finalize_twopart)
    "pi": ("u2", 1e-4, 65535), "amount": ("i2", 0.1, -32768), "samount": ("i2", 0.1, -32768)}


def enc(a, name):
    dt, sc, fill, *_ = ENC[name]; info = np.iinfo(np.dtype(dt))
    q = np.rint(np.asarray(a, float) / sc)
    q = np.where(np.isfinite(q), np.clip(q, info.min + (dt[0] == "i"), info.max - (dt[0] == "u")), fill)
    return q.astype(dt)


def value_enc(variable):
    return ("i2", 0.1 if variable == "precip" else 0.01, -32768)


# ----------------------------------------------------------------------------- per-day work
def transform(day_idx):
    """Value to interpolate at every station for this day + climatology at targets."""
    v = G["variable"]; mode = G["mode"]; vals = G["values"][day_idx].astype(float)
    if mode == "direct": return vals, None
    m = G["month"][day_idx] - 1; ndays = G["dim"][day_idx]
    if v == "precip":
        cd = np.maximum(G["NS"][:, m] / ndays, C.MIN_DAILY_NORMAL_MM)
        return np.minimum(vals / cd, C.MAX_PRECIP_RATIO), (m, ndays)
    return vals - G["NS"][:, m], (m, ndays)


def back(est, kvar, clim_target, m, ndays):
    v = G["variable"]; mode = G["mode"]; ks = np.sqrt(kvar)
    if mode == "direct":
        return (np.maximum(est, 0.0) if v == "precip" else est), ks
    if v == "precip":
        cd = np.maximum(clim_target / ndays, C.MIN_DAILY_NORMAL_MM); return np.maximum(est, 0.0) * cd, ks * cd
    return est + clim_target, ks


def twopart_day(day_idx):
    """Two-part precipitation for one day (pass 1): occurrence index pi and conditional amount at every target.
    Occurrence: anomalies o - p_station (o = 1 if P >= OCC_THRESHOLD) of ALL active stations are kriged (local variogram,
    50/7 stations as for the ratio method) and added to the target's frequency normal p(x): pi = p(x) + o_hat.
    Amount: intensity ratios r = P / i_station of the WET stations only (i = M / (p * days)) are kriged with their own
    local variogram and multiplied by the target's intensity normal i(x). The wet/dry decision (pi above a
    threshold that reproduces p(x) over the calibration period) is taken in finalize_twopart."""
    t0 = time.time(); v = G["variable"]; date = G["dates"][day_idx]; tag = f"{v}|twopart|{date:%Y-%m-%d}"
    m = G["month"][day_idx] - 1; nd = G["dim"][day_idx]
    vals = G["values"][day_idx].astype(float); fin = np.isfinite(vals)
    act = np.flatnonzero(fin & G["usable"]); res = dict(day_idx=day_idx, date=f"{date:%Y-%m-%d}", n_active=len(act))
    if len(act) < C.N_VARIO + 1: res["skipped"] = True; return res
    ps = G["PSF"][:, m]; i_st = G["NS"][:, m] / (np.maximum(ps, C.P_FLOOR) * nd)
    wet = fin & (vals >= C.OCC_THRESHOLD)
    y_occ = np.where(fin, wet.astype(float) - ps, np.nan)
    y_amt = np.where(wet, np.minimum(vals / np.maximum(i_st, C.MIN_DAILY_NORMAL_MM), C.MAX_PRECIP_RATIO), np.nan)
    act_w = act[wet[act] & np.isfinite(y_amt[act])]
    res.update(n_wet=len(act_w), wet_share=len(act_w) / len(act))
    g_occ = K.global_variogram(G["st_xy"], G["sdmat"], act, y_occ, tag + "|occ")
    g_amt = K.global_variogram(G["st_xy"], G["sdmat"], act_w, y_amt, tag + "|amt") if len(act_w) >= 10 else None
    def two(txy, tz, sidx, p_t, i_t):
        o1 = K.estimate_targets(txy, tz, sidx, act=act, st_xy=G["st_xy"], st_elev=G["st_elev"], sdmat=G["sdmat"], y=y_occ,
                                variable="occ", lapse_mode="none", day_tag=tag + "|occ", global_params=g_occ)
        pi = np.clip(p_t + o1["est"], 0.0, 1.0)
        if len(act_w) >= C.N_KRIG + 1:
            o2 = K.estimate_targets(txy, tz, sidx, act=act_w, st_xy=G["st_xy"], st_elev=G["st_elev"], sdmat=G["sdmat"], y=y_amt,
                                    variable="amt", lapse_mode="none", day_tag=tag + "|amt", global_params=g_amt)
            amt = np.maximum(o2["est"], 0.0) * i_t; samt = np.sqrt(o2["kvar"]) * i_t
        else:
            o2 = None; r0 = float(np.median(y_amt[act_w])) if len(act_w) else 0.0; amt = r0 * i_t; samt = np.zeros_like(i_t)
        return o1, o2, pi, amt, samt
    if G["targets"] is not None:
        cell = G["targets"]["cell"]; p_t = G["PGF"][m][cell]; i_t = G["NG"][m][cell] / (np.maximum(p_t, C.P_FLOOR) * nd)
        o1, o2, pi, amt, samt = two(G["targets"]["xy"], G["targets"]["elev"], np.full(len(cell), -1), p_t, i_t)
        g = {k: enc_tp(x, k) for k, x in (("pi", pi), ("amount", amt), ("samount", samt))}
        if G["diagnostics"]:
            tot = o1["nug"] + o1["psill"]
            fl = o1["flags"] | (G["targets"]["elev_filled"].astype(np.uint8) * QC_BITS["elevation_filled"])
            g.update(D1=enc(o1["D1"], "D1"), D7=enc(o1["D7"], "D7"), D50=enc(o1["D50"], "D50"), range=enc(o1["range"], "range"),
                     nugget_fraction=enc(o1["nug"] / np.where(tot > 0, tot, np.nan), "nugget_fraction"),
                     log10_sill=enc(np.log10(np.where(tot > 0, tot, np.nan)), "log10_sill"),
                     n_negative_weights=enc(o1["n_neg"], "n_negative_weights"), n_eff=enc(o1["neff"], "n_eff"),
                     R_pos_half=enc(o1["rpos12"], "R_pos_half"), R_pos_2_3=enc(o1["rpos23"], "R_pos_2_3"), R_pos_1_3=enc(o1["rpos13"], "R_pos_1_3"),
                     R_abs90=enc(o1["rabs90"], "R_abs90"), qc_flags=fl, wet_probability=enc(pi, "wet_probability"))
            if o2 is not None:
                t2 = o2["nug"] + o2["psill"]
                g.update(range_amount=enc(o2["range"], "range_amount"), nugget_fraction_amount=enc(o2["nug"] / np.where(t2 > 0, t2, np.nan), "nugget_fraction_amount"),
                         log10_sill_amount=enc(np.log10(np.where(t2 > 0, t2, np.nan)), "log10_sill_amount"), R_pos_half_amount=enc(o2["rpos12"], "R_pos_half_amount"))
            else:
                for k_ in ("range_amount", "nugget_fraction_amount", "log10_sill_amount", "R_pos_half_amount"): g[k_] = enc(np.full(len(cell), np.nan), k_)
        res["grid"] = g
        q = lambda a: np.nanpercentile(a, [10, 50, 90]) if np.isfinite(a).any() else [np.nan] * 3
        for k_ in ("D1", "D7", "D50", "rpos23", "rpos12", "rpos13", "rabs90", "neff", "range"):
            res.update({f"{k_}_p{p_}": x for p_, x in zip((10, 50, 90), q(o1[k_]))})
        res.update(pi_mean=float(np.nanmean(pi)), amount_mean=float(np.nanmean(amt)))
    if G["points"] is not None:
        pp = G["points"]; p_t = pp["pf"][:, m]; i_t = pp["clim"][:, m] / (np.maximum(p_t, C.P_FLOOR) * nd)
        o1, o2, pi, amt, samt = two(pp["xy"], pp["elev"], np.full(len(pp["xy"]), -1), p_t, i_t)
        res["points"] = pd.DataFrame(dict(name=pp["name"], date=res["date"], month=m + 1, pi=pi, amount=amt, samount=samt, p_target=p_t,
                                          D1=o1["D1"], D7=o1["D7"], qc_flags=o1["flags"]))
    if G["loo"]:
        if G["eval_excluded"]: tgt = np.flatnonzero(fin & ~G["usable"]); sidx = np.full(len(tgt), -1)
        else: tgt = act; sidx = act
        p_t = G["PLF"][tgt, m]; i_t = G["NL"][tgt, m] / (np.maximum(p_t, C.P_FLOOR) * nd)
        o1, o2, pi, amt, samt = two(G["st_xy"][tgt], G["st_elev"][tgt], sidx, p_t, i_t)
        res["loo"] = pd.DataFrame(dict(date=res["date"], station=G["ids"][tgt], month=np.int8(m + 1), obs=vals[tgt].astype(np.float32),
                                       pi=pi.astype(np.float32), amount=amt.astype(np.float32), samount=samt.astype(np.float32), p_target=p_t.astype(np.float32),
                                       D1=o1["D1"].astype(np.float32), D7=o1["D7"].astype(np.float32),
                                       qc_flags=(o1["flags"] | np.where(G["ns_fallback"][tgt, m], QC_BITS["station_normal_fallback"], 0)).astype(np.uint8)))
    res["seconds"] = time.time() - t0
    return res


def enc_tp(a, name):
    dt, sc, fill = TP[name]; info = np.iinfo(np.dtype(dt)); q = np.rint(np.asarray(a, float) / sc)
    return np.where(np.isfinite(q), np.clip(q, info.min + (dt[0] == "i"), info.max - (dt[0] == "u")), fill).astype(dt)


def dec_tp(q, name):
    dt, sc, fill = TP[name]; return np.where(q == fill, np.nan, q.astype(float) * sc)


def process_day(day_idx):
    if G.get("twopart"): return twopart_day(day_idx)
    t0 = time.time(); v = G["variable"]; date = G["dates"][day_idx]; tag = f"{v}|{G['mode']}|{date:%Y-%m-%d}"
    y, mm = transform(day_idx)
    act = np.flatnonzero(np.isfinite(y) & G["usable"])
    res = dict(day_idx=day_idx, date=f"{date:%Y-%m-%d}", n_active=len(act))
    if len(act) < C.N_VARIO + 1:          # the 50-station neighbourhood needs at least 51 active stations
        res["skipped"] = True; return res
    lapse = G["lapse_mode"]
    gpar = K.global_variogram(G["st_xy"], G["sdmat"], act, y if lapse == "none" else y - np.polyval(np.polyfit(G["st_elev"][act], y[act], 1), G["st_elev"]), tag)
    est_args = dict(act=act, st_xy=G["st_xy"], st_elev=G["st_elev"], sdmat=G["sdmat"], y=y, variable=v,
                    lapse_mode=lapse, day_tag=tag, global_params=gpar,
                    lapse_extrap_m=C.ANOM_LAPSE_MAX_EXTRAP_M if (G["mode"] == "anomaly" and lapse != "none") else None)
    m, ndays = mm if mm else (G["month"][day_idx] - 1, G["dim"][day_idx])
    # --- grid / sample targets
    if G["targets"] is not None:
        txy, tz = G["targets"]["xy"], G["targets"]["elev"]
        o = K.estimate_targets(txy, tz, np.full(len(txy), -1), **est_args)
        clim = G["NG"][m][G["targets"]["cell"]] if G["mode"] == "anomaly" else None
        val, ks = back(o["est"], o["kvar"], clim, m, ndays)
        fl = o["flags"] | (G["targets"]["elev_filled"].astype(np.uint8) * QC_BITS["elevation_filled"])
        tot = o["nug"] + o["psill"]
        ve = value_enc(v)
        res["grid"] = dict(value=np.where(np.isfinite(val), np.rint(val / ve[1]), ve[2]).astype(ve[0]),
                           kstd=np.where(np.isfinite(ks), np.rint(ks / ve[1]), ve[2]).astype(ve[0]))
        if G["diagnostics"]:
            res["grid"].update(D1=enc(o["D1"], "D1"), D7=enc(o["D7"], "D7"), D50=enc(o["D50"], "D50"), range=enc(o["range"], "range"),
                               nugget_fraction=enc(o["nug"] / np.where(tot > 0, tot, np.nan), "nugget_fraction"),
                               log10_sill=enc(np.log10(np.where(tot > 0, tot, np.nan)), "log10_sill"),
                               lapse=enc(o["b1"] * 1000, "lapse"), lapse_r2=enc(o["r2"], "lapse_r2"),
                               n_negative_weights=enc(o["n_neg"], "n_negative_weights"), n_eff=enc(o["neff"], "n_eff"),
                               R_pos_half=enc(o["rpos12"], "R_pos_half"),
                               R_pos_2_3=enc(o["rpos23"], "R_pos_2_3"), R_pos_1_3=enc(o["rpos13"], "R_pos_1_3"), R_abs90=enc(o["rabs90"], "R_abs90"), qc_flags=fl)
        q = lambda a: np.nanpercentile(a, [10, 50, 90]) if np.isfinite(a).any() else [np.nan] * 3
        for k_ in ("D1", "D7", "D50", "rpos23", "rpos12", "rpos13", "rabs90", "neff", "range"):
            res.update({f"{k_}_p{p}": x for p, x in zip((10, 50, 90), q(o[k_]))})
        nz = ~(o["flags"] & 1).astype(bool)
        res.update(pct_zero_shortcut=100 * (~nz).mean(), pct_negative_weight=100 * np.nanmean(o["n_neg"][nz] > 0) if nz.any() else np.nan,
                   pct_oscillatory=100 * ((o["flags"] & 8) > 0).mean(), pct_global_vario_fallback=100 * ((o["flags"] & 32) > 0).mean(),
                   pct_lapse_fallback=100 * ((o["flags"] & 2) > 0).mean(), pct_positive_lapse=100 * ((o["flags"] & 4) > 0).mean(),
                   field_mean=float(np.nanmean(val)), field_kstd_mean=float(np.nanmean(ks)))
    # --- points
    if G["points"] is not None:
        pp = G["points"]; o = K.estimate_targets(pp["xy"], pp["elev"], np.full(len(pp["xy"]), -1), **est_args)
        val, ks = back(o["est"], o["kvar"], pp["clim"][:, m] if G["mode"] == "anomaly" else None, m, ndays)
        res["points"] = pd.DataFrame(dict(name=pp["name"], date=res["date"], value=val, kriging_std=ks, D1=o["D1"], D7=o["D7"], qc_flags=o["flags"]))
    # --- leave-one-out at every active station
    if G["loo"]:
        if G["eval_excluded"]:        # predict the withheld stations from the retained network
            tgt = np.flatnonzero(np.isfinite(G["values"][day_idx]) & ~G["usable"]); sidx = np.full(len(tgt), -1)
        else:                         # classic leave-one-out over the active network
            tgt = act; sidx = act
        o = K.estimate_targets(G["st_xy"][tgt], G["st_elev"][tgt], sidx, **est_args)
        clim = G["NL"][tgt, m] if G["mode"] == "anomaly" else None
        val, ks = back(o["est"], o["kvar"], clim, m, ndays)
        obs = G["values"][day_idx][tgt].astype(float); act = tgt
        res["loo"] = pd.DataFrame(dict(date=res["date"], station=G["ids"][act], obs=obs.astype(np.float32), pred=val.astype(np.float32),
                                       kriging_std=ks.astype(np.float32), D1=o["D1"].astype(np.float32), D7=o["D7"].astype(np.float32),
                                       qc_flags=(o["flags"] | np.where(G["ns_fallback"][act, m], QC_BITS["station_normal_fallback"], 0)).astype(np.uint8)))
        e = val - obs; ok = np.isfinite(e)
        res.update(loo_n=int(ok.sum()), loo_bias=float(np.mean(e[ok])), loo_mae=float(np.mean(np.abs(e[ok]))), loo_rmse=float(np.sqrt(np.mean(e[ok] ** 2))),
                   loo_r=float(np.corrcoef(val[ok], obs[ok])[0, 1]) if ok.sum() > 3 and np.std(obs[ok]) > 0 else np.nan)
        if v == "precip":
            w = ok & (obs >= C.WET_THRESHOLD)
            res.update(loo_wet_n=int(w.sum()), loo_wet_mae=float(np.mean(np.abs(e[w]))) if w.any() else np.nan,
                       loo_wet_pbias=float(100 * e[w].sum() / obs[w].sum()) if w.any() and obs[w].sum() > 0 else np.nan)
    res["seconds"] = time.time() - t0
    return res


def init(shared): G.update(shared)


# ----------------------------------------------------------------------------- writing
def write_year(run_dir, variable, year, days, grid, tsel, arrays, diagnostics, value_file=True):
    import netCDF4
    ny, nx = int(grid["ny"]), int(grid["nx"]); iy, ix = grid["iy"][tsel], grid["ix"][tsel]; nd = len(days)
    def put(path, fields, attrs):
        tmp = path.with_suffix(".nc.part")
        with netCDF4.Dataset(tmp, "w", format="NETCDF4") as ds:
            ds.createDimension("time", nd); ds.createDimension("y", ny); ds.createDimension("x", nx)
            t = ds.createVariable("time", "i4", ("time",)); t.units = "days since 1991-01-01"; t.calendar = "standard"
            t[:] = (pd.DatetimeIndex(days) - pd.Timestamp("1991-01-01")).days
            ds.createVariable("y", "f8", ("y",))[:] = grid["y"]; ds.createVariable("x", "f8", ("x",))[:] = grid["x"]
            crs = ds.createVariable("crs", "i4"); crs.spatial_ref = "EPSG:3035"; crs.epsg_code = "EPSG:3035"
            for name, (dt, sc, fill, units, lname, data) in fields.items():
                var = ds.createVariable(name, dt, ("time", "y", "x"), zlib=True, complevel=6, shuffle=True, chunksizes=(1, ny, nx), fill_value=fill)
                var.set_auto_maskandscale(False)
                if sc != 1: var.scale_factor = sc
                var.units = units; var.long_name = lname; var.grid_mapping = "crs"
                block = np.full((nd, ny, nx), fill, dtype=dt); block[:, iy, ix] = data; var[:] = block
            for k, v in attrs.items(): setattr(ds, k, v)
            ds.complete = "yes"
        tmp.rename(path)
    ve = value_enc(variable)
    vname = "precipitation" if variable == "precip" else "temperature"; vunits = "mm day-1" if variable == "precip" else "degC"
    common = dict(source="Po River District 1 km daily %s, method v3 (%s)" % (vname, G["mode"]), crs="EPSG:3035",
                  method=G["method_text"], created=pd.Timestamp.now().isoformat(timespec="seconds"))
    if value_file:
        put(run_dir / f"{variable}_{year}.nc", {vname: (ve[0], ve[1], ve[2], vunits, f"daily {vname}", arrays["value"]),
                                                "kriging_std": (ve[0], ve[1], ve[2], vunits, "ordinary-kriging standard deviation (back-transformed)" if not G.get("twopart")
                                                                else "standard deviation of the two-part estimate (occurrence x amount)", arrays["kstd"])}, common)
    if diagnostics:
        f = {n: (ENC[n][0], ENC[n][1], ENC[n][2], ENC[n][3], ENC[n][4], arrays[n]) for n in ENC if n in arrays and not (variable == "precip" and n.startswith("lapse"))}
        put(run_dir / f"diagnostics_{variable}_{year}.nc", f, common)


def write_diagnostics(run_dir, variable, year, days, grid, tsel, arrays):
    """Diagnostics file only (two-part pass 1); same layout as write_year."""
    keep = {k: x for k, x in arrays.items() if k in ENC and not k.startswith("lapse")}
    write_year(run_dir, variable, year, days, grid, tsel, keep, True, value_file=False)


def _thresholds(pi, k):
    """pi (n_days, n_targets) with NaN; k (n_targets,) wet days wanted -> threshold so that about k days have pi >= t."""
    srt = -np.sort(-np.where(np.isfinite(pi), pi, -np.inf), axis=0)          # descending, NaN last
    n = np.isfinite(pi).sum(0); k = np.clip(np.rint(k).astype(int), 0, n)
    t = np.full(pi.shape[1], np.inf); pos = k > 0
    t[pos] = srt[k[pos] - 1, np.flatnonzero(pos)]
    return np.where(t <= 0, 1e-9, t)       # never call a cell wet with pi = 0


def _calibrate(df, key, pcol="p_target"):
    """Long table (date, month, key, pi, p_target): wet flag such that the number of wet days per key and calendar month
    equals sum(p_target) (monthly mode), or per key over all days (pooled mode, for short runs)."""
    df = df.copy()
    nmin = df.groupby("month").date.nunique().reindex(range(1, 13)).fillna(0).min()
    groups = ["month", key] if nmin >= C.CALIB_MIN_DAYS else [key]
    pi = df.pi.where(np.isfinite(df.pi)); gb = pi.groupby([df[g] for g in groups])
    rank = gb.rank(ascending=False, method="first")
    k = df[pcol].where(pi.notna()).groupby([df[g] for g in groups]).transform("sum").round()
    df["wet"] = (rank <= k) & (pi > 0)
    return df, ("monthly" if len(groups) == 2 else "pooled")


def _mix_std(pi, amt, samt):
    """standard deviation of a hurdle (Bernoulli(pi) x amount) estimate"""
    pi = np.clip(pi, 0, 1); return np.sqrt(np.maximum(pi * samt ** 2 + pi * (1 - pi) * amt ** 2, 0.0))


def finalize_twopart(run_dir, a, shared, grid, tsel, years):
    """Pass 2: occurrence calibration and final files (same names and layout as the ratio method)."""
    TPD = run_dir / "_twopart"; v = a.variable; ve = value_enc(v); t0 = time.time(); info = {}
    # ---- full grid
    gfiles = [TPD / f"grid_{y}.npz" for y in years if (TPD / f"grid_{y}.npz").exists()]
    if gfiles and not a.grid_sample:
        Z = [np.load(f) for f in gfiles]; days = pd.DatetimeIndex(np.concatenate([z["days"] for z in Z]))
        PI = np.concatenate([dec_tp(z["pi"], "pi").astype(np.float32) for z in Z]); mon = days.month.to_numpy() - 1
        cell = shared["targets"]["cell"]; PGF = shared["PGF"][:, cell]; WET = np.zeros(PI.shape, bool)
        monthly = np.bincount(mon, minlength=12).min() >= C.CALIB_MIN_DAYS
        if monthly:
            for m in range(12):
                r = np.flatnonzero(mon == m); t = _thresholds(PI[r], PGF[m] * np.isfinite(PI[r]).sum(0)); WET[r] = PI[r] >= t[None]
        else:
            k = np.nansum(np.where(np.isfinite(PI), PGF[mon], 0.0), 0); t = _thresholds(PI, k); WET = PI >= t[None]
        info["grid_calibration"] = "monthly" if monthly else "pooled"
        off = 0
        for f, z in zip(gfiles, Z):
            n = len(z["days"]); sl = slice(off, off + n); off += n
            amt = dec_tp(z["amount"], "amount"); samt = dec_tp(z["samount"], "samount"); pi = PI[sl].astype(float)
            val = np.where(WET[sl], amt, 0.0); val = np.where(np.isfinite(pi), val, np.nan); ks = _mix_std(pi, amt, samt)
            arr = dict(value=np.where(np.isfinite(val), np.rint(val / ve[1]), ve[2]).astype(ve[0]), kstd=np.where(np.isfinite(ks), np.rint(ks / ve[1]), ve[2]).astype(ve[0]))
            year = int(f.stem.split("_")[1]); write_year(run_dir, v, year, pd.DatetimeIndex(z["days"]), grid, tsel, arr, False)
        info["grid_wet_fraction"] = float(WET.mean()); del PI, WET
    # ---- grid sample (network thinning)
    sfiles = sorted(TPD.glob("grid_sample_*.parquet"))
    if sfiles:
        S = pd.concat([pd.read_parquet(f) for f in sfiles], ignore_index=True); S["date"] = S.date.astype(str)
        S["month"] = pd.to_datetime(S.date).dt.month; S["pi"] = dec_tp(S.pi.to_numpy(), "pi")
        S["p_target"] = shared["PGF"][S.month.to_numpy() - 1, S.cell.to_numpy()]
        S, info["sample_calibration"] = _calibrate(S, "cell")
        amt = dec_tp(S.amount.to_numpy(), "amount"); samt = dec_tp(S.samount.to_numpy(), "samount")
        val = np.where(S.wet, amt, 0.0); ks = _mix_std(S.pi.to_numpy(), amt, samt)
        S["value"] = np.where(np.isfinite(val), np.rint(val / ve[1]), ve[2]).astype(ve[0]); S["kstd"] = np.where(np.isfinite(ks), np.rint(ks / ve[1]), ve[2]).astype(ve[0])
        for y, g in S.groupby(S.date.str[:4]): g[["date", "cell", "value", "kstd"]].to_parquet(run_dir / f"grid_sample_{v}_{y}.parquet", index=False)
    # ---- station leave-one-out
    lfiles = sorted(TPD.glob("loo_*.parquet")); dsum = {}
    if lfiles:
        L = pd.concat([pd.read_parquet(f) for f in lfiles], ignore_index=True)
        L, info["loo_calibration"] = _calibrate(L, "station")
        L["pred"] = np.where(L.wet, L.amount, 0.0).astype(np.float32); L["kriging_std"] = _mix_std(L.pi.to_numpy(), L.amount.to_numpy(), L.samount.to_numpy()).astype(np.float32)
        L["wet_pred"] = L.wet.astype(np.uint8)
        for y, g in L.groupby(L.date.str[:4]):
            g[["date", "station", "obs", "pred", "kriging_std", "D1", "D7", "qc_flags", "pi", "amount", "p_target", "wet_pred"]].to_parquet(run_dir / f"loo_{v}_{y}.parquet", index=False)
        def day_metrics(g):
            e = (g.pred - g.obs).to_numpy(float); o = g.obs.to_numpy(float); p_ = g.pred.to_numpy(float); ok = np.isfinite(e); w = ok & (o >= C.WET_THRESHOLD)
            return pd.Series(dict(loo_n=int(ok.sum()), loo_bias=np.mean(e[ok]) if ok.any() else np.nan, loo_mae=np.mean(np.abs(e[ok])) if ok.any() else np.nan,
                                  loo_rmse=np.sqrt(np.mean(e[ok] ** 2)) if ok.any() else np.nan,
                                  loo_r=np.corrcoef(p_[ok], o[ok])[0, 1] if ok.sum() > 3 and np.std(o[ok]) > 0 and np.std(p_[ok]) > 0 else np.nan,
                                  loo_wet_n=int(w.sum()), loo_wet_mae=np.mean(np.abs(e[w])) if w.any() else np.nan,
                                  loo_wet_pbias=100 * e[w].sum() / o[w].sum() if w.any() and o[w].sum() > 0 else np.nan))
        dsum = L.groupby("date").apply(day_metrics)
        info["loo_wet_fraction_pred"] = float((L.pred >= C.OCC_THRESHOLD).mean()); info["loo_wet_fraction_obs"] = float((L.obs >= C.OCC_THRESHOLD).mean())
        info["loo_mae"] = float((L.pred - L.obs).abs().mean())
    # ---- points
    pfiles = sorted(TPD.glob("points_*.csv"))
    if pfiles:
        Pt = pd.concat([pd.read_csv(f) for f in pfiles], ignore_index=True); Pt, info["points_calibration"] = _calibrate(Pt, "name")
        Pt["value"] = np.where(Pt.wet, Pt.amount, 0.0); Pt["kriging_std"] = _mix_std(Pt.pi.to_numpy(), Pt.amount.to_numpy(), Pt.samount.to_numpy())
        for y, g in Pt.groupby(Pt.date.str[:4]):
            g[["name", "date", "value", "kriging_std", "D1", "D7", "qc_flags", "pi", "amount"]].sort_values(["name", "date"]).to_csv(run_dir / f"points_{v}_{y}.csv", index=False)
    # ---- daily summaries
    for y in years:
        f = TPD / f"daily_summary_{y}.csv"
        if not f.exists(): continue
        d = pd.read_csv(f)
        if len(dsum): d = d.drop(columns=[c for c in dsum.columns if c in d.columns]).merge(dsum, left_on="date", right_index=True, how="left")
        d.to_csv(run_dir / f"daily_summary_{v}_{y}.csv", index=False)
    (run_dir / "twopart_finalize.json").write_text(json.dumps(info, indent=1, default=float))
    print(f"two-part finalize ({(time.time()-t0)/60:.1f} min): " + ", ".join(f"{k}={x if isinstance(x, str) else round(x, 4)}" for k, x in info.items()), flush=True)


# ----------------------------------------------------------------------------- main
def build_shared(a):
    v = a.variable
    grid = P.load_grid()
    excl = [l.strip() for l in open(a.exclude_stations)] if a.exclude_stations else []
    st = P.attach_station_predictors(P.load_station_data(v), grid)
    S = len(st["ids"]); usable = ~np.isin(st["ids"], excl)
    dates = st["dates"]; sel = (dates >= pd.Timestamp(C.PERIOD[0])) & (dates <= pd.Timestamp(C.PERIOD[1]))
    shared = dict(variable=v, mode=a.mode, dates=dates, values=st["values"], ids=st["ids"], st_xy=st["xy"], st_elev=st["elev"],
                  sdmat=P.station_distance_km(st["xy"]), usable=usable, month=dates.month.to_numpy(), dim=dates.days_in_month.to_numpy().astype(float),
                  lapse_mode=("local50" if (v == "temp" and (a.mode == "direct" or a.anom_lapse == "local50")) else "none"),
                  diagnostics=not a.no_diagnostics, loo=not a.no_loo, ns_fallback=np.zeros((S, 12), bool),
                  eval_excluded=bool(a.eval_excluded))
    if a.mode == "anomaly":
        nd = pd.read_csv(C.OUT_ROOT / "normals" / f"station_normals_loo_{v}.csv", dtype={"station": str})
        NS = nd.pivot(index="station", columns="month", values="normal_station").reindex(st["ids"]).to_numpy(float, copy=True)
        NL = nd.pivot(index="station", columns="month", values="normal_grid_loo").reindex(st["ids"]).to_numpy(float, copy=True)
        fb = ~np.isfinite(NS); NS[fb] = NL[fb]
        import xarray as xr
        with xr.open_dataset(C.OUT_ROOT / "normals" / f"normals_{v}.nc") as ds:
            NG = ds["normal"].values[:, grid["iy"], grid["ix"]].astype(float)
        shared.update(NS=NS, NL=NL, NG=NG, ns_fallback=fb)
        if v == "precip" and C.PRECIP_METHOD == "twopart":
            wf = pd.read_csv(C.OUT_ROOT / "normals" / "station_wetfreq_loo_precip.csv", dtype={"station": str})
            PSF = wf.pivot(index="station", columns="month", values="wetfreq_station").reindex(st["ids"]).to_numpy(float, copy=True)
            PLF = wf.pivot(index="station", columns="month", values="wetfreq_grid_loo").reindex(st["ids"]).to_numpy(float, copy=True)
            fbf = ~np.isfinite(PSF); PSF[fbf] = PLF[fbf]
            with xr.open_dataset(C.OUT_ROOT / "normals" / "wetfreq_precip.nc") as ds:
                PGF = ds["wetfreq"].values[:, grid["iy"], grid["ix"]].astype(float)
            shared.update(twopart=True, PSF=PSF, PLF=PLF, PGF=PGF, ns_fallback=fb | fbf)
    # targets
    N = len(grid["xy"])
    if a.loo_only: tsel = None
    elif a.grid_sample: tsel = np.sort(np.random.default_rng(20261008).choice(N, a.grid_sample, replace=False))
    else: tsel = np.arange(N)
    shared["targets"] = None if tsel is None else dict(xy=grid["xy"][tsel], elev=grid["elev"][tsel], cell=tsel, elev_filled=grid["elev_filled"][tsel])
    shared["points"] = None
    if a.points:
        p = pd.read_csv(a.points)
        if {"lon", "lat"} <= set(p.columns):
            from pyproj import Transformer
            px, py = Transformer.from_crs("EPSG:4326", C.METRIC_CRS, always_xy=True).transform(p.lon.values, p.lat.values)
        else: px, py = p.x.values, p.y.values
        pxy = np.c_[px, py]; pr = P.sample_predictors_at(grid, pxy)
        pel = p["elevation"].to_numpy(float) if "elevation" in p else pr["dem"]
        clim = None
        if a.mode == "anomaly":
            from scipy.spatial import cKDTree
            _, ci = cKDTree(grid["xy"]).query(pxy); clim = shared["NG"][:, ci].T
        shared["points"] = dict(name=p["name"].astype(str).to_numpy(), xy=pxy, elev=pel, clim=clim,
                                pf=shared["PGF"][:, ci].T if shared.get("twopart") else None)
    shared["method_text"] = (f"local spherical variogram on the {C.N_VARIO} nearest active stations, PSO {C.PSO_PARTICLES}x{C.PSO_ITERATIONS}, "
                             f"ordinary kriging from the {C.N_KRIG} nearest; " + ("per-cell daily lapse on the same stations; " if shared["lapse_mode"] == "local50" else "")
                             + ("anomaly method on 1991-2020 PRISM-type monthly normals" if a.mode == "anomaly" else "direct interpolation"))
    if shared.get("twopart"):
        shared["method_text"] = (f"two-part anomaly method: occurrence = kriged occurrence anomalies (wet >= {C.OCC_THRESHOLD} mm) added to the "
                                 f"1991-2020 wet-day frequency normal and calibrated to it; amount = kriged wet-day intensity ratios times the intensity "
                                 f"normal M/(p days); local spherical variograms on the {C.N_VARIO} nearest (wet) stations, PSO {C.PSO_PARTICLES}x{C.PSO_ITERATIONS}, "
                                 f"ordinary kriging from the {C.N_KRIG} nearest")
    return shared, grid, tsel, sel


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variable", choices=["precip", "temp"], required=True)
    ap.add_argument("--mode", choices=["direct", "anomaly"], required=True)
    ap.add_argument("--anom-lapse", choices=["local50", "none"], default="local50", help="temperature, anomaly mode: local lapse on anomalies")
    ap.add_argument("--years", default=None, help="e.g. 1991-2020 or 2020")
    ap.add_argument("--dates-csv", default=None, help="restrict to dates listed in a CSV column Date/date")
    ap.add_argument("--loo-only", action="store_true", help="no grid output; only station LOO (fast method comparison)")
    ap.add_argument("--grid-sample", type=int, default=0, help="predict only N random grid cells (written to parquet)")
    ap.add_argument("--points", default=None, help="CSV name,lon,lat[,elevation] or name,x,y for point time series")
    ap.add_argument("--exclude-stations", default=None, help="file with station IDs to withhold from interpolation")
    ap.add_argument("--eval-excluded", action="store_true", help="evaluate at the --exclude-stations instead of LOO (network-thinning test)")
    ap.add_argument("--no-diagnostics", action="store_true"); ap.add_argument("--no-loo", action="store_true")
    ap.add_argument("--workers", type=int, default=12); ap.add_argument("--tag", default=None)
    a = ap.parse_args()
    run = a.tag or f"{a.variable}_{a.mode}" + ("_anomlapse-none" if a.variable == "temp" and a.mode == "anomaly" and a.anom_lapse == "none" else "")
    run_dir = C.OUT_ROOT / run; run_dir.mkdir(parents=True, exist_ok=True)
    shared, grid, tsel, sel = build_shared(a)
    dates = shared["dates"]; day_idx = np.flatnonzero(sel)
    if a.years:
        y0, y1 = (a.years.split("-") + [a.years])[:2]; yy = dates.year.to_numpy(); day_idx = day_idx[(yy[day_idx] >= int(y0)) & (yy[day_idx] <= int(y1))]
    if a.dates_csv:
        dc = pd.read_csv(a.dates_csv); col = "Date" if "Date" in dc else "date"; want = set(pd.to_datetime(dc[col]).dt.normalize())
        day_idx = np.array([i for i in day_idx if dates[i].normalize() in want])
    (run_dir / "run_config.json").write_text(json.dumps(dict(vars(a), method=shared["method_text"], n_days=int(len(day_idx)),
        n_stations=int(shared["usable"].sum()), out=str(run_dir), parameters={k: getattr(C, k) for k in dir(C) if k.isupper() and not isinstance(getattr(C, k), (dict, Path))}), indent=2, default=str))
    years = sorted(set(dates.year[day_idx]))
    print(f"run={run} kernels={K.KERNELS} variable={a.variable} mode={a.mode} lapse={shared['lapse_mode']} days={len(day_idx)} years={years[0]}-{years[-1]} "
          f"stations={int(shared['usable'].sum())} targets={'none' if tsel is None else len(tsel)} workers={a.workers}", flush=True)
    t_all = time.time()
    with mp.get_context("fork").Pool(a.workers, initializer=init, initargs=(shared,)) as pool:
        G.update(shared)
        TPD = run_dir / "_twopart"; tp = bool(shared.get("twopart"))
        if tp: TPD.mkdir(exist_ok=True)
        for year in years:
            final = run_dir / f"{a.variable}_{year}.nc"; loo_f = run_dir / f"loo_{a.variable}_{year}.parquet"
            done = (final.exists() or tsel is None or a.grid_sample) and (loo_f.exists() or not shared["loo"]) and (run_dir / f"daily_summary_{a.variable}_{year}.csv").exists()
            if tp:
                done = ((TPD / f"grid_{year}.npz").exists() or tsel is None or a.grid_sample) and ((TPD / f"loo_{year}.parquet").exists() or not shared["loo"]) \
                       and (TPD / f"daily_summary_{year}.csv").exists()
            if done: print(f"{year}: complete, skipped", flush=True); continue
            di = day_idx[dates.year[day_idx] == year]; days = dates[di]; pos = {d: k for k, d in enumerate(di)}
            arrays = None; loos, rows, pts, samples = [], [], [], []; t0 = time.time()
            for r in pool.imap_unordered(process_day, di, chunksize=1):
                if "grid" in r:
                    g = r.pop("grid")
                    if a.grid_sample:
                        keys = ("pi", "amount", "samount") if tp else ("value", "kstd")
                        samples.append(pd.DataFrame(dict(date=r["date"], cell=tsel, **{k: g[k] for k in keys})))
                    else:
                        if arrays is None:
                            fills = {k: (value_enc(a.variable)[2] if k in ("value", "kstd") else TP[k][2] if k in TP else ENC[k][2]) for k in g}
                            arrays = {k: np.full((len(di),) + x.shape, fills[k], x.dtype) for k, x in g.items()}
                        for k, x in g.items(): arrays[k][pos[r["day_idx"]]] = x
                if "loo" in r: loos.append(r.pop("loo"))
                if "points" in r: pts.append(r.pop("points"))
                rows.append(r)
            summ = pd.DataFrame(rows).sort_values("date")
            if tp:                      # pass 1 of the two-part method: keep pi / amount until the occurrence calibration
                if arrays is not None:
                    np.savez_compressed(TPD / f"grid_{year}.npz", days=np.array([f"{d:%Y-%m-%d}" for d in days]), **{k: arrays[k] for k in TP})
                    if shared["diagnostics"]: write_diagnostics(run_dir, a.variable, year, days, grid, tsel, arrays)
                if samples: pd.concat(samples).to_parquet(TPD / f"grid_sample_{year}.parquet", index=False)
                if loos: pd.concat(loos).to_parquet(TPD / f"loo_{year}.parquet", index=False)
                if pts: pd.concat(pts).to_csv(TPD / f"points_{year}.csv", index=False)
                summ.to_csv(TPD / f"daily_summary_{year}.csv", index=False)
                print(f"{year}: {len(di)} days in {(time.time()-t0)/60:.1f} min | mean worker s/day {summ.seconds.mean():.1f} | pass 1", flush=True)
                continue
            if arrays is not None:      # days skipped for too few stations keep the fill value
                write_year(run_dir, a.variable, year, days, grid, tsel, arrays, shared["diagnostics"])
            if samples: pd.concat(samples).to_parquet(run_dir / f"grid_sample_{a.variable}_{year}.parquet", index=False)
            if loos: pd.concat(loos).to_parquet(loo_f, index=False)
            if pts: pd.concat(pts).sort_values(["name", "date"]).to_csv(run_dir / f"points_{a.variable}_{year}.csv", index=False)
            summ.to_csv(run_dir / f"daily_summary_{a.variable}_{year}.csv", index=False)
            l = summ.get("loo_mae", pd.Series(dtype=float))
            print(f"{year}: {len(di)} days in {(time.time()-t0)/60:.1f} min | mean worker s/day {summ.seconds.mean():.1f} | LOO MAE {l.mean():.3f}", flush=True)
    if shared.get("twopart"): finalize_twopart(run_dir, a, shared, grid, tsel, years)
    print(f"DONE in {(time.time()-t_all)/60:.1f} min -> {run_dir}")


if __name__ == "__main__":
    main()
