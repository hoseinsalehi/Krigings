# -*- coding: utf-8 -*-
"""Vectorised kernels: local empirical variogram, batched PSO fit of a spherical model,
ordinary kriging with variance and weight diagnostics, per-target daily lapse.

All functions work on chunks of B targets at once. Numerics follow the 360-day
diagnostic code (binning, PSO update rule, bounds); only the model is fixed to spherical.
"""
import hashlib
import sys
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).resolve().parent))
import numpy as np
from scipy.spatial import cKDTree
import os
import po_config as C


# ----------------------------------------------------------------------------- model
def spherical_core(r, out=None):
    rr = np.minimum(r, 1.0)
    if out is None: return 1.5 * rr - 0.5 * rr ** 3
    np.multiply(rr, rr, out=out); np.multiply(out, rr, out=out); np.multiply(out, -0.5, out=out); np.add(out, 1.5 * rr, out=out)
    return out


def gamma_sph(h, nug, c, a):
    """gamma(h) for h (B,...) and per-row parameters (B,). gamma(0) = 0."""
    tail = (1,) * (h.ndim - 1)
    g = nug.reshape((-1,) + tail) + c.reshape((-1,) + tail) * spherical_core(h / np.maximum(a, 1e-12).reshape((-1,) + tail))
    return np.where(h < 1e-9, 0.0, g)


# ----------------------------------------------------------------------------- empirical variogram
def empirical_batch_np(Y, D, n_lags=C.N_LAGS, min_pairs=C.MIN_PAIRS_PER_LAG):
    """Y (B,n) values, D (B,n,n) km. Equal-width bins from 0 to the row's maximum pair distance."""
    B, n = Y.shape; iu0, iu1 = np.triu_indices(n, 1)
    d = D[:, iu0, iu1].astype(np.float64); sv = 0.5 * (Y[:, iu0] - Y[:, iu1]) ** 2
    ok = np.isfinite(d) & np.isfinite(sv) & (d > 0); npairs = ok.sum(1)
    md = np.where(ok, d, -np.inf).max(1); md = np.where(np.isfinite(md) & (md > 0), md, 1.0)
    edges = np.linspace(np.zeros(B), md * (1 + 1e-10), n_lags + 1, axis=-1)
    b = np.minimum((d[:, :, None] >= edges[:, None, 1:]).sum(-1), n_lags - 1)
    flat = (np.broadcast_to(np.arange(B)[:, None], d.shape) * n_lags + b)[ok]
    cnt = np.bincount(flat, minlength=B * n_lags).reshape(B, n_lags).astype(float)
    lag = np.bincount(flat, weights=d[ok], minlength=B * n_lags).reshape(B, n_lags)
    gam = np.bincount(flat, weights=sv[ok], minlength=B * n_lags).reshape(B, n_lags)
    mask = cnt >= min_pairs
    with np.errstate(invalid="ignore", divide="ignore"):
        lag = np.where(mask, lag / cnt, np.nan); gam = np.where(mask, gam / cnt, np.nan)
    cnt = np.where(mask, cnt, 0.0)
    gmax = np.nanmax(np.where(mask, gam, np.nan), 1) if mask.any() else np.zeros(B)
    valid = (npairs >= max(20, min_pairs * 3)) & (mask.sum(1) >= 3) & (np.nan_to_num(gmax) > 1e-12)
    return lag, gam, cnt, mask, valid


# ----------------------------------------------------------------------------- PSO (spherical)
def seed_for(*parts):
    h = hashlib.blake2b(digest_size=8)
    for p in parts: h.update(str(p).encode())
    return int.from_bytes(h.digest(), "little") & 0xFFFFFFFF


def pso_spherical_np(lag, gam, cnt, mask, seed):
    """Batched PSO, same update rule and bounds as the 360-day diagnostic.
    Returns params (B,3) = nugget, partial sill, range [km]; weighted SSE (B,); range_bound_hit (B,)."""
    B, L = lag.shape; P = C.PSO_PARTICLES; rng = np.random.default_rng(seed)
    with np.errstate(invalid="ignore"):
        maxg = np.maximum(np.nanmax(np.where(mask, gam, np.nan), 1), 1e-8)
        minh = np.maximum(np.nanmin(np.where(mask, lag, np.nan), 1), 1e-4)
        maxh = np.maximum(np.nanmax(np.where(mask, lag, np.nan), 1), minh * 1.01)
    lo = np.stack([np.zeros(B), np.full(B, 1e-12), np.maximum(0.25 * minh, 1e-4)], 1)[:, None, :]
    hi = np.stack([1.5 * maxg, 3 * maxg, 2 * maxh], 1)[:, None, :]
    pos = lo + rng.random((B, P, 3)) * (hi - lo); vel = rng.normal(0, .06, (B, P, 3)) * (hi - lo)
    if C.FIT_WEIGHTING == "pair-count": ww = cnt / np.maximum(cnt.sum(1, keepdims=True), 1)
    else: ww = mask / np.maximum(mask.sum(1, keepdims=True), 1)
    ww = ww[:, None, :]; g0 = np.where(mask, gam, 0.0)[:, None, :]; h = np.where(mask, lag, 1.0)[:, None, :]
    buf = np.empty((B, P, L))
    def obj(X):
        spherical_core(h / np.maximum(X[..., 2, None], 1e-12), buf)
        np.multiply(buf, X[..., 1, None], out=buf); np.add(buf, X[..., 0, None], out=buf)
        np.subtract(buf, g0, out=buf); np.square(buf, out=buf); np.multiply(buf, ww, out=buf)
        return buf.sum(-1)
    sc = obj(pos); pb = pos.copy(); pbs = sc.copy(); gi = sc.argmin(1)
    gb = np.take_along_axis(pos, gi[:, None, None], 1)[:, 0]; gbs = np.take_along_axis(sc, gi[:, None], 1)[:, 0]
    for _ in range(C.PSO_ITERATIONS):
        r1 = rng.random(pos.shape); r2 = rng.random(pos.shape)
        vel = C.PSO_INERTIA * vel + C.PSO_C1 * r1 * (pb - pos) + C.PSO_C2 * r2 * (gb[:, None, :] - pos)
        pos = np.clip(pos + vel, lo, hi); sc = obj(pos)
        imp = sc < pbs; pb[imp] = pos[imp]; pbs[imp] = sc[imp]
        gi = pbs.argmin(1); cand = np.take_along_axis(pbs, gi[:, None], 1)[:, 0]; better = cand < gbs
        gbs = np.where(better, cand, gbs); gb = np.where(better[:, None], np.take_along_axis(pb, gi[:, None, None], 1)[:, 0], gb)
    bound_hit = gb[:, 2] >= 0.975 * hi[:, 0, 2]
    return gb, gbs, bound_hit


# ----------------------------------------------------------------------------- lapse
def lapse_batch(E, T):
    """Row-wise OLS T = b0 + b1 E on finite pairs. Returns b0, b1, r2, span, n."""
    ok = np.isfinite(E) & np.isfinite(T); n = ok.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(ok, E, 0.); t = np.where(ok, T, 0.); zm = z.sum(1) / n; tm = t.sum(1) / n
        dz = np.where(ok, E - zm[:, None], 0.); dt = np.where(ok, T - tm[:, None], 0.)
        szz = (dz * dz).sum(1); szt = (dz * dt).sum(1); stt = (dt * dt).sum(1)
        b1 = szt / szz; b0 = tm - b1 * zm; r2 = np.where(stt > 0, 1 - (stt - b1 * szt) / stt, np.nan)
    span = np.where(ok, E, -np.inf).max(1) - np.where(ok, E, np.inf).min(1)
    good = (n >= 3) & (szz > 1e-8)
    nan = np.nan
    return np.where(good, b0, nan), np.where(good, b1, nan), np.where(good, r2, nan), np.where(good, span, nan), n


# ----------------------------------------------------------------------------- kriging
def ok_solve(xyk, dk, nug, c, a):
    """Ordinary kriging with k neighbours. xyk (B,k,2) m, dk (B,k) km.
    Returns weights (B,k), OK variance (B,), ok (B,)."""
    B, k, _ = xyk.shape
    dd = np.sqrt(((xyk[:, :, None, :] - xyk[:, None, :, :]) ** 2).sum(-1)) / 1000.0
    A = np.ones((B, k + 1, k + 1)); A[:, :k, :k] = gamma_sph(dd, nug, c, a); A[:, k, k] = 0.0
    g0 = gamma_sph(dk, nug, c, a); rhs = np.ones((B, k + 1)); rhs[:, :k] = g0
    w = np.full((B, k), np.nan); mu = np.full(B, np.nan); ok = np.zeros(B, bool)
    try:
        s = np.linalg.solve(A, rhs[..., None])[..., 0]; w[:] = s[:, :k]; mu[:] = s[:, k]; ok[:] = np.isfinite(s).all(1)
    except np.linalg.LinAlgError:
        for i in range(B):
            try: s = np.linalg.solve(A[i], rhs[i]); w[i] = s[:k]; mu[i] = s[k]; ok[i] = np.isfinite(s).all()
            except np.linalg.LinAlgError: pass
    var = np.maximum((w * g0).sum(1) + mu, 0.0)
    return w, var, ok


def weight_diagnostics(w, d):
    """Support metrics of the kriging weights; d (B,k) km sorted ascending."""
    aw = np.abs(w); s = aw.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = aw / s[:, None]; neff = 1.0 / (p * p).sum(1)
        pos = np.where(w > 0, w, np.nan); mx = np.nanmax(pos, 1)
        rpos = np.max(np.where(w >= 0.5 * mx[:, None], d, -np.inf), 1)
        rpos23 = np.max(np.where(w >= (2 / 3) * mx[:, None], d, -np.inf), 1)
        rpos13 = np.max(np.where(w >= (1 / 3) * mx[:, None], d, -np.inf), 1)
        cum = np.cumsum(p, 1); i90 = np.argmax(cum >= 0.9, 1); rabs90 = np.take_along_axis(d, i90[:, None], 1)[:, 0]
    return dict(n_neg=(w < -1e-10).sum(1), neff=neff, rpos12=rpos, rpos23=rpos23, rpos13=rpos13, rabs90=rabs90, sum_abs_w=s)


# ----------------------------------------------------------------------------- one day, many targets
def estimate_targets(txy, telev, self_idx, act, st_xy, st_elev, sdmat, y, variable, lapse_mode, day_tag,
                     global_params=None, lapse_extrap_m=None):
    """Kriging estimates at targets for one day.

    txy (M,2), telev (M,), self_idx (M,) station index of the target or -1 (LOO exclusion),
    act (A,) indices of active stations, y (S,) value to krige at every station (raw, or anomaly),
    lapse_mode 'local50' | 'none'. Returns dict of (M,) arrays in the kriged quantity's units
    (trend added back if a lapse is used).
    """
    M = len(txy); A = len(act); kq = min(C.N_VARIO + 1, A)
    tree = cKDTree(st_xy[act]); dm, il = tree.query(txy, k=kq)
    nb_all = act[il]
    # drop the target's own station (LOO): stable sort moves it to the end, then keep N_VARIO
    order = np.argsort(nb_all == self_idx[:, None], axis=1, kind="stable")[:, :C.N_VARIO]
    nb = np.take_along_axis(nb_all, order, 1); dkm = np.take_along_axis(dm, order, 1) / 1000.0
    out = {k: np.full(M, np.nan) for k in ("est", "kvar", "nug", "psill", "range", "sse", "b1", "r2", "span",
                                            "n_neg", "neff", "rpos12", "rpos23", "rpos13", "rabs90", "sum_abs_w")}
    out.update(D1=dkm[:, 0], D7=dkm[:, C.N_KRIG - 1], D50=dkm[:, -1])
    flags = np.zeros(M, np.uint8)            # bit0 zero-shortcut, bit1 lapse fallback, bit2 inversion, bit3 oscillatory,
                                              # bit4 range bound hit, bit5 global-variogram fallback
    trend = np.zeros(M)
    for s0 in range(0, M, C.CHUNK):
        sl = slice(s0, min(M, s0 + C.CHUNK)); n_ = nb[sl]; B = n_.shape[0]
        Y = y[n_].astype(float)
        if lapse_mode == "local50":
            b0, b1, r2, span, nz = lapse_batch(st_elev[n_], Y)
            fb = ~np.isfinite(b1) | (span < C.LAPSE_MIN_SPAN_M) | (nz < C.LAPSE_MIN_STATIONS)
            if fb.any():
                k2 = min(C.LAPSE_FALLBACK_N + 1, A); _, il2 = tree.query(txy[sl][fb], k=k2); nb2 = act[il2]
                si = self_idx[sl][fb]; Yv = np.where(nb2 == si[:, None], np.nan, y[nb2])
                b0f, b1f, r2f, spf, _ = lapse_batch(st_elev[nb2], Yv)
                b0[fb], b1[fb], r2[fb], span[fb] = b0f, b1f, r2f, spf
            b0 = np.nan_to_num(b0); b1 = np.nan_to_num(b1)
            Y = Y - (b0[:, None] + b1[:, None] * st_elev[n_]); ze = telev[sl]
            if lapse_extrap_m is not None:          # no daily-lapse extrapolation far beyond the local station range
                zs = np.where(np.isfinite(Y), st_elev[n_], np.nan)
                ze = np.clip(ze, np.nanmin(zs, 1) - lapse_extrap_m, np.nanmax(zs, 1) + lapse_extrap_m)
            trend[sl] = b0 + b1 * ze
            out["b1"][sl], out["r2"][sl], out["span"][sl] = b1, r2, span
            flags[sl] |= (fb.astype(np.uint8) << 1) | ((b1 > 0).astype(np.uint8) << 2)
        zero = (Y[:, :C.N_KRIG] <= C.ZERO_TOL).all(1) if variable == "precip" and lapse_mode == "none" else np.zeros(B, bool)
        need = ~zero
        params = np.tile([0.0, 1.0, 100.0], (B, 1)); sse = np.full(B, np.nan); bh = np.zeros(B, bool); vfb = np.zeros(B, bool)
        if need.any():
            q = np.flatnonzero(need); nq = n_[q]
            lag, gam, cnt, msk, valid = empirical_batch(Y[q], sdmat[nq[:, :, None], nq[:, None, :]])
            if not valid.all():                       # neutralise invalid rows before PSO
                bad = ~valid; msk[bad] = False; msk[bad, :3] = True; lag[bad, :3] = [1., 2., 3.]; gam[bad, :3] = 1.; cnt[bad, :3] = 1.
            p_, s_, b_ = pso_spherical(lag, gam, cnt, msk, seed_for(day_tag, s0))
            if global_params is not None and (~valid).any():
                p_[~valid] = global_params; vfb[q[~valid]] = True
            params[q], sse[q], bh[q] = p_, s_, b_
        w, var, ok = ok_solve(st_xy[n_[:, :C.N_KRIG]], dkm[sl][:, :C.N_KRIG], params[:, 0], params[:, 1], params[:, 2])
        est = (w * Y[:, :C.N_KRIG]).sum(1)
        est = np.where(zero, 0.0, est); var = np.where(zero, 0.0, var); ok |= zero
        wd = weight_diagnostics(w, dkm[sl][:, :C.N_KRIG])
        out["est"][sl] = np.where(ok, est, np.nan) + trend[sl]; out["kvar"][sl] = np.where(ok, var, np.nan)
        pz = np.where(zero[:, None], np.nan, params)
        out["nug"][sl], out["psill"][sl], out["range"][sl], out["sse"][sl] = pz[:, 0], pz[:, 1], pz[:, 2], sse
        for k_ in ("n_neg", "neff", "rpos12", "rpos23", "rpos13", "rabs90", "sum_abs_w"): out[k_][sl] = np.where(zero, np.nan, wd[k_])
        flags[sl] |= zero.astype(np.uint8) | ((wd["sum_abs_w"] > C.OSC_SUM_ABS_W) & ~zero).astype(np.uint8) << 3 \
                     | (bh.astype(np.uint8) << 4) | (vfb.astype(np.uint8) << 5)
    out["flags"] = flags
    return out


def global_variogram(st_xy, sdmat, act, yv, day_tag):
    """One spherical variogram on all active stations (fallback for invalid local variograms)."""
    D = sdmat[np.ix_(act, act)][None].astype(float); Y = yv[act][None].astype(float)
    lag, gam, cnt, msk, valid = empirical_batch(Y, D)
    if not valid[0]: return None
    p, _, _ = pso_spherical(lag, gam, cnt, msk, seed_for(day_tag, "global"))
    return p[0]


# ----------------------------------------------------------------------------- kernel dispatch
# numba versions are bit-identical to the numpy ones (checked by test_numba.py) and ~3x faster.
# Set PO_NUMBA=0 to force the numpy kernels.
try:
    if os.environ.get("PO_NUMBA", "1") == "0": raise ImportError
    import po_numba as _NB
    empirical_batch, pso_spherical, KERNELS = _NB.empirical_batch, _NB.pso_spherical, "numba"
except ImportError:
    empirical_batch, pso_spherical, KERNELS = empirical_batch_np, pso_spherical_np, "numpy"
