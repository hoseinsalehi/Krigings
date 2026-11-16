# -*- coding: utf-8 -*-
"""numba-compiled versions of the two hot kernels of po_core (93 % of run time).

They reproduce po_core.empirical_batch and po_core.pso_spherical operation by operation
(same random stream from numpy's Generator, same IEEE operation order, no fastmath), so the
results are expected to be bit-identical; test_numba_equivalence() checks this on real data.
"""
import math
import numpy as np
import numba as nb
import sys
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).resolve().parent))

import po_config as C


# ----------------------------------------------------------------------------- empirical variogram
@nb.njit(cache=True, parallel=False)
def _empirical(Y, D, n_lags, min_pairs):
    B, n = Y.shape
    lag = np.zeros((B, n_lags)); gam = np.zeros((B, n_lags)); cnt = np.zeros((B, n_lags))
    npairs = np.zeros(B, np.int64)
    for b in range(B):
        # maximum valid pair distance
        md = -np.inf
        for i in range(n):
            yi = Y[b, i]
            for j in range(i + 1, n):
                d = np.float64(D[b, i, j]); sv = 0.5 * (yi - Y[b, j]) ** 2
                if math.isfinite(d) and math.isfinite(sv) and d > 0:
                    if d > md: md = d
        if not (math.isfinite(md) and md > 0): md = 1.0
        stop = md * (1 + 1e-10); step = stop / n_lags
        # edges as np.linspace(0, stop, n_lags+1): k*step, last edge = stop exactly
        for i in range(n):
            yi = Y[b, i]
            for j in range(i + 1, n):
                d = np.float64(D[b, i, j]); sv = 0.5 * (yi - Y[b, j]) ** 2
                if not (math.isfinite(d) and math.isfinite(sv) and d > 0): continue
                npairs[b] += 1
                k = 0
                for e in range(1, n_lags + 1):
                    edge = stop if e == n_lags else e * step
                    if d >= edge: k += 1
                if k > n_lags - 1: k = n_lags - 1
                cnt[b, k] += 1.0; lag[b, k] += d; gam[b, k] += sv
    return lag, gam, cnt, npairs


def empirical_batch(Y, D, n_lags=C.N_LAGS, min_pairs=C.MIN_PAIRS_PER_LAG):
    lag, gam, cnt, npairs = _empirical(np.ascontiguousarray(Y, np.float64), np.ascontiguousarray(D), n_lags, min_pairs)
    mask = cnt >= min_pairs
    with np.errstate(invalid="ignore", divide="ignore"):
        lag = np.where(mask, lag / cnt, np.nan); gam = np.where(mask, gam / cnt, np.nan)
    cnt = np.where(mask, cnt, 0.0)
    B = Y.shape[0]
    gmax = np.nanmax(np.where(mask, gam, np.nan), 1) if mask.any() else np.zeros(B)
    valid = (npairs >= max(20, min_pairs * 3)) & (mask.sum(1) >= 3) & (np.nan_to_num(gmax) > 1e-12)
    return lag, gam, cnt, mask, valid


# ----------------------------------------------------------------------------- PSO
@nb.njit(cache=True, inline="always")
def _pairwise_sum12(a):
    # numpy's pairwise summation for a contiguous reduction of length n (8-way unrolled block)
    n = a.shape[0]
    if n < 8:
        s = 0.0
        for i in range(n): s += a[i]
        return s
    r0, r1, r2, r3, r4, r5, r6, r7 = a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]
    i = 8
    while i < n - (n % 8):
        r0 += a[i]; r1 += a[i + 1]; r2 += a[i + 2]; r3 += a[i + 3]; r4 += a[i + 4]; r5 += a[i + 5]; r6 += a[i + 6]; r7 += a[i + 7]
        i += 8
    s = ((r0 + r1) + (r2 + r3)) + ((r4 + r5) + (r6 + r7))
    while i < n:
        s += a[i]; i += 1
    return s


@nb.njit(cache=True, inline="always")
def _obj(x0, x1, x2, h, g0, ww, buf):
    L = h.shape[0]; a = max(x2, 1e-12)
    for l in range(L):
        r = h[l] / a; rr = r if r < 1.0 else 1.0
        v = rr * rr; v = v * rr; v = v * -0.5; v = v + 1.5 * rr
        v = v * x1; v = v + x0; v = v - g0[l]; v = v * v; v = v * ww[l]
        buf[l] = v
    return _pairwise_sum12(buf)


@nb.njit(cache=True)
def _pso(lo, hi, h, g0, ww, U0, N0, R, inertia, c1, c2):
    B, P, _ = U0.shape; it = R.shape[0]; L = h.shape[1]
    gb_out = np.empty((B, 3)); gbs_out = np.empty(B)
    pos = np.empty((P, 3)); vel = np.empty((P, 3)); pb = np.empty((P, 3)); pbs = np.empty(P); sc = np.empty(P); buf = np.empty(L)
    for b in range(B):
        for p in range(P):
            for k in range(3):
                span = hi[b, k] - lo[b, k]
                pos[p, k] = lo[b, k] + U0[b, p, k] * span
                vel[p, k] = N0[b, p, k] * span
        for p in range(P):
            sc[p] = _obj(pos[p, 0], pos[p, 1], pos[p, 2], h[b], g0[b], ww[b], buf)
            pbs[p] = sc[p]
            for k in range(3): pb[p, k] = pos[p, k]
        gi = 0
        for p in range(1, P):
            if sc[p] < sc[gi]: gi = p
        g = np.empty(3)
        for k in range(3): g[k] = pos[gi, k]
        gbs = sc[gi]
        for t in range(it):
            for p in range(P):
                for k in range(3):
                    v = inertia * vel[p, k]
                    v = v + (c1 * R[t, 0, b, p, k]) * (pb[p, k] - pos[p, k])
                    v = v + (c2 * R[t, 1, b, p, k]) * (g[k] - pos[p, k])
                    vel[p, k] = v
                    x = pos[p, k] + v
                    if x < lo[b, k]: x = lo[b, k]
                    if x > hi[b, k]: x = hi[b, k]
                    pos[p, k] = x
            for p in range(P):
                sc[p] = _obj(pos[p, 0], pos[p, 1], pos[p, 2], h[b], g0[b], ww[b], buf)
                if sc[p] < pbs[p]:
                    pbs[p] = sc[p]
                    for k in range(3): pb[p, k] = pos[p, k]
            gi = 0
            for p in range(1, P):
                if pbs[p] < pbs[gi]: gi = p
            if pbs[gi] < gbs:
                gbs = pbs[gi]
                for k in range(3): g[k] = pb[gi, k]
        for k in range(3): gb_out[b, k] = g[k]
        gbs_out[b] = gbs
    return gb_out, gbs_out


def pso_spherical(lag, gam, cnt, mask, seed):
    B, L = lag.shape; P = C.PSO_PARTICLES; rng = np.random.default_rng(seed)
    with np.errstate(invalid="ignore"):
        maxg = np.maximum(np.nanmax(np.where(mask, gam, np.nan), 1), 1e-8)
        minh = np.maximum(np.nanmin(np.where(mask, lag, np.nan), 1), 1e-4)
        maxh = np.maximum(np.nanmax(np.where(mask, lag, np.nan), 1), minh * 1.01)
    lo = np.stack([np.zeros(B), np.full(B, 1e-12), np.maximum(0.25 * minh, 1e-4)], 1)
    hi = np.stack([1.5 * maxg, 3 * maxg, 2 * maxh], 1)
    # identical random stream to po_core.pso_spherical: random, normal, then (r1, r2) per iteration
    U0 = rng.random((B, P, 3)); N0 = rng.normal(0, .06, (B, P, 3))
    R = rng.random((C.PSO_ITERATIONS, 2, B, P, 3))
    if C.FIT_WEIGHTING == "pair-count": ww = cnt / np.maximum(cnt.sum(1, keepdims=True), 1)
    else: ww = mask / np.maximum(mask.sum(1, keepdims=True), 1)
    g0 = np.where(mask, gam, 0.0); h = np.where(mask, lag, 1.0)
    gb, gbs = _pso(lo, hi, np.ascontiguousarray(h), np.ascontiguousarray(g0), np.ascontiguousarray(ww), U0, N0, R,
                   C.PSO_INERTIA, C.PSO_C1, C.PSO_C2)
    return gb, gbs, gb[:, 2] >= 0.975 * hi[:, 2]
