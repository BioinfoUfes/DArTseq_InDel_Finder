"""Local prediction of the melting curve (helicity vs temperature) for HRM.

Zimm-Bragg/Ising helix-coil model (transfer matrix, exactly solvable, O(N)) with
nearest-neighbor parameters (SantaLucia 1998) and salt correction. It reproduces the SAME
kind of output as uMelt (helicity and the derivative -dθ/dT) and Tm = the derivative peak,
WITHOUT relying on any external site. It is a cooperative approximation (it does not
include loop entropy à la Poland-Scheraga), honest for screening and comparing alleles;
bench validation and, if desired, a spot check in uMelt still apply.
"""
import math

R = 1.98720425864083e-3  # kcal/mol/K

# SantaLucia (1998) unified NN: (dH kcal/mol, dS cal/mol/K)
NN = {
    'AA': (-7.6, -21.3), 'TT': (-7.6, -21.3),
    'AT': (-7.2, -20.4), 'TA': (-7.2, -21.3),
    'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
    'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
    'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0),
    'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
    'CG': (-10.6, -27.2),
    'GC': (-9.8, -24.4),
    'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9),
}


def na_equivalent_mM(mono_mM, mg_mM, dntp_mM):
    mg_free = max(mg_mM - dntp_mM, 0.0)
    return mono_mM + 120.0 * math.sqrt(mg_free) if mg_free > 0 else mono_mM


def _stacks(seq):
    seq = seq.upper().replace('-', '')
    out = []
    for i in range(len(seq) - 1):
        d = seq[i:i + 2]
        if d in NN:
            out.append(NN[d])
        else:
            out.append((-8.0, -22.0))  # fallback for an ambiguous base
    return out


def _lnZ(stacks, T, na_M, sigma, lam=1.0, open_mask=None):
    """ln of the partition function (Ising/Zimm-Bragg) with scaling (overflow-safe).
    open_mask: set of stack indices forced to COIL (e.g. the heteroduplex loop/bulge)."""
    salt = 0.368 * math.log(na_M)      # per-NN salt correction (cal/mol/K)
    v_h, v_c = 0.0, 1.0                 # initial state: coil
    lnscale = 0.0
    for i, (dH, dS) in enumerate(stacks):
        if open_mask is not None and i in open_mask:
            s = 0.0                     # unpaired region (bulge): never helical
        else:
            dS_c = dS + salt
            dG = dH - T * dS_c / 1000.0
            s = lam * math.exp(-dG / (R * T))
        n_h = v_h * s + v_c * sigma * s
        n_c = v_h + v_c
        tot = n_h + n_c
        v_h, v_c = n_h / tot, n_c / tot
        lnscale += math.log(tot)
    return lnscale + math.log(v_h + v_c)


def _theta_curve(stacks, na_M, sigma, tmin, tmax, step, open_mask=None):
    n = len(stacks); eps = 1e-3
    temps, theta = [], []
    for i in range(int((tmax - tmin) / step) + 1):
        Tc = tmin + i * step; Tk = Tc + 273.15
        lnZ1 = _lnZ(stacks, Tk, na_M, sigma, 1.0, open_mask)
        lnZ2 = _lnZ(stacks, Tk, na_M, sigma, 1.0 + eps, open_mask)
        th = (lnZ2 - lnZ1) / (n * math.log(1.0 + eps))
        temps.append(Tc); theta.append(min(max(th, 0.0), 1.0))
    return temps, theta


def _deriv_tm(temps, theta, step):
    deriv = [0.0] * len(theta)
    for i in range(1, len(theta) - 1):
        deriv[i] = -(theta[i + 1] - theta[i - 1]) / (2 * step)
    Tm = None
    if len(deriv) > 2:
        k = max(range(1, len(deriv) - 1), key=lambda j: deriv[j])
        Tm = round(temps[k], 2)
    return deriv, Tm


def helicity(seq, mono_mM=20.0, mg_mM=3.0, dntp_mM=0.0, sigma=1e-4,
             tmin=45.0, tmax=98.0, step=0.5, open_mask=None):
    """Returns (temps[C], theta[0..1], deriv[-dtheta/dT], Tm[C])."""
    stacks = _stacks(seq)
    if len(stacks) < 5:
        return [], [], [], None
    na_M = max(na_equivalent_mM(mono_mM, mg_mM, dntp_mM), 1.0) / 1000.0
    temps, theta = _theta_curve(stacks, na_M, sigma, tmin, tmax, step, open_mask)
    deriv, Tm = _deriv_tm(temps, theta, step)
    return temps, theta, deriv, Tm


def het_curve(seq_long, indel_start, indel_len, mono_mM=20.0, mg_mM=3.0, dntp_mM=0.0,
              sigma=1e-4, tmin=45.0, tmax=98.0, step=0.5):
    """Curve of a HETEROZYGOUS sample (long + short alleles).
    Mixture: 25% long homoduplex + 25% short homoduplex + 50% heteroduplex
    (long x short strands with an unpaired loop of 'indel_len' bases -> melts earlier)."""
    na_M = max(na_equivalent_mM(mono_mM, mg_mM, dntp_mM), 1.0) / 1000.0
    long_stacks = _stacks(seq_long)
    seq_short = seq_long[:indel_start] + seq_long[indel_start + indel_len:]
    short_stacks = _stacks(seq_short)
    # heteroduplex: the long backbone with the InDel region forced open (bulge)
    open_mask = set(range(max(0, indel_start - 1), min(len(long_stacks), indel_start + indel_len + 1)))
    tL, thL = _theta_curve(long_stacks, na_M, sigma, tmin, tmax, step)
    tS, thS = _theta_curve(short_stacks, na_M, sigma, tmin, tmax, step)
    tH, thH = _theta_curve(long_stacks, na_M, sigma, tmin, tmax, step, open_mask)
    # align on a common grid
    import numpy as np
    grid = np.arange(tmin, tmax, step)
    iL = np.interp(grid, tL, thL); iS = np.interp(grid, tS, thS); iH = np.interp(grid, tH, thH)
    theta = (0.25 * iL + 0.25 * iS + 0.5 * iH)
    temps = list(grid)
    deriv, Tm = _deriv_tm(temps, list(theta), step)
    return temps, list(theta), deriv, Tm


def allele_delta(seq_long, seq_short, **kw):
    tL, thL, dL, TmL = helicity(seq_long, **kw)
    tS, thS, dS, TmS = helicity(seq_short, **kw)
    dtm = round(abs(TmL - TmS), 2) if (TmL is not None and TmS is not None) else None
    return dict(temps_long=tL, theta_long=thL, deriv_long=dL, tm_long=TmL,
                temps_short=tS, theta_short=thS, deriv_short=dS, tm_short=TmS, dtm=dtm)
