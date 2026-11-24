"""RG-C Method 1, Eqs. 7, 8, 10 and ND3 Table II, October 16, 2026.

Rates are selected counts / (Q+ + Q-). Coefficients include thermal corrections.
MT in the note is He here. F (foils only) in the note is ET here.
"""
import numpy as np

COEFFICIENTS = dict(A=-0.765506, B=-0.000757116, C=-0.233737,
                    D=-0.170112, E=-0.863533, F=0.0336457, G=0.441609)
SOURCE = 'https://clasweb.jlab.org/wiki/images/7/72/Dilution_Factor_Determination_for_RG_C_10-16-26.pdf'


def divide(a, b):
    a, b = np.broadcast_arrays(a, b)
    out = np.full(a.shape, np.nan)
    return np.divide(a, b, out=out, where=np.isfinite(b) & (b != 0))


def transfer_ratio(sp23):
    """Charge-normalized CD2/CH2 ratio, per bin, as proposed in Sec. III.3."""
    return divide(sp23['CD2'], sp23['CH2'])


def method1(rates, ratio=None):
    """Direct Method 1. Sp23 uses measured CD2, other periods use ratio*CH2.

    Packing fraction is a cross-check only. It does not enter the DF estimator.
    Return unconstrained estimates, leaving nonphysical results visible.
    """
    c = COEFFICIENTS
    cd = rates['CD2'] if ratio is None else ratio*rates['CH2']
    signal = cd + c['A']*rates['C'] + c['B']*rates['He'] + c['C']*rates['ET']
    denominator = cd + c['D']*rates['C'] + c['E']*rates['He'] + c['F']*rates['ET']
    amplitude = divide(rates['ND3'] - rates['He'], denominator)
    return dict(dilution=divide(amplitude*signal, rates['ND3']),
                packing_fraction=(5.86/5)*c['G']*amplitude,
                packing_fraction_bath=c['G']*amplitude,
                denominator=denominator, signal=signal, effective_cd2_rate=cd)
