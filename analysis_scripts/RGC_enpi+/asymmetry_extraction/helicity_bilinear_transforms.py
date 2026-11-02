"""Definitions and exact linear transformations for RGC helicity bilinears.

The source structure-function order is the order used by
``extract_structure_function_ratios.py``:

    (u1, u2, lu1, ul1, ul2, ll0, ll1)

The transformed quantities are normalized helicity-amplitude bilinears in the
notation of equations (27)--(35) of ``formulas-SF-vs-HSF.pdf``.  ``C`` below
is the common two-body phase-space factor and the denominator is the
unseparated F_UU = F_UU,T + epsilon F_UU,L used by the RGC analysis.
"""

from __future__ import annotations

import numpy as np


SOURCE_PARAMETERS = ("u1", "u2", "lu1", "ul1", "ul2", "ll0", "ll1")

HELICITY_PARAMETERS = (
    "im_T0p_T1p",
    "im_T0m_T1m",
    "im_Tm1p_T1p",
    "delta_T1_strength",
    "re_T0p_T1p",
    "re_T0m_T1m",
    "re_Tm1p_T1p",
)

SQRT2 = np.sqrt(2.0)

# h = TRANSFORM @ f, with
# f = (FUUcosphi/FUU, FUUcos2phi/FUU, FLUsinphi/FUU,
#      FULsinphi/FUU, FULsin2phi/FUU, FLL/FUU, FLLcosphi/FUU).
TRANSFORM = np.asarray(
    [
        [0.0, 0.0, -1.0 / SQRT2, -1.0 / SQRT2, 0.0, 0.0, 0.0],
        [0.0, 0.0, -1.0 / SQRT2, +1.0 / SQRT2, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0],
        [-1.0 / SQRT2, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 / SQRT2],
        [-1.0 / SQRT2, 0.0, 0.0, 0.0, 0.0, 0.0, +1.0 / SQRT2],
        [0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ],
    dtype=float,
)


HELICITY_LABELS = {
    "im_T0p_T1p": (
        r"$\frac{\mathcal{C}}{F_{UU}}\,\mathrm{Im}\sum_{\sigma}"
        r"T_{0,+}^{0\sigma}(T_{+1,+}^{0\sigma})^{*}$"
    ),
    "im_T0m_T1m": (
        r"$\frac{\mathcal{C}}{F_{UU}}\,\mathrm{Im}\sum_{\sigma}"
        r"T_{0,-}^{0\sigma}(T_{+1,-}^{0\sigma})^{*}$"
    ),
    "im_Tm1p_T1p": (
        r"$\frac{\mathcal{C}}{F_{UU}}\,\mathrm{Im}\sum_{\sigma}"
        r"T_{-1,+}^{0\sigma}(T_{+1,+}^{0\sigma})^{*}$"
    ),
    "delta_T1_strength": (
        r"$\frac{\mathcal{C}}{F_{UU}}\sum_{\sigma}"
        r"\left(|T_{+1,+}^{0\sigma}|^{2}-|T_{+1,-}^{0\sigma}|^{2}\right)$"
    ),
    "re_T0p_T1p": (
        r"$\frac{\mathcal{C}}{F_{UU}}\,\mathrm{Re}\sum_{\sigma}"
        r"T_{0,+}^{0\sigma}(T_{+1,+}^{0\sigma})^{*}$"
    ),
    "re_T0m_T1m": (
        r"$\frac{\mathcal{C}}{F_{UU}}\,\mathrm{Re}\sum_{\sigma}"
        r"T_{0,-}^{0\sigma}(T_{+1,-}^{0\sigma})^{*}$"
    ),
    "re_Tm1p_T1p": (
        r"$\frac{\mathcal{C}}{F_{UU}}\,\mathrm{Re}\sum_{\sigma}"
        r"T_{-1,+}^{0\sigma}(T_{+1,+}^{0\sigma})^{*}$"
    ),
}


HELICITY_DESCRIPTIONS = {
    "im_T0p_T1p": (
        "Imaginary longitudinal/transverse interference for initial-proton "
        "helicity +1/2, summed over recoil-neutron helicity."
    ),
    "im_T0m_T1m": (
        "Imaginary longitudinal/transverse interference for initial-proton "
        "helicity -1/2, summed over recoil-neutron helicity."
    ),
    "im_Tm1p_T1p": (
        "Imaginary interference of photon helicities -1 and +1 for "
        "initial-proton helicity +1/2."
    ),
    "delta_T1_strength": (
        "Difference of the photon-helicity +1 strengths for initial-proton "
        "helicities +1/2 and -1/2."
    ),
    "re_T0p_T1p": (
        "Real longitudinal/transverse interference for initial-proton "
        "helicity +1/2, summed over recoil-neutron helicity."
    ),
    "re_T0m_T1m": (
        "Real longitudinal/transverse interference for initial-proton "
        "helicity -1/2, summed over recoil-neutron helicity."
    ),
    "re_Tm1p_T1p": (
        "Real interference of photon helicities -1 and +1 for "
        "initial-proton helicity +1/2."
    ),
}


PUBLISHED_FIVE_PARAMETERS = (
    "im_T0p_T1p",
    "im_T0m_T1m",
    "im_Tm1p_T1p",
    "delta_T1_strength",
    "re_LT_difference",
)

PUBLISHED_FIVE_LABELS = {
    **{key: HELICITY_LABELS[key] for key in HELICITY_PARAMETERS[:4]},
    "re_LT_difference": (
        r"$\frac{\mathcal{C}}{F_{UU}}\,\mathrm{Re}\sum_{\sigma}\left["
        r"T_{0,+}^{0\sigma}(T_{+1,+}^{0\sigma})^{*}-"
        r"T_{0,-}^{0\sigma}(T_{+1,-}^{0\sigma})^{*}\right]$"
    ),
}


def transform_values(values: np.ndarray) -> np.ndarray:
    """Transform one or more source vectors; the final axis must have length 7."""
    values = np.asarray(values, dtype=float)
    if values.shape[-1] != len(SOURCE_PARAMETERS):
        raise ValueError("The final values axis must have length 7.")
    return np.einsum("ij,...j->...i", TRANSFORM, values)


def transform_covariance(covariance: np.ndarray) -> np.ndarray:
    """Apply the exact covariance transformation V_H = M V_F M^T."""
    covariance = np.asarray(covariance, dtype=float)
    if covariance.shape != (7, 7):
        raise ValueError("Source covariance must have shape (7, 7).")
    return TRANSFORM @ covariance @ TRANSFORM.T


def published_five_from_ratios(values: np.ndarray) -> np.ndarray:
    """Transform (lu1, ul1, ul2, ll0, ll1) into five publishable bilinears."""
    values = np.asarray(values, dtype=float)
    if values.shape[-1] != 5:
        raise ValueError("The final values axis must have length 5.")
    lu1, ul1, ul2, ll0, ll1 = np.moveaxis(values, -1, 0)
    return np.stack(
        (
            -(lu1 + ul1) / SQRT2,
            -(lu1 - ul1) / SQRT2,
            -ul2,
            2.0 * ll0,
            -SQRT2 * ll1,
        ),
        axis=-1,
    )


def published_five_diagonal_uncertainties(errors: np.ndarray) -> np.ndarray:
    """Propagate five marginal errors assuming their covariance is unavailable."""
    errors = np.asarray(errors, dtype=float)
    if errors.shape[-1] != 5:
        raise ValueError("The final error axis must have length 5.")
    lu1, ul1, ul2, ll0, ll1 = np.moveaxis(errors, -1, 0)
    return np.stack(
        (
            np.hypot(lu1, ul1) / SQRT2,
            np.hypot(lu1, ul1) / SQRT2,
            np.abs(ul2),
            2.0 * np.abs(ll0),
            SQRT2 * np.abs(ll1),
        ),
        axis=-1,
    )
