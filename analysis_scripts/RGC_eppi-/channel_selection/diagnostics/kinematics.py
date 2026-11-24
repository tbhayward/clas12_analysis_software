"""Vectorized exclusive pion diagnostics. Four-vectors have order (E, px, py, pz)."""
import numpy as np
MN, MP, MPI, ME, MD = .9395654133, .9382720813, .13957039, .00051099895, 1.87561294257
XB = np.array([.10, .25, .35, .45, .60])
TP = np.array([.05, .25, .45, .65, .85, 1.05, 1.25])
BRANCHES = ['beam_energy', 'x', 'tprime', 'Q2', 'W', 'y', 'Mx2',
            'fiducial_status', 'detector1', 'detector2'] + [
                f'{p}_{v}' for p in ('e', 'p1', 'p2') for v in ('p', 'theta', 'phi')]
# key: axis label, lower edge, upper edge, reference position
METRICS = {
 'mx2_n': (r'$M_{X,n}^2$ (GeV$^2$)', -2, 2, 0),
 'mx2_d': (r'$M_{X,d}^2$ (GeV$^2$)', -1, 4, MP**2),
 'delta_e': (r'$E_\pi-E_{\pi,\mathrm{pred}}$ (GeV)', -1.5, 1.5, 0),
 'delta_p': (r'$|p_\pi|-|p_{\pi,\mathrm{pred}}|$ (GeV)', -1.5, 1.5, 0),
 'delta_px': (r'$\Delta p_x$ (GeV)', -1, 1, 0),
 'delta_py': (r'$\Delta p_y$ (GeV)', -1, 1, 0),
 'delta_pz': (r'$\Delta p_z$ (GeV)', -1.5, 1.5, 0),
 'delta_theta': (r'$\theta_\pi-\theta_{\pi,\mathrm{pred}}$ (deg)', -25, 25, 0),
 'delta_phi': (r'Wrapped $\varphi_\pi-\varphi_{\pi,\mathrm{pred}}$ (deg)', -40, 40, 0),
 'cone': (r'Pion opening angle to prediction (deg)', 0, 40, 0),
 'pred_mass': (r'$p_{\pi,\mathrm{pred}}^2-m_\pi^2$ (GeV$^2$)', -2, 2, 0),
 'ptmiss': (r'$p_{T,\mathrm{miss}}$ (GeV)', 0, 1.5, 0),
}
CORRELATIONS = [('mx2_n', 'delta_e'), ('delta_e', 'ptmiss'), ('delta_p', 'delta_theta')]

def track(a, prefix, mass):
    p, t, f = (a[f'{prefix}_{v}'] for v in ('p', 'theta', 'phi'))
    return np.column_stack((np.sqrt(p*p+mass*mass), p*np.sin(t)*np.cos(f),
                            p*np.sin(t)*np.sin(f), p*np.cos(t)))

def mass2(v):
    return v[:, 0]**2 - np.sum(v[:, 1:]**2, axis=1)

def calculate(a):
    e, p, pi = track(a, 'e', ME), track(a, 'p1', MP), track(a, 'p2', MPI)
    beam = np.zeros_like(e)
    beam[:, 0] = a['beam_energy']
    beam[:, 3] = np.sqrt(beam[:, 0]**2-ME**2)
    pred = beam-e-p
    pred[:, 0] += MN
    residual = pi-pred
    missing_d = -residual.copy()
    missing_d[:, 0] += MD-MN
    pmag = np.linalg.norm(pred[:, 1:], axis=1)
    pimag = np.linalg.norm(pi[:, 1:], axis=1)
    angle_valid = (pmag > 1e-12) & (pimag > 1e-12)
    with np.errstate(divide='ignore', invalid='ignore'):
        theta = np.arccos(np.clip(pred[:, 3]/pmag, -1, 1))
        phi = a['p2_phi']-np.arctan2(pred[:, 2], pred[:, 1])
        cone = np.arccos(np.clip(np.sum(pi[:, 1:]*pred[:, 1:], axis=1)/(pmag*pimag), -1, 1))
    phi = np.where((np.hypot(pred[:, 1], pred[:, 2]) > 1e-12) &
                   (np.hypot(pi[:, 1], pi[:, 2]) > 1e-12), np.arctan2(np.sin(phi), np.cos(phi)), np.nan)
    return dict(mx2_n=mass2(residual), mx2_d=mass2(missing_d), delta_e=residual[:, 0],
                delta_p=pimag-pmag, delta_px=residual[:, 1], delta_py=residual[:, 2],
                delta_pz=residual[:, 3], delta_theta=np.where(angle_valid, np.degrees(a['p2_theta']-theta), np.nan),
                delta_phi=np.degrees(phi), cone=np.where(angle_valid, np.degrees(cone), np.nan),
                pred_mass=mass2(pred)-MPI**2, ptmiss=np.hypot(residual[:, 1], residual[:, 2]))

def bin_indices(x, tp):
    """Same half-open bin edges as pi+ (including exclusion of the final upper edge)."""
    xi, ti = np.searchsorted(XB, x, side='right')-1, np.searchsorted(TP, tp, side='right')-1
    valid = np.isfinite(x) & np.isfinite(tp) & (xi >= 0) & (xi < 4) & (ti >= 0) & (ti < 6)
    return np.where(valid, xi*6+ti, -1)
