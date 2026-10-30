#!/usr/bin/env python3
"""Fa18 Out reconstructed |t_gamma| for FD versus FT photons.

Read DataEvents colXX branches produced by rebuild_pilot/convert_data_stream.cpp.
Angles in the ROOT tree are already in radians; no conversion is applied.
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import uproot

M_E = 0.00051099895  # GeV
BEAM_E = 10.604  # GeV, Fa18 Out
DEFAULT_DATA = Path('/volatile/clas12/thayward/dvcs_pass2_full_rebuild/data/fa18_out/dvcs_data/parts')
DEFAULT_OUTPUT = Path('/u/home/thayward/t_distribution.png')
BRANCHES = ['col09', 'col10', 'col11', 'col12', 'col18', 'col19', 'col20']


def calculate_abs_t(arrays):
    """Calculate -(k-k'-p_gamma)^2 using electron and photon four-vectors."""
    ep = np.asarray(arrays['col10'], dtype=np.float64)
    et = np.asarray(arrays['col11'], dtype=np.float64)  # radians in ROOT
    ef = np.asarray(arrays['col12'], dtype=np.float64)  # radians in ROOT
    gp = np.asarray(arrays['col18'], dtype=np.float64)
    gt = np.asarray(arrays['col19'], dtype=np.float64)  # radians in ROOT
    gf = np.asarray(arrays['col20'], dtype=np.float64)  # radians in ROOT
    beam_p = np.sqrt(BEAM_E**2 - M_E**2)
    electron_e = np.sqrt(ep**2 + M_E**2)
    q_e = BEAM_E - electron_e
    q_x = -ep * np.sin(et) * np.cos(ef)
    q_y = -ep * np.sin(et) * np.sin(ef)
    q_z = beam_p - ep * np.cos(et)
    dx = q_x - gp * np.sin(gt) * np.cos(gf)
    dy = q_y - gp * np.sin(gt) * np.sin(gf)
    dz = q_z - gp * np.cos(gt)
    t = (q_e - gp)**2 - dx**2 - dy**2 - dz**2
    return -t


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=DEFAULT_DATA)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--t-max', type=float, default=2.0)
    parser.add_argument('--bins', type=int, default=100)
    parser.add_argument('--step-size', default='100 MB')
    args = parser.parse_args()
    files = sorted(args.input_dir.glob('*.root'))
    # Exclude unfinished files, and require a successful marker if markers exist.
    files = [p for p in files if '.partial.' not in p.name and p.with_name(p.name + '.ok').exists()]
    if not files:
        raise SystemExit(f'No validated ROOT parts found in {args.input_dir}')
    edges = np.linspace(0.0, args.t_max, args.bins + 1)
    hist_fd = np.zeros(args.bins, dtype=np.int64)
    hist_ft = np.zeros(args.bins, dtype=np.int64)
    total_fd = total_ft = 0
    for index, path in enumerate(files, 1):
        with uproot.open(path) as root:
            if 'DataEvents' not in root:
                print(f'[SKIP] Missing DataEvents: {path}')
                continue
            #endif
            tree = root['DataEvents']
            for chunk in tree.iterate(BRANCHES, library='np', step_size=args.step_size):
                t_abs = calculate_abs_t(chunk)
                detector = np.asarray(chunk['col09'], dtype=np.int32)
                good = np.isfinite(t_abs) & (t_abs >= 0)
                fd = good & (detector == 1)
                ft = good & (detector == 0)
                hist_fd += np.histogram(t_abs[fd], bins=edges)[0]
                hist_ft += np.histogram(t_abs[ft], bins=edges)[0]
                total_fd += int(np.count_nonzero(fd))
                total_ft += int(np.count_nonzero(ft))
            #endfor
        print(f'[{index}/{len(files)}] {path.name}', flush=True)
    #endfor
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.stairs(hist_fd, edges, color='#0072B2', linewidth=1.8, label=f'FD photon (N={total_fd:,})')
    ax.stairs(hist_ft, edges, color='#D55E00', linewidth=1.8, label=f'FT photon (N={total_ft:,})')
    ax.set_xlabel(r'$|t_{\gamma}|$ (GeV$^2$)')
    ax.set_ylabel('Events / bin')
    ax.set_title(r'Fa18 Out (10.604 GeV): $|t_{\gamma}|$ from $q-p_{\gamma}$')
    ax.set_xlim(0, args.t_max)
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)
    print(f'[SAVED] {args.output} (FD={total_fd:,}, FT={total_ft:,})')


if __name__ == '__main__':
    main()
