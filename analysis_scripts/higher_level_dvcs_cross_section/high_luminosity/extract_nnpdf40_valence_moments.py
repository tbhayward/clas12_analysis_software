#!/usr/bin/env python3
"""Extract valence A20 moments from the official NNPDF4.0 LHAPDF grid archive.

No LHAPDF Python installation is required. The script reads the .tar.gz directly.
It deliberately refuses Q values outside the tabulated grid rather than silently
extrapolating PDFs.
"""
from __future__ import annotations
import argparse, io, math, re, tarfile
from pathlib import Path
import numpy as np

DEFAULT_ARCHIVE = Path(__file__).resolve().parent / "import" / "NNPDF40_nnlo_as_01180.tar.gz"
DEFAULT_OUT = Path(__file__).resolve().parent / "output" / "stage5_ji" / "tables"


def _blocks(raw: str):
    # LHAPDF6 .dat files: metadata, then one or more --- delimited grid blocks.
    parts = raw.split("---")
    for part in parts[1:]:
        lines = [x.strip() for x in part.splitlines() if x.strip()]
        if len(lines) < 4:
            continue
        try:
            x = np.fromstring(lines[0], sep=" ")
            q = np.fromstring(lines[1], sep=" ")
            flav = np.fromstring(lines[2], sep=" ", dtype=int)
        except Exception:
            continue
        if x.size < 2 or q.size < 2 or flav.size < 2:
            continue
        vals = np.fromstring("\n".join(lines[3:]), sep=" ")
        need = x.size * q.size * flav.size
        if vals.size < need:
            continue
        vals = vals[:need].reshape(x.size, q.size, flav.size)  # Q varies fastest
        yield x, q, flav, vals


def _read_raw(raw: str, name: str, qtarget: float):
    grids = list(_blocks(raw))
    if not grids:
        raise RuntimeError(f"No LHAPDF grid found in {name}")
    # Select a subgrid that actually contains Q. Never extrapolate.
    candidates = [g for g in grids if g[1][0] <= qtarget <= g[1][-1]]
    if not candidates:
        qmin = min(g[1][0] for g in grids); qmax = max(g[1][-1] for g in grids)
        raise ValueError(f"Requested Q={qtarget:.8g} GeV is outside tabulated range [{qmin:.8g}, {qmax:.8g}] GeV")
    x, q, flav, vals = candidates[0]
    # LHAPDF grids store x*f(x,Q). Interpolate in log(Q), as appropriate for grid lookup.
    lq = np.log(q); lt = math.log(qtarget)
    j = int(np.searchsorted(lq, lt))
    if j == 0:
        xf = vals[:, 0, :]
    elif j == len(q):
        xf = vals[:, -1, :]
    else:
        w = (lt-lq[j-1])/(lq[j]-lq[j-1])
        xf = (1-w)*vals[:, j-1, :] + w*vals[:, j, :]
    idx = {int(f): i for i, f in enumerate(flav)}
    for f in (-2, -1, 1, 2):
        if f not in idx:
            raise RuntimeError(f"Required flavor {f} absent from {name}")
    xuv = xf[:, idx[2]] - xf[:, idx[-2]]
    xdv = xf[:, idx[1]] - xf[:, idx[-1]]
    # A20 = integral dx [x q_v]. Valence number = integral dx [x q_v]/x.
    a_u = float(np.trapezoid(xuv, x))
    a_d = float(np.trapezoid(xdv, x))
    n_u = float(np.trapezoid(xuv/x, x))
    n_d = float(np.trapezoid(xdv/x, x))
    return a_u, a_d, n_u, n_d, float(x[0]), float(x[-1]), float(q[0]), float(q[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    ap.add_argument("--q2", type=float, default=3.0, help="Q^2 in GeV^2 (default 3.0; must lie inside the grid)")
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    qtarget = math.sqrt(args.q2)
    if args.q2 <= 0: raise SystemExit("--q2 must be positive")
    if not args.archive.exists(): raise SystemExit(f"Missing archive: {args.archive}")
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    probe = None
    with tarfile.open(args.archive, "r:gz") as tf:
        for member in tf:
            if not member.isfile() or not re.search(r"_\d{4}\.dat$", member.name) or member.name.endswith("_0000.dat"):
                continue
            raw = tf.extractfile(member).read().decode("utf-8", errors="replace")
            try:
                result = _read_raw(raw, member.name, qtarget)
            except ValueError as e:
                raise SystemExit(f"ERROR: {e}\nThis script intentionally refuses PDF extrapolation. Choose a Q^2 inside the official grid.")
            if probe is None:
                probe = result
            a_u,a_d,n_u,n_d,*_ = result
            rid = int(re.search(r"_(\d{4})\.dat$", member.name).group(1))
            rows.append((rid, a_u,a_d,n_u,n_d))
    rows.sort(key=lambda r: r[0])
    if not rows:
        raise SystemExit("No replica members found")

    arr = np.asarray(rows, float)
    au, ad = arr[:,1], arr[:,2]
    cov = np.cov(np.vstack([au,ad]), ddof=1)
    corr = cov[0,1]/math.sqrt(cov[0,0]*cov[1,1])
    nu, nd = arr[:,3], arr[:,4]

    tag = f"Q2_{args.q2:g}".replace(".","p")
    repfile = args.outdir / f"nnpdf40_valence_A20_replicas_{tag}.csv"
    np.savetxt(repfile, arr, delimiter=",", header="replica,A20_uv,A20_dv,valence_number_u,valence_number_d", comments="", fmt=["%d","%.10g","%.10g","%.10g","%.10g"])
    summary = args.outdir / f"nnpdf40_valence_A20_summary_{tag}.csv"
    with summary.open("w") as f:
        f.write("quantity,value\n")
        for key,val in [
            ("pdf_set","NNPDF40_nnlo_as_01180"),("Q2_GeV2",args.q2),("Q_GeV",qtarget),("n_replicas",len(rows)),
            ("A20_uv_mean",au.mean()),("A20_uv_std",au.std(ddof=1)),("A20_dv_mean",ad.mean()),("A20_dv_std",ad.std(ddof=1)),
            ("cov_uv_dv",cov[0,1]),("corr_uv_dv",corr),
            ("valence_u_mean",nu.mean()),("valence_u_max_abs_dev_from_2",np.max(np.abs(nu-2))),
            ("valence_d_mean",nd.mean()),("valence_d_max_abs_dev_from_1",np.max(np.abs(nd-1))),
            ("x_min",probe[4]),("x_max",probe[5]),("selected_subgrid_Q_min_GeV",probe[6]),("selected_subgrid_Q_max_GeV",probe[7])]:
            f.write(f"{key},{val}\n")

    print("NNPDF4.0 valence-moment extraction")
    print(f"  set: NNPDF40_nnlo_as_01180")
    print(f"  Q^2 = {args.q2:g} GeV^2  (Q = {qtarget:.6f} GeV)")
    print(f"  replicas: {len(rows)}")
    print(f"  A20_uv = {au.mean():.6f} +/- {au.std(ddof=1):.6f}")
    print(f"  A20_dv = {ad.mean():.6f} +/- {ad.std(ddof=1):.6f}")
    print(f"  corr(A20_uv,A20_dv) = {corr:.6f}")
    print(f"  valence-number check: u_v = {nu.mean():.6f} (target 2), d_v = {nd.mean():.6f} (target 1)")
    print(f"  wrote: {repfile}")
    print(f"  wrote: {summary}")

if __name__ == "__main__": main()
