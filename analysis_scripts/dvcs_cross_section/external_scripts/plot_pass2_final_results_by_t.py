#!/usr/bin/env python3
"""Pass-2 DVCS final-results comparisons and kinematic diagnostics.

The script writes separate 10.2/10.6-GeV and BSA/cross-section trees.  For
each observable it produces the original Pass-2 value canvases, Pass-2 pulls
relative to KM15, and pull summaries grouped by |t|, Q2, and xB.  When
published Pass-1 inputs are available it also writes parallel value and pull
canvases, direct percentage-difference canvases and distributions, and marks
the published Pass-1 correlated-normalization band separately from the
point-to-point errors.  CSV-level kinematic and correction profiles are
written on the same |t|-grouped Q2/xB layout.  BH and KM15 curves remain on
all value canvases.

Place this file alongside plot_pass2_vs_pass1_model_comparison.py.  Published
Pass-1 inputs are auto-discovered in common import/imports directories or may
be supplied explicitly with --pass1-xs and --pass1-bsa-text.

Example:
  python3 plot_pass2_final_results_by_t.py imports/dvcs_pass2_analysis.csv \
    --output-dir output/final_results_by_t --phi-dense 24

Requires the existing dvcsgen and km15_cli.py installations used by the
comparison code.  A persistent JSON model cache avoids repeated evaluation.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import json
import math
import os
from pathlib import Path
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

# Avoid rewriting tracked .pyc files on the ifarm when importing local helpers.
sys.dont_write_bytecode = True
try:
    import plot_pass2_vs_pass1_model_comparison as base
except ModuleNotFoundError:
    helper_dir = Path(__file__).resolve().parent / "external_scripts"
    if not helper_dir.is_dir():
        raise
    sys.path.insert(0, str(helper_dir))
    import plot_pass2_vs_pass1_model_comparison as base


PASS1_NORMALIZATION_FRACTION = 0.31

# These are bin-level means/corrections already stored in the final CSV.  They
# are profiles, not substitutes for event-level distributions from the ROOT
# trees.  Keeping the list here makes the diagnostic output self-documenting.
DIAGNOSTIC_SPECS = (
    ("xBavg", r"Mean $x_B$", "mean_xB"),
    ("Q2avg", r"Mean $Q^2$ (GeV$^2$)", "mean_Q2"),
    ("t_abs_avg", r"Mean $|t|$ (GeV$^2$)", "mean_t"),
    ("e_theta", r"Mean $\theta_e$ (deg)", "mean_e_theta"),
    ("p_theta", r"Mean $\theta_p$ (deg)", "mean_p_theta"),
    ("g_theta", r"Mean $\theta_\gamma$ (deg)", "mean_g_theta"),
    ("Frad", r"$F_{\rm rad}$", "Frad"),
    ("Fbin", r"$F_{\rm bin}$", "Fbin"),
    ("phase_space_fraction", r"Physical / rectangular bin volume", "phase_space_fraction"),
)


def finite(value):
    return math.isfinite(value)


def energy_tag(energy):
    return f"{energy:.1f}_GeV"


def observable_tag(kind):
    return "sigma" if kind == "xs" else "BSA"


def model_cache_key(k, args, kind):
    energy = args.e_beam if args.e_beam is not None else args._active_energy
    parts = [
        k.xb_min, k.xb_max, k.q2_min, k.q2_max, k.t_min, k.t_max,
        energy, args.phi_dense, args.dvcsgen_dir, args.km15_cli,
        args.python_km15, kind,
    ]
    if getattr(args, "reference_kinematics", "means") == "means":
        parts += ["measured-means-v1", args._model_references[(energy,k,kind)]]
    return json.dumps(parts)


def compute_cell_models(job):
    """Compute the missing model observables for one kinematic cell."""
    if len(job) == 4:
        k, cfg, kinds, references = job
        result = {}
        requests = []
        for kind in kinds:
            for x,q,t,phi,e in references[kind]:
                for h in ([0] if kind == "xs" else [+1,-1]):requests.append((x,q,t,phi,h))
        predictions = iter(base.km15_xs_many(requests,cfg))
        for kind in kinds:
            points = references[kind]
            curve = {"phi": [p[3] for p in points], "km15": [], "bh": [], "reference": "measured_means"}
            for x,q,t,phi,e in points:
                if kind == "xs":
                    curve["km15"].append(next(predictions))
                    curve["bh"].append(base.bh_xs(x,q,t,phi,cfg))
                else:
                    plus,minus=next(predictions),next(predictions);denominator=plus+minus
                    curve["km15"].append((plus-minus)/denominator if finite(denominator) and denominator!=0 else float("nan"))
                    curve["bh"].append(0.)
            result[kind]=curve
        return k,result
    k, cfg, kinds = job
    phi = [360.0 * i / (cfg.phi_dense - 1) for i in range(cfg.phi_dense)]
    xb = 0.5 * (k.xb_min + k.xb_max)
    q2 = 0.5 * (k.q2_min + k.q2_max)
    t = 0.5 * (k.t_min + k.t_max)
    result = {}
    if "xs" in kinds:
        result["xs"] = {"phi": phi, "km15": [], "bh": []}
    if "bsa" in kinds:
        result["bsa"] = {"phi": phi, "km15": [], "bh": [0.0] * len(phi)}
    requests = [(xb,q2,t,angle,h) for angle in phi
                for h in (([0] if "xs" in kinds else []) +
                          ([+1,-1] if "bsa" in kinds else []))]
    predictions = iter(base.km15_xs_many(requests, cfg))
    for angle in phi:
        if "xs" in kinds:
            result["xs"]["bh"].append(base.bh_xs(xb, q2, t, angle, cfg))
            result["xs"]["km15"].append(
                next(predictions)
            )
        if "bsa" in kinds:
            plus = next(predictions)
            minus = next(predictions)
            denominator = plus + minus
            result["bsa"]["km15"].append(
                (plus - minus) / denominator
                if finite(denominator) and denominator != 0 else float("nan")
            )
    return k, result


def save_cache_atomic(path, cache):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(cache, separators=(",", ":")))
    os.replace(temporary, path)


def dvcs_tmin(q2, xb, proton_mass=0.9382720813):
    """Exact negative t_min for gamma* p -> gamma p."""
    eps2 = 4.0 * proton_mass**2 * xb**2 / q2
    return -q2 * (
        2.0 * (1.0 - xb) * (1.0 - math.sqrt(1.0 + eps2)) + eps2
    ) / (4.0 * xb * (1.0 - xb) + eps2)


def point_in_key(xb, q2, t_value, key, tolerance=1.0e-9):
    return (
        key.xb_min - tolerance <= xb <= key.xb_max + tolerance
        and key.q2_min - tolerance <= q2 <= key.q2_max + tolerance
        and key.t_min - tolerance <= t_value <= key.t_max + tolerance
    )


def best_key_for_point(keys, xb, q2, t_value):
    candidates = [k for k in keys if point_in_key(xb, q2, t_value, k)]
    if not candidates:
        return None
    return min(
        candidates,
        key=lambda k: (
            ((xb - 0.5 * (k.xb_min + k.xb_max)) / max(k.xb_max - k.xb_min, 1e-12)) ** 2
            + ((q2 - 0.5 * (k.q2_min + k.q2_max)) / max(k.q2_max - k.q2_min, 1e-12)) ** 2
            + ((t_value - 0.5 * (k.t_min + k.t_max)) / max(k.t_max - k.t_min, 1e-12)) ** 2
        ),
    )


def resolve_optional_input(explicit, candidates, label):
    if explicit is not None:
        path = explicit.expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"{label} does not exist: {path}")
        return path
    for candidate in candidates:
        candidate = candidate.expanduser().resolve()
        if candidate.exists():
            print(f"[pass1] Auto-detected {label}: {candidate}", flush=True)
            return candidate
    print(f"[pass1] No {label} found; the corresponding Pass-1 outputs will be skipped.", flush=True)
    return None


def load_published_pass1_xs(path, keys):
    """Load E214M1 or a loader-compatible CSV into the current Pass-2 cells."""
    out = {k: [] for k in keys}
    if path.suffix.lower() == ".csv":
        loaded = base.load_pass1_csv(
            path=path, no_point_to_point_systematics=False, print_columns=False
        )
        for key, points in loaded.items():
            if key not in out:
                continue
            for point in points:
                err = 0.5 * (abs(point.err_low) + abs(point.err_high))
                if finite(point.phi) and finite(point.xs) and finite(err) and err > 0:
                    out[key].append((point.phi % 360.0, point.xs, err))
        return {k: sorted(v) for k, v in out.items() if v}

    parsed = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.strip().replace("\t", " ").split()
            if len(fields) < 8:
                continue
            try:
                _, xb, q2, t_abs, phi, xs, stat, syst_total = map(float, fields[:8])
            except ValueError:
                continue
            if not all(finite(v) for v in (xb, q2, t_abs, phi, xs, stat, syst_total)):
                continue
            if xs <= 0 or stat < 0 or syst_total < 0:
                continue
            normalization = PASS1_NORMALIZATION_FRACTION * abs(xs)
            ptp = math.sqrt(max(syst_total**2 - normalization**2, 0.0))
            error = math.hypot(stat, ptp)
            key = best_key_for_point(keys, xb, q2, t_abs)
            if key is not None and error > 0:
                out[key].append((phi % 360.0, xs, error))
                parsed.append(1)
    print(
        f"[pass1] Loaded {len(parsed)} E214M1 cross-section points in "
        f"{sum(bool(v) for v in out.values())} Pass-2 cells.", flush=True
    )
    return {k: sorted(v) for k, v in out.items() if v}


def load_pass1_bsa(path, keys, energy):
    """Load the published seven-column RGA BSA table and assign by t'."""
    out = {k: [] for k in keys}
    kept = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.strip().split()
            if not fields or fields[0].startswith("#") or len(fields) < 7:
                continue
            try:
                phi_rad, q2, xb, t_negative, ebeam, asymmetry, error = map(float, fields[:7])
            except ValueError:
                continue
            values = (phi_rad, q2, xb, t_negative, ebeam, asymmetry, error)
            if not all(finite(v) for v in values) or error <= 0:
                continue
            if not math.isclose(ebeam, energy, abs_tol=0.05):
                continue
            t_prime = dvcs_tmin(q2, xb) - t_negative
            key = best_key_for_point(keys, xb, q2, t_prime)
            if key is None:
                continue
            out[key].append((math.degrees(phi_rad) % 360.0, asymmetry, abs(error)))
            kept += 1
    print(
        f"[pass1] Loaded {kept} BSA points at {energy:g} GeV in "
        f"{sum(bool(v) for v in out.values())} Pass-2 cells.", flush=True
    )
    return {k: sorted(v) for k, v in out.items() if v}


def interpolate_model(curves, model, phi):
    xx = np.asarray(curves["phi"], dtype=float)
    yy = np.asarray(curves[model], dtype=float)
    valid = np.isfinite(xx) & np.isfinite(yy)
    exact = valid & np.isclose(xx, phi % 360.0, rtol=0., atol=1e-9)
    if np.any(exact):return float(yy[np.flatnonzero(exact)[0]])
    if curves.get("reference") == "measured_means" or np.count_nonzero(valid) < 2:
        return float("nan")
    xx, yy = xx[valid], yy[valid]
    order = np.argsort(xx)
    return float(np.interp(phi % 360.0, xx[order], yy[order]))


def km15_pulls(data, curves):
    pulls = []
    for phi, value, error in data:
        model = interpolate_model(curves, "km15", phi)
        if finite(model) and finite(error) and error > 0:
            pulls.append((phi, (value - model) / error))
    return pulls


def angular_distance(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def nearest_phi_pairs(pass2_data, pass1_data, tolerance):
    """Return one-to-one nearest-phi Pass-2/Pass-1 point pairs."""
    available = list(range(len(pass1_data)))
    pairs = []
    for phi2, value2, error2 in sorted(pass2_data):
        if not available:
            break
        index = min(available, key=lambda i: angular_distance(phi2, pass1_data[i][0]))
        phi1, value1, error1 = pass1_data[index]
        if angular_distance(phi2, phi1) > tolerance:
            continue
        pairs.append((phi2, value2, error2, phi1, value1, error1))
        available.remove(index)
    return pairs


def pass1_pulls(pass2_data, pass1_data, tolerance):
    """One-to-one nearest-phi pulls with point-to-point errors in quadrature."""
    pulls = []
    for phi2, value2, error2, _, value1, error1 in nearest_phi_pairs(
        pass2_data, pass1_data, tolerance
    ):
        denominator = math.hypot(error2, error1)
        if finite(denominator) and denominator > 0:
            pulls.append((phi2, (value2 - value1) / denominator))
    return pulls


def pass1_percentage_differences(pass2_data, pass1_data, tolerance):
    """Return 100*(Pass-2/Pass-1 - 1) with point-to-point uncertainty.

    The published Pass-1 correlated normalization is deliberately excluded
    from each error bar and is displayed as a common shaded band instead.
    """
    differences = []
    for phi2, value2, error2, _, value1, error1 in nearest_phi_pairs(
        pass2_data, pass1_data, tolerance
    ):
        if not all(finite(v) for v in (value1, value2, error1, error2)) or value1 <= 0:
            continue
        ratio = value2 / value1
        relative_variance = 0.0
        if value2 != 0:
            relative_variance += (error2 / value2) ** 2
        relative_variance += (error1 / value1) ** 2
        error = 100.0 * abs(ratio) * math.sqrt(relative_variance)
        difference = 100.0 * (ratio - 1.0)
        if finite(difference) and finite(error):
            differences.append((phi2, difference, error))
    return differences


def key_range(key, axis):
    if axis == "t":
        return key.t_min, key.t_max
    if axis == "q2":
        return key.q2_min, key.q2_max
    if axis == "xb":
        return key.xb_min, key.xb_max
    raise ValueError(f"Unknown grouping axis: {axis}")


def group_values_by_axis(per_key, axis):
    grouped = {}
    for key, values in per_key.items():
        grouped.setdefault(key_range(key, axis), []).extend(values)
    return grouped


def axis_title(axis, limits):
    lo, hi = limits
    if axis == "t":
        return rf"$|t|=[{lo:.3f},{hi:.3f}]\ \mathrm{{GeV}}^2$"
    if axis == "q2":
        return rf"$Q^2=[{lo:.3f},{hi:.3f}]\ \mathrm{{GeV}}^2$"
    return rf"$x_B=[{lo:.3f},{hi:.3f}]$"


def axis_stem(axis, index, limits):
    lo, hi = limits
    return f"{axis}bin_{index:02d}_{lo:.3f}_{hi:.3f}"


def layout_for_group(keys, panels, kind):
    xbins = sorted(set((k.xb_min, k.xb_max) for k in keys))
    qbins = sorted(set((k.q2_min, k.q2_max) for k in keys))
    occupied = {
        (len(qbins) - 1 - iq, ix)
        for iq, qb in enumerate(qbins)
        for ix, xb in enumerate(xbins)
        if any(
            (k.xb_min, k.xb_max) == xb
            and (k.q2_min, k.q2_max) == qb
            and (
                panels[k].get(kind, [])
                if kind in ("xs", "bsa")
                else panels[k].get("diagnostics", {}).get(kind, [])
            )
            for k in keys
        )
    }
    bottom_by_col = {
        ix: max(row for row, col in occupied if col == ix)
        for ix in range(len(xbins))
        if any(col == ix for row, col in occupied)
    }
    left_by_row = {
        row: min(col for rr, col in occupied if rr == row)
        for row in range(len(qbins))
        if any(rr == row for rr, col in occupied)
    }
    return xbins, qbins, occupied, bottom_by_col, left_by_row


def panel_key(keys, xb, qb):
    matches = [
        k for k in keys
        if (k.xb_min, k.xb_max) == xb and (k.q2_min, k.q2_max) == qb
    ]
    return matches[0] if matches else None


def format_canvas_axes(ax, xb, qb, show_x, show_y, kind, pull=False):
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.grid(alpha=0.18)
    ax.text(
        0.035, 0.96,
        rf"$x_B=[{xb[0]:.3f},{xb[1]:.3f}]$" + "\n"
        + rf"$Q^2=[{qb[0]:.2f},{qb[1]:.2f}]\ \mathrm{{GeV}}^2$",
        transform=ax.transAxes, va="top", ha="left", fontsize=8,
        bbox={"facecolor": "white", "alpha": 0.78, "edgecolor": "none", "pad": 2},
    )
    ax.tick_params(axis="x", labelbottom=show_x, bottom=show_x)
    ax.tick_params(axis="y", labelleft=show_y, left=show_y)
    if show_x:
        ax.set_xlabel(r"$\phi$ (deg)")
    if show_y:
        if pull:
            ax.set_ylabel("Pull")
        else:
            ax.set_ylabel(r"$d^4\sigma$ (analysis units)" if kind == "xs" else r"$A_{LU}$")


def finish_endpoint_ticks(axs, bottom_by_col, left_by_row):
    for ix, row in bottom_by_col.items():
        ax = axs[row, ix]
        ax.tick_params(axis="x", labelbottom=True, bottom=True)
        labels = ["0" if ix == left_by_row[row] else "", "90", "180", "270", "360"]
        ax.set_xticks([0, 90, 180, 270, 360])
        ax.set_xticklabels(labels)
        tick_labels = ax.get_xticklabels()
        tick_labels[0].set_horizontalalignment("left")
        tick_labels[-1].set_horizontalalignment("right")


def set_cross_section_row_limits(axs, keys, panels, cache, args, kind,
                                 xbins, qbins, occupied, left_by_row,
                                 pass1=None):
    for row in range(len(qbins)):
        positives = []
        for col in range(len(xbins)):
            if (row, col) not in occupied:
                continue
            qb = qbins[len(qbins) - 1 - row]
            key = panel_key(keys, xbins[col], qb)
            if key is None or not panels[key][kind]:
                continue
            for _, value, error in panels[key][kind]:
                positives.extend(v for v in (value, value - error, value + error) if finite(v) and v > 0)
            if pass1:
                for _, value, error in pass1.get(key, []):
                    positives.extend(v for v in (value, value - error, value + error) if finite(v) and v > 0)
            curves = cache[model_cache_key(key, args, kind)]
            for model in ("km15", "bh"):
                positives.extend(float(v) for v in curves[model] if finite(float(v)) and float(v) > 0)
        if positives:
            lo, hi = min(positives), max(positives)
            padding = max(0.10, 0.06 * math.log10(hi / lo)) if hi > lo else 0.15
            limits = (10 ** (math.log10(lo) - padding), 10 ** (math.log10(hi) + padding))
            for col in range(len(xbins)):
                if (row, col) in occupied:
                    axs[row, col].set_ylim(*limits)
                    if col != left_by_row[row]:
                        axs[row, col].tick_params(
                            axis="y", which="both", left=False, right=False,
                            labelleft=False, labelright=False,
                        )
                        axs[row, col].yaxis.get_offset_text().set_visible(False)


def draw_value_canvas(keys, trange, kind, energy, panels, pass1, cache, args, output):
    xbins, qbins, occupied, bottom_by_col, left_by_row = layout_for_group(keys, panels, kind)
    figure_width=max(8.0,4.0*len(xbins))
    figure_height=max(4.8,2.85*len(qbins))
    fig, axs = plt.subplots(
        len(qbins), len(xbins), figsize=(figure_width,figure_height),
        squeeze=False, sharex=True, sharey=(kind == "bsa"),
    )
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.10, top=1.0-1.15/figure_height, hspace=0.0, wspace=0.0)
    data_label = "Pass-2 + published Pass-1" if pass1 is not None else "Pass-2"
    fig.suptitle(
        rf"CLAS12 {data_label} $ep\to ep\gamma$ — "
        + ("cross section" if kind == "xs" else "beam-spin asymmetry") + "\n"
        + rf"$|t|\in[{trange[0]:.3f},\,{trange[1]:.3f}]\ \mathrm{{GeV}}^2$, "
        + rf"$E_{{beam}}={energy:g}$ GeV"
        + ("\nTheory and data at measured mean kinematics" if getattr(args,"reference_kinematics","means")=="means" else "\nTheory at bin midpoints"),
        fontsize=13.5, y=1.0-.04/figure_height,
    )
    for iq, qb in enumerate(qbins):
        for ix, xb in enumerate(xbins):
            row = len(qbins) - 1 - iq
            ax = axs[row, ix]
            key = panel_key(keys, xb, qb)
            if key is None or not panels[key][kind]:
                ax.axis("off")
                continue
            phi, values, errors = np.asarray(sorted(panels[key][kind])).T
            ax.errorbar(
                phi, values, yerr=errors, fmt="o", ms=3.5, lw=0.9,
                capsize=1.5, color="black", label="CLAS12 Pass-2", zorder=4,
            )
            if pass1 is not None and pass1.get(key):
                pphi, pvalues, perrors = np.asarray(sorted(pass1[key])).T
                ax.errorbar(
                    pphi, pvalues, yerr=perrors, fmt="s", ms=3.2, lw=0.9,
                    capsize=1.5, color="#2ca02c", mfc="white",
                    label="CLAS12 published Pass-1", zorder=3,
                )
            curves = cache[model_cache_key(key, args, kind)]
            for model, color, name, linestyle in [
                ("km15", "#1f77b4", "KM15", "-"),
                ("bh", "#ff7f0e", "BH", "--"),
            ]:
                yy = np.asarray(curves[model], dtype=float)
                xx = np.asarray(curves["phi"], dtype=float)
                valid = np.isfinite(yy) & np.isfinite(xx)
                if np.any(valid):
                    order = np.argsort(xx[valid])
                    # Connect exact reference-point predictions for readability.
                    # Pulls still use the original model values at each mean.
                    ax.plot(xx[valid][order], yy[valid][order], linestyle,
                            lw=1.65, marker="." if valid.sum()==1 else None,
                            color=color, label=name, zorder=2)
            show_x = row == bottom_by_col.get(ix)
            show_y = ix == left_by_row.get(row)
            format_canvas_axes(ax, xb, qb, show_x, show_y, kind)
            if kind == "xs":
                ax.set_yscale("log")
            else:
                ax.set_ylim(-0.65, 0.65)
    if kind == "xs":
        set_cross_section_row_limits(
            axs, keys, panels, cache, args, kind, xbins, qbins,
            occupied, left_by_row, pass1,
        )
    finish_endpoint_ticks(axs, bottom_by_col, left_by_row)
    handles = [
        Line2D([], [], color="black", marker="o", ls="none", markersize=5, label="CLAS12 Pass-2"),
    ]
    if pass1 is not None:
        handles.append(Line2D([], [], color="#2ca02c", marker="s", mfc="white", ls="none", markersize=5, label="CLAS12 published Pass-1"))
    handles.extend([
        Line2D([], [], color="#1f77b4", lw=2, ls="-", label="KM15 at means" if args.reference_kinematics=="means" else "KM15"),
        Line2D([], [], color="#ff7f0e", lw=2, ls="--", label="BH at means" if args.reference_kinematics=="means" else "BH"),
    ])
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.0-.80/figure_height), ncol=len(handles), frameon=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=args.dpi)
    plt.close(fig)
    print(f"[write] {output}", flush=True)


def positive_finite(value):
    number=float(value)
    if not finite(number) or number<=0:
        raise argparse.ArgumentTypeError("must be a positive finite number")
    return number


def pull_display_points(pulls, limit):
    """Clip display only; preserve full pulls for summaries and distributions."""
    values=np.asarray(pulls,dtype=float)
    outside=np.abs(values)>limit
    display=np.where(outside,np.sign(values)*limit*.97,values)
    return display,outside


def draw_pull_canvas(keys, trange, kind, energy, panels, cache, args, output,
                     comparison, pass1=None):
    xbins, qbins, occupied, bottom_by_col, left_by_row = layout_for_group(keys, panels, kind)
    figure_height=max(4.8,2.85*len(qbins))
    fig, axs = plt.subplots(
        len(qbins), len(xbins), figsize=(max(8.,4.*len(xbins)),figure_height),
        squeeze=False, sharex=True, sharey=True,
    )
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.10, top=1.-1.05/figure_height, hspace=0.0, wspace=0.0)
    reference = "KM15" if comparison == "km15" else "published Pass-1"
    fig.suptitle(
        rf"CLAS12 Pass-2 pulls relative to {reference} — "
        + ("cross section" if kind == "xs" else "beam-spin asymmetry") + "\n"
        + rf"$|t|\in[{trange[0]:.3f},\,{trange[1]:.3f}]\ \mathrm{{GeV}}^2$, "
        + rf"$E_{{beam}}={energy:g}$ GeV",
        fontsize=13.5, y=1.-.04/figure_height,
    )
    all_pulls = []
    per_key = {}
    for key in keys:
        if comparison == "km15":
            values = km15_pulls(panels[key][kind], cache[model_cache_key(key, args, kind)])
        else:
            values = pass1_pulls(
                panels[key][kind], (pass1 or {}).get(key, []), args.pass1_phi_tolerance
            )
        per_key[key] = values
        all_pulls.extend(pull for _, pull in values)
    finite_pulls = [p for p in all_pulls if finite(p)]
    ylimit = getattr(args,"pull_limit",10.0)
    for iq, qb in enumerate(qbins):
        for ix, xb in enumerate(xbins):
            row = len(qbins) - 1 - iq
            ax = axs[row, ix]
            key = panel_key(keys, xb, qb)
            if key is None or not panels[key][kind]:
                ax.axis("off")
                continue
            values = per_key.get(key, [])
            if values:
                phi, pulls = np.asarray(values).T
                display,outside=pull_display_points(pulls,ylimit)
                ax.plot(phi[~outside],display[~outside],"o",ms=3.8,color="#7f3c8d",zorder=3)
                for sign,marker in ((1,"^"),(-1,"v")):
                    mask=outside & (np.sign(pulls)==sign)
                    ax.plot(phi[mask],display[mask],marker,ms=5,color="#c44e52",zorder=4)
                if outside.any():
                    ax.text(.96,.94,f"{outside.sum()} outside ±{ylimit:g}\nmax |pull|={max(abs(pulls)):g}",
                            transform=ax.transAxes,ha="right",va="top",fontsize=7,color="#a03030",
                            bbox=dict(facecolor="white",edgecolor="none",alpha=.9,pad=1))
            ax.axhline(0.0, color="black", lw=1.0, alpha=0.7)
            ax.axhline(+1.0, color="0.5", lw=0.8, ls="--", alpha=0.55)
            ax.axhline(-1.0, color="0.5", lw=0.8, ls="--", alpha=0.55)
            ax.set_ylim(-ylimit, ylimit)
            show_x = row == bottom_by_col.get(ix)
            show_y = ix == left_by_row.get(row)
            format_canvas_axes(ax, xb, qb, show_x, show_y, kind, pull=True)
            ax.text(
                0.96, 0.06, f"N={len(values)}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8,
            )
    finish_endpoint_ticks(axs, bottom_by_col, left_by_row)
    fig.legend(
        handles=[Line2D([], [], color="#7f3c8d", marker="o", ls="none", markersize=5, label=f"(Pass-2 − {reference}) / combined uncertainty" if comparison == "pass1" else "(Pass-2 − KM15) / Pass-2 uncertainty")],
        loc="upper center", bbox_to_anchor=(0.5, 1.-.63/figure_height), frameon=False,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=args.dpi)
    plt.close(fig)
    print(f"[write] {output}", flush=True)
    return finite_pulls


def draw_percentage_canvas(keys, trange, kind, energy, panels, pass1, args, output):
    """Draw direct Pass-2/Pass-1 percentage differences for one |t| canvas."""
    xbins, qbins, occupied, bottom_by_col, left_by_row = layout_for_group(
        keys, panels, kind
    )
    figure_height=max(4.8,2.85*len(qbins))
    fig, axs = plt.subplots(
        len(qbins), len(xbins), figsize=(max(8.,4.*len(xbins)),figure_height),
        squeeze=False, sharex=True, sharey=True,
    )
    fig.subplots_adjust(
        left=max(.075,.95/max(8.,4.*len(xbins))), right=0.985, bottom=0.10, top=1.-1.05/figure_height, hspace=0.0, wspace=0.0
    )
    fig.suptitle(
        r"CLAS12 $100\,(\mathrm{Pass\!\!-\!2}/\mathrm{Pass\!\!-\!1}-1)$ — "
        + ("cross section" if kind == "xs" else "beam-spin asymmetry") + "\n"
        + rf"$|t|\in[{trange[0]:.3f},\,{trange[1]:.3f}]\ \mathrm{{GeV}}^2$, "
        + rf"$E_{{beam}}={energy:g}$ GeV",
        fontsize=13.5, y=1.-.04/figure_height,
    )
    per_key = {
        key: pass1_percentage_differences(
            panels[key][kind], pass1.get(key, []), args.pass1_phi_tolerance
        )
        for key in keys
    }
    all_values = [value for values in per_key.values() for _, value, _ in values]
    # Fixed display range; full values remain in percentage summaries.
    ylimit = 100.0
    for iq, qb in enumerate(qbins):
        for ix, xb in enumerate(xbins):
            row = len(qbins) - 1 - iq
            ax = axs[row, ix]
            key = panel_key(keys, xb, qb)
            if key is None or not panels[key][kind]:
                ax.axis("off")
                continue
            values = per_key.get(key, [])
            ax.axhspan(
                -100.0 * PASS1_NORMALIZATION_FRACTION,
                +100.0 * PASS1_NORMALIZATION_FRACTION,
                color="#2ca02c", alpha=0.12, zorder=0,
            )
            if values:
                phi, differences, errors = np.asarray(values).T
                ax.errorbar(
                    phi, differences, yerr=errors, fmt="o", ms=3.6, lw=0.8,
                    capsize=1.4, color="#8c2d04", zorder=3,
                )
            ax.axhline(0.0, color="black", lw=1.0, alpha=0.75)
            ax.set_ylim(-ylimit, ylimit)
            show_x = row == bottom_by_col.get(ix)
            show_y = ix == left_by_row.get(row)
            format_canvas_axes(ax, xb, qb, show_x, show_y, kind)
            if show_y:
                ax.set_ylabel(r"$100\,(\mathrm{P2}/\mathrm{P1}-1)$ (%)")
            ax.text(
                0.96, 0.06, f"N={len(values)}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8,
            )
    finish_endpoint_ticks(axs, bottom_by_col, left_by_row)
    fig.legend(
        handles=[
            Line2D([], [], color="#8c2d04", marker="o", ls="none", markersize=5,
                   label="Pass-2 / Pass-1 (point-to-point errors)"),
            Line2D([], [], color="#2ca02c", lw=8, alpha=0.20,
                   label="Pass-1 ±31% correlated normalization"),
        ],
        loc="upper center", bbox_to_anchor=(0.5, 1.-.63/figure_height), ncol=2, frameon=False,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=args.dpi)
    plt.close(fig)
    print(f"[write] {output}", flush=True)
    return [value for value in all_values if finite(value)]


def draw_diagnostic_canvas(keys, trange, energy, panels, spec, args, output):
    """Draw a stored CSV bin-mean/correction profile on the standard layout."""
    name, ylabel, _ = spec
    xbins, qbins, occupied, bottom_by_col, left_by_row = layout_for_group(
        keys, panels, name
    )
    if not occupied:
        return
    fig, axs = plt.subplots(
        len(qbins), len(xbins), figsize=(4.0 * len(xbins), 2.85 * len(qbins)),
        squeeze=False, sharex=True, sharey=True,
    )
    fig.subplots_adjust(
        left=0.075, right=0.985, bottom=0.10, top=0.87, hspace=0.0, wspace=0.0
    )
    fig.suptitle(
        f"CLAS12 Pass-2 CSV profile — {ylabel}\n"
        + rf"$|t|\in[{trange[0]:.3f},\,{trange[1]:.3f}]\ \mathrm{{GeV}}^2$, "
        + rf"$E_{{beam}}={energy:g}$ GeV",
        fontsize=15, y=0.985,
    )
    for iq, qb in enumerate(qbins):
        for ix, xb in enumerate(xbins):
            row = len(qbins) - 1 - iq
            ax = axs[row, ix]
            key = panel_key(keys, xb, qb)
            values = (
                panels.get(key, {}).get("diagnostics", {}).get(name, [])
                if key is not None else []
            )
            if not values:
                ax.axis("off")
                continue
            phi, profile = np.asarray(sorted(values)).T
            ax.plot(phi, profile, "o-", ms=3.2, lw=0.9, color="#1b9e77")
            show_x = row == bottom_by_col.get(ix)
            show_y = ix == left_by_row.get(row)
            format_canvas_axes(ax, xb, qb, show_x, show_y, "xs")
            if show_y:
                ax.set_ylabel(ylabel)
    finish_endpoint_ticks(axs, bottom_by_col, left_by_row)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=args.dpi)
    plt.close(fig)
    print(f"[write] {output}", flush=True)


def draw_pull_histogram(pulls, title, output, dpi):
    values = np.asarray([v for v in pulls if finite(v)], dtype=float)
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    if values.size:
        spread = max(4.0, float(np.max(np.abs(values))) * 1.08)
        bins = np.linspace(-spread, spread, max(13, min(41, int(math.sqrt(values.size) * 3) + 1)))
        counts, edges, _ = ax.hist(values, bins=bins, color="#7f3c8d", alpha=0.70, edgecolor="white", label="Pulls")
        mu = float(np.mean(values))
        sigma = float(np.std(values, ddof=0))
        if values.size >= 2 and sigma > 0 and finite(sigma):
            xx = np.linspace(edges[0], edges[-1], 500)
            bin_width = float(np.mean(np.diff(edges)))
            yy = values.size * bin_width * np.exp(-0.5 * ((xx - mu) / sigma) ** 2) / (sigma * math.sqrt(2.0 * math.pi))
            ax.plot(xx, yy, color="#d62728", lw=2.0, label="Gaussian fit")
        ax.text(
            0.97, 0.95,
            rf"$N={values.size}$" + "\n" + rf"$\mu={mu:+.3f}$" + "\n" + rf"$\sigma={sigma:.3f}$",
            transform=ax.transAxes, ha="right", va="top",
            bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "0.75"},
        )
    else:
        ax.text(0.5, 0.5, "No matched pulls", transform=ax.transAxes, ha="center", va="center")
    ax.axvline(0.0, color="black", lw=1.0, alpha=0.7)
    ax.set_xlabel("Pull")
    ax.set_ylabel("Entries")
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.18)
    if values.size:
        ax.legend(frameon=False)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    print(f"[write] {output}", flush=True)


def draw_percentage_histogram(differences, title, output, dpi):
    """Summarize direct percentage differences without imposing a Gaussian."""
    values = np.asarray([v for v in differences if finite(v)], dtype=float)
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    normalization = 100.0 * PASS1_NORMALIZATION_FRACTION
    ax.axvspan(-normalization, normalization, color="#2ca02c", alpha=0.13,
               label="Pass-1 ±31% correlated normalization")
    if values.size:
        spread = max(50.0, float(np.max(np.abs(values))) * 1.08)
        bins = np.linspace(
            -spread, spread,
            max(13, min(41, int(math.sqrt(values.size) * 3) + 1)),
        )
        ax.hist(values, bins=bins, color="#8c2d04", alpha=0.72,
                edgecolor="white", label="Matched points")
        median = float(np.median(values))
        low, high = (float(v) for v in np.percentile(values, [16.0, 84.0]))
        within = float(np.mean(np.abs(values) <= normalization))
        ax.axvline(median, color="#d95f02", lw=2.0, label="Median")
        ax.text(
            0.97, 0.95,
            rf"$N={values.size}$" + "\n"
            + rf"median $={median:+.1f}\%$" + "\n"
            + rf"16–84% $=[{low:+.1f},{high:+.1f}]\%$" + "\n"
            + rf"within ±31% $={100.0 * within:.1f}\%$",
            transform=ax.transAxes, ha="right", va="top",
            bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "0.75"},
        )
    else:
        ax.text(0.5, 0.5, "No matched points", transform=ax.transAxes,
                ha="center", va="center")
    ax.axvline(0.0, color="black", lw=1.0, alpha=0.75)
    ax.set_xlabel(r"$100\,(\mathrm{Pass\!\!-\!2}/\mathrm{Pass\!\!-\!1}-1)$ (%)")
    ax.set_ylabel("Entries")
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.18)
    if values.size:
        ax.legend(frameon=False)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    print(f"[write] {output}", flush=True)


def csv_numeric_value(row, column):
    """Read either a scalar CSV cell or the central value of a tuple cell."""
    if column is None:
        return float("nan")
    raw = row.get(column, "")
    parsed = base.parse_tuple_cell(raw)
    if parsed.ok and finite(parsed.value):
        return parsed.value
    return base.to_float(raw, default=float("nan"))


def diagnostic_columns(fields, energy):
    suffix = "10.6 GeV" if math.isclose(energy, 10.6) else "Sp19 Inb"
    energy_label = "10.6 GeV" if math.isclose(energy, 10.6) else "10.2 GeV"
    columns = {}
    for name, _, _ in DIAGNOSTIC_SPECS:
        if name == "phase_space_fraction":
            columns[name] = (
                base.find_column(fields, [f"bin_volume, {energy_label}"], required=False),
                base.find_column(fields, [f"cubic bin_volume, {energy_label}"], required=False),
            )
        elif name in ("Frad", "Fbin"):
            columns[name] = base.find_column(
                fields, [f"{name}, {energy_label}"], required=False
            )
        else:
            columns[name] = base.find_column(
                fields, [f"{name}, {suffix}"], required=False
            )
    return columns


def parse_pass2(args):
    with args.csv.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        cols = base.detect_common_columns(fields, "Pass-2")
        xs_col = args.xs_column or base.find_column(
            fields, base.PASS2_CENTRAL_XS_CANDIDATES, required=False
        )
        xs_10p2_candidates = [
            "normed cross sections, ep->epg, exp, Sp19 Inb, unpol",
            "cross sections, ep->epg, exp, Sp19 Inb, unpol",
            "normed cross sections, ep->epg, exp, 10.2 GeV, unpol",
            "cross sections, ep->epg, exp, 10.2 GeV, unpol",
        ]
        xs_10p2 = args.xs_column_10p2 or base.find_column(
            fields, xs_10p2_candidates, required=False
        )
        if xs_10p2 is not None and xs_10p2 not in fields:
            raise ValueError(f"10.2 GeV XS column not found: {xs_10p2}")
        if xs_col is None and (args.e_beam is None or args.e_beam == 10.6):
            raise ValueError("No 10.6 GeV cross-section column found")
        if xs_10p2 is None and (args.e_beam is None or args.e_beam == 10.2):
            print(
                "[warning] No 10.2 GeV cross-section column detected; 10.2-GeV "
                "BSA will still be plotted. Use --xs-column-10p2 to specify it.",
                flush=True,
            )
        bsa_col = base.find_column(
            fields, [base.PASS2_BSA_TARGET_COLUMNS["10.6 GeV"]], required=False
        )
        bsa_10p2 = base.find_column(
            fields, [base.PASS2_BSA_TARGET_COLUMNS["Sp19 Inb"]], required=False
        )
        xs_ptp_col = base.find_column(
            fields, [base.PASS2_POINT_TO_POINT_TOTAL_COLUMN], required=False
        )
        bsa_ptp_col = base.find_column(
            fields,
            [
                "BSA, counts, 10.6 GeV, point-to-point total sys",
                "BSA, counts, 10.6 GeV, point-to-point sys",
                "BSA, counts, 10.6 GeV, total point-to-point sys",
            ],
            required=False,
        )
        if args.bsa_total_ptp and bsa_ptp_col is None:
            raise ValueError("No recognized BSA point-to-point systematic column")
        panels_by_energy = {
            energy: {} for energy in (10.2, 10.6)
            if args.e_beam is None or args.e_beam == energy
        }
        diagnostics_by_energy = {
            energy: diagnostic_columns(fields, energy) for energy in panels_by_energy
        }
        for row in reader:
            if not base.row_is_valid(row, cols["valid"]):
                continue
            key = base.make_bin_key(row, cols)
            phi = base.row_phi(row, cols)
            if not finite(phi):
                continue
            for energy, panels in panels_by_energy.items():
                item = panels.setdefault(
                    key, {"xs": [], "bsa": [], "diagnostics": {}, "references": {"xs": [], "bsa": []}}
                )
                plot_phi = phi
                reference = None
                if getattr(args,"reference_kinematics","means") == "means":
                    label = "10.6 GeV" if energy == 10.6 else "Sp19 Inb"
                    coordinates = tuple(base.to_float(row.get(p+", "+label,""),default=float("nan")) for p in ("xBavg","Q2avg","t_abs_avg","phiavg"))
                    reference = (*coordinates[:3],coordinates[3] % 360.,10.6 if energy==10.6 else 10.1998)
                    plot_phi = coordinates[3]

                for name, _, _ in DIAGNOSTIC_SPECS:
                    column = diagnostics_by_energy[energy].get(name)
                    if name == "phase_space_fraction":
                        physical = csv_numeric_value(row, column[0]) if column else float("nan")
                        rectangular = csv_numeric_value(row, column[1]) if column else float("nan")
                        diagnostic = (
                            physical / rectangular
                            if finite(physical) and finite(rectangular) and rectangular > 0
                            else float("nan")
                        )
                    else:
                        diagnostic = csv_numeric_value(row, column)
                    if finite(diagnostic):
                        item["diagnostics"].setdefault(name, []).append(
                            (phi % 360.0, diagnostic)
                        )
                selected_xs = xs_col if energy == 10.6 else xs_10p2
                selected_bsa = bsa_col if energy == 10.6 else bsa_10p2
                rowspecs = [
                    ("xs", selected_xs, xs_ptp_col if energy == 10.6 else None, args.xs_total_ptp),
                    ("bsa", selected_bsa, bsa_ptp_col if energy == 10.6 else None, args.bsa_total_ptp),
                ]
                for observable, column, systematics_column, include_systematics in rowspecs:
                    if column is None:
                        continue
                    value = base.parse_tuple_cell(row.get(column, ""))
                    if not value.ok or not finite(value.value) or not finite(value.stat) or value.stat <= 0:
                        continue
                    if observable == "xs" and value.value <= 0:
                        continue
                    error = value.stat
                    if include_systematics and systematics_column:
                        systematic = base.to_float(row.get(systematics_column, ""), default=0.0)
                        if finite(systematic):
                            error = math.hypot(error, systematic)
                    if reference is not None:
                        if not all(finite(v) for v in reference):raise ValueError(f"Missing measured means for {key}, {energy:g} GeV")
                        item["references"][observable].append(reference)
                    item[observable].append((plot_phi % 360.0, value.value, error))
    selected = {
        energy: {k: v for k, v in panels.items() if v["xs"] or v["bsa"]}
        for energy, panels in panels_by_energy.items()
    }
    args._model_references = {(e,k,kind): tuple(sorted(v["references"][kind],key=lambda p:p[3]))
        for e,panels in selected.items() for k,v in panels.items() for kind in ("xs","bsa")}
    return selected


def make_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="Pass-2 extracted 4D CSV")
    parser.add_argument("--output-dir", type=Path, default=Path("output/final_results_by_t"))
    parser.add_argument("--pass1-xs", type=Path, default=None, help="Published E214M1 text table or loader-compatible Pass-1 CSV")
    parser.add_argument("--pass1-bsa-text", type=Path, default=None, help="Published seven-column RGA Pass-1 BSA table")
    parser.add_argument("--no-pass1", action="store_true", help="Do not auto-detect or plot any Pass-1 data")
    parser.add_argument("--pass1-phi-tolerance", type=float, default=2.0, help="Maximum circular phi separation in degrees for Pass-2/Pass-1 pulls (default: 2)")
    parser.add_argument("--reference-kinematics",choices=("midpoints","means"),default="means",
                        help="Theory/data reference: measured row means (default, exact model points) or legacy bin midpoints")
    parser.add_argument("--pull-limit",type=positive_finite,default=10.0,help="Symmetric pull-panel range (default: 10); outliers marked at boundary; statistics unchanged")
    parser.add_argument("--phi-dense", type=int, default=24, help="Phi points per model curve (default: 24)")
    parser.add_argument("--e-beam", type=float, choices=(10.2, 10.6), default=None, help="Restrict to one beam energy (default: both)")
    parser.add_argument("--xs-column-10p2", default=None, help="Override the Sp19 10.2-GeV cross-section column")
    parser.add_argument("--dvcsgen-dir", default=base.DEFAULT_DVCSGEN_DIR)
    parser.add_argument("--km15-cli", default=base.DEFAULT_KM15_CLI)
    parser.add_argument("--python-km15", default=sys.executable)
    parser.add_argument("--xs-column", default=None)
    parser.add_argument("--xs-total-ptp", action="store_true", help="XS error bars/pulls: stat and point-to-point systematics in quadrature")
    parser.add_argument("--bsa-total-ptp", action="store_true", help="BSA error bars/pulls: stat and point-to-point systematics in quadrature")
    parser.add_argument("--recompute-models", action="store_true")
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--workers", type=int, default=10, help="Parallel model workers (1–10; default 10)")
    return parser


def main():
    args = make_parser().parse_args()
    if not 1 <= args.workers <= 10:
        raise SystemExit("--workers must be between 1 and 10")
    if args.phi_dense < 2:
        raise SystemExit("--phi-dense must be >= 2")
    if args.pass1_phi_tolerance < 0:
        raise SystemExit("--pass1-phi-tolerance must be nonnegative")
    args.csv = args.csv.expanduser().resolve()
    if not args.csv.exists():
        raise FileNotFoundError(f"Pass-2 CSV does not exist: {args.csv}")
    args.output_dir = args.output_dir.expanduser().resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    panels_by_energy = parse_pass2(args)
    if not any(panels_by_energy.values()):
        raise RuntimeError("No valid Pass-2 results found in CSV")

    script_dir = Path(__file__).resolve().parent
    search_roots = [script_dir, script_dir.parent, args.csv.parent, args.csv.parent.parent]
    pass1_xs_path = None
    pass1_bsa_path = None
    if not args.no_pass1:
        pass1_xs_path = resolve_optional_input(
            args.pass1_xs,
            [root / folder / "clasdb_E214M1.txt" for root in search_roots for folder in ("import", "imports")],
            "published Pass-1 E214M1 cross-section table",
        )
        pass1_bsa_path = resolve_optional_input(
            args.pass1_bsa_text,
            [root / folder / name for root in search_roots for folder in ("import", "imports") for name in ("rga_pass1_bsa.txt", "RGA_pass1_BSA.txt")],
            "published Pass-1 BSA table",
        )

    for energy, panels in panels_by_energy.items():
        if not panels:
            continue
        args._active_energy = energy
        energy_dir = args.output_dir / energy_tag(energy)
        energy_dir.mkdir(parents=True, exist_ok=True)
        cfg = base.ModelConfig(
            e_beam=(10.1998 if energy==10.2 and args.reference_kinematics=="means" else energy), dvcsgen_dir=args.dvcsgen_dir,
            km15_cli=args.km15_cli, py_km15=args.python_km15,
            phi_dense=args.phi_dense, use_cache=True,
            allow_missing_models=False, verbose_worker_models=False,
        )
        cache_path = energy_dir / f"model_curves_cache_{energy:g}GeV.json"
        cache = (
            json.loads(cache_path.read_text())
            if cache_path.exists() and not args.recompute_models else {}
        )
        groups = {}
        for key in panels:
            groups.setdefault((key.t_min, key.t_max), []).append(key)
        print(
            f"[results] {energy:g} GeV: {len(panels)} populated cells; "
            f"{len(groups)} |t| bins", flush=True,
        )
        jobs = []
        for key in panels:
            kinds = tuple(
                kind for kind in ("xs", "bsa")
                if panels[key][kind] and model_cache_key(key, args, kind) not in cache
            )
            if kinds:
                if args.reference_kinematics == "means":
                    jobs.append((key,cfg,kinds,{kind:args._model_references[(energy,key,kind)] for kind in kinds}))
                else:
                    jobs.append((key, cfg, kinds))
        if jobs:
            print(f"[models] Computing {len(jobs)} cells with up to {args.workers} workers", flush=True)
            with ProcessPoolExecutor(max_workers=min(args.workers, len(jobs))) as pool:
                futures = [pool.submit(compute_cell_models, job) for job in jobs]
                try:
                    for index, future in enumerate(as_completed(futures), 1):
                        key, curves = future.result()
                        for kind, entry in curves.items():
                            cache[model_cache_key(key, args, kind)] = entry
                        if index % 8 == 0 or index == len(jobs):
                            save_cache_atomic(cache_path, cache)
                        print(f"[models] {index}/{len(jobs)} cells complete", flush=True)
                except BaseException:
                    save_cache_atomic(cache_path, cache)
                    raise
        else:
            print("[models] All required curves loaded from cache", flush=True)

        keys = list(panels)
        pass1_by_kind = {"xs": {}, "bsa": {}}
        # Published E214M1 is the 10.6-GeV cross-section dataset.
        if pass1_xs_path is not None and math.isclose(energy, 10.6):
            pass1_by_kind["xs"] = load_published_pass1_xs(pass1_xs_path, keys)
        if pass1_bsa_path is not None:
            pass1_by_kind["bsa"] = load_pass1_bsa(pass1_bsa_path, keys, energy)

        # These profiles expose where each phi-bin population sits inside the
        # nominal 4D bin and how the model/phase-space corrections vary.  They
        # use only values already stored in the final CSV.
        for it, (trange, group_keys) in enumerate(sorted(groups.items()), 1):
            stem = f"tbin_{it:02d}_{trange[0]:.3f}_{trange[1]:.3f}"
            for spec in DIAGNOSTIC_SPECS:
                name, _, file_tag = spec
                if not any(
                    panels[key].get("diagnostics", {}).get(name)
                    for key in group_keys
                ):
                    continue
                draw_diagnostic_canvas(
                    group_keys, trange, energy, panels, spec, args,
                    energy_dir / "kinematic_profiles" / "by_t" / file_tag
                    / f"{file_tag}_{energy:g}GeV_{stem}.png",
                )

        for kind in ("xs", "bsa"):
            if not any(panels[key][kind] for key in panels):
                continue
            root = energy_dir / observable_tag(kind)
            km15_overall = []
            pass1_overall = []
            percentage_overall = []
            km15_by_key = {
                key: [pull for _, pull in km15_pulls(
                    panels[key][kind], cache[model_cache_key(key, args, kind)]
                )]
                for key in panels if panels[key][kind]
            }
            p1 = pass1_by_kind[kind]
            pass1_by_key = {
                key: [pull for _, pull in pass1_pulls(
                    panels[key][kind], p1.get(key, []), args.pass1_phi_tolerance
                )]
                for key in panels if panels[key][kind]
            } if p1 else {}
            percentage_by_key = {
                key: [difference for _, difference, _ in pass1_percentage_differences(
                    panels[key][kind], p1.get(key, []), args.pass1_phi_tolerance
                )]
                for key in panels if panels[key][kind]
            } if p1 else {}
            for it, (trange, group_keys) in enumerate(sorted(groups.items()), 1):
                group_keys = [key for key in group_keys if panels[key][kind]]
                if not group_keys:
                    continue
                stem = f"tbin_{it:02d}_{trange[0]:.3f}_{trange[1]:.3f}"
                km15_root = root / "pass2_vs_KM15"
                draw_value_canvas(
                    group_keys, trange, kind, energy, panels, None, cache, args,
                    km15_root / "values" / f"pass2_{energy:g}GeV_{kind}_{stem}.png",
                )
                pulls = draw_pull_canvas(
                    group_keys, trange, kind, energy, panels, cache, args,
                    km15_root / "pulls" / f"pass2_vs_KM15_{energy:g}GeV_{kind}_{stem}.png",
                    comparison="km15",
                )
                km15_overall.extend(pulls)
                draw_pull_histogram(
                    pulls,
                    f"Pass-2 vs KM15 pulls — {energy:g} GeV {observable_tag(kind)}, "
                    + rf"$|t|=[{trange[0]:.3f},{trange[1]:.3f}]$",
                    km15_root / "pull_distributions" / f"pull_distribution_{stem}.png",
                    args.dpi,
                )

                if p1 and any(p1.get(key) for key in group_keys):
                    pass1_root = root / "pass2_vs_pass1"
                    draw_value_canvas(
                        group_keys, trange, kind, energy, panels, p1, cache, args,
                        pass1_root / "values" / f"pass2_pass1_{energy:g}GeV_{kind}_{stem}.png",
                    )
                    p1_pulls = draw_pull_canvas(
                        group_keys, trange, kind, energy, panels, cache, args,
                        pass1_root / "pulls" / f"pass2_vs_pass1_{energy:g}GeV_{kind}_{stem}.png",
                        comparison="pass1", pass1=p1,
                    )
                    pass1_overall.extend(p1_pulls)
                    draw_pull_histogram(
                        p1_pulls,
                        f"Pass-2 vs Pass-1 pulls — {energy:g} GeV {observable_tag(kind)}, "
                        + rf"$|t|=[{trange[0]:.3f},{trange[1]:.3f}]$",
                        pass1_root / "pull_distributions" / f"pull_distribution_{stem}.png",
                        args.dpi,
                    )
                    percentage_values = draw_percentage_canvas(
                        group_keys, trange, kind, energy, panels, p1, args,
                        pass1_root / "percentage_differences" / "canvases_by_t"
                        / f"pass2_over_pass1_percent_{energy:g}GeV_{kind}_{stem}.png",
                    )
                    percentage_overall.extend(percentage_values)
                    draw_percentage_histogram(
                        percentage_values,
                        f"Pass-2 / Pass-1 percentage differences — {energy:g} GeV "
                        f"{observable_tag(kind)}, "
                        + rf"$|t|=[{trange[0]:.3f},{trange[1]:.3f}]$",
                        pass1_root / "percentage_differences" / "distributions" / "by_t"
                        / f"percentage_difference_{stem}.png",
                        args.dpi,
                    )

            # The original summaries were organized only by |t|.  Add the
            # orthogonal Q2 and xB groupings without changing the legacy paths.
            for axis in ("q2", "xb"):
                for index, (limits, values) in enumerate(
                    sorted(group_values_by_axis(km15_by_key, axis).items()), 1
                ):
                    draw_pull_histogram(
                        values,
                        f"Pass-2 vs KM15 pulls — {energy:g} GeV {observable_tag(kind)}, "
                        + axis_title(axis, limits),
                        root / "pass2_vs_KM15" / "pull_distributions" / f"by_{axis}"
                        / f"pull_distribution_{axis_stem(axis, index, limits)}.png",
                        args.dpi,
                    )
                if pass1_by_key:
                    for index, (limits, values) in enumerate(
                        sorted(group_values_by_axis(pass1_by_key, axis).items()), 1
                    ):
                        draw_pull_histogram(
                            values,
                            f"Pass-2 vs Pass-1 pulls — {energy:g} GeV {observable_tag(kind)}, "
                            + axis_title(axis, limits),
                            root / "pass2_vs_pass1" / "pull_distributions" / f"by_{axis}"
                            / f"pull_distribution_{axis_stem(axis, index, limits)}.png",
                            args.dpi,
                        )
                if percentage_by_key:
                    for index, (limits, values) in enumerate(
                        sorted(group_values_by_axis(percentage_by_key, axis).items()), 1
                    ):
                        draw_percentage_histogram(
                            values,
                            f"Pass-2 / Pass-1 percentage differences — {energy:g} GeV "
                            f"{observable_tag(kind)}, " + axis_title(axis, limits),
                            root / "pass2_vs_pass1" / "percentage_differences"
                            / "distributions" / f"by_{axis}"
                            / f"percentage_difference_{axis_stem(axis, index, limits)}.png",
                            args.dpi,
                        )
            draw_pull_histogram(
                km15_overall,
                f"Overall Pass-2 vs KM15 pulls — {energy:g} GeV {observable_tag(kind)}",
                root / "pass2_vs_KM15" / "pull_distributions" / "pull_distribution_overall.png",
                args.dpi,
            )
            if pass1_overall:
                draw_pull_histogram(
                    pass1_overall,
                    f"Overall Pass-2 vs Pass-1 pulls — {energy:g} GeV {observable_tag(kind)}",
                    root / "pass2_vs_pass1" / "pull_distributions" / "pull_distribution_overall.png",
                    args.dpi,
                )
                draw_percentage_histogram(
                    percentage_overall,
                    f"Overall Pass-2 / Pass-1 percentage differences — {energy:g} GeV "
                    f"{observable_tag(kind)}",
                    root / "pass2_vs_pass1" / "percentage_differences"
                    / "distributions" / "percentage_difference_overall.png",
                    args.dpi,
                )
            elif pass1_by_kind[kind]:
                warnings.warn(
                    f"Pass-1 {observable_tag(kind)} data were loaded at {energy:g} GeV, "
                    "but no phi-matched Pass-2/Pass-1 pulls were found."
                )
    print("[done] All requested PNG canvases and pull summaries generated.", flush=True)


if __name__ == "__main__":
    main()
