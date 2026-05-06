#ifndef PI0_CORRECTED_COUNTS_H
#define PI0_CORRECTED_COUNTS_H

#include <string>

/**
 * Update the pass-2 CSV with pi0-corrected DVCS signal yields and
 * produce per-period yield plots (pos/neg helicities versus phi).
 *
 * CSV behavior:
 *   For each row, each period, and each helicity {unpol,pos,neg}:
 *     - Sum the normalized normalized raw yields over all three topologies.
 *     - Read the contamination ratio triple
 *           "contamination ratio, <period>" = (c, c_stat, c_sys)
 *       If this cell is empty, assume c = 0 and c_stat = c_sys = 0.
 *     - Compute the signal yield S and its statistical error:
 *           N_norm = sum_topologies normalized raw yield
 *           S      = (1 - c) * N_norm
 *           Var(S) = (1 - c)^2 * N_norm + N_norm^2 * c_stat^2
 *     - Write the result to the column
 *           "signal yield, ep->epg, exp, <period>, <hel>"
 *       as a triple "(S, S_stat, S_sys)" with S_sys = 0 for now.
 *
 * Plot behavior:
 *   - For each period (Fa18 Inb, Fa18 Out, Sp18 Inb, Sp18 Out, Sp19 Inb),
 *     build a grid in xB-Q2-t, and in each cell plot:
 *       - + helicity (red open circles) and
 *       - - helicity (blue solid circles)
 *     vs phi (degrees), using the pi0-corrected signal yields and their
 *     propagated statistical errors.
 *   - Plots are saved to:
 *         <out_root_dir>/signal_yield_plots/<PeriodDir>/
 *         plot_signal_yield_<PeriodDir>_xB_<idx>.png
 *     where PeriodDir is the canonical directory name (e.g. Fa18_Inb)
 *     and idx = round(1000 * xBmin).
 *
 * Returns:
 *   true on success, false on failure (with diagnostics to stderr).
 */
bool update_pi0_corrected_counts_csv(const std::string& csv_path,
                                     const std::string& out_root_dir);

#endif // PI0_CORRECTED_COUNTS_H