#!/usr/bin/env python3
"""
Validate the AAOgen/VPK pi0 structure-function model against the CLAS6
Bedlinskiy et al. 5.75-GeV data archived in HEPData (record ins1294143).

Physics purpose
---------------
The VPK parameters used by AAOgen were fitted to the CLAS 6-GeV pi0
structure-function data.  This script is therefore primarily a fit-domain
closure/reproduction test, with special attention to Q2 ~ 1-3 GeV^2 where
the RGA pi0 normalization has most of its statistical weight.

The script:
  * downloads the official HEPData YAML archive if it is not already cached;
  * parses all Bedlinskiy structure-function tables;
  * evaluates the same VPK XSIGMA_U, XSIGMA_TT and XSIGMA_LT functions used
    by validate_aao_pi0_model.py at each published point;
  * writes point-level and setting-level data/VPK comparisons;
  * makes PNG diagnostics, including a combined CLAS6 + Hall-A Q2 plot when
    the Hall-A output already exists.

Run from external_pi0_normalization:
    python3 validate_clas6_bedlinskiy_vpk.py

Requirements:
    numpy pandas matplotlib pyyaml
"""

from __future__ import annotations

import io
import math
import re
import tarfile
import zipfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import matplotlib.pyplot as plt

# Reuse the already-validated literal transcription of dvmpx.F.
from validate_aao_pi0_model import (
    epsilon, tminq, xcheck_kine, xsigma_t, xsigma_l, xsigma_tt, xsigma_lt
)

HERE = Path(__file__).resolve().parent
IMPORT = HERE / "import"
OUTPUT = HERE / "output"
PNG = OUTPUT / "png"

HEPDATA_URLS = [
    "https://www.hepdata.net/download/submission/ins1294143/1/yaml",
    "https://www.hepdata.net/download/submission/64122/1/yaml",
]
CACHE_CANDIDATES = [
    IMPORT / "HEPData-ins1294143-v1.zip",
    IMPORT / "bedlinskiy2014_hepdata_yaml.zip",
    IMPORT / "bedlinskiy2014_hepdata_yaml.tar.gz",
]
EBEAM_GEV = 5.75


def _load_hepdata_archive() -> tuple[bytes, Path]:
    """Use a local HEPData archive first; download only as a fallback."""
    IMPORT.mkdir(parents=True, exist_ok=True)

    for path in CACHE_CANDIDATES:
        if path.exists():
            print(f"[Bedlinskiy] using local HEPData archive:\n  {path}")
            return path.read_bytes(), path

    errors = []
    for url in HEPDATA_URLS:
        try:
            print(f"[Bedlinskiy] downloading official HEPData archive:\n  {url}")
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                blob = r.read()
            if len(blob) < 1000:
                raise RuntimeError(f"download suspiciously small ({len(blob)} bytes)")

            # HEPData's YAML endpoint currently returns ZIP. Detect rather than
            # trusting the URL or Content-Type.
            if zipfile.is_zipfile(io.BytesIO(blob)):
                path = IMPORT / "HEPData-ins1294143-v1.zip"
            elif tarfile.is_tarfile(fileobj := io.BytesIO(blob)):
                path = IMPORT / "bedlinskiy2014_hepdata_yaml.tar.gz"
            else:
                raise RuntimeError("download is neither ZIP nor tar archive")
            path.write_bytes(blob)
            return blob, path
        except Exception as exc:
            errors.append(f"{url}: {exc}")

    expected = "\n".join(f"  {p}" for p in CACHE_CANDIDATES[:2])
    raise RuntimeError(
        "Could not download HEPData. Download the YAML archive for HEPData "
        "record ins1294143 and save it as either:\n"
        f"{expected}\n\n" + "\n".join(errors)
    )


def _yaml_documents_from_archive(blob: bytes):
    """Yield YAML documents from either the ZIP or tar archive HEPData supplies."""
    bio = io.BytesIO(blob)

    if zipfile.is_zipfile(bio):
        bio.seek(0)
        with zipfile.ZipFile(bio) as zf:
            for name in zf.namelist():
                if name.endswith("/") or not name.lower().endswith((".yaml", ".yml")):
                    continue
                text = zf.read(name).decode("utf-8")
                # submission.yaml is a multi-document YAML stream; table
                # files are normally single-document. load_all handles both.
                for i, doc in enumerate(yaml.safe_load_all(text)):
                    if doc is not None:
                        yield (name if i == 0 else f"{name}#doc{i+1}"), doc
        return

    bio.seek(0)
    try:
        with tarfile.open(fileobj=bio, mode="r:*") as tf:
            for member in tf.getmembers():
                if not member.isfile() or not member.name.lower().endswith((".yaml", ".yml")):
                    continue
                f = tf.extractfile(member)
                if f is None:
                    continue
                text = f.read().decode("utf-8")
                for i, doc in enumerate(yaml.safe_load_all(text)):
                    if doc is not None:
                        yield (member.name if i == 0 else f"{member.name}#doc{i+1}"), doc
        return
    except tarfile.TarError as exc:
        raise RuntimeError("HEPData input is neither a valid ZIP nor tar archive.") from exc


def _number(x):
    if x is None:
        return np.nan
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip().replace("−", "-")
    try:
        return float(s)
    except ValueError:
        return np.nan


def _central(v):
    """Central value for HEPData scalar or low/high bin."""
    if "value" in v:
        return _number(v["value"])
    lo, hi = _number(v.get("low")), _number(v.get("high"))
    return 0.5 * (lo + hi) if np.isfinite(lo) and np.isfinite(hi) else np.nan


def _sym_error(errors):
    """Quadrature of all published symmetric/asymmetric errors for one value."""
    sq = 0.0
    labels = []
    for e in errors or []:
        labels.append(str(e.get("label", "")))
        if "symerror" in e:
            z = _number(str(e["symerror"]).replace("%", ""))
            if "%" in str(e["symerror"]):
                # Caller cannot convert percent without the central value.
                return np.nan, ";".join(labels), True
            if np.isfinite(z):
                sq += z*z
        elif "asymerror" in e:
            plus = abs(_number(e["asymerror"].get("plus")))
            minus = abs(_number(e["asymerror"].get("minus")))
            z = 0.5*(plus+minus)
            if np.isfinite(z):
                sq += z*z
    return math.sqrt(sq), ";".join(labels), False


def _qualifier_map(dep):
    q = {}
    for item in dep.get("qualifiers", []) or []:
        q[str(item.get("name", "")).upper()] = item.get("value")
    return q


def _extract_range(text, names):
    """Extract central or midpoint from strings such as '1.14 - 1.16'."""
    if text is None:
        return np.nan
    s = str(text).replace("−", "-")
    nums = re.findall(r"(?<![A-Za-z])[-+]?\d+(?:\.\d+)?", s)
    vals = [float(z) for z in nums]
    if len(vals) >= 2:
        return 0.5*(vals[0]+vals[1])
    if vals:
        return vals[0]
    return np.nan


def _setting_from_description(desc):
    # HEPData descriptions use e.g.:
    # "The structure functions for Q**2 = 1.14 - 1.16 GeV**2 and
    #  XB = 0.131 - 0.133 as functions of t."
    s = (desc or "").replace("−", "-")
    mq = re.search(r"Q\*?\*?2\s*=\s*([0-9.]+)(?:\s*-\s*([0-9.]+))?", s, re.I)
    mx = re.search(r"X\s*B|XB", s, re.I)
    q2 = np.nan
    xb = np.nan
    if mq:
        a = float(mq.group(1)); b = float(mq.group(2)) if mq.group(2) else a
        q2 = 0.5*(a+b)
    if mx:
        tail = s[mx.start():]
        m = re.search(r"(?:X\s*B|XB)\s*=\s*([0-9.]+)(?:\s*-\s*([0-9.]+))?", tail, re.I)
        if m:
            a = float(m.group(1)); b = float(m.group(2)) if m.group(2) else a
            xb = 0.5*(a+b)
    return q2, xb


def _classify_observable(name):
    s = str(name).upper().replace(" ", "")
    # Be deliberately conservative: only consume the three published
    # structure functions, never an unrelated dependent variable.
    if ("SIG(T)" in s or "SIGMA_T" in s) and ("SIG(L)" in s or "SIGMA_L" in s):
        return "U"
    if "TT" in s:
        return "TT"
    if "LT" in s or "TL" in s:
        return "LT"
    return None


def parse_hepdata(blob: bytes) -> pd.DataFrame:
    """
    Parse the actual Bedlinskiy HEPData schema.

    Each table carries Q2, xB and -t as three point-by-point independent
    variables.  They are NOT table-description metadata.  This matters because
    Q2/xB vary slightly across the t points within a nominal setting.
    """
    rows = []
    n_tables = 0

    for name, doc in _yaml_documents_from_archive(blob):
        if not isinstance(doc, dict) or "dependent_variables" not in doc:
            continue

        indeps = doc.get("independent_variables", [])
        if not indeps:
            continue

        # Map independent-variable arrays by their actual HEPData headers.
        ivars = {}
        for iv in indeps:
            key = str(iv.get("header", {}).get("name", "")).upper().replace(" ", "")
            ivars[key] = iv

        qvar = next((v for k,v in ivars.items()
                     if k in ("Q**2","Q^2","Q2")), None)
        xbvar = next((v for k,v in ivars.items()
                      if k in ("XB","X_B")), None)
        tvar = next((v for k,v in ivars.items()
                     if k in ("-T","-T'","ABS(T)","|T|","T")), None)

        if qvar is None or xbvar is None or tvar is None:
            # submission.yaml and any non-data documents land here.
            continue

        qvals = qvar.get("values", [])
        xbvals = xbvar.get("values", [])
        tvals = tvar.get("values", [])
        n = len(tvals)
        if len(qvals) != n or len(xbvals) != n:
            raise RuntimeError(
                f"{name}: independent-variable lengths differ: "
                f"Q2={len(qvals)}, xB={len(xbvals)}, t={len(tvals)}"
            )

        n_tables += 1
        tname = str(tvar.get("header", {}).get("name", ""))
        tunit = str(tvar.get("header", {}).get("units", ""))

        for dep in doc.get("dependent_variables", []):
            obs_name = dep.get("header", {}).get("name", "")
            obs_unit = str(dep.get("header", {}).get("units", ""))
            tag = _classify_observable(obs_name)
            if tag is None:
                continue

            dvals = dep.get("values", [])
            if len(dvals) != n:
                raise RuntimeError(
                    f"{name} {obs_name}: dependent length {len(dvals)} != {n}"
                )

            for qv, xbv, tv, dv in zip(qvals, xbvals, tvals, dvals):
                q2 = _central(qv)
                xb = _central(xbv)
                mt = _central(tv)
                val = _number(dv.get("value"))
                if not all(np.isfinite(z) for z in (q2, xb, mt, val)):
                    continue

                err, labels, percent = _sym_error(dv.get("errors", []))
                if percent:
                    sq = 0.0
                    for e in dv.get("errors", []) or []:
                        if "symerror" in e:
                            raw = str(e["symerror"])
                            z = _number(raw.replace("%", ""))
                            z = abs(val)*z/100.0 if "%" in raw else z
                            if np.isfinite(z):
                                sq += z*z
                        elif "asymerror" in e:
                            plus = abs(_number(e["asymerror"].get("plus")))
                            minus = abs(_number(e["asymerror"].get("minus")))
                            z = 0.5*(plus+minus)
                            if np.isfinite(z):
                                sq += z*z
                    err = math.sqrt(sq)

                # HEPData explicitly labels this coordinate "-t"; VPK expects
                # signed t.  Preserve the published positive -t separately.
                signed_t = -abs(mt) if "-T" in tname.upper() else (mt if mt < 0 else -mt)

                rows.append({
                    "source_table": name,
                    "Q2_GeV2": q2,
                    "xB": xb,
                    "Ebeam_GeV": EBEAM_GEV,
                    "minus_t_GeV2": abs(mt),
                    "t_GeV2": signed_t,
                    "t_header": tname,
                    "t_units": tunit,
                    "observable": tag,
                    "observable_header": obs_name,
                    "observable_units": obs_unit,
                    "data_value": val,
                    "data_totalerr": err,
                    "error_labels": labels,
                })

    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("Parsed HEPData archive but found no Bedlinskiy structure functions.")

    print(f"[Bedlinskiy] parsed {n_tables} data tables, {len(df)} structure-function rows")
    print(f"[Bedlinskiy] Q2 range: {df.Q2_GeV2.min():.3f}--{df.Q2_GeV2.max():.3f} GeV^2")
    print(f"[Bedlinskiy] xB range: {df.xB.min():.3f}--{df.xB.max():.3f}")
    print(f"[Bedlinskiy] -t range: {df.minus_t_GeV2.min():.3f}--{df.minus_t_GeV2.max():.3f} GeV^2")
    return df


def evaluate(df):
    out = df.copy()
    model = []
    valid = []
    epsv = []
    tminv = []
    for r in out.itertuples(index=False):
        t, x, q2, e = float(r.t_GeV2), float(r.xB), float(r.Q2_GeV2), float(r.Ebeam_GeV)
        eps = epsilon(x, q2, e)
        tag = r.observable
        if tag == "U":
            m = xsigma_t(t,x,q2,e) + eps*xsigma_l(t,x,q2,e)
        elif tag == "TT":
            m = xsigma_tt(t,x,q2,e)
        elif tag == "LT":
            m = xsigma_lt(t,x,q2,e)
        else:
            m = np.nan
        model.append(m)
        valid.append(xcheck_kine(t,x,q2,e))
        epsv.append(eps)
        tminv.append(tminq(q2,x))
    out["vpk_value"] = model
    out["vpk_kine_valid"] = valid
    out["epsilon_vpk"] = epsv
    out["minus_tmin_vpk_GeV2"] = tminv
    out["data_over_vpk"] = np.where(np.asarray(model) != 0, out.data_value/np.asarray(model), np.nan)
    out["data_over_vpk_err"] = np.where(np.asarray(model) != 0, out.data_totalerr/np.abs(np.asarray(model)), np.nan)
    return out


def weighted_setting_summary(points):
    u = points[(points.observable=="U") & points.vpk_kine_valid &
               np.isfinite(points.data_over_vpk) & (points.data_over_vpk_err>0)].copy()
    rows=[]
    for table,g in u.groupby("source_table",sort=True):
        w=1/g.data_over_vpk_err.to_numpy(float)**2
        ratio=float(np.sum(w*g.data_over_vpk)/np.sum(w))
        err=float(math.sqrt(1/np.sum(w)))
        rows.append({
            "source_table": table,
            "Q2_GeV2": float(np.average(g.Q2_GeV2, weights=w)),
            "xB": float(np.average(g.xB, weights=w)),
            "n_t_points": len(g),
            "weighted_data_over_vpk_U": ratio,
            "weighted_data_over_vpk_U_err": err,
            "mean_minus_t_GeV2": float(np.average(g.minus_t_GeV2, weights=w)),
            "Q2_min_GeV2": float(g.Q2_GeV2.min()),
            "Q2_max_GeV2": float(g.Q2_GeV2.max()),
            "xB_min": float(g.xB.min()),
            "xB_max": float(g.xB.max()),
        })
    return pd.DataFrame(rows).sort_values("Q2_GeV2").reset_index(drop=True)


def make_plots(points, summary):
    PNG.mkdir(parents=True,exist_ok=True)

    # CLAS6-only setting summary: this is the key low-Q2 closure plot.
    fig,ax=plt.subplots(figsize=(8.4,5.6))
    sc=ax.scatter(summary.Q2_GeV2,summary.weighted_data_over_vpk_U,
                  c=summary.xB,s=45,zorder=3)
    ax.errorbar(summary.Q2_GeV2,summary.weighted_data_over_vpk_U,
                yerr=summary.weighted_data_over_vpk_U_err,fmt="none",capsize=2,zorder=2)
    ax.axhline(1,linewidth=1)
    cb=fig.colorbar(sc,ax=ax); cb.set_label(r"$x_B$")
    ax.set(xlabel=r"$Q^2$ (GeV$^2$)",ylabel=r"CLAS6 / VPK",
           title=r"Bedlinskiy CLAS6 $\sigma_U$: VPK fit-domain closure")
    ax.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(PNG/"bedlinskiy2014_data_over_vpk_vs_Q2.png",dpi=220)
    plt.close(fig)

    # Point-level residual versus -t, colored by Q2.
    u=points[(points.observable=="U") & points.vpk_kine_valid].copy()
    fig,ax=plt.subplots(figsize=(8.4,5.6))
    sc=ax.scatter(u.minus_t_GeV2,u.data_over_vpk,c=u.Q2_GeV2,s=24)
    ax.axhline(1,linewidth=1)
    cb=fig.colorbar(sc,ax=ax); cb.set_label(r"$Q^2$ (GeV$^2$)")
    ax.set(xlabel=r"$-t$ (GeV$^2$)",ylabel=r"CLAS6 / VPK",
           title=r"Bedlinskiy CLAS6 $\sigma_U$: point-level VPK residuals")
    ax.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(PNG/"bedlinskiy2014_data_over_vpk_vs_minus_t.png",dpi=220)
    plt.close(fig)

    # Combined CLAS6 + Hall-A setting-level diagnostic, if Hall-A has been run.
    hall = OUTPUT/"dlamini2021_vpk_setting_summary.csv"
    if hall.exists():
        h=pd.read_csv(hall)
        fig,ax=plt.subplots(figsize=(8.8,5.8))
        ax.errorbar(summary.Q2_GeV2,summary.weighted_data_over_vpk_U,
                    yerr=summary.weighted_data_over_vpk_U_err,fmt="o",capsize=2,
                    label="CLAS6 Bedlinskiy 2014")
        ax.errorbar(h.Q2_GeV2,h.weighted_data_over_vpk_U,
                    yerr=h.weighted_data_over_vpk_U_err,fmt="s",capsize=2,
                    label="Hall A Dlamini 2021")
        ax.axhline(1,linewidth=1)
        ax.set(xlabel=r"$Q^2$ (GeV$^2$)",ylabel="Data / VPK",
               title=r"VPK $\pi^0$ model: fit domain and higher-$Q^2$ validation")
        ax.grid(alpha=.25); ax.legend()
        fig.tight_layout()
        fig.savefig(PNG/"clas6_halla_data_over_vpk_vs_Q2.png",dpi=220)
        plt.close(fig)


def main():
    OUTPUT.mkdir(parents=True,exist_ok=True)
    blob, archive_path = _load_hepdata_archive()
    raw=parse_hepdata(blob)
    points=evaluate(raw)

    invalid=int((~points.vpk_kine_valid).sum())
    if invalid:
        print(f"[WARNING] {invalid}/{len(points)} parsed rows fail VPK kinematic checks.")
        print("          Inspect these rows before interpreting the comparison.")

    # Unit guard. We expect the published structure functions and VPK output
    # to both be nb/GeV^2. Abort rather than silently invent a conversion.
    units=sorted(set(str(x) for x in points.observable_units.dropna().unique()))
    print("[Bedlinskiy] published dependent-variable units:", units)
    if units and not all(("NB" in u.upper() and "GEV" in u.upper()) for u in units):
        raise RuntimeError(
            "Unexpected HEPData structure-function units. No automatic conversion "
            "was applied. Inspect observable_units in the parsed output."
        )

    summary=weighted_setting_summary(points)
    points.to_csv(OUTPUT/"bedlinskiy2014_vpk_point_comparison.csv",index=False)
    summary.to_csv(OUTPUT/"bedlinskiy2014_vpk_setting_summary.csv",index=False)
    make_plots(points,summary)

    print("\n=== Bedlinskiy 2014 CLAS6 sigma_U / VPK ===")
    print(summary.to_string(index=False,formatters={
        "Q2_GeV2":"{:.3f}".format,
        "xB":"{:.4f}".format,
        "weighted_data_over_vpk_U":"{:.4f}".format,
        "weighted_data_over_vpk_U_err":"{:.4f}".format,
    }))
    low=summary[(summary.Q2_GeV2>=1.5)&(summary.Q2_GeV2<=2.5)]
    if len(low):
        w=1/low.weighted_data_over_vpk_U_err.to_numpy(float)**2
        r=np.sum(w*low.weighted_data_over_vpk_U)/np.sum(w)
        er=math.sqrt(1/np.sum(w))
        print(f"\nQ2=1.5--2.5 GeV^2 setting-weighted closure: data/VPK = {r:.4f} +/- {er:.4f}")
        print("(Diagnostic only: setting systematics may be correlated.)")
    print(f"\nWrote:\n  {OUTPUT/'bedlinskiy2014_vpk_point_comparison.csv'}")
    print(f"  {OUTPUT/'bedlinskiy2014_vpk_setting_summary.csv'}")
    print(f"  {PNG}/")


if __name__ == "__main__":
    main()
