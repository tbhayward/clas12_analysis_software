#!/usr/bin/env python3
"""
evaluate_pi0_gk_partons.py

Self-contained PARTONS GK pi0 evaluator for the high_luminosity study.

This is NEW code.  It does not import or call the older DVCS/BH scripts.
It only reuses the already-proven farm configuration:
  /work/clas12/thayward/partons/partons-example
  ./bin/PARTONS_example
  partons_v4.sif (historical filename; installed executable may be PARTONS 5)

First goal:
  validate the DVMP XML/module chain with one common Stage-3 point and retain
  complete stdout/stderr.  Once the probe succeeds, --all evaluates a compact
  phi grid at both RGA and RGK beam energies for every common point.

Why total cross sections first?
  PARTONS exposes DVMPCrossSectionUUUMinus publicly, whereas the GK06
  CrossSectionT/L/TT/LT methods are private.  We therefore validate the public
  observable and its conventions before deriving reduced partial cross sections.
  No preliminary CLAS12 central values are used.

Expected repo location:
  .../analysis_scripts/higher_level_dvcs_cross_section/high_luminosity/
"""

from __future__ import annotations
import argparse, math, os, re, shutil, subprocess, sys, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

DEFAULT_PROJECT = Path("/work/clas12/thayward/partons/partons-example")
DEFAULT_EXECUTABLE = "./bin/PARTONS_example"
DEFAULT_WORKERS = 4
DEFAULT_CHUNK_SIZE = 30
# PARTONS MesonType::fromString() expects its canonical string representation.
# The enum integer (PI0=9) is NOT the XML representation; "PI0" also maps to
# UNDEFINED in the installed build.  PARTONS canonical type strings are lower-case.
DEFAULT_MESON_VALUE = "pi0"
# Four points determine const + cos(phi) + cos(2phi), with redundancy.
DEFAULT_PHI_DEG = (0.0, 60.0, 120.0, 180.0, 240.0, 300.0)
PROBE_PHI_DEG = (0.0, 90.0, 180.0)
PROTON_MASS_GEV = 0.9382720813

RESULT_RE = re.compile(
    r"Result:\s*([+-]?(?:[\d,]+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*\[([^\]]+)\]"
)

def parse_args():
    here=Path(__file__).resolve().parent
    p=argparse.ArgumentParser()
    p.add_argument("--stage3",type=Path,default=here/"output"/"pi0_gk_stage3")
    p.add_argument("--output",type=Path,default=here/"output"/"pi0_gk_stage3"/"partons_gk")
    p.add_argument("--project",type=Path,default=DEFAULT_PROJECT)
    p.add_argument("--executable",default=DEFAULT_EXECUTABLE)
    p.add_argument("--sif",type=Path,default=None,
                   help="PARTONS SIF. If omitted, search sensible high_luminosity/project/work locations.")
    p.add_argument("--workers",type=int,default=DEFAULT_WORKERS)
    p.add_argument("--chunk-size",type=int,default=DEFAULT_CHUNK_SIZE)
    p.add_argument("--all",action="store_true",
                   help="After preflight, evaluate all 154 points. Default is a one-point probe.")
    p.add_argument("--force",action="store_true")
    p.add_argument("--meson-value",default=DEFAULT_MESON_VALUE,
                   help="MesonType string passed to PARTONS XML (default: pi0).")
    p.add_argument("--rga-input",type=Path,default=None,
                   help="RGA combined_reduced_cross_sections.csv. Auto-discovered when omitted.")
    p.add_argument("--rgk-input",type=Path,default=None,
                   help="RGK rgk6535_reduced_cross_sections.csv. Auto-discovered when omitted.")
    p.add_argument("--dry-run",action="store_true",
                   help="Write XML/maps but do not invoke PARTONS.")
    return p.parse_args()

def find_sif(explicit: Path|None, here: Path, project: Path) -> Path:
    if explicit is not None:
        return explicit.expanduser().resolve()
    candidates=[
        here/"partons_v4.sif",
        project/"partons_v4.sif",
        project.parent/"partons_v4.sif",
        # Established container used by the existing PARTONS work on ifarm.
        Path("/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/external_scripts/partons_v4.sif"),
        Path("/work/clas12/thayward/partons/partons_v4.sif"),
        Path.cwd()/"partons_v4.sif",
    ]
    for p in candidates:
        if p.exists(): return p.resolve()
    raise FileNotFoundError(
        "Could not locate partons_v4.sif. Pass --sif /full/path/partons_v4.sif. "
        "Searched:\n  " + "\n  ".join(str(x) for x in candidates)
    )

def preflight(sif:Path,project:Path,executable:str):
    exe=project/executable.removeprefix("./")
    print("[PARTONS DVMP preflight]")
    print(f"  apptainer : {shutil.which('apptainer') or 'NOT FOUND'}")
    print(f"  SIF       : {sif}")
    print(f"  project   : {project}")
    print(f"  executable: {exe}")
    bad=[]
    if shutil.which("apptainer") is None: bad.append("apptainer not on PATH")
    if not sif.exists(): bad.append(f"SIF missing: {sif}")
    if not project.exists(): bad.append(f"project missing: {project}")
    if not exe.exists(): bad.append(f"executable missing: {exe}")
    elif not os.access(exe,os.X_OK): bad.append(f"executable not executable: {exe}")
    tmp=project/"bin"/"tmp"
    if not tmp.exists(): bad.append(f"logger directory missing: {tmp}")
    elif not os.access(tmp,os.W_OK): bad.append(f"logger directory not writable: {tmp}")
    if bad: raise RuntimeError("Preflight failed:\n  - "+"\n  - ".join(bad))

def module_xml():
    # Public PARTONS DVMP chain:
    # DVMPCrossSectionUUUMinus -> DVMPProcessGK06 -> DVMPCFFGK06 -> GPDGK19.
    # Explicit scales/xi converter and GK alpha_s follow the module dependencies
    # documented for the DVMP channel.
    return """<module type="DVMPObservableModule" name="DVMPCrossSectionUUUMinus">
<module type="DVMPProcessModule" name="DVMPProcessGK06">
<module type="DVMPScalesModule" name="DVMPScalesQ2Multiplier">
<param name="lambda" value="1." />
</module>
<module type="DVMPXiConverterModule" name="DVMPXiConverterXBToXi"></module>
<module type="DVMPConvolCoeffFunctionModule" name="DVMPCFFGK06">
<param name="qcd_order_type" value="LO" />
<module type="RunningAlphaStrongModule" name="RunningAlphaStrongGK"></module>
<module type="GPDModule" name="GPDGK19"></module>
</module>
</module>
</module>"""

def task_xml(row, meson_value=DEFAULT_MESON_VALUE):
    # PARTONS prints DVMPObservableKinematic::phi in degrees.  Pass the
    # requested Trento azimuth directly rather than converting it to radians.
    phi_deg=float(row.phi_deg)
    return f"""<task service="DVMPObservableService" method="computeSingleKinematic" storeInDB="0">
<kinematics type="DVMPObservableKinematic">
<param name="xB" value="{float(row.xB):.15g}" />
<param name="t" value="{-abs(float(row.minus_t)):.15g}" />
<param name="Q2" value="{float(row.Q2):.15g}" />
<param name="E" value="{float(row.E):.15g}" />
<param name="phi" value="{phi_deg:.15g}" />
<param name="meson" value="{meson_value}" />
</kinematics>
<computation_configuration>
{module_xml()}
</computation_configuration>
</task>
<task service="DVMPObservableService" method="printResults"></task>"""

def _find_campaign_input(explicit:Path|None, stage3:Path, filename:str) -> Path:
    """Locate the native Rosenbluth campaign CSV without silently using common-bin centers."""
    if explicit is not None:
        path=explicit.expanduser().resolve()
        if not path.exists(): raise FileNotFoundError(path)
        return path

    roots=[stage3,stage3.parent,Path(__file__).resolve().parent]
    candidates=[]
    for root in roots:
        candidates += [
            root/filename,
            root/"rga_10604"/filename,
            root/"rgk_6535"/filename,
            root/"inputs"/"rga_10604"/filename,
            root/"inputs"/"rgk_6535"/filename,
            root/"fa18_rosenbluth_inputs"/"rga_10604"/filename,
            root/"fa18_rosenbluth_inputs"/"rgk_6535"/filename,
        ]
    for path in candidates:
        if path.exists(): return path.resolve()

    # A bounded recursive search is useful because exported input packages carry
    # timestamped top-level directory names.
    seen=set()
    for root in roots:
        if not root.exists(): continue
        for path in root.glob(f"**/{filename}"):
            rp=path.resolve()
            if rp not in seen:
                seen.add(rp)
                return rp
    raise FileNotFoundError(
        f"Could not locate native campaign input {filename}. "
        f"Pass --{'rga-input' if filename.startswith('combined') else 'rgk-input'} /full/path/{filename}."
    )


def _native_campaign_bins(path:Path, campaign:str) -> pd.DataFrame:
    """Collapse phi rows to one native flux-coordinate record per (iq2,ixb,it)."""
    d=pd.read_csv(path)
    req=["iq2","ixb","it","Q2_flux_coordinate_GeV2","xB_flux_coordinate",
         "minus_t_center_GeV2","beam_energy_GeV","virtual_photon_epsilon"]
    miss=[c for c in req if c not in d]
    if miss: raise RuntimeError(f"{path} missing required columns {miss}")

    keys=["iq2","ixb","it"]
    vals=req[3:]
    # Every phi bin belonging to one hadronic bin must carry identical model
    # coordinates.  Check this explicitly before taking the first row.
    nunique=d.groupby(keys,dropna=False)[vals].nunique(dropna=False)
    bad=nunique[(nunique>1).any(axis=1)]
    if len(bad):
        raise RuntimeError(f"{path}: {len(bad)} hadronic bins have inconsistent native coordinates across phi")
    out=d.groupby(keys,as_index=False,dropna=False)[vals].first()
    out["campaign"]=campaign
    return out


def load_queries(stage3:Path, rga_input:Path|None=None, rgk_input:Path|None=None):
    qpath=stage3/"model_queries"/"01_common_gk_query_points.csv"
    if not qpath.exists(): raise FileNotFoundError(qpath)
    q=pd.read_csv(qpath)
    req=["point_id","Q2_common_GeV2","xB_common","minus_t_common_GeV2",
         "iq2_rga","ixb_rga","it_rga","iq2_rgk","ixb_rgk","it_rgk"]
    miss=[c for c in req if c not in q]
    if miss: raise RuntimeError(f"Stage-3 query file missing {miss}")

    rga_path=_find_campaign_input(rga_input,stage3,"combined_reduced_cross_sections.csv")
    rgk_path=_find_campaign_input(rgk_input,stage3,"rgk6535_reduced_cross_sections.csv")
    print("[native campaign coordinates]")
    print(f"  RGA input : {rga_path}")
    print(f"  RGK input : {rgk_path}")

    for camp,path in (("rga",rga_path),("rgk",rgk_path)):
        native=_native_campaign_bins(path,camp)
        native=native.rename(columns={
            "iq2":f"iq2_{camp}","ixb":f"ixb_{camp}","it":f"it_{camp}",
            "Q2_flux_coordinate_GeV2":f"Q2_{camp}_GeV2",
            "xB_flux_coordinate":f"xB_{camp}",
            "minus_t_center_GeV2":f"minus_t_{camp}_GeV2",
            "beam_energy_GeV":f"E_{camp}_GeV",
            "virtual_photon_epsilon":f"epsilon_{camp}_native",
        }).drop(columns="campaign")
        keys=[f"iq2_{camp}",f"ixb_{camp}",f"it_{camp}"]
        q=q.merge(native,on=keys,how="left",validate="many_to_one")
        missing=q[f"Q2_{camp}_GeV2"].isna()
        if missing.any():
            ids=", ".join(q.loc[missing,"point_id"].astype(str).head(8))
            raise RuntimeError(f"Could not match {missing.sum()} Stage-3 points to native {camp.upper()} bins; first: {ids}")

        # The epsilon already carried by Stage 3 should originate from the same
        # native coordinates.  Preserve both and demand agreement.
        old=f"epsilon_{camp}"
        new=f"epsilon_{camp}_native"
        if old in q:
            diff=np.abs(q[old].astype(float)-q[new].astype(float))
            if np.nanmax(diff)>1e-10:
                raise RuntimeError(f"Stage-3/native {camp.upper()} epsilon mismatch: max |delta|={np.nanmax(diff):.3e}")

    return q


def expanded_points(q, phi_grid=DEFAULT_PHI_DEG):
    rows=[]
    skipped=[]
    # Evaluate each campaign at its OWN measured flux-weighted (Q2,xB)
    # coordinate and beam energy.  Common coordinates remain metadata for the
    # later explicit model-based translation to the Rosenbluth comparison point.
    for r in q.itertuples(index=False):
        for camp in ("rga","rgk"):
            E=float(getattr(r,f"E_{camp}_GeV"))
            Q2=float(getattr(r,f"Q2_{camp}_GeV2"))
            xB=float(getattr(r,f"xB_{camp}"))
            mt=float(getattr(r,f"minus_t_{camp}_GeV2"))
            eps=float(getattr(r,f"epsilon_{camp}_native"))
            y=Q2/(2.0*PROTON_MASS_GEV*E*xB)
            base=dict(
                point_id=r.point_id,campaign=camp,E=E,y=y,epsilon=eps,
                Q2=Q2,xB=xB,minus_t=mt,
                Q2_common=float(r.Q2_common_GeV2),xB_common=float(r.xB_common),
                minus_t_common=float(r.minus_t_common_GeV2),
                delta_Q2_from_common=Q2-float(r.Q2_common_GeV2),
                delta_xB_from_common=xB-float(r.xB_common),
                delta_minus_t_from_common=mt-float(r.minus_t_common_GeV2),
            )
            if not (0.0 < y < 1.0):
                skipped.append({**base,"reason":"native campaign flux coordinate has y outside (0,1)"})
                continue
            for phi in phi_grid:
                rows.append({**base,"phi_deg":float(phi)})
    out=pd.DataFrame(rows)
    if len(out): out.insert(0,"eval_id",[f"E{i:06d}" for i in range(len(out))])
    return out,pd.DataFrame(skipped)

def write_chunks(points,out,chunk_size,meson_value=DEFAULT_MESON_VALUE):
    xml_dir=out/"xml"; xml_dir.mkdir(parents=True,exist_ok=True)
    jobs=[]
    for ic,start in enumerate(range(0,len(points),chunk_size)):
        g=points.iloc[start:start+chunk_size].copy()
        xml=xml_dir/f"gk_pi0_{ic:05d}.xml"
        mp=xml_dir/f"gk_pi0_{ic:05d}_point_map.csv"
        body="\n".join(task_xml(r, meson_value) for r in g.itertuples(index=False))
        xml.write_text(
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            '<scenario date="2026-09-24" description="high-luminosity pi0 GK projection">\n'
            +body+'\n</scenario>\n')
        g.to_csv(mp,index=False)
        jobs.append(dict(chunk=ic,xml=str(xml),map=str(mp),n=len(g)))
    return jobs

def parse_stdout(txt):
    vals=[]
    for v,u in RESULT_RE.findall(txt):
        vals.append((float(v.replace(",","")),u.strip()))
    return vals

def run_job(job,sif,project,executable,out,force=False,dry=False):
    ic=job["chunk"]; logs=out/"logs"; logs.mkdir(parents=True,exist_ok=True)
    stdout=logs/f"gk_pi0_{ic:05d}.stdout.txt"
    stderr=logs/f"gk_pi0_{ic:05d}.stderr.txt"
    status=logs/f"gk_pi0_{ic:05d}.status.txt"
    if not force and stdout.exists() and status.exists():
        try:
            rc=int(status.read_text().strip())
            vals=parse_stdout(stdout.read_text(errors="replace"))
            cached_txt=stdout.read_text(errors="replace")
            cached_err=stderr.read_text(errors="replace") if stderr.exists() else ""
            logger_error=("[ERROR]" in cached_txt) or ("[ERROR]" in cached_err)
            if rc==0 and not logger_error and len(vals)==job["n"] and all(np.isfinite(v) for v,_ in vals):
                return dict(chunk=ic,returncode=0,process_returncode=0,
                            logger_error=False,reused=True,nresults=len(vals),
                            expected=job["n"])
        except Exception: pass
    if dry:
        return dict(chunk=ic,returncode=0,reused=False,nresults=0,dry_run=True)

    xml=Path(job["xml"]).resolve()
    binds=[project.resolve(),xml.parent]
    cmd=["apptainer","exec"]
    for b in binds: cmd += ["--bind",f"{b}:{b}"]
    cmd += ["--pwd",str(project.resolve()),str(sif),executable,str(xml)]
    env=os.environ.copy()
    for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
        env[k]="1"
    t=time.perf_counter()
    p=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env)
    stdout.write_text(p.stdout); stderr.write_text(p.stderr); status.write_text(str(p.returncode)+"\n")
    vals=parse_stdout(p.stdout)
    logger_error = ("[ERROR]" in p.stdout) or ("[ERROR]" in p.stderr)
    effective_rc = p.returncode if p.returncode != 0 else (90 if logger_error else 0)
    return dict(chunk=ic,returncode=effective_rc,process_returncode=p.returncode,
                logger_error=logger_error,reused=False,nresults=len(vals),
                expected=job["n"],seconds=time.perf_counter()-t)

def collect(jobs,out):
    pieces=[]
    units=set()
    for j in jobs:
        mp=pd.read_csv(j["map"])
        stdout=out/"logs"/f"gk_pi0_{j['chunk']:05d}.stdout.txt"
        if not stdout.exists(): continue
        vals=parse_stdout(stdout.read_text(errors="replace"))
        if len(vals)!=len(mp):
            raise RuntimeError(f"Chunk {j['chunk']}: parsed {len(vals)} results for {len(mp)} points")
        mp["partons_value"]=[x[0] for x in vals]
        mp["partons_unit"]=[x[1] for x in vals]
        units.update(mp.partons_unit.unique())
        pieces.append(mp)
    if not pieces: return pd.DataFrame(),units
    return pd.concat(pieces,ignore_index=True),units

def harmonic_closure(df):
    """Fit PARTONS observable vs phi at fixed point/beam: a0+a1 cos(phi)+a2 cos(2phi)."""
    rows=[]
    for (pid,camp),g in df.groupby(["point_id","campaign"]):
        ph=np.deg2rad(g.phi_deg.to_numpy(float))
        A=np.column_stack([np.ones(len(g)),np.cos(ph),np.cos(2*ph)])
        y=g.partons_value.to_numpy(float)
        b,*_=np.linalg.lstsq(A,y,rcond=None)
        pred=A@b
        scale=max(np.max(np.abs(y)),1e-300)
        rows.append(dict(point_id=pid,campaign=camp,E=float(g.E.iloc[0]),
            Q2=float(g.Q2.iloc[0]),xB=float(g.xB.iloc[0]),minus_t=float(g.minus_t.iloc[0]),
            coefficient_const=b[0],coefficient_cosphi=b[1],coefficient_cos2phi=b[2],
            max_fractional_harmonic_residual=float(np.max(np.abs(y-pred))/scale)))
    return pd.DataFrame(rows)

def main():
    a=parse_args(); here=Path(__file__).resolve().parent
    out=a.output.resolve(); out.mkdir(parents=True,exist_ok=True)
    project=a.project.expanduser().resolve()
    sif=find_sif(a.sif,here,project)
    preflight(sif,project,a.executable)

    q=load_queries(a.stage3.resolve(),a.rga_input,a.rgk_input)
    if not a.all:
        first=q.iloc[[0]]
        pts,skipped=expanded_points(first,PROBE_PHI_DEG)
        print(f"\n[probe] point {first.iloc[0].point_id}: evaluating {len(PROBE_PHI_DEG)} phi values for each physically allowed beam energy")
    else:
        pts,skipped=expanded_points(q,DEFAULT_PHI_DEG)
        print(f"\n[production] {len(q)} matched Rosenbluth points; {len(pts)} PARTONS evaluations at native campaign coordinates")

    pts.to_csv(out/"01_partons_evaluation_points.csv",index=False)
    skipped.to_csv(out/"01b_partons_skipped_unphysical.csv",index=False)
    if len(skipped):
        print(f"  prefilter : skipped {len(skipped)} point/campaign combination(s) with y outside (0,1)")
        for r in skipped.head(8).itertuples(index=False):
            print(f"              {r.point_id} {r.campaign}: E={r.E:.3f} GeV, Q2={r.Q2:.4g} GeV^2, xB={r.xB:.4g}, y={r.y:.4f}")
        if len(skipped)>8: print(f"              ... and {len(skipped)-8} more (see 01b_partons_skipped_unphysical.csv)")
    print(f"  meson XML : {a.meson_value!r} (MesonType::fromString representation)")
    jobs=write_chunks(pts,out,a.chunk_size,a.meson_value)
    pd.DataFrame(jobs).to_csv(out/"02_chunk_manifest.csv",index=False)

    if a.dry_run:
        print(f"[dry-run] wrote {len(jobs)} XML chunk(s) to {out/'xml'}")
        return

    results=[]
    if a.workers<=1 or len(jobs)==1:
        for j in jobs:
            results.append(run_job(j,sif,project,a.executable,out,a.force,False))
    else:
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            fut=[ex.submit(run_job,j,sif,project,a.executable,out,a.force,False) for j in jobs]
            for f in as_completed(fut): results.append(f.result())
    stat=pd.DataFrame(results).sort_values("chunk")
    stat.to_csv(out/"03_chunk_status.csv",index=False)
    print(stat.to_string(index=False))

    bad=stat[(stat.returncode!=0)|(stat.nresults!=stat.get("expected",stat.nresults))]
    if len(bad):
        print(f"\nERROR: {len(bad)} chunk(s) failed or returned the wrong result count.")
        first_chunk=int(bad.iloc[0].chunk)
        so=out/"logs"/f"gk_pi0_{first_chunk:05d}.stdout.txt"
        se=out/"logs"/f"gk_pi0_{first_chunk:05d}.stderr.txt"
        combined="\n".join([
            so.read_text(errors="replace") if so.exists() else "",
            se.read_text(errors="replace") if se.exists() else "",
        ])
        err_lines=[re.sub(r"\x1b\[[0-9;]*m", "", x) for x in combined.splitlines()
                   if "[ERROR]" in x or "Exception" in x or "UNDEFINED" in x]
        if err_lines:
            print("First PARTONS diagnostic:")
            for line in err_lines[:8]: print("  "+line)
        if "MesonType::getPossibleGPDTypes" in combined and "UNDEFINED" in combined:
            print(f"\nMeson serialization {a.meson_value!r} still maps to UNDEFINED in this build.")
            print("Do not run --all. Try a build-local canonical MesonType string only after")
            print("checking MesonType::toString()/fromString() in the installed PARTONS source.")
        if "QCD order: UNDEFINED not implemented" in combined:
            print("\nDVMPCFFGK06 still sees an undefined QCD order despite the explicit LO setting.")
            print("Do not run --all; inspect how this PARTONS build serializes qcd_order_type.")
        print(f"Full logs: {out/'logs'}")
        sys.exit(2)

    raw,units=collect(jobs,out)
    raw.to_csv(out/"04_partons_raw_phi_grid.csv",index=False)
    harm=harmonic_closure(raw)
    harm.to_csv(out/"05_partons_phi_harmonic_closure.csv",index=False)

    print(f"\nParsed {len(raw)} finite PARTONS results.")
    print("Units reported by PARTONS:",", ".join(sorted(units)))
    print("Max phi-harmonic closure residual:",
          f"{harm.max_fractional_harmonic_residual.max():.3e}")
    if not a.all:
        print("\nProbe succeeded. Review the reported units and harmonic closure.")
        print("Native RGA/RGK coordinates were used. Review the reported units and harmonic closure before --all.")
    else:
        print(f"\nWrote complete public-observable GK grid to {out}")
        print("This file deliberately does NOT yet label the harmonic coefficients")
        print("as sigma_T/L/TT/LT; that conversion will be made only after the")
        print("PARTONS electroproduction/reduced-cross-section normalization is verified.")

if __name__=="__main__":
    main()
