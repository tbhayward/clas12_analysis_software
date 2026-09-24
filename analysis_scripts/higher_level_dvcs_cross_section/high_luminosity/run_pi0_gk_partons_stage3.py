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
# Four points determine const + cos(phi) + cos(2phi), with redundancy.
DEFAULT_PHI_DEG = (0.0, 60.0, 120.0, 180.0, 240.0, 300.0)

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
<module type="RunningAlphaStrongModule" name="RunningAlphaStrongGK"></module>
<module type="GPDModule" name="GPDGK19"></module>
</module>
</module>
</module>"""

def task_xml(row):
    # PARTONS DVMP phi is radians in the Trento convention.
    phi_rad=math.radians(float(row.phi_deg))
    return f"""<task service="DVMPObservableService" method="computeSingleKinematic" storeInDB="0">
<kinematics type="DVMPObservableKinematic">
<param name="xB" value="{float(row.xB):.15g}" />
<param name="t" value="{-abs(float(row.minus_t)):.15g}" />
<param name="Q2" value="{float(row.Q2):.15g}" />
<param name="E" value="{float(row.E):.15g}" />
<param name="phi" value="{phi_rad:.15g}" />
<param name="meson" value="9" />
</kinematics>
<computation_configuration>
{module_xml()}
</computation_configuration>
</task>
<task service="DVMPObservableService" method="printResults"></task>"""

def load_queries(stage3:Path):
    qpath=stage3/"model_queries"/"01_common_gk_query_points.csv"
    if not qpath.exists(): raise FileNotFoundError(qpath)
    q=pd.read_csv(qpath)
    req=["point_id","Q2_common_GeV2","xB_common","minus_t_common_GeV2"]
    miss=[c for c in req if c not in q]
    if miss: raise RuntimeError(f"Stage-3 query file missing {miss}")
    return q

def expanded_points(q):
    rows=[]
    # Beam energies are the actual supplied campaign energies.
    for r in q.itertuples(index=False):
        for camp,E in (("rga",10.604),("rgk",6.535)):
            for phi in DEFAULT_PHI_DEG:
                rows.append(dict(point_id=r.point_id,campaign=camp,E=E,
                    Q2=r.Q2_common_GeV2,xB=r.xB_common,
                    minus_t=r.minus_t_common_GeV2,phi_deg=phi))
    out=pd.DataFrame(rows)
    out.insert(0,"eval_id",[f"E{i:06d}" for i in range(len(out))])
    return out

def write_chunks(points,out,chunk_size):
    xml_dir=out/"xml"; xml_dir.mkdir(parents=True,exist_ok=True)
    jobs=[]
    for ic,start in enumerate(range(0,len(points),chunk_size)):
        g=points.iloc[start:start+chunk_size].copy()
        xml=xml_dir/f"gk_pi0_{ic:05d}.xml"
        mp=xml_dir/f"gk_pi0_{ic:05d}_point_map.csv"
        body="\n".join(task_xml(r) for r in g.itertuples(index=False))
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

    q=load_queries(a.stage3.resolve())
    pts=expanded_points(q)
    if not a.all:
        first=q.iloc[[0]]
        pts=expanded_points(first)
        print(f"\n[probe] point {first.iloc[0].point_id}: evaluating 6 phi values at each of 2 beam energies")
    else:
        print(f"\n[production] {len(q)} common points x 2 energies x 6 phi = {len(pts)} PARTONS evaluations")

    pts.to_csv(out/"01_partons_evaluation_points.csv",index=False)
    jobs=write_chunks(pts,out,a.chunk_size)
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
        print(f"Inspect {out/'logs'} and send me the first stdout/stderr pair.")
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
        print("If they are sensible, rerun with --all.")
    else:
        print(f"\nWrote complete public-observable GK grid to {out}")
        print("This file deliberately does NOT yet label the harmonic coefficients")
        print("as sigma_T/L/TT/LT; that conversion will be made only after the")
        print("PARTONS electroproduction/reduced-cross-section normalization is verified.")

if __name__=="__main__":
    main()
