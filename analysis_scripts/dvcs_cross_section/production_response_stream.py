"""Bounded ROOT streaming around the persistent native response evaluator."""
from __future__ import annotations
import json,os,struct,subprocess,uuid
from pathlib import Path
import numpy as np
import uproot

def metadata(name):
    return name in ('runnum','evnum','event_sequence','helicity','weight','fiducial_status','rec_kinematics_valid') or name.startswith(('detector','num_','vx_','vy_','vz_','gen_','mc_'))

def read_exact(stream,n):
    result=bytearray(n);view=memoryview(result);offset=0
    while offset<n:
        count=stream.readinto(view[offset:])
        if not count:raise RuntimeError('Native particle-response evaluator closed its output unexpectedly')
        offset+=count
    return result

def stream_file(job,worker,log,batch_size=16384,progress=None):
    if not 1<=batch_size<=65536:raise ValueError('Response batch_size outside 1–65536')
    # Disable uproot's nested pools: file workers are the sole concurrency budget.
    executor=uproot.source.futures.TrivialExecutor();summary=dict(source=job['source'],entries=0,written=0,unresolved=0,nonphysical=0,note_rows=0,note_probability=job['note_probability'])
    protocol=Path(job['target']).with_name('protocol.'+uuid.uuid4().hex+'.json');process=None
    try:
        with uproot.open(job['source'],decompression_executor=executor,interpretation_executor=executor,num_workers=1) as source:
            tree=source['PhysicsEvents'];empty=tree.arrays(entry_stop=0,library='np');empty=dict(empty) if isinstance(empty,dict) else {n:empty[n] for n in empty.dtype.names}
            if 'response_applied' in empty:raise ValueError('Refusing already modified source')
            for n,v in empty.items():
                if v.ndim!=1 or v.dtype.kind not in 'fiub':raise ValueError('Response requires scalar numeric branch: '+n)
            # Integer identifiers/weights never pass through floating-point storage.
            physics=[n for n in empty if not metadata(n) or n=='runnum' or n.startswith('detector')]
            config=dict(job,columns=physics);protocol.write_text(json.dumps(config));process=subprocess.Popen([str(worker),str(protocol)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=log)
            with uproot.recreate(job['target'],compression=uproot.LZ4(4)) as output,Path(job['notes']).open('wb') as notes:
                types={n:v.dtype for n,v in empty.items()};types.update(response_applied=np.dtype('int32'),response_valid=np.dtype('int32'),response_unresolved=np.dtype('int32'));result=output.mktree('PhysicsEvents',types,title=tree.title)
                for chunk in tree.iterate(step_size=batch_size,library='np'):
                    chunk=dict(chunk) if isinstance(chunk,dict) else {n:chunk[n] for n in chunk.dtype.names};n=len(next(iter(chunk.values())));matrix=np.ascontiguousarray(np.column_stack([chunk[k] for k in physics]),dtype='<f8');process.stdin.write(struct.pack('<Q',n));process.stdin.write(matrix.data);process.stdin.flush()
                    transformed=np.frombuffer(read_exact(process.stdout,matrix.nbytes),dtype='<f8').reshape(matrix.shape);flags=np.frombuffer(read_exact(process.stdout,n*3*4),dtype='<i4').reshape(n,3);count=struct.unpack('<Q',read_exact(process.stdout,8))[0]
                    if count>n:raise ValueError('Invalid diagnostic sample count from native evaluator')
                    sample=read_exact(process.stdout,count*22*31*8);notes.write(sample)
                    changed=flags[:,0]!=0;invalid=(flags[:,1]==0)&any(x!=0 for x in job['parameters'])
                    for j,name in enumerate(physics):
                        if chunk[name].dtype.kind=='f' and not metadata(name):
                            keep=changed|invalid if name in ('x','Q2') else changed
                            if np.any(keep):chunk[name]=chunk[name].copy();chunk[name][keep]=transformed[keep,j]
                    for j,name in enumerate(('response_applied','response_valid','response_unresolved')):chunk[name]=flags[:,j]
                    result.extend(chunk);summary['entries']+=n;summary['written']+=n;summary['unresolved']+=int(flags[:,2].sum());summary['nonphysical']+=int((flags[:,1]==0).sum());summary['note_rows']+=count
                    if progress:progress(summary['entries'],tree.num_entries)
                process.stdin.close()
                if process.wait()!=0:raise RuntimeError('Native particle-response evaluator failed; see writer.log')
            if summary['entries']!=tree.num_entries:raise ValueError('Response changed input entry count')
        summary['note_shape']=[summary['note_rows'],22,31];Path(job['summary']).write_text(json.dumps(summary,indent=2));return summary
    finally:
        if process and process.poll() is None:process.terminate();process.wait()
        if process and process.stdout:process.stdout.close()
        if process and process.stdin and not process.stdin.closed:process.stdin.close()
        protocol.unlink(missing_ok=True)
