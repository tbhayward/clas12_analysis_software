"""Generated phase space shared by numerical CSV corrections.
Matches kinematic_bin_integration.h and passes_bin_volume_phase_space_cuts.
Detector momentum/angle selections belong to the acceptance, not this volume.
"""
from dataclasses import dataclass
from functools import lru_cache
import json,math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

MP=.9382720813
YMIN=.19
YMAX=.85
ENERGIES={10.6:10.60,10.2:10.1998}

def energy_for_label(e):return ENERGIES[float(e)]
def minus_t_min(x,q):
    eps=4*MP*MP*x*x/q
    return q*(2*(1-x)*(1-math.sqrt(1+eps))+eps)/(4*x*(1-x)+eps)

@dataclass(frozen=True)
class Domain:
    apply_w: bool=True
    w_min: float=2.
    y_min: float=YMIN
    y_max: float=YMAX

    @classmethod
    def from_json(cls,path):
        d=json.loads(Path(path).read_text())
        return cls(bool(d['enable_analysis_phase_space_cuts']),float(d['analysis_W_min_GeV']),
                   float(d.get('generated_y_min',YMIN)),float(d.get('generated_y_max',YMAX)))

    def q_range(self,x,box,e):
        if not 0<x<1:return 0.,0.
        lo=max(box[2],2*MP*e*x*self.y_min);hi=min(box[3],2*MP*e*x*self.y_max)
        if self.apply_w:lo=max(lo,(self.w_min*self.w_min-MP*MP)*x/(1-x))
        return lo,hi

    def contains(self,point):
        x,q,t,phi,e=point
        if not all(math.isfinite(v) for v in point) or not(0<x<1 and q>0):return False
        y=q/(2*MP*x*e)
        return (self.y_min<y<self.y_max and (not self.apply_w or MP*MP+q*(1/x-1)>self.w_min**2)
                and t>=minus_t_min(x,q))

    def x_breaks(self,box,e):
        x0,x1,q0,q1,t0,t1=box;v=[x0,x1]
        def add(x):
            if math.isfinite(x) and x0<x<x1:v.append(x)
        for y in (self.y_min,self.y_max):
            slope=2*MP*e*y
            for q in (q0,q1):add(q/slope)
            if self.apply_w:add(1-(self.w_min**2-MP*MP)/slope)
        if self.apply_w:
            for q in (q0,q1):add(q/(q+self.w_min**2-MP*MP))
        base=sorted(set(v));extra=[]
        for a,b in zip(base[:-1],base[1:]):
            if b-a<1e-13:continue
            for side in (0,1):
                for threshold in (t0,t1):
                    def f(x):
                        q=self.q_range(x,box,e)[side]
                        return minus_t_min(x,q)-threshold if q>0 else -threshold
                    left=a+1e-12*(b-a);right=b-1e-12*(b-a)
                    if f(left)*f(right)<0:extra.append(brentq(f,left,right,xtol=1e-13))
        return sorted(set(base+extra))

    def q_segments(self,x,box,e):
        lo,hi=self.q_range(x,box,e)
        if hi<=lo:return []
        cuts=[lo,hi]
        # t_min is monotonic in Q² at fixed x. Split its crossings instead
        # of sampling a rectangular grid then throwing away forbidden nodes.
        for threshold in (box[4],box[5]):
            f=lambda q:minus_t_min(x,q)-threshold
            if f(lo)*f(hi)<0:cuts.append(brentq(f,lo,hi,xtol=1e-12))
        cuts=sorted(set(cuts))
        return [(a,b) for a,b in zip(cuts[:-1],cuts[1:]) if minus_t_min(x,(a+b)/2)<box[5]]

    @lru_cache(maxsize=2048)
    def fraction(self,box,e):
        x0,x1,q0,q1,t0,t1=box
        if not(x1>x0 and q1>q0 and t1>t0):raise ValueError('Invalid bin edges')
        def fx(x):
            total=0.
            for lo,hi in self.q_segments(x,box,e):
                f=lambda q:max(0.,t1-max(t0,minus_t_min(x,q)))/(t1-t0)
                value,error=quad(f,lo,hi,epsabs=1e-11*(hi-lo),epsrel=1e-10,limit=150)
                if error>max(1e-10*(hi-lo),1e-8*abs(value)):raise ArithmeticError('Q² volume integration did not converge')
                total+=value
            return total/(q1-q0)
        cuts=self.x_breaks(box,e);total=0.
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            value,error=quad(fx,lo,hi,epsabs=1e-11*(hi-lo),epsrel=1e-10,limit=200)
            if error>max(1e-9*(hi-lo),1e-8*abs(value)):raise ArithmeticError('xB volume integration did not converge')
            total+=value
        return min(1.,max(0.,total/(x1-x0)))

    def volume(self,box,phi_edges,e):
        rectangle=(box[1]-box[0])*(box[3]-box[2])*(box[5]-box[4])*math.radians(phi_edges[1]-phi_edges[0])
        return rectangle*self.fraction(tuple(box),e),rectangle

    def quadrature(self,box,phi_edges,e,order):
        nodes,weights=gauss(order);points=[];wts=[]
        def mapped(a,b):return (a+b)/2+(b-a)/2*nodes,(b-a)/2*weights
        cuts=self.x_breaks(box,e)
        for xl,xh in zip(cuts[:-1],cuts[1:]):
            xs,xw=mapped(xl,xh)
            for x,wx in zip(xs,xw):
                for ql,qh in self.q_segments(float(x),box,e):
                    qs,qw=mapped(ql,qh)
                    for q,wq in zip(qs,qw):
                        tl=max(box[4],minus_t_min(float(x),float(q)));th=box[5]
                        if th<=tl:continue
                        # Squared distance from the lower t boundary smooths
                        # sqrt(t-t_min) terms without approximating the model.
                        u=(nodes+1)/2;ts=tl+(th-tl)*u*u;tw=(th-tl)*u*weights
                        ps,pw=mapped(*phi_edges)
                        for t,wt in zip(ts,tw):
                            for phi,wp in zip(ps,pw):
                                points.append((float(x),float(q),float(t),float(phi),e))
                                wts.append(float(wx*wq*wt*wp*math.pi/180))
        return points,np.asarray(wts)

@lru_cache(maxsize=32)
def gauss(order):return np.polynomial.legendre.leggauss(order)

def box_from_row(row):return tuple(float(row[k]) for k in ('xBmin','xBmax','Q2min','Q2max','t_abs_min','t_abs_max'))

def find_config(csv_path,explicit=None):
    if explicit:
        p=Path(explicit).resolve()
        if not p.is_file():raise FileNotFoundError(p)
        return p
    for parent in Path(csv_path).resolve().parents:
        for relative in ('jsons/global_cuts_config.json','output/jsons/global_cuts_config.json','global_cuts_config.json'):
            p=parent/relative
            if p.is_file():return p
    raise FileNotFoundError('Saved cuts configuration not found; supply --cuts-json')
