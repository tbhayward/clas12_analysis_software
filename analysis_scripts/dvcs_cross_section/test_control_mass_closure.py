#!/usr/bin/env python3
import numpy as np
import study_control_mass_closure as study
import study_unordered_control_response as control

a=np.zeros((4,29));a[:,20]=3;a[:,21]=.1;a[:,22]=.3;a[:,24]=2;a[:,25]=.2;a[:,26]=1.;p=np.array([3*np.sin(.1)*np.cos(.3),3*np.sin(.1)*np.sin(.3),3*np.cos(.1)]);q=np.array([2*np.sin(.2)*np.cos(1),2*np.sin(.2)*np.sin(1),2*np.cos(.2)]);expected=np.sqrt(25-(p+q)@(p+q));a[:,28]=expected
a[1,28]-=.001;a[2,20]=np.nan;a[3,28]-=.002
mass,valid=study.rebuilt_mass(a);assert np.allclose(mass[valid],expected,atol=1e-13)
m,s,r,d=study.metrics(a,np.array([True,True,True,False]));assert m['selected']==3 and m['resolved']==2 and m['unresolved']==1 and np.allclose(d,[0,1],atol=1e-9)
a[:,20:23],a[:,24:27]=a[:,24:27].copy(),a[:,20:23].copy();other,v=study.rebuilt_mass(a);assert np.allclose(other[v],mass[v])
assert control.ACTIVE_GROUPS==(0,3,4,5)
print('PASS: independent invariant-mass formula, zero closure, known mass offset, matched-event selection, unresolved accounting, photon-exchange symmetry and excluded FD/FT groups')
