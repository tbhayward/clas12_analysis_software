import numpy as np
from study_production_pi0_audit import calculate,triple,sums,TOPS
assert np.array_equal(triple(''),np.zeros(3))
a=calculate(np.array([10,1,0]),np.array([20,2,0]),np.array([40,4,0]),np.array([100,10,0]));assert np.isclose(a[0],.05) and np.isclose(a[1],.01)
assert np.array_equal(calculate(*[np.zeros(3)]*4),np.zeros(3))
row={f'normalized raw yield, ep->epg, {t}, exp, Fa18 Inb, unpol':'(10,2,0)' for t in TOPS}
assert np.allclose(sums(row,'ep->epg','Fa18 Inb',True),[30,np.sqrt(12),0])
print('PASS: production triples, topology sums, zero-input behavior and statistical propagation')

for text in [' (1, 2, 3) ', '"(1,2,3)"', '(1e3,2e-2,0)', '1,2,3']:
    assert len(triple(text))==3
for text in ['(1oops,2,3)','(1,2)','(nan,2,3)','(1,-2,0)']:
    try:triple(text)
    except ValueError:pass
    else:raise AssertionError('Accepted malformed yield '+text)
print('PASS: production quoted/scientific triples; strict rejection of malformed/nonfinite/negative uncertainty')
