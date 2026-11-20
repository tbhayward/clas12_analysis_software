import unittest
from unittest.mock import patch
import audit_cosine_minima as a

class AuditTests(unittest.TestCase):
    def test_starts_do_not_use_injection(self):
        n={'values':{k:0.1 for k in a.p.PHYSICS_PARAMETERS}}
        x=a.starts_for(n,{},3,17)
        self.assertEqual(x,a.starts_for(n,{},3,17))
        self.assertNotEqual(x,a.starts_for(n,{},3,18))
        self.assertEqual(x[1],{k:0. for k in a.p.PHYSICS_PARAMETERS})

    def test_lower_profile_is_promoted_and_rechecked(self):
        fit=lambda v:{'minimum_nll':v,'values':{},'valid':True}
        sequence=[(fit(10),[]),(fit(9),[]),(fit(11),[]),(fit(8),[]),(fit(9),[]),(fit(10),[])]
        with patch.object(a.p,'make_bin_nll',return_value=(lambda **kw:20,{})), patch.object(a,'starts_for',return_value=[]), patch.object(a,'search',side_effect=sequence), patch.object(a.study,'evaluate_recovery',return_value={'usable':True}):
            f,c,_=a.refit(None,None,None,6,{}, {},{'u1':0.,'ll1':0.},{},1,1,100)
        self.assertEqual(f['minimum_nll'],8)
        self.assertEqual(c['profile_promotions'],1)
        self.assertFalse(c['unresolved_lower_profile'])
        self.assertEqual(c['profiles']['u1']['two_delta_nll'],2)

    def test_persistent_lower_profile_is_excluded(self):
        fit=lambda v:{'minimum_nll':v,'values':{},'valid':True}
        sequence=[(fit(10),[])]
        for cycle in range(3):
            sequence.extend([(fit(7-cycle),[]),(fit(20),[])])
            if cycle<2:sequence.append((fit(9-cycle),[]))
        with patch.object(a.p,'make_bin_nll',return_value=(lambda **kw:20,{})), patch.object(a,'starts_for',return_value=[]), patch.object(a,'search',side_effect=sequence), patch.object(a.study,'evaluate_recovery',return_value={'usable':True}):
            _,c,_=a.refit(None,None,None,6,{}, {},{'u1':0.,'ll1':0.},{},1,1,100)
        self.assertFalse(c['usable'])
        self.assertTrue(c['unresolved_lower_profile'])

if __name__=='__main__':unittest.main()
