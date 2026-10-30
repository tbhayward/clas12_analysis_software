"""Authoritative positional schema from processing_dvcs.groovy and processing_exclusive_pi0.groovy."""
COMMON = '''fiducial_status num_pos num_neg num_neutrals runnum evnum helicity detector1 detector2 e_p e_theta e_phi vz_e p1_p p1_theta p1_phi vz_p1 p2_p p2_theta p2_phi vz_p2 open_angle_ep open_angle_ep1 open_angle_ep2 open_angle_p1p2 Q2 W Mx2 Mx2_1 Mx2_2 x t t1 t2 tmin y z z1 z2 Mh xF xF1 xF2 pT pT1 pT2 pTpT xi xi1 xi2 eta eta1 eta2 Delta_eta eta1_gN eta2_gN phi1 phi2 Delta_phi phih phiR theta DepA DepB DepC DepV DepW'''.split()
DVCS = COMMON + 'Emiss2 theta_gamma_gamma pTmiss'.split()
PI0 = COMMON + 'Mh_gammagamma detector_gamma1 detector_gamma2 open_angle_egamma1 open_angle_egamma2 gamma_phi1 gamma_phi2 Emiss2 theta_pi0_pi0 pTmiss'.split()
assert len(DVCS) == 70 and len(PI0) == 77
SCHEMAS = {'dvcs_data': DVCS, 'pi0_data': PI0, 'pi0_as_dvcs': DVCS}
