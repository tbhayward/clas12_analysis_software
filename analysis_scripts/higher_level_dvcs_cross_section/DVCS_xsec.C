
double PI = TMath::Pi();

double alpha = 1./137.036;
double alpha3 = TMath::Power(alpha,3.);
double hbarc2 = 0.38938;//GeV^{2} mbarn

double m = 0.000511;
double M = 0.93827;
double muP = 2.79285;

double GetF1(double T);
double GetF2(double T);
double GetGMP(double tau);
double GetGEP(double tau);

double GetImH(double xi, double t);
double GetImHt(double xi, double t);
double GetImE(double xi, double t);
double GetImEt(double xi, double t);

double GetReH(double xi, double t);
double GetReHt(double xi, double t);
double GetReE(double xi, double t);
double GetReEt(double xi, double t);

bool hasH = false;
bool hasHt = false;
bool hasE = false;
bool hasEt = false;

////////////////////////////////////////////////////////////////////////////////////////////////////////////

class BMK_DVCS{
	public:
		double q_beam, L_beam, L_target;
		double EB, xB, Q2, t, phi, theta_Tpol, phi_Tpol;//primary variables ; EB on fixed target ; Trento convention 

		double xi, nu, y, eps, eps2, phi_BMK, t_min, K2, K, J, Ktild2, Ktilda;//secondary variables ;  phi_BMK = pi - phi_Trento
		double Jacob;//Jacobian from (xB,y) to (xB,Q2)
		double F1, F2, FF_comb1, FF_comb2, FF_comb3;
		double ImH, ImHt, ImE, ImEt;
		double ReH, ReHt, ReE, ReEt;

		bool VERB;

		BMK_DVCS(double rq_beam, double rL_beam, double rL_target, double rEB, double rxB, double rQ2, double rt, double rphi, double rtheta_Tpol=0, double rphi_Tpol=0);
		void setSecondaryVars(void);
		void setPrimaryVars(double rq_beam, double rL_beam, double rL_target, double rEB, double rxB, double rQ2, double rt, double rphi, double rtheta_Tpol=0, double rphi_Tpol=0);

		double CrossSection(void);
		double TPolCrossSection(void);
		double BSA(void);
		double pBSA(void);
		double TLSA(void);
		double TLLSA(void);
		double TTSAx(void);
		double TTSAy(void);
		double TTSSAx(void);
		double TTSSAy(void);
		double BCA(void);
		double BCSA(void);
		double BC0SA(void);
		
		double BCLA(void);
		double BCLLA(void);
		double BCTxA(void);
		double BCTyA(void);

		double T2(void);
		double BH2(void);
		double DVCS2(void);
		double BHDVCS(void);

		double c0_BH(void);
		double c1_BH(void);
		double c2_BH(void);
		double c0_BH_LP(void);
		double c1_BH_LP(void);
		double c0_BH_TP(void);
		double c1_BH_TP(void);
		double s1_BH_TP(void);
		double BHP1(void);
		double BHP2(void);

		double c0_I(void);
		double c1_I(void);
		double s1_I(void);
		double c0_I_LP(void);
		double c1_I_LP(void);
		double s1_I_LP(void);
		double c0_I_TP(void);
		double c1_I_TP(void);
		double s1_I_TP(void);
		
		double c0_DVCS(void);
		double c0_DVCS_LP(void);
		double c0_DVCS_TP(void);
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////

BMK_DVCS::BMK_DVCS(double rq_beam, double rL_beam, double rL_target, double rEB, double rxB, double rQ2, double rt, double rphi, double rtheta_Tpol, double rphi_Tpol){
	this->q_beam=rq_beam;
	this->L_beam = rL_beam;
	this->L_target = rL_target;
	this->EB=rEB;
	this->xB=rxB;
	this->Q2=rQ2;
	this->t=-TMath::Abs(rt);
	this->phi=rphi;
	this->theta_Tpol=PI*rtheta_Tpol/180;
	this->phi_Tpol=PI*rphi_Tpol/180;
	VERB = false;
	this->setSecondaryVars();
}

void BMK_DVCS::setPrimaryVars(double rq_beam, double rL_beam, double rL_target, double rEB, double rxB, double rQ2, double rt, double rphi, double rtheta_Tpol, double rphi_Tpol){
	this->q_beam=rq_beam;
	this->L_beam = rL_beam;
	this->L_target = rL_target;
	this->EB=rEB;
	this->xB=rxB;
	this->Q2=rQ2;
	this->t=-TMath::Abs(rt);
	this->phi=rphi;
	this->theta_Tpol=PI*rtheta_Tpol/180;;
	this->phi_Tpol=PI*rphi_Tpol/180;
	this->setSecondaryVars();
}

void BMK_DVCS::setSecondaryVars(void){
	xi = xB * (1 + 0.5*t/Q2) / ( 2-xB + xB*t/Q2 );
	phi_BMK = PI * ( 1 - phi/180);
	nu = Q2/(2.*M*xB);
	Jacob = y/Q2;
	y = nu / EB;
	eps = 2*xB*M / TMath::Sqrt(Q2);
	eps2 = eps*eps;

	t_min = -Q2 * ( 2.*(1.-xB)*(1.-TMath::Sqrt(1.+eps2))+eps2 ) / (4.*xB*(1.-xB)+eps2);
	if(t_min < t )K2=0;
	else K2 = -(t-t_min)/Q2 * (1.-xB) * (1.-y-0.25*y*y*eps2) * ( TMath::Sqrt(1.+eps2) + (4.*xB*(1.-xB)+eps2)/(4.*(1.-xB)) * (t-t_min)/Q2 );
	K = TMath::Sqrt(K2);
	J = (1.-y-0.5*y*eps2)*(1.+t/Q2)-(1.-xB)*(2.-y)*t/Q2;

	if(t_min < t )Ktild2 = 0;
	else Ktild2 = (t_min-t) * ( (1-xB) * (1+eps2) + (t_min-t) * (eps2 +4*(1-xB)*xB)/(4*Q2) );
	Ktilda = TMath::Sqrt(Ktild2);
	
	F1 = GetF1(t);
	F2 = GetF2(t);
	FF_comb1 = F1*F1 - t*F2*F2/(4*M*M);
	FF_comb2 = TMath::Power(F1+F2,2);
	FF_comb3 = F1 + t*F2/(4*M*M);

	ImH  = GetImH(  xi , t);
	ImHt = GetImHt( xi , t);
	ImE  = GetImE(  xi , t);
	ImEt = GetImEt( xi , t);

	ReH  = GetReH(  xi , t);
	ReHt = GetReHt( xi , t);
	ReE  = GetReE(  xi , t);
	ReEt = GetReEt( xi , t);

	if(VERB){
		cout << "primary kine --- " << EB << " " << xB << " " << Q2 << " " << t << " " << phi << endl;
		cout << "             ---  phi_BMK=" << phi_BMK << " , y=" << y << " , eps2=" << eps2 << " , F1=" << F1 << " , F2=" << F2 << endl;
		cout << "             --- tmin=" << t_min << " , t=" << t << " , J=" << J << " , K2=" << K2 << endl;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

double BMK_DVCS::CrossSection(){
	if(VERB)cout << "prefactor=" << 1e9*hbarc2*2.*PI*Jacob * alpha3 * xB * y / (16*PI*PI * Q2 * TMath::Sqrt( 1 + eps2 ) ) << endl;
	//                              e_phi
	//                mbarn 
	//      to nb                                 | (22) from 012108 = (1.1) from 1005.5209 x Jacobian
	return  1e6   *   hbarc2       * 2*PI  *Jacob * alpha3 * xB * y / (16*PI*PI * Q2 * TMath::Sqrt( 1 + eps2 ) ) * T2() ;
}

double BMK_DVCS::TPolCrossSection(){
	//       (1) from 1212.6674 is differential in phi_Tpol
	return  1e6   *   hbarc2 * alpha3 * xB * y*y / (16*PI*PI * Q2*Q2 * TMath::Sqrt( 1 + eps2 ) ) * T2() ;
}

double BMK_DVCS::BSA(void){
	q_beam = -1;
	L_target = 0;
	L_beam = 1;
	double xsec1 = CrossSection();
	L_beam = -1;
	double xsec2 = CrossSection();
	return (xsec1-xsec2)/(xsec1+xsec2);
}

double BMK_DVCS::pBSA(void){
	q_beam =  1;
	L_target = 0;
	L_beam = 1;
	double xsec1 = CrossSection();
	L_beam = -1;
	double xsec2 = CrossSection();
	return (xsec1-xsec2)/(xsec1+xsec2);
}

double BMK_DVCS::TLSA(void){
	q_beam = -1;
	L_beam = 0;
	L_target = 1;
	theta_Tpol = 0;
	double xsec1 = CrossSection();
	L_target = -1;
	double xsec2 = CrossSection();
	return (xsec1-xsec2)/(xsec1+xsec2);
}

double BMK_DVCS::TLLSA(void){
	q_beam = -1;
	L_beam   =  1;
	L_target =  1;
	theta_Tpol = 0;
	double xsec1 = CrossSection();
	L_beam   = -1;
	L_target = -1;
	double xsec2 = CrossSection();
	L_beam   = -1;
	L_target =  1;
	double xsec3 = CrossSection();
	L_beam   =  1;
	L_target = -1;
	double xsec4 = CrossSection();
	return (xsec1+xsec2-xsec3-xsec4) / (xsec1+xsec2+xsec3+xsec4);
}

double BMK_DVCS::TTSAx(void){
	q_beam = -1;
	L_beam = 0;
	L_target =  1;
	theta_Tpol = PI/2;
	phi_Tpol = 0;
	double xsec1 = CrossSection();
	L_target = -1;
	double xsec2 = CrossSection();
	return (xsec1-xsec2)/(xsec1+xsec2);
}

double BMK_DVCS::TTSAy(void){
	q_beam = -1;
	L_beam = 0;
	L_target =  1;
	theta_Tpol = PI/2;
	phi_Tpol = PI/2;
	double xsec1 = CrossSection();
	L_target = -1;
	double xsec2 = CrossSection();
	return (xsec1-xsec2)/(xsec1+xsec2);
}

double BMK_DVCS::TTSSAx(void){
	q_beam = -1;
	L_beam = 1;
	L_target =  1;
	theta_Tpol = PI/2;
	phi_Tpol = 0;
	double xsec1 = CrossSection();
	L_target = -1;
	double xsec2 = CrossSection();
	L_beam = -1;
	L_target =  1;
	double xsec3 = CrossSection();
	L_target = -1;
	double xsec4 = CrossSection();
	return (xsec1+xsec4-xsec2-xsec3)/(xsec1+xsec2+xsec3+xsec4);
	//return (xsec1-xsec4-xsec2+xsec3)/(xsec1+xsec2+xsec3+xsec4);// is wrong
}

double BMK_DVCS::TTSSAy(void){
	q_beam = -1;
	L_beam = 1;
	L_target =  1;
	theta_Tpol = PI/2;
	phi_Tpol = PI/2;
	double xsec1 = CrossSection();
	L_target = -1;
	double xsec2 = CrossSection();
	L_beam = -1;
	L_target =  1;
	double xsec3 = CrossSection();
	L_target = -1;
	double xsec4 = CrossSection();
	return (xsec1+xsec4-xsec2-xsec3)/(xsec1+xsec2+xsec3+xsec4);
	//return (xsec1-xsec4-xsec2+xsec3)/(xsec1+xsec2+xsec3+xsec4);// is wrong
}

double BMK_DVCS::BCA(void){
	L_beam = 0;
	L_target = 0;
	q_beam = 1;
	double xsec1 = CrossSection();
	q_beam = -1;
	double xsec2 = CrossSection();
	return (xsec1-xsec2)/(xsec1+xsec2);
}

double BMK_DVCS::BCSA(void){
	L_target = 0;
	q_beam = 1;
	L_beam = 1;
	double xsec1 = CrossSection();
	L_beam = -1;
	double xsec2 = CrossSection();
	q_beam = -1;
	L_beam =  1;
	double xsec3 = CrossSection();
	L_beam = -1;
	double xsec4 = CrossSection();
	return ( xsec1-xsec2 -xsec3+xsec4 )/(xsec1+xsec2+xsec3+xsec4);
}

double BMK_DVCS::BC0SA(void){
	L_target = 0;
	q_beam = 1;
	L_beam = 1;
	double xsec1 = CrossSection();
	L_beam = -1;
	double xsec2 = CrossSection();
	q_beam = -1;
	L_beam =  1;
	double xsec3 = CrossSection();
	L_beam = -1;
	double xsec4 = CrossSection();
	return ( xsec1-xsec2 +xsec3-xsec4 )/(xsec1+xsec2+xsec3+xsec4);
}

// HERE
double BMK_DVCS::BCLA(void){
	q_beam = 1;
	L_beam = 0;
	L_target = 1;
	double xsec1 = CrossSection();
	L_target = -1;
	double xsec2 = CrossSection();
	q_beam = -1;
	L_target = 1;
	double xsec3 = CrossSection();
	L_target = -1;
	double xsec4 = CrossSection();
	return ( xsec1-xsec2 -xsec3-+xsec4 )/(xsec1+xsec2+xsec3+xsec4);
}

double BMK_DVCS::BCLLA(void){
	q_beam = 1;
	L_beam = 1;
	L_target = 1;
	double xsec1 = CrossSection();
	L_target = -1;
	double xsec2 = CrossSection();
	q_beam = -1;
	L_target = 1;
	double xsec3 = CrossSection();
	L_target = -1;
	double xsec4 = CrossSection();
	
	q_beam = 1;
	L_beam = -1;
	L_target = 1;
	double xsec5 = CrossSection();
	L_target = -1;
	double xsec6 = CrossSection();
	q_beam = -1;
	L_target = 1;
	double xsec7 = CrossSection();
	L_target = -1;
	double xsec8 = CrossSection();
	return ( xsec1-xsec2 +xsec3-xsec4 -xsec5+xsec6 -xsec7+xsec8 )/(xsec1+xsec2+xsec3+xsec4+xsec5+xsec6+xsec7+xsec8);
	return 0;
}

double BMK_DVCS::BCTxA(void){
	q_beam = 1;
	L_beam = 0;
	theta_Tpol = PI/2;
	phi_Tpol = 0;
	double xsec1 = CrossSection();
	L_target = -1;
	double xsec2 = CrossSection();
	q_beam = -1;
	L_target =  1;
	double xsec3 = CrossSection();
	L_target = -1;
	double xsec4 = CrossSection();
	return (xsec1+xsec4-xsec2-xsec3)/(xsec1+xsec2+xsec3+xsec4);
}

double BMK_DVCS::BCTyA(void){
	q_beam = 1;
	L_beam = 0;
	theta_Tpol = PI/2;
	phi_Tpol = PI/2;
	double xsec1 = CrossSection();
	L_target = -1;
	double xsec2 = CrossSection();
	q_beam = -1;
	L_target =  1;
	double xsec3 = CrossSection();
	L_target = -1;
	double xsec4 = CrossSection();
	return (xsec1+xsec4-xsec2-xsec3)/(xsec1+xsec2+xsec3+xsec4);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

double BMK_DVCS::T2(void){
	// (2.34)  from 1005.5209 and (27) from 012108
	return BH2() + DVCS2() - q_beam * BHDVCS();
}

double BMK_DVCS::BH2(void){
	double DENOM = xB*xB *y*y *t * BHP1() * BHP2() * TMath::Power(1+eps2 , 2) ;
	double res_c0_BH = c0_BH();
	if( L_target != 0 ){
		res_c0_BH += L_target * ( TMath::Cos(theta_Tpol) * c0_BH_LP() + TMath::Sin(theta_Tpol) * c0_BH_TP() );
	}
	double res_c1_BH = c1_BH();
	if( L_target != 0 ){
		res_c1_BH += L_target * ( TMath::Cos(theta_Tpol) * c1_BH_LP() + TMath::Sin(theta_Tpol) * c1_BH_TP() );
	}
	double res_c2_BH = c2_BH();
	double res_s1_BH = 0;
	if( L_target != 0 ){
		res_s1_BH = L_target * TMath::Sin(theta_Tpol) * s1_BH_TP() ;
	}
	double HARMONICS = res_c0_BH + res_c1_BH*TMath::Cos(phi_BMK) + res_c2_BH*TMath::Cos(2*phi_BMK) + res_s1_BH*TMath::Sin(phi_BMK);
	if(VERB)cout << "c0_BH=" << res_c0_BH << " , c1_BH=" << res_c1_BH << " , c2_BH=" << res_c2_BH << " , numerator=" << HARMONICS << " , denom=" << DENOM << endl;
	return HARMONICS / DENOM;
}

double BMK_DVCS::DVCS2(void){
	double DENOM = y*y * Q2;
	double res_c0_DVCS = c0_DVCS();
	if( L_target != 0 ){
		res_c0_DVCS += L_target  * ( TMath::Cos(theta_Tpol) * c0_DVCS_LP() + TMath::Sin(theta_Tpol) * c0_DVCS_TP() );
	}
	//only twist 2 coefficients
	return res_c0_DVCS / DENOM;
}

double BMK_DVCS::BHDVCS(void){
	// (2.34) from 1005.5209
	double DENOM = xB *y*y*y *t * BHP1() * BHP2();
	double res_c0_I = c0_I();
	if( L_target != 0 ){
		res_c0_I +=  L_target * ( TMath::Cos(theta_Tpol) * c0_I_LP() + TMath::Sin(theta_Tpol) * c0_I_TP() );
	}
	double res_c1_I = c1_I();
	if( L_target != 0 ){
		res_c1_I +=  L_target * ( TMath::Cos(theta_Tpol) * c1_I_LP() + TMath::Sin(theta_Tpol) * c1_I_TP() );
	}
	double res_s1_I = s1_I();
	if( L_target != 0 ){
		res_s1_I +=  L_target * ( TMath::Cos(theta_Tpol) * s1_I_LP() + TMath::Sin(theta_Tpol) * s1_I_TP() );
	}
	//only twist 2 coefficients
	double HARMONICS = res_c0_I + res_c1_I*TMath::Cos(phi_BMK) + res_s1_I*TMath::Sin(phi_BMK);
	return HARMONICS / DENOM;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

double BMK_DVCS::c0_BH(void){
	// (35) from 0112108
	double c0_BH_unp = 0;
	c0_BH_unp  = (2+eps2)*( 4*xB*xB*M*M*TMath::Power(1+t/Q2,2)/t + 4*(1-xB)*(1+xB*t/Q2) )*FF_comb1;
	c0_BH_unp += 4*xB*xB*( xB+(1-xB+0.5*eps2)*TMath::Power(1-t/Q2,2)-xB*(1-2*xB)*t*t/(Q2*Q2) )*FF_comb2;
	c0_BH_unp *= TMath::Power(2-y,2.);
	c0_BH_unp += 8*K2*( (2+3*eps2)*Q2*FF_comb1/t + 2*xB*xB*FF_comb2  );
	c0_BH_unp += 8*(1+eps2)*(1-y-eps2*y*y/4)*(  2*eps2*(1-t/(4*M*M))*FF_comb1 - xB*xB*TMath::Power(1-t/Q2,2)*FF_comb2 );
	return c0_BH_unp;
}

double BMK_DVCS::c0_BH_LP(void){
	// (38) from 0112108
	double c0_BH_LP = 0;
	c0_BH_LP  = ( 1 - (1-xB)*t/Q2 ) * ( xB*xB*M*M/t * TMath::Power(1+t/Q2,2) + (1-xB)*(1+xB*t/Q2)  ) * FF_comb3;
	c0_BH_LP += 0.5*(0.5*xB*(1-t/Q2)-0.25*t/(M*M))*( 2-xB-2*TMath::Power(1-xB,2)*t/Q2 + eps2*(1-t/Q2) - xB*(1-2*xB)*t*t/(Q2*Q2) ) * (F1+F2);
	c0_BH_LP  = 8 * L_beam * xB * (2-y) * y * TMath::Sqrt(1+eps2)/( 1-0.25*t/(M*M) ) * (F1+F2) * c0_BH_LP;
	return c0_BH_LP;
}

double BMK_DVCS::c0_BH_TP(void){
	// (40) from 0112108
	double c0_BH_TP = 0;
	c0_BH_TP = xB*xB*xB *M*M/Q2 * ( 1 - t/Q2 ) * (F1+F2) + ( 1-(1-xB)*t/Q2 ) * ( xB*xB*M*M/t *(1-t/Q2)*F1 + xB*F2*0.5 ) ;
	c0_BH_TP = -8 * L_beam * TMath::Cos(phi_Tpol) * (2-y) * y * TMath::Sqrt(Q2)/M *TMath::Sqrt(1+eps2) * K / TMath::Sqrt(1-y-eps2*y*y*0.25) * (F1+F2) * c0_BH_TP;
	return c0_BH_TP;
}

double BMK_DVCS::c1_BH(void){
	// (36) from 0112108
	double c1_BH_unp = 0;
	c1_BH_unp = 8*K*(2.-y)*( (4*xB*xB*M*M/t-2*xB-eps2)*FF_comb1 + 2*xB*xB*(1-(1-2*xB)*t/Q2)*FF_comb2 );
	return c1_BH_unp;
}

double BMK_DVCS::c1_BH_LP(void){
	// (39) from 0112108
	double c1_BH_LP = 0;
	c1_BH_LP  = ( 1+xB-(3-2*xB)*(1+xB*t/Q2)-4*xB*xB*M*M/t*(1+t*t/(Q2*Q2)) ) * FF_comb3;
	c1_BH_LP += ( 0.5*t/(M*M) - xB*(1-t/Q2) ) * (1-xB+xB*t/Q2) * (F1+F2);
	c1_BH_LP  = -8 * L_beam * xB * y * K * TMath::Sqrt(1+eps2)/( 1-0.25*t/(M*M) )  * (F1+F2) * c1_BH_LP;
	return c1_BH_LP;
}

double BMK_DVCS::c1_BH_TP(void){
	// (41) from 0112108
	double c1_BH_TP = 0;
	c1_BH_TP = 2*K2*Q2/( t*(1-y-eps2*y*y*0.25) ) * ( xB*(1-t/Q2)*F1 +0.25*t/(M*M)*F2 ) + (1+eps2)*xB*(1-t/Q2)*FF_comb3;
	c1_BH_TP = -16 * L_beam * TMath::Cos(phi_Tpol) * xB * y * TMath::Sqrt(1-y-eps2*y*y*0.25) * M/TMath::Sqrt(Q2) * TMath::Sqrt(1+eps2) * (F1+F2) * c1_BH_TP;
	return c1_BH_TP;
}

double BMK_DVCS::s1_BH_TP(void){
	// (42) from 0112108
	double s1_BH_TP = 0;
	s1_BH_TP = -16 * L_beam * TMath::Sin(phi_Tpol) * y * xB*xB * TMath::Sqrt(1-y-0.25*eps2*y*y) * M/TMath::Sqrt(Q2) * TMath::Power(1+eps2,1.5) *(1-t/Q2)*(F1+F2)*FF_comb3;
	return s1_BH_TP;
}

double BMK_DVCS::c2_BH(void){
	// (37) from 0112108
	double c2_BH_unp = 0;
	c2_BH_unp = 8*xB*xB*K2*(4*M*M*FF_comb1/t+2*FF_comb2);
	return c2_BH_unp;
}

double BMK_DVCS::BHP1(void){
	double result = - ( J + 2*TMath::Sqrt(K2) * TMath::Cos(phi_BMK) ) / ( y*(1.+eps2) );
	//cout << "IN BHP1" << J << " " << K2 << " , P1=" << result << endl;
	return result;//- ( J + 2*TMath::Sqrt(K2) * TMath::Cos(phi_BMK) ) / ( y*(1.+eps2) ); 
}
double BMK_DVCS::BHP2(void){return 1 + t/Q2 - BHP1(); }

////////////////////////////////////////////////////////////////////////////////////////////////////////////

double BMK_DVCS::c0_I(void){
	double C_pp0_unp0 = 0;
	double C_pp0_unpV = 0;
	double C_pp0_unpA = 0;
	double CFF_unp0 = 0;
	double CFF_unpV = 0;
	double CFF_unpA = 0;
	
	//unpolarized (A.1) from 1005.5209
	// (B.1) from 1212.6674
	C_pp0_unp0 = Ktild2/Q2 * TMath::Power(2-y,2) / TMath::Sqrt(1+eps2);
	C_pp0_unp0 += t/Q2 * (1-y-0.25*eps2*y*y)*(2-xB)*( 1 + ( 2*xB*(2-xB+0.5*(TMath::Sqrt(1+eps2)-1) + 0.5*eps2/xB )*t/Q2 +eps2 )/((2-xB)*(1+TMath::Sqrt(1+eps2))) );
	C_pp0_unp0 *= -4*(2-y)*(1+TMath::Sqrt(1+eps2)) / TMath::Power(1+eps2,2) ;

	C_pp0_unpV = (1-y-0.25*y*y*eps2)* 0.5*(1+TMath::Sqrt(1+eps2))*( 1 + t/Q2 ) * ( 1 + ( TMath::Sqrt(1+eps2) - 1 + 2*xB)/(1+TMath::Sqrt(1+eps2)) * t/Q2 );
	C_pp0_unpV += TMath::Power(2-y,2)*Ktild2 / ( TMath::Sqrt(1+eps2)*Q2 );
	C_pp0_unpV *= 8*(2-y)*xB*t / ( TMath::Power(1+eps2,2)*Q2 );

	C_pp0_unpA = 0.5*(1+TMath::Sqrt(1+eps2) ) * ( 1 + TMath::Sqrt(1+eps2) - xB + ( TMath::Sqrt(1+eps2) - 1 + xB*(3+TMath::Sqrt(1+eps2)-2*xB)/(1+TMath::Sqrt(1+eps2) ) ) * t/Q2 ) -2*Ktild2/Q2;
	C_pp0_unpA = 8*(2-y)/TMath::Power(1+eps2,2) * t/Q2 * ( TMath::Power(2-y,2) * Ktild2/(TMath::Sqrt(1+eps2)*Q2) * 0.5*(1+TMath::Sqrt(1+eps2)-2*xB) + (1-y-0.25*eps2*y*y) * C_pp0_unpA );

	// (83) (84) (85) from 1212.6674
	CFF_unp0 = F1 * ReH - t/(4*M*M) * F2 * ReE + xB * ( F1 +F2 ) * ReHt / ( 2-xB + xB*t/Q2 ) ;
	CFF_unpV = xB * (F1+F2) * ( ReH + ReE )  / ( 2-xB + xB*t/Q2 );
	CFF_unpA = xB * (F1+F2) *      ReHt      / ( 2-xB + xB*t/Q2 );

	return C_pp0_unp0 * CFF_unp0 + C_pp0_unpV * CFF_unpV + C_pp0_unpA * CFF_unpA;
}

double BMK_DVCS::c0_I_LP(void){
	double C_pp0_LP0 = 0;
	double C_pp0_LPV = 0;
	double C_pp0_LPA = 0;
	double CFF_LP0 = 0;
	double CFF_LPV = 0;
	double CFF_LPA = 0;

	//polarized (A.5) from 1005.5209
	C_pp0_LP0 = TMath::Power(2-y,2)*Ktild2/Q2 + (1-y-0.25*eps2*y*y)*(xB*t/Q2-(1-t/Q2)*eps2*0.5)*( 1 + ( TMath::Sqrt(1+eps2) -1+2*xB )/(1+TMath::Sqrt(1+eps2))*t/Q2 );
	C_pp0_LP0 = -4*L_beam*y*(1+TMath::Sqrt(1+eps2))/TMath::Power(1+eps2,2.5)*C_pp0_LP0;

	C_pp0_LPV = (2-xB+3*eps2/2) * (1 + (4*(1-xB)*xB+eps2)/(4-2*xB+3*eps2) *t/Q2 ) * ( 1 + (TMath::Sqrt(1+eps2)-1+2*xB)/(1+TMath::Sqrt(1+eps2)) *t/Q2 );
	C_pp0_LPV = TMath::Power(2-y,2) * ( 1 + TMath::Sqrt(1+eps2) + 2*xB )/(1+TMath::Sqrt(1+eps2)) * Ktild2/Q2 + (1-y-0.25*eps2*y*y)*C_pp0_LPV;
	C_pp0_LPV = 4*L_beam*y*(1+TMath::Sqrt(1+eps2))/TMath::Power(1+eps2,2.5) * t/Q2 * C_pp0_LPV;

	C_pp0_LPA = (1+TMath::Sqrt(1+eps2)) * (1 - (1-2*xB)*t/Q2 ) * ( 1 + (TMath::Sqrt(1+eps2)-1+2*xB)/(1+TMath::Sqrt(1+eps2)) *t/Q2 );
	C_pp0_LPA = 4*L_beam*y/TMath::Power(1+eps2,2.5) * xB*t/Q2 * ( 2*TMath::Power(2-y,2)*Ktild2/Q2 + (1-y-0.25*eps2*y*y)*C_pp0_LPA );

	// (86) (87) (88) from 1212.6674
	CFF_LP0  = xB / (2-xB+xB*t/Q2) * (F1+F2) * ( ReH + 0.5*xB*(1-t/Q2) * ReE -(1-2*xB)*t*ReHt/Q2 -0.25*t*ReEt/(M*M) );
	CFF_LP0 += 2/ (2-xB+xB*t/Q2) * F1 * ( ( (1-xB)*(1+xB*t/Q2) + 0.5*xB +xB*xB*M*M*(3+t/Q2)/Q2 ) * ReHt + 0.5*xB *( 0.25*t/(M*M) - 0.5*xB*(1-t/Q2) ) * ReEt );

	CFF_LPV = xB / (2-xB+xB*t/Q2) * (F1+F2) * ( ReH + 0.5*xB*(1-t/Q2) * ReE );

	CFF_LPA = xB / (2-xB+xB*t/Q2) * (F1+F2) * ( ReHt + 2*xB*M*M*ReHt/Q2 + xB*ReEt*0.5 );
	
	return C_pp0_LP0 * CFF_LP0 + C_pp0_LPV * CFF_LPV + C_pp0_LPA * CFF_LPA;
}

double BMK_DVCS::c0_I_TP(void){

	// (71) from 0112108
	double CI_TPM = ( xB*xB*F1 - (1-xB)*t*F2/(M*M) )/(2-xB) * ImH + ( 0.25*t/(M*M)*( (2-xB)*F1 + xB*xB*F2/(2-xB) ) + xB*xB*F1/(2-xB)  ) * ImE - xB*xB/(2-xB) * (F1+F2) * (ImHt + 0.25*t/(M*M)*ImEt);
	double CI_TPP = (F1+F2)*( xB*xB/(2-xB)*(ReH+0.5*xB*ReE) + 0.25*xB*t/(M*M)*ReE ) - xB*xB*F1/(2-xB) *(ReHt+0.5*xB*ReEt) + 0.25*t/(M*M)*( 4*(1-xB)*F2*ReHt/(2-xB) - ( xB*F1+xB*xB*F2/(2-xB) )*ReEt );

	// (61) from 0112108
	double res_c0_I_TP = (2-y)*TMath::Sin(phi_Tpol) * TMath::Power(2-y,2)/(1-y) * CI_TPM - L_beam*y*TMath::Cos(phi_Tpol) * ( TMath::Power(2-y,2)/(1-y) + 2 ) * CI_TPP;
	res_c0_I_TP = 8*M*TMath::Sqrt(1-y)*K/TMath::Sqrt(Q2) * res_c0_I_TP;
	return res_c0_I_TP;
}

double BMK_DVCS::c1_I(void){
	double C_pp1_unp0 = 0;
	double C_pp1_unpV = 0;
	double C_pp1_unpA = 0;
	double CFF_unp0 = 0;
	double CFF_unpV = 0;
	double CFF_unpA = 0;

	//unpolarized (A.1) from 1005.5209
	C_pp1_unp0 = 1 - (1-3*xB)*t/Q2 + ( 1-TMath::Sqrt(1+eps2)+3*eps2 )/( 1+TMath::Sqrt(1+eps2)-eps2 ) * xB*t/Q2;
	C_pp1_unp0 = -4*K*( 2 - 2*y + y*y + 0.5*eps2*y*y ) * ( 1 + TMath::Sqrt(1+eps2) - eps2 )/TMath::Power(1+eps2,2.5) * C_pp1_unp0;
	C_pp1_unp0 = -16*K * (1-y-0.25*y*y*eps2)/TMath::Power(1+eps2,2.5) * ( ( 1 + (1-xB) * 0.5*(TMath::Sqrt(1+eps2)-1)/xB ) * xB*t/Q2 - 0.75*eps2 ) + C_pp1_unp0 ; 

	C_pp1_unpV = 16*K*xB*t / (Q2*TMath::Power(1+eps2,2.5)) * ( TMath::Power(2-y,2)*(1-(1-2*xB)*t/Q2) + (1-y-0.25*eps2*y*y) * 0.5*(1+TMath::Sqrt(1+eps2)-2*xB)*(t-t_min)/Q2 );

	C_pp1_unpA = -TMath::Power(2-y,2) * (1-0.5*xB + 0.25*(1+TMath::Sqrt(1+eps2)-2*xB)*(1-t/Q2) + (4*xB*(1-xB)+eps2)* 0.5/TMath::Sqrt(1+eps2) *(t-t_min)/Q2 ) ;
	C_pp1_unpA = -16*K  *t / (Q2*TMath::Power(1+eps2,2)) * ( (1-y-0.25*eps2*y*y) * ( 1-(1-2*xB)*t/Q2 + (4*xB*(1-xB)+eps2)*0.25/TMath::Sqrt(1+eps2) *(t-t_min)/Q2 ) + C_pp1_unpA );

	CFF_unp0 = F1 * ReH - t/(4*M*M) * F2 * ReE + xB * ( F1 +F2 ) * ReHt / ( 2-xB + xB*t/Q2 ) ;
	CFF_unpV = xB * (F1+F2) * ( ReH + ReE )  / ( 2-xB + xB*t/Q2 );
	CFF_unpA = xB * (F1+F2) *      ReHt      / ( 2-xB + xB*t/Q2 );

	return C_pp1_unp0 * CFF_unp0 + C_pp1_unpV * CFF_unpV + C_pp1_unpA * CFF_unpA;
}

double BMK_DVCS::c1_I_LP(void){
	double C_pp1_LP0 = 0;
	double C_pp1_LPV = 0;
	double C_pp1_LPA = 0;
	double CFF_LP0 = 0;
	double CFF_LPV = 0;
	double CFF_LPA = 0;

	//polarized (A.5) from 1005.5209
	C_pp1_LP0 = 1 - (1 - 2*xB*(2+TMath::Sqrt(1+eps2))/(1+TMath::Sqrt(1+eps2)-eps2) ) * t/Q2;
	C_pp1_LP0 = -4*L_beam*K*y*(2-y)/TMath::Power(1+eps2,2.5) * (1+TMath::Sqrt(1+eps2)-eps2) * C_pp1_LP0;
	//C_pp1_LP0 = -4*L_target*K*y*(2-y)/TMath::Power(1+eps2,2.5) * (1+TMath::Sqrt(1+eps2)-eps2) * C_pp1_LP0;

	C_pp1_LPV = 1 - ( 1 + (1-eps2)/TMath::Sqrt(1+eps2) -2*xB*(1 + 4*(1-xB)/TMath::Sqrt(1+eps2)) ) / ( 2*(TMath::Sqrt(1+eps2)+2*(1-xB)) ) * (t-t_min)/Q2;
	C_pp1_LPV = 8*L_beam*K*y*(2-y)/TMath::Power(1+eps2,2) * (TMath::Sqrt(1+eps2) + 2*(1-xB) ) * t/Q2 * C_pp1_LPV;
	//C_pp1_LPV = 8*L_target*K*y*(2-y)/TMath::Power(1+eps2,2) * (TMath::Sqrt(1+eps2) + 2*(1-xB) ) * t/Q2 * C_pp1_LPV;

	C_pp1_LPA = 16*L_beam*K*y*(2-y)/TMath::Power(1+eps2,2.5) * xB*t/Q2 * (1-(1-2*xB)*t/Q2);
	//C_pp1_LPA = 16*L_target*K*y*(2-y)/TMath::Power(1+eps2,2.5) * xB*t/Q2 * (1-(1-2*xB)*t/Q2);

	CFF_LP0 = xB / (2-xB+xB*t/Q2) * (F1+F2) * ( ReH + 0.5*xB*(1-t/Q2) * ReE ) + ( 1 + M*M/Q2 * xB*xB*(3+t/Q2)/(2-xB+xB*t/Q2) ) * F1 * ReHt;
	CFF_LP0 += -t/Q2 * (2*xB*(1-2*xB))/(2-xB+xB*t/Q2) * F2*ReHt - xB/(2-xB+xB*t/Q2) * ( 0.5*xB*(1-t/Q2)*F1 + 0.25*t*F2/(M*M) ) * ReEt;

	CFF_LPV = xB / (2-xB+xB*t/Q2) * (F1+F2) * ( ReH + 0.5*xB*(1-t/Q2) * ReE );

	CFF_LPA = xB / (2-xB+xB*t/Q2) * (F1+F2) * ( ReHt + 2*xB*M*M*ReHt/Q2 + xB*ReEt*0.5 );

	return C_pp1_LP0 * CFF_LP0 + C_pp1_LPV * CFF_LPV + C_pp1_LPA * CFF_LPA;
}

double BMK_DVCS::c1_I_TP(void){
	// (71) from 0112108
	double CI_TPM = ( xB*xB*F1 - (1-xB)*t*F2/(M*M) )/(2-xB) * ImH + ( 0.25*t/(M*M)*( (2-xB)*F1 + xB*xB*F2/(2-xB) ) + xB*xB*F1/(2-xB)  ) * ImE - xB*xB/(2-xB) * (F1+F2) * (ImHt + 0.25*t/(M*M)*ImEt);//IM
	double CI_TPP = (F1+F2)*( xB*xB/(2-xB)*(ReH+0.5*xB*ReE) + 0.25*xB*t/(M*M)*ReE ) - xB*xB*F1/(2-xB) *(ReHt+0.5*xB*ReEt) + 0.25*t/(M*M)*( 4*(1-xB)*F2*ReHt/(2-xB) - ( xB*F1+xB*xB*F2/(2-xB) )*ReEt );//RE

	// (62) from 0112108
	double res_c1_I_TP = TMath::Sin(phi_Tpol) * (2-2*y+y*y)*CI_TPM - L_beam * y * (2-y) * TMath::Cos(phi_Tpol) * CI_TPP;
	res_c1_I_TP = 8*M*TMath::Sqrt(1-y)/TMath::Sqrt(Q2) * res_c1_I_TP;
	return res_c1_I_TP;	
}

double BMK_DVCS::s1_I(void){
	double S_pp1_unp0 =0;
	double S_pp1_unpV =0;
	double S_pp1_unpA =0;
	double CFF_unp0 =0;
	double CFF_unpV =0;
	double CFF_unpA =0;

	//unpolarized (A.1) from 1005.5209
	S_pp1_unp0 =  8*L_beam*K*(2-y)*y/(1+eps2) * ( 1 + ( 1-xB + 0.5*( TMath::Sqrt(1+eps2)-1) ) / (1+eps2) * (t-t_min)/Q2 );

	S_pp1_unpV = -8*L_beam*K*(2-y)*y/TMath::Power(1+eps2,2) * xB*t/Q2 * ( TMath::Sqrt(1+eps2)-1 + (1+TMath::Sqrt(1+eps2)-2*xB) * t/Q2 );

	S_pp1_unpA =  8*L_beam*K*(2-y)*y/(1+eps2) * t/Q2 * ( 1 - (1-2*xB)*(1+TMath::Sqrt(1+eps2)-2*xB)*0.5/TMath::Sqrt(1+eps2) * (t-t_min)/Q2 ) ;

	CFF_unp0 = F1 * ImH - t/(4*M*M) * F2 * ImE + xB * ( F1 +F2 ) * ImHt / ( 2-xB + xB*t/Q2 ) ;
	CFF_unpV = xB * (F1+F2) * ( ImH + ImE )  / ( 2-xB + xB*t/Q2 );
	CFF_unpA = xB * (F1+F2) *      ImHt      / ( 2-xB + xB*t/Q2 );
	
	return S_pp1_unp0 * CFF_unp0 + S_pp1_unpV * CFF_unpV + S_pp1_unpA * CFF_unpA;
}

double BMK_DVCS::s1_I_LP(void){
	double S_pp1_LP0 =0;
	double S_pp1_LPV =0;
	double S_pp1_LPA =0;
	double CFF_LP0 =0;
	double CFF_LPV =0;
	double CFF_LPA =0;

	S_pp1_LP0 = 4*K*(2-2*y+y*y+0.5*eps2*y*y)/TMath::Power(1+eps2,3)*(1+TMath::Sqrt(1+eps2)) * ( 2*TMath::Sqrt(1+eps2)-1 +(1+TMath::Sqrt(1+eps2)-2*xB)/(1+TMath::Sqrt(1+eps2)) * t/Q2 );
	S_pp1_LP0 += 8*K*(2-2*y+0.5*eps2*y*y)/TMath::Power(1+eps2,3)*( 3*eps2/2 +(1-TMath::Sqrt(1+eps2)-0.5*eps2-xB*(3-TMath::Sqrt(1+eps2))) * t/Q2);

	S_pp1_LPV = 1 - ( 1 - TMath::Sqrt(1+eps2) + 0.5*eps2 -2*xB*( 3*(1-xB)-TMath::Sqrt(1+eps2) ) ) / ( 4 - xB*(TMath::Sqrt(1+eps2)+3) +2.5*eps2 ) * t/Q2; 
	S_pp1_LPV *= 32*K*(1-y+0.25*eps2*y*y)/TMath::Power(1+eps2,3) * (1 - 0.25*(3+TMath::Sqrt(1+eps2)) * xB +5*eps2/8) *t/Q2;
	S_pp1_LPV += 8*K*(2-2*y+y*y+0.5*eps2*y*y)/TMath::Power(1+eps2,2) * t/Q2 * (1 - (1-2*xB)*(1+TMath::Sqrt(1+eps2)-2*xB)/(2*(1+eps2)) * (t-t_min)/Q2 );

	S_pp1_LPA  = -8*K*(2-2*y+y*y+0.5*eps2*y*y)/TMath::Power(1+eps2,3) * xB*t/Q2 * (TMath::Sqrt(1+eps2)-1+(1+TMath::Sqrt(1+eps2)-2*xB)*t/Q2);
	S_pp1_LPA +=  8*K*(1-y+0.25*eps2*y*y)/TMath::Power(1+eps2,3) * (3+TMath::Sqrt(1+eps2)) * xB*t/Q2 * (1 - (3-TMath::Sqrt(1+eps2)-6*xB)/(3+TMath::Sqrt(1+eps2)) *t/Q2 );

	CFF_LP0 = xB / (2-xB+xB*t/Q2) * (F1+F2) * ( ImH + 0.5*xB*(1-t/Q2) * ImE ) + ( 1 + M*M/Q2 * xB*xB*(3+t/Q2)/(2-xB+xB*t/Q2) ) * F1 * ImHt;
	CFF_LP0 += -t/Q2 * (2*xB*(1-2*xB))/(2-xB+xB*t/Q2) * F2*ImHt - xB/(2-xB+xB*t/Q2) * ( 0.5*xB*(1-t/Q2)*F1 + 0.25*t*F2/(M*M) ) * ImEt;

	CFF_LPV = xB / (2-xB+xB*t/Q2) * (F1+F2) * ( ImH + 0.5*xB*(1-t/Q2) * ImE );

	CFF_LPA = xB / (2-xB+xB*t/Q2) * (F1+F2) * ( ImHt + 2*xB*M*M*ImHt/Q2 + xB*ImEt*0.5 );

       	return S_pp1_LP0 * CFF_LP0 + S_pp1_LPV * CFF_LPV + S_pp1_LPA * CFF_LPA;
}
		
double BMK_DVCS::s1_I_TP(void){
	// (71) from 0112108
	double CI_TPM = ( xB*xB*F1 - (1-xB)*t*F2/(M*M) )/(2-xB) * ReH + ( 0.25*t/(M*M)*( (2-xB)*F1 + xB*xB*F2/(2-xB) ) + xB*xB*F1/(2-xB)  ) * ReE - xB*xB/(2-xB) * (F1+F2) * (ReHt + 0.25*t/(M*M)*ReEt);//RE
	double CI_TPP = (F1+F2)*( xB*xB/(2-xB)*(ImH+0.5*xB*ImE) + 0.25*xB*t/(M*M)*ImE ) - xB*xB*F1/(2-xB) *(ImHt+0.5*xB*ImEt) + 0.25*t/(M*M)*( 4*(1-xB)*F2*ImHt/(2-xB) - ( xB*F1+xB*xB*F2/(2-xB) )*ImEt );//IM

	// (62) from 0112108
	double res_s1_I_TP = TMath::Cos(phi_Tpol) * (2-2*y+y*y)*CI_TPP + L_beam * y * (2-y) * TMath::Sin(phi_Tpol) * CI_TPM;
	res_s1_I_TP = 8*M*TMath::Sqrt(1-y)/TMath::Sqrt(Q2) * res_s1_I_TP;
	return res_s1_I_TP;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

double BMK_DVCS::c0_DVCS(void){
	double C_0_unp = 0;
	double C_DVCS_unp = 0;

	// (2.18) from 1005.5209
	// (37) from 1212.6674
	C_0_unp = 2 * (2 - 2*y + y*y + 0.5*eps2*y*y ) / ( 1 + eps2 );


	// 48 from 1212.6674
	C_DVCS_unp = 4*(1-xB)*(1+xB*t/Q2)/TMath::Power(2-xB+xB*t/Q2,2) * ( ReH*ReH + ImH*ImH + ReHt*ReHt + ImHt*ImHt);
	C_DVCS_unp += (2+t/Q2)*eps2 /TMath::Power(2-xB+xB*t/Q2,2) * ( ReHt*ReHt + ImHt*ImHt ) - 0.25*t/(M*M) * ( ReE*ReE + ImE*ImE );
	C_DVCS_unp += -xB*xB/TMath::Power(2-xB+xB*t/Q2,2) * ( TMath::Power(1+t/Q2,2) * (2*ReH*ImE+2*ReE*ImH+ReE*ReE+ImE*ImE) + (2*ReHt*ImEt+2*ReEt*ImHt) + 0.25*t/(M*M) * ( ReE*ReE + ImE*ImE ) );

	return C_0_unp * C_DVCS_unp;
}

double BMK_DVCS::c0_DVCS_LP(void){
	double C_0_LP = 0;
	double C_DVCS_LP = 0;

	// (2.20) from 1005.5209
	// (40) from 1212.6674
	C_0_LP = 2 * L_beam * y * (2-y) / TMath::Sqrt(1+eps2);

	// (49) from 1212.6674
	C_DVCS_LP = ( 4*(1-xB)*(1+xB*t/Q2) + 2*( 1-xB+0.5*(1+t/Q2) )*eps2 )/TMath::Power(2-xB+xB*t/Q2,2) * ( 2*ReH*ImHt + 2*ReHt*ImH );
	C_DVCS_LP += -( xB*xB*(1+xB*t/Q2-(1-xB)*t/Q2)*t/Q2  )/TMath::Power(2-xB+xB*t/Q2,2) * ( 2*ReH*ImEt + 2*ReEt*ImH + 2*ReHt*ImE + 2*ReE*ImHt );
	C_DVCS_LP += -0.5*( 4*xB*(1-xB)*(1+xB*t/Q2)*t/Q2 + xB*TMath::Power(1+t/Q2,2)*eps2 )/TMath::Power(2-xB+xB*t/Q2,2) * ( 2*ReHt*ImE + 2*ReE*ImHt );
	C_DVCS_LP += -xB*( 0.5*xB*xB*TMath::Power(1+t/Q2,2)/TMath::Power(2-xB+xB*t/Q2,2) + 0.25*t/(M*M) )/TMath::Power(2-xB+xB*t/Q2,2) * ( 2*ReE*ImEt + 2*ReEt*ImE );

	return C_0_LP * C_DVCS_LP;
}

double BMK_DVCS::c0_DVCS_TP(void){
	double C_0_TPP = 0;
	double C_0_TPM = 0;
	double C_DVCS_TPP = 0;
	double C_DVCS_TPM = 0;

	// (43) from 1212.6674
	C_0_TPP =  (2-y)/(1+eps2) * Ktilda/M * L_beam * TMath::Cos(phi_Tpol);

	C_0_TPM = -(2-y)/(1+eps2) * Ktilda/M * TMath::Sin(phi_Tpol) * (2-2*y+y*y+0.5*eps2*y*y)/(2-y);

	// (50) from 1212.6674
	C_DVCS_TPP = xB*(2*ReH*ImEt+2*ReEt*ImH) + 4*xB*(1-2*xB)*M*M/Q2*(2*ReH*ImHt+2*ReHt*ImH) - (2-xB+xB*t/Q2+0.5*(3+t/Q2)*eps2)*(2*ReHt*ImE+2*ReE*ImHt) + 0.5*xB*xB*(1-t/Q2)*(2*ReE*ImEt+2*ReEt*ImE) ;
	C_DVCS_TPP = C_DVCS_TPP * 2/TMath::Power(2-xB+xB*t/Q2,2);

	// (51) from 1212.6674
	C_DVCS_TPM = ( 2 * ( 2*ImH*ReE - 2*ReH*ImE ) - 2*xB*( 2*ImHt*ReEt - 2*ReHt*ImEt ) )/TMath::Power(2-xB+xB*t/Q2,2);

	return C_0_TPP * C_DVCS_TPP + C_0_TPM * C_DVCS_TPM;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

double GetF1(double T){
	T = TMath::Abs(T);
	double tau = T/(4*M*M);
	return (GetGEP(tau)+tau*GetGMP(tau))/(1+tau);
}

double GetF2(double T){
	T = TMath::Abs(T);
	double tau = T/(4*M*M);
	return (GetGMP(tau)-GetGEP(tau))/(1+tau);
}

double GetGMP(double tau){
	double res;
	double Q2 = 4*M*M*tau;
	double tcut = 4*0.13957*0.13957;
	double t0 = -0.7;
	double z = ( TMath::Sqrt(tcut+Q2) - TMath::Sqrt(tcut-t0) ) / ( TMath::Sqrt(tcut+Q2) + TMath::Sqrt(tcut-t0) ) ;
	res = 0;
	res += 0.264142994136;
	res += -1.095306122120*z;
	res += 1.218553781780 * TMath::Power(z,2);
	res += 0.661136493537 * TMath::Power(z,3);
	res += -1.405678925030 * TMath::Power(z,4);
	res += -1.356418438880 * TMath::Power(z,5);
	res += 1.447029155340 * TMath::Power(z,6);
	res += 4.235669735900 * TMath::Power(z,7);
	res += -5.334045653410 * TMath::Power(z,8);
	res += -2.916300520960 * TMath::Power(z,9);
	res += 8.707403067570 * TMath::Power(z,10);
	res += -5.706999943750 * TMath::Power(z,11);
	res += 1.280814375890 * TMath::Power(z,12);
	return muP*res;
}

double GetGEP(double tau){
	double res;
	double Q2 = 4*M*M*tau;
	double tcut = 4*0.13957*0.13957;
	double t0 = -0.7;
	double z = ( TMath::Sqrt(tcut+Q2) - TMath::Sqrt(tcut-t0) ) / ( TMath::Sqrt(tcut+Q2) + TMath::Sqrt(tcut-t0) ) ;
	res = 0;
	res += 0.239163298067;
	res += -1.109858574410 * z;
	res += 1.444380813060 * TMath::Power(z,2);
	res += 0.479569465603 * TMath::Power(z,3);
	res += -2.286894741870 * TMath::Power(z,4);
	res += 1.126632984980 * TMath::Power(z,5);
	res += 1.250619843540 * TMath::Power(z,6);
	res += -3.631020471590 * TMath::Power(z,7);
	res += 4.082217023790 * TMath::Power(z,8);
	res += 0.504097346499 * TMath::Power(z,9);
	res += -5.085120460510 * TMath::Power(z,10);
	res += 3.967742543950 * TMath::Power(z,11);
	res += -0.981529071103 * TMath::Power(z,12);
	return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

double renormImag = 1;

double GetImH(double xi, double t){
	//return 0;
	double res = 0;
	if(hasH){
		res=(1-xi)/t;
		// VALENCE
		double r = 0.9;
		double alpha = 0.43 + 0.85*t;
		double n = 1.35;
		double b = 0.4;
		double Mm2 = 0.64;
		double P = 1;
		res = TMath::Pi()*5.0/9.0 * n * r /(1+xi) * TMath::Power( 2*xi/(1+xi),-alpha ) * TMath::Power( (1-xi)/(1+xi) , b ) * TMath::Power( 1 - (1-xi)/(1+xi) *t/Mm2 , -P );
		// SEA
		if(false){
			r = 1;
			alpha = 1.13 + 0.15*t;
			n = 1.5;
			b = 0.4;
			Mm2 = 0.5;
			P = 2;
			res += TMath::Pi()*2.0/9.0 * n * r /(1+xi) * TMath::Power( 2*xi/(1+xi),-alpha ) * TMath::Power( (1-xi)/(1+xi) , b ) * TMath::Power( 1 - (1-xi)/(1+xi) *t/Mm2 , -P );
		}
		res *= 2.0;// correction from VGG
	}
	return res * renormImag;
}

double GetImHt(double xi, double t){
	//return 0;
	double res = 0;
	if(hasHt){
		res=0.1*(1-xi)/t;
		double r = 7;
		double alpha = 0.43 + 0.85*t;
		double n = 0.6;
		double b = 2;
		double Mm2 = 0.8;
		double P = 1;
		res = TMath::Pi()*5.0/9.0 * n * r /(1+xi) * TMath::Power( 2*xi/(1+xi),-alpha ) * TMath::Power( (1-xi)/(1+xi) , b ) * TMath::Power( 1 - (1-xi)/(1+xi) *t/Mm2 , -P );
		res *= 0.4;// correction from VGG
	}
	return res * renormImag;
}

double GetImE(double xi, double t){
	//return 0;
	// from HERMES CFF paper 1301.1230 argue ImE = 0.5 x ImH
	double res = 0;
	if(hasE){
		res=1.5*(1-xi)/t;
		// VALENCE
		double r = 0.9;
		double alpha = 0.43 + 0.85*t;
		double n = 1.35;
		double b = 0.4;
		double Mm2 = 0.64;
		double P = 1;
		res = TMath::Pi()*5.0/9.0 * n * r /(1+xi) * TMath::Power( 2*xi/(1+xi),-alpha ) * TMath::Power( (1-xi)/(1+xi) , b ) * TMath::Power( 1 - (1-xi)/(1+xi) *t/Mm2 , -P );
		if(false){
			// SEA
			r = 1;
			alpha = 1.13 + 0.15*t;
			n = 1.5;
			b = 0.4;
			Mm2 = 0.5;
			P = 2;
			res += TMath::Pi()*2.0/9.0 * n * r /(1+xi) * TMath::Power( 2*xi/(1+xi),-alpha ) * TMath::Power( (1-xi)/(1+xi) , b ) * TMath::Power( 1 - (1-xi)/(1+xi) *t/Mm2 , -P );
		}
	}
	return res * renormImag;
}

double GetImEt(double xi, double t){
	return 0;
}

double renormReal = 1.0; //test Vadim comparison

double GetReH(double xi, double t){
	double res = 0;
	if(hasH){
		// in theory PV \int ImH, to be checked
		res  = -12*xi*TMath::Power(1-xi,2)*TMath::Sqrt(TMath::Abs(t))/TMath::Power(1-t/0.7,2);
		// D-term
		res += -3*TMath::Power(1-xi,4)/TMath::Power(1-t/1.1,2); 
	}
	return res/( 1 + TMath::Power( t/0.8 , 4) );
	return res*renormReal/( 1 + TMath::Power( t , 4) );
}

double GetReHt(double xi, double t){
	double res = 0;
	if(hasHt)res=-12*xi*TMath::Power(1-xi,2) /TMath::Power(1-t/1.5,2);
	return res*renormReal;
}

double GetReE(double xi, double t){
	double res = 0;
	if(hasE){
		res = -7 * xi*TMath::Power(1-xi,2) *TMath::Sqrt(TMath::Abs(t))/TMath::Power(1-t/0.7,2);
		res += -3*TMath::Power(1-xi,2)/TMath::Power(1-t/1.2,2);

	}
	return res*renormReal / ( 1 + TMath::Power( t , 4) );
	return res*renormReal;
}

double GetReEt(double xi, double t){
	//return 0;
	double res = 0;
	if(hasEt)res=10/t * 1/( 1 + TMath::Power( 3*xi , 4) );
	return res*renormReal;
}

