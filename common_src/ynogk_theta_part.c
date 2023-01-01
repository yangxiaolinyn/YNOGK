/*
 * common_src/yongk_theta_part.c
 *
 * This module aim on solve cubic and quartic polynomial equations.
 * One can use these subroutine root3 and root4 to find roots of cubic and 
 * quartic equations respectively. 
 *
 * Author       Yang Xiao-lin
 *
 * Date         2022-12-15
 *
 * City         Kun-ming Yunnan Provice, China
 *
 * 2022-12-15   Stating the writting of the module.
 */
 

#include "ynogk_theta_part.h"


/*****************************************************************************\
!*
!*  Private Variables for this module.
!*
\*****************************************************************************/
static double _Complex dd[4];
static unsigned int del;
static int t1, t2;


static int index_p4[5], cases_int;
static double rff_p = zero, integ[5], integ4[5], integ04[5], integ14[5];
//static int N_temp, tt1, tt2;

static double a4, b4;
static double a_spin, a2, lambda, lam2, q, muobs, sinobs, scal;


static int mobseqmtp;
static double mu_tp1, mu_tp12, mu_tp2;
static double b0, b1, b2, b3, g2, g3;

static double tobs, tp2;
static double Wmum, Wmup, tminus, tplus;


static double PI0, PI01;

static double half_period_wp;
static double period_wp;

static double h;

static double pp, p1, p2, p_mt1_mt2;
static double PI1_p, PI2_p;
static double PI1_phi, PI2_phi;
static double PI1_time, PI2_time;
/*****************************************************************************/
static void Get_mu_t1_t2_New( double p, double f12342, int mobseqmtp, 
	double p_mt1_mt2, double muobs, double mu_tp1, 
	double mu_tp2, int *t1, int *t2 );
/*****************************************************************************/

/*
!*
!*
!*     PURPOSE:  Computes \mu part of integrals in coordinates \phi, t and affine parameter \sigma,
!*               expressed by equation (71) and (72) in Yang & Wang (2012).    
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f12342---------p_\theta, which is the \theta component of four momentum of a photon 
!*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
!*               f12343---------p_\phi, which is the \phi component of four momentum of a photon 
!*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  phyt-----------value of integral \phi_\theta expressed by equation (72) in
!*                              Yang & Wang (2012).  
!*               timet----------value of integral \t_\theta expressed by equation (71) in
!*                              Yang & Wang (2012). And \sigma_\theta=time_\theta.
!*               mucos----------value of function \mu(p).
!*               t1,t2----------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively.
!*               sign_pth-------sign of \theta component of 4-momentum of the photon.         
!*     ROUTINES CALLED: mutp, root3, weierstrass_int_J3 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: 
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-15.
!*
!*/
int Integration_Theta_part( ptcl *pt, double p, double *phit, double *timet, 
		double *mucos, double *sign_pth, int *tm1, int *tm2 )
{
	double pp_phi, pp_t;
	double p1_phi, p1_t;
	double p2_phi, p2_t;

	double tmu, u;
	static unsigned int count_num = 1;

	*tm1 = 0;
	*tm2 = 0;
	a_spin = pt->a_spin;
	a2 = pt->a2;
	lambda = pt->lambda;
	lam2 = pt->lam2;
	q = pt->q;
	muobs = pt->muobs;
	sinobs = pt->sinobs;
	scal = pt->scal;

	if ( pt->f1234[3] == zero && pt->f1234[2] == zero && fabs( muobs ) == one) {
		*sign_pth = zero;
		*mucos = sign( one, muobs );
		*timet = zero;             // this is because that mu==1 for ever
		*phit = zero;              // this is because that mu==1 for ever,this
		return 0;                    // because that Theta_mu=-a^2(1-mu^2), 
					   // so,mu must =+1 or -1 for ever.
	}
	if( muobs == zero && ( fabs( lambda ) < fabs( a_spin ) ) && q == zero ) {
		*sign_pth = zero;
		*timet = zero;
		*phit = zero;
		*mucos = zero;
		return 0;
	}

	mutp( pt );
	mu_tp1 = pt->mu_tp1;
	mu_tp12 = sq(mu_tp1);
	mu_tp2 = pt->mu_tp2;
	mobseqmtp = pt->mobseqmtp;

	if ( mu_tp1 == zero ) {
		// photons are confined in the equatorial plane, 
		// so the integrations about \theta are valished.
                *sign_pth = zero;
                *timet = zero;
                *phit = zero;
                *mucos = zero;
                count_num = count_num + 1;
		*tm1 = 0;
		*tm2 = 0;
                return 0;
	}

	if ( a_spin == zero ) {
		*timet = zero;
		//phyt_schwarz(p, pt->f1234[3], pt->f1234[2], lambda, q,
		//	sinobs, muobs, scal, phyt, mucos, t1, t2, sign_pth );
		count_num = count_num + 1;
		return 0;
	}
	a4 = zero;
	b4 = one;


	// equations (26)-(29) in Yang & Wang (2013). 
	b0 = - four * a2 * pt->mutp3 + two * pt->mu_tp1 * (a2 - lam2 - q);
	b1 = - two * a2 * pt->mutp2 + one / three * (a2 - lam2 - q);
	b2 = - four / three * a2 * pt->mu_tp1;
	b3 = - a2;
		
	/* equation (31) in Yang & Wang (2012). */
	g2 = three / four * ( sq(b1) - b0*b2 );
	g3 = one / 16.0 * ( three * b0 * b1 * b2 - two * sq3(b1) - sq(b0) * b3 );

	root3( zero, -g2/four, -g3/four, dd, &del );
 
	// equation (30) in Yang & Wang (2012). 
	if ( muobs != mu_tp1)
		tobs = b0 / four / ( muobs - mu_tp1 ) + b1 / four;
	else
		tobs = infinity;


	if ( mu_tp1 - one != zero) {
		// equation (64) in Yang & Wang (2012). 
		Wmum = b0 / ( eight * sq( - one - mu_tp1 ) );
		Wmup = b0 / ( eight * sq( one - mu_tp1 ) );
		tminus = b0 / four / ( - one - mu_tp1 ) + b1 / four;
		tplus = b0 / four / ( one - mu_tp1 ) + b1 / four;
	}
	index_p4[1] = 0;
	cases_int = 1;
	weierstrass_int_J3( tobs, infinity, dd, del, a4, b4, 
				index_p4, rff_p, integ04, cases_int );
	PI0 = integ04[1];
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	tp2 = b0 / four / ( mu_tp2 - mu_tp1 ) + b1 / four;
	weierstrass_int_J3( tp2, infinity, dd, del, a4, b4,
				index_p4, rff_p, integ14, cases_int );
	half_period_wp = integ14[1]; 
	period_wp = two * half_period_wp;
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	// equation (34) in Yang & Wang (2013).    
	if ( pt->f1234[2] < zero )
		PI01 = - PI0;
	else
		PI01 = PI0;

	tmu = weierstrassP( p + PI01, g2, g3, dd, del );
	// equation (32) in Yang & Wang (2013). 
	*mucos = mu_tp1 + b0 / ( four * tmu - b1 );
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	! Now we get the value of parameter sign_pth
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	u = p+PI01;
	if ( u <= zero )
		*sign_pth = - one;
	else {
		u = fmod(u, period_wp);
		if ( u <= half_period_wp )
			*sign_pth = one;
		else
			*sign_pth = - one;
	}

	return 0;

	/***********************************************************/
	h = - b1 / four;
	/* to get number of turning points of t1 and t2. */
	weierstrass_int_J3( tobs, tmu, dd, del, a4, b4,
			index_p4, rff_p, integ4, cases_int );
	pp = integ4[1];

	//equation (51) in Yang & Wang (2012).        
	p_mt1_mt2 = half_period_wp;
	PI1_p = PI0;
	PI2_p = p_mt1_mt2 - PI0;
	p1 = PI0 - pp;
	p2 = p_mt1_mt2 - p1;
	PI1_phi = zero;
	PI2_phi = zero;
	PI1_time = zero;
	PI2_time = zero;

	index_p4[1] = -1;
	index_p4[2] = -2;
	index_p4[3] = 0;
	index_p4[4] = -4;

	Get_mu_t1_t2_New( p, pt->f1234[2], mobseqmtp, p_mt1_mt2, 
		muobs, mu_tp1, mu_tp2, &t1, &t2 );


	/*************integration for pp part*********************************/
	if ( lambda != zero ) {
		cases_int = 2;
		weierstrass_int_J3( tobs, tmu, dd, del, -tplus, 
				b4, index_p4, fabs(pp), integ4, cases_int );
		weierstrass_int_J3( tobs, tmu, dd, del, -tminus,
				b4, index_p4, fabs(pp), integ14, cases_int );
		// equation (72) in Yang & Wang (2012). 
		pp_phi = lambda * ( pp / ( one - mu_tp12 ) + 
				integ4[2] * Wmup - integ14[2] * Wmum );
	} else 
		pp_phi = zero;

	cases_int = 4;
	weierstrass_int_J3( tobs, tmu, dd, del, h, b4,
			index_p4, fabs(pp), integ, cases_int );

	// equation (71) in Yang & Wang (2012). 
	pp_t = a2 * ( pp * mu_tp12 + integ[2] * mu_tp1 * 
			b0 / two + integ[4] * sq(b0) / sixteen );

	/*************integration for p1 part*********************************/
	if ( t1 == 0 ) {
		p1_phi = zero;
		p1_t = zero;
	} else {
		if ( lambda != zero ) {
			if ( PI1_phi == zero ) {
				cases_int = 2;
				weierstrass_int_J3( tobs, infinity, dd, del,
					-tplus, b4, index_p4, PI0, integ4, cases_int );
				weierstrass_int_J3( tobs,infinity,dd,del,
					-tminus, b4, index_p4, PI0, integ14, cases_int );
				// equation (72) in Yang & Wang (2012). 
				PI1_phi = lambda * ( PI0 / ( one - mu_tp12 ) +
					integ4[2] * Wmup - integ14[2] * Wmum );
			}
			// equation (51) in Yang & Wang (2012). 
			p1_phi = PI1_phi - pp_phi;
		} else 
			p1_phi = zero;


		if ( PI1_time == zero) {
			cases_int = 4;
			weierstrass_int_J3( tobs, infinity, dd, del,
				h, b4, index_p4, PI0, integ, cases_int );
			// equation (62) in Yang & Wang (2012). 
			PI1_time = a2 * ( PI0 * mu_tp12 + integ[2] *
				mu_tp1 * b0 / two + integ[4] * b0 * b0 / sixteen );
		}
		// equation (51) in Yang & Wang (2012). 
		p1_t = PI1_time - pp_t;
	}
	/*************integration for p2 part*********************************/
	if ( t2 == 0) {
                p2_phi = zero;
                p2_t = zero;
	} else {
		if ( lambda != zero ) {
			if ( PI2_phi == zero) {
				cases_int = 2;        
				weierstrass_int_J3( tp2, tobs, dd, del,
					-tplus, b4, index_p4, PI2_p, integ4, cases_int );
				weierstrass_int_J3( tp2, tobs, dd, del,
					-tminus, b4, index_p4, PI2_p, integ14, cases_int );
				// equation (72) in Yang & Wang (2013). 
		                PI2_phi = lambda * ( PI2_p / ( one - mu_tp12 ) +
					integ4[2] * Wmup - integ14[2] * Wmum );
			}
			// equation (51) in Yang & Wang (2013). 
			p2_phi = PI2_phi + pp_phi;
		} else
			p2_phi = zero;


		if ( PI2_time == zero) {
			cases_int = 4;
			weierstrass_int_J3( tp2, tobs, dd, del, h, b4, 
				index_p4, PI2_p, integ, cases_int );
			// equation (71) in Yang & Wang (2012).  
			PI2_time = a2 * ( PI2_p * mu_tp12 + integ[2] *
				mu_tp1 * b0 / two + integ[4] * b0 * b0 / sixteen );
		}
		// equation (51) in Yang & Wang (2012).      
		p2_t = PI2_time + pp_t;
	}


	/*********************************************************************/
	// equation (52) in Yang & Wang (2012). 
	if ( mobseqmtp ) {
		if ( muobs == mu_tp1 ) {
			*phit = -pp_phi + two * ( t1 * p1_phi + t2 * p2_phi );
			*timet = -pp_t + two * ( t1 * p1_t + t2 * p2_t );
		} else {
			*phit = pp_phi + two * ( t1 * p1_phi + t2 * p2_phi );
			*timet = pp_t + two * ( t1 * p1_t + t2 * p2_t );
		}
	} else {
		if ( pt->f1234[2] < zero ) {
			*phit = pp_phi + two * ( t1 * p1_phi + t2 * p2_phi );
			*timet = pp_t + two * ( t1 * p1_t + t2 * p2_t );
		}
		if ( pt->f1234[2] > zero ) {
			*phit = -pp_phi + two * ( t1 * p1_phi + t2 * p2_phi );
			*timet = -pp_t + two * ( t1 * p1_t + t2 * p2_t );
		}
	}
}

/*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-15.
!*
!*/
int Integration_Theta_part_Settings( ptcl *pt )
{
	a_spin = pt->a_spin;
	a2 = pt->a2;
	lambda = pt->lambda;
	lam2 = pt->lam2;
	q = pt->q;
	muobs = pt->muobs;
	sinobs = pt->sinobs;
	scal = pt->scal;

	mutp( pt );
	mu_tp1 = pt->mu_tp1;
	mu_tp12 = sq(mu_tp1);
	mu_tp2 = pt->mu_tp2;
	mobseqmtp = pt->mobseqmtp;

	if ( mu_tp1 == zero ) {
		// photons are confined in the equatorial plane, 
		// so the integrations about \theta are valished.
                /**sign_pth = zero;
                *timet = zero;
                *phit = zero;
                *mucos = zero;
                count_num = count_num + 1;
		*tm1 = 0;
		*tm2 = 0;*/
                return 0;
	}

	if ( a_spin == zero ) {
		//*timet = zero;
		//phyt_schwarz(p, pt->f1234[3], pt->f1234[2], lambda, q,
		//	sinobs, muobs, scal, phyt, mucos, t1, t2, sign_pth );
		//count_num = count_num + 1;
		return 0;
	}
	a4 = zero;
	b4 = one;


	// equations (26)-(29) in Yang & Wang (2013). 
	b0 = - four * a2 * pt->mutp3 + two * pt->mu_tp1 * (a2 - lam2 - q);
	b1 = - two * a2 * pt->mutp2 + one / three * (a2 - lam2 - q);
	b2 = - four / three * a2 * pt->mu_tp1;
	b3 = - a2;
		
	/* equation (31) in Yang & Wang (2012). */
	g2 = three / four * ( sq(b1) - b0*b2 );
	g3 = one / 16.0 * ( three * b0 * b1 * b2 - two * sq3(b1) - sq(b0) * b3 );

	root3( zero, -g2/four, -g3/four, dd, &del );
 
	// equation (30) in Yang & Wang (2012). 
	if ( muobs != mu_tp1)
		tobs = b0 / four / ( muobs - mu_tp1 ) + b1 / four;
	else
		tobs = infinity;


	if ( mu_tp1 != one ) {
		// equation (64) in Yang & Wang (2012). 
		Wmum = b0 / ( eight * sq( - one - mu_tp1 ) );
		Wmup = b0 / ( eight * sq( one - mu_tp1 ) );
		tminus = b0 / four / ( - one - mu_tp1 ) + b1 / four;
		tplus = b0 / four / ( one - mu_tp1 ) + b1 / four;
	}
	index_p4[1] = 0;
	cases_int = 1;
	weierstrass_int_J3( tobs, infinity, dd, del, a4, b4, 
				index_p4, rff_p, integ04, cases_int );
	PI0 = integ04[1];
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	tp2 = b0 / four / ( mu_tp2 - mu_tp1 ) + b1 / four;
	weierstrass_int_J3( tp2, infinity, dd, del, a4, b4,
				index_p4, rff_p, integ14, cases_int );
	half_period_wp = integ14[1]; 
	period_wp = two * half_period_wp;
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	// equation (34) in Yang & Wang (2013).    
	if ( pt->f1234[2] < zero )
		PI01 = - PI0;
	else
		PI01 = PI0;
 

	/***********************************************************/
	h = - b1 / four; 

	//equation (51) in Yang & Wang (2012).        
	p_mt1_mt2 = half_period_wp;
	PI1_p = PI0;
	PI2_p = p_mt1_mt2 - PI0;

	PI1_phi = zero;
	PI2_phi = zero;
	PI1_time = zero;
	PI2_time = zero;

	/*index_p4[1] = -1;
	index_p4[2] = -2;
	index_p4[3] = 0;
	index_p4[4] = -4;*/

	//Get_mu_t1_t2_New( p, pt->f1234[2], mobseqmtp, p_mt1_mt2, 
	//	muobs, mu_tp1, mu_tp2, &t1, &t2 );
	return 0;
}

 
 





/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!* To determine the number of time_0 N_t1, N_t2 that the particle meets the
!* two turn points mu_tp1, mu_tp2 
!*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
static void Get_mu_t1_t2_New( double p, double f12342, int mobseqmtp, double p_mt1_mt2, 
		double muobs, double mu_tp1, double mu_tp2, int *t1, int *t2 )
{
	double p1_temp, p2_temp;
	int N_temp;
	if ( mobseqmtp ) {
		p1_temp = zero;
		p2_temp = p_mt1_mt2;
		*t1 = 0;
		*t2 = 0;
		N_temp = 0;
		if ( muobs == mu_tp1 ) {
			do {
				if ( p1_temp <= p && p <= p2_temp )
					break;
				else {
					N_temp = N_temp + 1;
					p1_temp = p2_temp;
					p2_temp = p2_temp + p_mt1_mt2;
					*t1 = *t2;
					*t2 = N_temp - *t1;
				}
			} while (true);
		} else if ( muobs == mu_tp2 ) {
			do {
				if ( p1_temp <= p && p <= p2_temp )
					break;
				else {
					N_temp = N_temp + 1;
					p1_temp = p2_temp;
					p2_temp = p2_temp + p_mt1_mt2;
					*t2 = *t1;
					*t1 = N_temp - *t2;
				}
			} while (true);
		}
	} else {
		p1_temp = zero;
		*t1 = 0;
		*t2 = 0;
		N_temp = 0;
		if ( f12342 > zero ) {
			p2_temp = PI2_p;
			do {
				if ( p1_temp <= p && p <= p2_temp )
					break;
				else {
					N_temp = N_temp + 1;
					p1_temp = p2_temp;
					p2_temp = p2_temp + p_mt1_mt2;
					*t1 = *t2;
					*t2 = N_temp - *t1;
				}
			} while (true);
		}
		if ( f12342 < zero ) {
			p2_temp = PI1_p;
			do {
				if ( p1_temp <= p && p <= p2_temp )
					break;
				else {
                                      N_temp = N_temp + 1;
                                      p1_temp = p2_temp;
                                      p2_temp = p2_temp + p_mt1_mt2;
                                      *t2 = *t1;
                                      *t1 = N_temp - *t2;
				}
			} while (true);
		}
	}
}




/*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-15.
!*
!*/
int Get_Integrations_of_Theta_part( ptcl *pt, double p, double *phit, double *timet, 
		double *mucos, double *sign_pth, int *tm1, int *tm2 )
{
	double pp_phi, pp_t;
	double p1_phi, p1_t;
	double p2_phi, p2_t;
	double tmu, u;

	*tm1 = 0;
	*tm2 = 0;
	if ( pt->f1234[3] == zero && pt->f1234[2] == zero && fabs( muobs ) == one) {
		*sign_pth = zero;
		*mucos = sign( one, muobs );
		*timet = zero;             // this is because that mu==1 for ever
		*phit = zero;              // this is because that mu==1 for ever,this
		return 0;                    // because that Theta_mu=-a^2(1-mu^2), 
					   // so,mu must =+1 or -1 for ever.
	}
	if( muobs == zero && ( fabs( lambda ) < fabs( a_spin ) ) && q == zero ) {
		*sign_pth = zero;
		*timet = zero;
		*phit = zero;
		*mucos = zero;
		return 0;
	}
  
	if ( fabs( mu_tp1 - zero ) < 1.e-11 ) {
		// photons are confined in the equatorial plane, 
		// so the integrations about \theta are valished.
                *sign_pth = zero;
                *timet = zero;
                *phit = zero;
                *mucos = zero;
		*tm1 = 0;
		*tm2 = 0;
                return 0;
	}

	if ( a_spin == zero ) {
		*timet = zero;
		//phyt_schwarz(p, pt->f1234[3], pt->f1234[2], lambda, q,
		//	sinobs, muobs, scal, phyt, mucos, t1, t2, sign_pth );
		return 0;
	}
 
	tmu = weierstrassP( p + PI01, g2, g3, dd, del );
	// equation (32) in Yang & Wang (2013). 
	*mucos = mu_tp1 + b0 / ( four * tmu - b1 );
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	! Now we get the value of parameter sign_pth
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	u = p+PI01;
	if ( u <= zero )
		*sign_pth = - one;
	else {
		u = fmod(u, period_wp);
		if ( u <= half_period_wp )
			*sign_pth = one;
		else
			*sign_pth = - one;
	}

	/***********************************************************/
	index_p4[1] = 0;
	cases_int = 1;
	weierstrass_int_J3( tobs, tmu, dd, del, a4, b4,
			index_p4, rff_p, integ4, cases_int );
	pp = integ4[1];
	p1 = PI0 - pp;
	p2 = p_mt1_mt2 - p1;

	Get_mu_t1_t2_New( p, pt->f1234[2], mobseqmtp, p_mt1_mt2, 
		muobs, mu_tp1, mu_tp2, &t1, &t2 );

	*tm1 = t1;
	*tm2 = t2;

	index_p4[1] = -1;
	index_p4[2] = -2;
	index_p4[3] = 0;
	index_p4[4] = -4;
 
	/*************integration for pp part*********************************/
	if ( lambda != zero ) {
		cases_int = 2;
		weierstrass_int_J3( tobs, tmu, dd, del, -tplus, 
				b4, index_p4, fabs(pp), integ4, cases_int );
		weierstrass_int_J3( tobs, tmu, dd, del, -tminus,
				b4, index_p4, fabs(pp), integ14, cases_int );
		// equation (72) in Yang & Wang (2012). 
		pp_phi = lambda * ( pp / ( one - mu_tp12 ) + 
				integ4[2] * Wmup - integ14[2] * Wmum );
	} else 
		pp_phi = zero;

	cases_int = 4;
	weierstrass_int_J3( tobs, tmu, dd, del, h, b4,
			index_p4, fabs(pp), integ, cases_int );

	// equation (71) in Yang & Wang (2012). 
	pp_t = a2 * ( pp * mu_tp12 + integ[2] * mu_tp1 * 
			b0 / two + integ[4] * sq(b0) / sixteen );
	/*printf( " pp_theta  = %f  %f  %f %f  %f \n ", mu_tp12, integ[2], integ[4], mu_tp1, b0 );
	printf( " pp_theta  = %f  %f  %f %f  %f \n ", tobs, tmu, pp, h, b4 );
	printf("%20.16f + %20.16fi \n", creal(dd[1]), cimag(dd[1]));
	printf("%20.16f + %20.16fi \n", creal(dd[2]), cimag(dd[2]));
	printf("%20.16f + %20.16fi \n", creal(dd[3]), cimag(dd[3]));
	printf( " muobs = %40.35f  %40.35f  %40.35f \n ", muobs, one, muobs - one );*/

	/*************integration for p1 part*********************************/
	if ( t1 == 0 ) {
		p1_phi = zero;
		p1_t = zero;
	} else {
		if ( lambda != zero ) {
			if ( PI1_phi == zero ) {
				cases_int = 2;
				weierstrass_int_J3( tobs, infinity, dd, del,
					-tplus, b4, index_p4, PI0, integ4, cases_int );
				weierstrass_int_J3( tobs,infinity,dd,del,
					-tminus, b4, index_p4, PI0, integ14, cases_int );
				// equation (72) in Yang & Wang (2012). 
				PI1_phi = lambda * ( PI0 / ( one - mu_tp12 ) +
					integ4[2] * Wmup - integ14[2] * Wmum );
			}
			// equation (51) in Yang & Wang (2012). 
			p1_phi = PI1_phi - pp_phi;
		} else 
			p1_phi = zero;


		if ( PI1_time == zero) {
			cases_int = 4;
			weierstrass_int_J3( tobs, infinity, dd, del,
				h, b4, index_p4, PI0, integ, cases_int );
			// equation (62) in Yang & Wang (2012). 
			PI1_time = a2 * ( PI0 * mu_tp12 + integ[2] *
				mu_tp1 * b0 / two + integ[4] * b0 * b0 / sixteen );
		}
		// equation (51) in Yang & Wang (2012). 
		p1_t = PI1_time - pp_t;
	}
	/*************integration for p2 part*********************************/
	if ( t2 == 0) {
                p2_phi = zero;
                p2_t = zero;
	} else {
		if ( lambda != zero ) {
			if ( PI2_phi == zero) {
				cases_int = 2;        
				weierstrass_int_J3( tp2, tobs, dd, del,
					-tplus, b4, index_p4, PI2_p, integ4, cases_int );
				weierstrass_int_J3( tp2, tobs, dd, del,
					-tminus, b4, index_p4, PI2_p, integ14, cases_int );
				// equation (72) in Yang & Wang (2013). 
		                PI2_phi = lambda * ( PI2_p / ( one - mu_tp12 ) +
					integ4[2] * Wmup - integ14[2] * Wmum );
			}
			// equation (51) in Yang & Wang (2013). 
			p2_phi = PI2_phi + pp_phi;
		} else
			p2_phi = zero;


		if ( PI2_time == zero) {
			cases_int = 4;
			weierstrass_int_J3( tp2, tobs, dd, del, h, b4, 
				index_p4, PI2_p, integ, cases_int );
			// equation (71) in Yang & Wang (2012).  
			PI2_time = a2 * ( PI2_p * mu_tp12 + integ[2] *
				mu_tp1 * b0 / two + integ[4] * b0 * b0 / sixteen );
		}
		// equation (51) in Yang & Wang (2012).      
		p2_t = PI2_time + pp_t;
	}


	/*********************************************************************/
	// equation (52) in Yang & Wang (2012). 
	if ( mobseqmtp ) {
		if ( muobs == mu_tp1 ) {
			*phit = -pp_phi + two * ( t1 * p1_phi + t2 * p2_phi );
			*timet = -pp_t + two * ( t1 * p1_t + t2 * p2_t );
			//printf( "time_theta1 = %f  %f %f \n ", pp_t, p1_t, p2_t );
		} else {
			*phit = pp_phi + two * ( t1 * p1_phi + t2 * p2_phi );
			*timet = pp_t + two * ( t1 * p1_t + t2 * p2_t );
			//printf( "time_theta2 = %f  %f %f \n ", pp_t, p1_t, p2_t );
		}
	} else {
		if ( pt->f1234[2] < zero ) {
			*phit = pp_phi + two * ( t1 * p1_phi + t2 * p2_phi );
			*timet = pp_t + two * ( t1 * p1_t + t2 * p2_t );
			//printf( "time_theta3 = %f  %f %f \n ", pp_t, p1_t, p2_t );
		}
		if ( pt->f1234[2] > zero ) {
			*phit = -pp_phi + two * ( t1 * p1_phi + t2 * p2_phi );
			*timet = -pp_t + two * ( t1 * p1_t + t2 * p2_t );
			//printf( "time_theta4 = %f  %f %f \n ", pp_t, p1_t, p2_t );
		}
	}
	return 0;
}





static double a_spin_1;
static double scal_1;
static double lambda_1;
static double q_1;
static double f12343_1;
static double f12342_1;
static double muobs_1;
static double sinobs_1;


/*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-17.
!*
!*/
void Get_Integrals_For_Theta_Part( ptcl *pt, double p, double *phit, double *timet, 
		double *mucos, double *sign_pth, int *tm1, int *tm2 )
{
	static unsigned int count_num = 1;

mus_restarting:
	//printf("here33ss\n");
	if ( count_num == 1 ) {

		Integration_Theta_part_Settings( pt );
		count_num += 1;

		a_spin_1 = a_spin;
		scal_1 = scal;
		lambda_1 = lambda;
		q_1 = q;
		f12343_1 = pt->f1234[3];
		f12342_1 = pt->f1234[2];
		muobs_1 = muobs;
		sinobs_1 = sinobs;
 
		Get_Integrations_of_Theta_part( pt, p, phit, timet, 
			mucos, sign_pth, tm1, tm2 ); 
    	} else {
		if ( f12343_1 == pt->f1234[3]
		&& f12342_1 == pt->f1234[2]
		&& lambda_1 == lambda
		&& q_1 == q
		&& sinobs_1 == sinobs
		&& muobs_1 == muobs
		&& a_spin_1 == a_spin
		&& scal_1 == scal ) {

			Get_Integrations_of_Theta_part( pt, p, phit, timet, 
				mucos, sign_pth, tm1, tm2 );
		} else {
			count_num = 1;
			goto mus_restarting;
		}

	}
}

 







