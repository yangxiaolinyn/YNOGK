/*
 * common_src/sphericalmotion.c
 *
!********************************************************************
!* This module aim on the calculation of the geodesic trajectory 
!* of a photon doing spherical motion.
!********************************************************************
 *
 * Author       Yang Xiao-lin
 *
 * Date         2023-01-03.
 *
 * City         Kun-ming Yunnan Provice, China
 *
 * 2022-12-11   Stating the writting of the module.
 */
 

#include "ynogkBL.h"


static double c_m, c_add, a_m, a_add;
static double rff_p;
static double PI0, PI01, p_mt1_mt2, PI2_p;
static double PI1_phi, PI2_phi, PI1_sig, PI2_sig, PI1_time, PI2_time;
static double tmu, h;
static double integ5[6], integ05[6], integ15[6];
static double integ[5];
static int index_p5[6], cases_int;

static double a4, b4;
static double b0, b1, b2, b3, g2, g3;
static double _Complex dd[4];
static unsigned int del;
static double tobs, tp1, tp2;

static double come, Delta, c_phi, c_time;
static double pp_phi;

static double sinmax, mumax;
static double Omega, c_tau;




static double radiusofsphericalmotion( double a, double c );
static double mveone_int( double x, double y, double z0 );
static void SPINZERO( ptcl *pt, double p, double kp, double kt, double theta_max, 
		double *phyt, double *timet, double *sigmat,
		double *mucos, int *t1, int *t2 );



static double kp_1;
static double kt_1;
static double lambda_1;
static double q_1;
static double muobs_1;
static double sinobs_1;
static double a_spin_1;
static double robs_1;
static double thetamax_1;


/*
!*
!*     PURPOSE:   This routine computs the four B_L coordiants r,\mu,\phi,t and affine 
!*                parameter \sigma of the spherical motion. 
!*     INPUTS:   p--------------the independent variable.
!*               kt-------------p_\theta, the \theta component of four momentum of the photon.
!*               kp-------------p_\phi,  the \phi component of four momentum of the photon.
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{ini}), muobs=cos(\theta_{ini}), where 
!*                              \theta_{ini} is the initial theta angle of the photon.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or the initial position of photon. 
!*               theta_max-------the maximum or minimum value of the \theta coordinate of the geodesic.
!*                              which also the \theta turning point of the spherical motion.       
!*     OUTPUTS:  phyt-----------the \phi coordinat of the photon.
!*               timet----------the t coordinats of the photon.
!*               sigmat---------the affine parameter \sigma.
!*               mucos----------the \theta coordinates of the photon, and mucos=cos(\theta).  
!*               t1,t2----------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively.      
!*     ROUTINES CALLED: metricg, circ_mb, weierstrass_int_J3, weierstrassP
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2013)  
!*     DATE WRITTEN:  5 Jan 2013
!*     REVISIONS: 
!*
!*
!*     C VERSION:  Yang Xiao-lin    2023-01-03.
!*
!*/
void SPHERICALMOTION_BL( ptcl *pt, double p, double theta_max, double kp, 
	double kt, double *phyt, double *timet, double *sigmat, 
	double *mucos, int *tt1, int *tt2 )
{
	double p1, p2, pp, p_mu = zero;
	double pp_sig, pp_t, p1_t, p1_sig, p1_phi;
	double p2_t, p2_sig, p2_phi;
	int t1, t2;
	 
	/*  kp_1,kt_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,robs_1,a4,b4,mu_tp1,mu_tp2,reals,&
            mobseqmtp,b0,b1,b2,b3,g2,g3,dd,del,PI0,c_m,c_add,a_m,a_add,tp2,tobs,h,p_mt1_mt2,&
            PI1_phi,PI2_phi,PI1_time,PI2_time,PI2_p,come,tp1,bigT,Delta,c_phi,c_time,&
            PI01,sinmax,mumax,theta_max_1,Omega,c_tau*/

	int static count_num = 1;
 
	t1 = 0;
	t2 = 0;
	if ( theta_max != 90.0 && theta_max != 180.0 ) {
		sinmax = sin( theta_max * dtor );
		mumax = cos( theta_max * dtor );
	} else {
		if ( theta_max == 90.0 ) {
			sinmax = one;
			mumax = zero;
		}
		if ( theta_max == 180.0 ) {
			sinmax = zero;
			mumax = one;
		}
	}
	pt->mu_tp1 = fabs( mumax );
	pt->mu_tp2 = - pt->mu_tp2;
	pt->mu_tp12 = pow( pt->mu_tp1, 2 );
	pt->mu_tp13 = pow( pt->mu_tp1, 3 );

	if ( pt->mu_tp1 == zero ) {
		/* photons are confined in the equatorial plane, 
		so the integrations about \theta are valished. */
		double Omega, c_tau;
		double somiga, expnu, exppsi, expmu1, expmu2;

		Omega = one / ( pow( pt->robs, (three/two) ) + pt->a_spin );
		*mucos = zero;
		*phyt = p;
		*timet = (*phyt) / Omega;
		metricgij( pt->robs, pt->muobs, pt->sinobs, pt->a_spin, 
			&somiga, &expnu, &exppsi, &expmu1, &expmu2 );
		c_tau = sqrt( -expnu * expnu + somiga*somiga*exppsi*exppsi
				- two*exppsi*somiga*Omega
				+ exppsi*exppsi*Omega*Omega );
		*sigmat = c_tau * (*timet);
		count_num = count_num + 1;
		return;
	}


	if ( pt->a_spin == zero ) {
		*timet = zero;
		SPINZERO( pt, p, kp, kt, theta_max, phyt, timet, sigmat, mucos, tt1, tt2);
		count_num ++;
		return;
	}

	a4 = zero;
	b4 = one;
	// equations (26)-(29) of Yang and Wang (2013).
	come = -one;
	b0 = four*pt->a2 * pt->mu_tp13 * come - two * pt->mu_tp1
			* ( pt->a2 * come + pt->lam2 + pt->q );
	b1 = two * pt->a2 * sq( pt->mu_tp1 ) * come
		- ( pt->a2 * come + pt->lam2 + pt->q ) / three;
	b2 = four / three * pt->a2 * pt->mu_tp1 * come;
	b3 = pt->a2 * come;

	// equations (31) of Yang and Wang (2013).
	g2 = three / four * ( b1 * b1 - b0*b2 );
	g3 = ( three * b0 * b1 * b2 - two * pow(b1, 3) - b0 * b0 * b3 ) / 16.0;  
	root3( zero, -g2/four, -g3/four, dd, &del );

	// equations (30) of Yang and Wang (2013).
	if ( pt->muobs != pt->mu_tp1 )
		tobs = b0 / four / ( pt->muobs - pt->mu_tp1 ) + b1 / four;
	else
		tobs = infinity;

	tp1 = infinity;
	tp2 = b0 / four / ( pt->mu_tp2 - pt->mu_tp1 ) + b1 / four;
	//equations (72)-(73) of Yang and Wang (2013).
	if ( pt->mu_tp1 - one != zero ) {
		c_m = b0 / ( four * sq( - one - pt->mu_tp1 ) );
		c_add = b0 / ( four * sq( one - pt->mu_tp1 ) );
		a_m = b0 / four / ( - one - pt->mu_tp1 ) + b1 / four;
		a_add = b0 / four / ( one - pt->mu_tp1 ) + b1 / four;
	}

	index_p5[1] = 0;
	cases_int = 1;
	// equations (53) of Yang and Wang (2013).
	weierstrass_int_J3( tobs, tp1, dd, del, a4, b4, 
		index_p5, rff_p, integ05, cases_int );
	PI0=integ05[1];
	if ( kt < zero )
              PI01 = - PI0;
	else
              PI01 = PI0;

	tmu = weierstrassP( p + PI01, g2, g3, dd, del );
	// equations (32) of Yang and Wang (2013).
	*mucos = pt->mu_tp1 + b0 / ( four * tmu - b1 );
	h = - b1 / four;
	// to get number of turn points of t1 and t2.
	/******************************************************************/
	// mu = mu_tp + b0 / ( four * tmu - b1 );
	weierstrass_int_J3(tp2,tp1,dd,del,a4,b4,index_p5,rff_p,integ15,cases_int);
	weierstrass_int_J3(tobs,tmu,dd,del,a4,b4,index_p5,rff_p,integ5,cases_int);
	// equations (53) of Yang and Wang (2013).	
	p_mt1_mt2 = integ15[1];
	PI2_p = p_mt1_mt2 - PI0;
	pp = integ5[1];
	p1 = PI0 - pp;
	// equations (53) of Yang and Wang (2013).
	p2 = p_mt1_mt2 - p1;
	PI1_phi = zero;
	PI2_phi = zero;
	PI1_sig = zero;
	PI2_sig = zero;
	PI1_time = zero;
	PI2_time = zero;
	// equations (54) of Yang and Wang (2013).
	for ( long int j = 0; j < 1.e5; j++ ) {
		for ( long int i = j; i <= j + 1; i++ ) {
			if ( pt->mobseqmtp ) {
				if( pt->muobs == pt->mu_tp1 ) {
					t1 = j;
					t2 = i;
					p_mu = -pp + two*( t1*p1 + t2*p2 );
				} else {
					t1 = i;
					t2 = j;
					p_mu = pp + two*( t1*p1 + t2*p2 );
				}
			} else {
				if (kt < zero) {
					t1 = i;
					t2 = j;
					p_mu = pp + two*( t1*p1 + t2*p2 );
				}
				if (kt > zero) {
					t1 = j;
					t2 = i;
					p_mu = -pp + two*( t1*p1 + t2*p2 );
				}
			}  
			if ( fabs( p - p_mu ) < 1.e-3 )
				goto gooutfor;
		}
	}

gooutfor:
	// equations (71)-(73) of Yang and Wang (2013).
	Delta = pt->r2 + pt->a2 - two * pt->robs;
	c_phi = pt->a_spin * ( pt->robs * (two) - ( pt->lambda * pt->a_spin ) ) / Delta;
	c_time = ( (two) * pow( pt->robs, three ) + ( two * pt->a_spin 
		* ( pt->a_spin - pt->lambda ) ) * pt->robs ) / Delta;
	index_p5[1] = -1;
	index_p5[2] = -2;
	index_p5[3] = 0;
	index_p5[4] = -4;
	index_p5[5] = 0;
	// *****pp part***************************************
	if ( pt->lambda != zero ) {
		cases_int = 2;
		weierstrass_int_J3(tobs,tmu,dd,del,-a_add,b4,index_p5,abs(pp),integ5,cases_int);
		weierstrass_int_J3(tobs,tmu,dd,del,-a_m,b4,index_p5,abs(pp),integ15,cases_int);
		// equations (21) (72) of Yang and Wang (2013).
		pp_phi = ( pp / ( one- pt->mu_tp12 ) + (integ5[2] * c_add
			- integ15[2] * c_m ) / two ) * pt->lambda + c_phi * pp;
	} else 
		pp_phi = c_phi * pp;

	cases_int = 4;
	weierstrass_int_J3(tobs,tmu,dd,del,h,b4,index_p5,abs(pp),integ,cases_int);
	// equations (20) (71) of Yang and Wang (2013).
	pp_sig = ( pp * pt->mu_tp12 + integ[2] * pt->mu_tp1 * b0 / two + integ[4] * b0 * b0
			/ sixteen ) * pt->a2 + pt->r2 * pp;
	// equations (22) of Yang and Wang (2013).
	pp_t = pp_sig + c_time * pp;


	// *****p1 part***************************************
	if ( t1 == 0 ) {
		p1_phi = zero;
		p1_sig = zero;
		p1_t = zero;
	} else {
		if ( pt->lambda != zero ) {
			if ( PI1_phi == zero ) {
				cases_int = 2;
				weierstrass_int_J3(tobs,infinity,dd,del,-a_add,
						b4,index_p5,PI0,integ5,cases_int);
				weierstrass_int_J3(tobs,infinity,dd,del,-a_m,
						b4,index_p5,PI0,integ15,cases_int);
				 //equations (21) (72) of Yang and Wang (2013).
				PI1_phi = ( PI0 / ( one - pt->mu_tp12 )
					+ ( integ5[2] * c_add - integ15[2] * c_m ) / two )
					* pt->lambda + c_phi * PI0;
			} 
			p1_phi = PI1_phi - pp_phi;
		} else {
			if ( PI1_phi == zero ) PI1_phi = c_phi * PI0;
			p1_phi = PI1_phi - pp_phi;
		}
		if ( PI1_time == zero || PI1_sig == zero) {
			cases_int = 4;
			weierstrass_int_J3(tobs,infinity,dd,del,h,b4,index_p5,PI0,integ,cases_int);
			// equations (20) (71) of Yang and Wang (2013).
			PI1_sig = ( PI0 * pt->mu_tp12 + integ[2] * pt->mu_tp1 * b0 / two
				+ integ[4] * b0 * b0 / sixteen ) * pt->a2 + pt->r2 * PI0;
			// equations (22) of Yang and Wang (2013).
			PI1_time = PI1_sig + c_time * PI0;
		}
		p1_sig = PI1_sig - pp_sig;
		p1_t = PI1_time - pp_t;
	}

	// *****p2 part***************************************
	if ( t2 == 0 ) {
		p2_sig = zero;
		p2_phi = zero;
		p2_t = zero;
	} else {
		if ( pt->lambda != zero ) {
			if ( PI2_phi == zero ) {
				cases_int = 2;	
				weierstrass_int_J3(tp2,tobs,dd,del,-a_add,b4,
						index_p5,PI2_p,integ5,cases_int);
				weierstrass_int_J3(tp2,tobs,dd,del,-a_m,b4,
						index_p5,PI2_p,integ15,cases_int);
				// equations (21) (72) of Yang and Wang (2013).
				PI2_phi = ( PI2_p / ( one - pt->mu_tp12 )
					+ ( integ5[2] * c_add - integ15[2] * c_m )
					/ two ) * pt->lambda + c_phi * PI2_p;
			}
			p2_phi = PI2_phi + pp_phi;
		} else {
			if ( PI2_phi == zero ) PI2_phi = c_phi * PI2_p;
			p2_phi = PI2_phi + pp_phi;
		}
                
		if ( PI2_time == zero ) {
			cases_int = 4;
			weierstrass_int_J3(tp2,tobs,dd,del,h,b4,index_p5,PI2_p,integ,cases_int);
			// equations (20) (71) of Yang and Wang (2013).
			PI2_sig = ( PI2_p * pt->mu_tp12 + integ[2] * pt->mu_tp1 * b0
				/ two + integ[4] * b0 * b0 / sixteen )
				* pt->a2 + pt->r2 * PI2_p;
		// equations (22) of Yang and Wang (2013).
			PI2_time = PI2_sig + c_time * PI2_p;
		}
		p2_sig = PI2_sig + pp_sig;
		p2_t = PI2_time + pp_t;
	}



	if ( pt->mobseqmtp ) {
		if ( pt->muobs == pt->mu_tp1 ) {
			// equations (52) of Yang and Wang (2013).
			*phyt = -pp_phi + two * ( t1*p1_phi + t2*p2_phi );
			*timet = -pp_t + two * ( t1*p1_t + t2*p2_t );
			*sigmat = -pp_sig + two * ( t1*p1_sig + t2*p2_sig );
		} else {
			// equations (52) of Yang and Wang (2013).
			*phyt = pp_phi + two*( t1*p1_phi + t2*p2_phi );
			*timet = pp_t + two*( t1*p1_t + t2*p2_t );
			*sigmat = pp_sig + two*( t1*p1_sig + t2*p2_sig );
		}
	} else {
		if ( kt < zero ) {
			// equations (52) of Yang and Wang (2013).
			*phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi);	
			*timet=pp_t+two*(t1*p1_t+t2*p2_t);
			*sigmat=pp_sig+two*(t1*p1_sig+t2*p2_sig);
		}
		if ( kt > zero) {
			// equations (52) of Yang and Wang (2013).		
			*phyt = -pp_phi + two*( t1*p1_phi + t2*p2_phi );
			*timet = -pp_t + two*( t1*p1_t + t2*p2_t );
			*sigmat = -pp_sig + two*( t1*p1_sig + t2*p2_sig );
		}
	}
	if ( pt->mu_tp1 == one ) *phyt = *phyt + ( t1 + t2) * pi;
	// phyt = fmod(phyt, twopi);
	// if ( phyt < zero ) phyt = phyt + twopi;
	count_num++;
}




/*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2023-01-04.
!*
!*/
static void Spherical_motion_Settings( ptcl *pt, double theta_max, double kt )
{
	if ( theta_max != 90.0 && theta_max != 180.0 ) {
		sinmax = sin( theta_max * dtor );
		mumax = cos( theta_max * dtor );
	} else {
		if ( theta_max == 90.0 ) {
			sinmax = one;
			mumax = zero;
		}
		if ( theta_max == 180.0 ) {
			sinmax = zero;
			mumax = one;
		}
	}
	pt->mu_tp1 = fabs( mumax );
	pt->mu_tp2 = - pt->mu_tp2;
	pt->mu_tp12 = pow( pt->mu_tp1, 2 );
	pt->mu_tp13 = pow( pt->mu_tp1, 3 );
 

	a4 = zero;
	b4 = one;
	// equations (26)-(29) of Yang and Wang (2013).
	come = -one;
	b0 = four*pt->a2 * pt->mu_tp13 * come - two * pt->mu_tp1
			* ( pt->a2 * come + pt->lam2 + pt->q );
	b1 = two * pt->a2 * sq( pt->mu_tp1 ) * come
		- ( pt->a2 * come + pt->lam2 + pt->q ) / three;
	b2 = four / three * pt->a2 * pt->mu_tp1 * come;
	b3 = pt->a2 * come;

	// equations (31) of Yang and Wang (2013).
	g2 = three / four * ( b1 * b1 - b0*b2 );
	g3 = ( three * b0 * b1 * b2 - two * pow(b1, 3) - b0 * b0 * b3 ) / 16.0;  
	root3( zero, -g2/four, -g3/four, dd, &del );

	// equations (30) of Yang and Wang (2013).
	if ( pt->muobs != pt->mu_tp1 )
		tobs = b0 / four / ( pt->muobs - pt->mu_tp1 ) + b1 / four;
	else
		tobs = infinity;

	tp1 = infinity;
	tp2 = b0 / four / ( pt->mu_tp2 - pt->mu_tp1 ) + b1 / four;
	//equations (72)-(73) of Yang and Wang (2013).
	if ( pt->mu_tp1 - one != zero ) {
		c_m = b0 / ( four * sq( - one - pt->mu_tp1 ) );
		c_add = b0 / ( four * sq( one - pt->mu_tp1 ) );
		a_m = b0 / four / ( - one - pt->mu_tp1 ) + b1 / four;
		a_add = b0 / four / ( one - pt->mu_tp1 ) + b1 / four;
	}

	index_p5[1] = 0;
	cases_int = 1;
	// equations (53) of Yang and Wang (2013).
	weierstrass_int_J3( tobs, tp1, dd, del, a4, b4, 
		index_p5, rff_p, integ05, cases_int );
	PI0 = integ05[1];

	//printf("111 = %f  %f  \n", tobs, tp1 );
	//printf("111 = %f + %f  \n", creal(dd[1]), cimag(dd[1]) );
	//printf("111 = %f + %f  \n", creal(dd[2]), cimag(dd[2]) );
	//printf("111 = %f + %f  \n", creal(dd[3]), cimag(dd[3]) );

	if ( kt < zero )
              PI01 = - PI0;
	else
              PI01 = PI0;


	h = - b1 / four;
	// to get number of turn points of t1 and t2.
	/******************************************************************/
	// mu = mu_tp + b0 / ( four * tmu - b1 );
	weierstrass_int_J3(tp2,tp1,dd,del,a4,b4,index_p5,rff_p,integ15,cases_int);
	// equations (53) of Yang and Wang (2013).	
	p_mt1_mt2 = integ15[1];
	PI2_p = p_mt1_mt2 - PI0;

	// equations (53) of Yang and Wang (2013).

	PI1_phi = zero;
	PI2_phi = zero;
	PI1_sig = zero;
	PI2_sig = zero;
	PI1_time = zero;
	PI2_time = zero; 
	// equations (71)-(73) of Yang and Wang (2013).
	Delta = pt->r2 + pt->a2 - two * pt->robs;
	c_phi = pt->a_spin * ( pt->robs * two - ( pt->lambda * pt->a_spin ) ) / Delta;
	c_time = ( two * pow( pt->robs, three ) + ( two * pt->a_spin 
		* ( pt->a_spin - pt->lambda ) ) * pt->robs ) / Delta;

	index_p5[1] = -1;
	index_p5[2] = -2;
	index_p5[3] = 0;
	index_p5[4] = -4;
	index_p5[5] = 0;
}




/*
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2023-01-05.
!*
!*
!*/
static void Get_Sphmot_Results( ptcl *pt, double p, double theta_max, double kp, double kt, double *phyt,
	double *timet, double *sigmat, double *mucos, int *tt1, int *tt2 )
{
	double p1, p2, pp, p_mu = zero;
	double pp_sig, pp_t, p1_t, p1_sig, p1_phi;
	double p2_t, p2_sig, p2_phi;
	int t1, t2;

	t1 = 0;
	t2 = 0;

	if ( pt->mu_tp1 == zero ) {
		/* photons are confined in the equatorial plane, 
		so the integrations about \theta are valished. */ 
		*mucos = zero;
		*phyt = p;
		*timet = (*phyt) / Omega;
		*sigmat = c_tau * (*timet);
		*tt1 = 0;
		*tt2 = 0;
		return;
	}


	if ( pt->a_spin == zero ) {
		*timet = zero;
		SPINZERO( pt, p, kp, kt, theta_max, phyt, timet, sigmat, mucos, tt1, tt2 );
		return;
	}


	tmu = weierstrassP( p + PI01, g2, g3, dd, del );
	// equations (32) of Yang and Wang (2013).
	*mucos = pt->mu_tp1 + b0 / ( four * tmu - b1 ); 

	/******************************************************************/
	index_p5[1] = 0;
	cases_int = 1;
	weierstrass_int_J3(tobs,tmu,dd,del,a4,b4,index_p5,rff_p,integ5,cases_int);
	// equations (53) of Yang and Wang (2013).
	pp = integ5[1];
	p1 = PI0 - pp;
	p2 = p_mt1_mt2 - p1;
	// equations (54) of Yang and Wang (2013).
	for ( long int j = 0; j < 1.e5; j++ ) {
		for ( long int i = j; i <= j + 1; i++ ) {
			if ( pt->mobseqmtp ) {
				if( pt->muobs == pt->mu_tp1 ) {
					t1 = j;
					t2 = i;
					p_mu = -pp + two*( t1*p1 + t2*p2 );
				} else {
					t1 = i;
					t2 = j;
					p_mu = pp + two*( t1*p1 + t2*p2 );
				}
			} else {
				if (kt < zero) {
					t1 = i;
					t2 = j;
					p_mu = pp + two*( t1*p1 + t2*p2 );
				}
				if (kt > zero) {
					t1 = j;
					t2 = i;
					p_mu = -pp + two*( t1*p1 + t2*p2 );
				}
			}  
			if ( fabs( p - p_mu ) < 1.e-3 )
				goto gooutfor;
		}
	}

gooutfor:
	index_p5[1] = -1;
	index_p5[2] = -2;
	index_p5[3] = 0;
	index_p5[4] = -4;
	index_p5[5] = 0;
	// *****pp part***************************************
	if ( pt->lambda != zero ) {
		cases_int = 2;
		weierstrass_int_J3(tobs,tmu,dd,del,-a_add,b4,index_p5, fabs(pp),integ5,cases_int);
		weierstrass_int_J3(tobs,tmu,dd,del,-a_m,b4,index_p5, fabs(pp),integ15,cases_int);
		// equations (21) (72) of Yang and Wang (2013).
		pp_phi = ( pp / ( one- pt->mu_tp12 ) + (integ5[2] * c_add
			- integ15[2] * c_m ) / two ) * pt->lambda + c_phi * pp;
	} else 
		pp_phi = c_phi * pp;

	cases_int = 4;
	weierstrass_int_J3(tobs,tmu,dd,del,h,b4,index_p5, fabs(pp),integ,cases_int);
	// equations (20) (71) of Yang and Wang (2013).
	pp_sig = ( pp * pt->mu_tp12 + integ[2] * pt->mu_tp1 * b0 / two + integ[4] * b0 * b0
			/ sixteen ) * pt->a2 + pt->r2 * pp;
	// equations (22) of Yang and Wang (2013).
	pp_t = pp_sig + c_time * pp;


	// *****p1 part***************************************
	if ( t1 == 0 ) {
		p1_phi = zero;
		p1_sig = zero;
		p1_t = zero;
	} else {
		if ( pt->lambda != zero ) {
			if ( PI1_phi == zero ) {
				cases_int = 2;
				weierstrass_int_J3(tobs,infinity,dd,del,-a_add,
						b4,index_p5,PI0,integ5,cases_int);
				weierstrass_int_J3(tobs,infinity,dd,del,-a_m,
						b4,index_p5,PI0,integ15,cases_int);
				 //equations (21) (72) of Yang and Wang (2013).
				PI1_phi = ( PI0 / ( one - pt->mu_tp12 )
					+ ( integ5[2] * c_add - integ15[2] * c_m ) / two )
					* pt->lambda + c_phi * PI0;
			} 
			p1_phi = PI1_phi - pp_phi;
		} else {
			if ( PI1_phi == zero ) PI1_phi = c_phi * PI0;
			p1_phi = PI1_phi - pp_phi;
		}
		if ( PI1_time == zero || PI1_sig == zero) {
			cases_int = 4;
			weierstrass_int_J3(tobs,infinity,dd,del,h,b4,index_p5,PI0,integ,cases_int);
			// equations (20) (71) of Yang and Wang (2013).
			PI1_sig = ( PI0 * pt->mu_tp12 + integ[2] * pt->mu_tp1 * b0 / two
				+ integ[4] * b0 * b0 / sixteen ) * pt->a2 + pt->r2 * PI0;
			// equations (22) of Yang and Wang (2013).
			PI1_time = PI1_sig + c_time * PI0;
		}
		p1_sig = PI1_sig - pp_sig;
		p1_t = PI1_time - pp_t;
	}

	// *****p2 part***************************************
	if ( t2 == 0 ) {
		p2_sig = zero;
		p2_phi = zero;
		p2_t = zero;
	} else {
		if ( pt->lambda != zero ) {
			if ( PI2_phi == zero ) {
				cases_int = 2;	
				weierstrass_int_J3(tp2,tobs,dd,del,-a_add,b4,
						index_p5,PI2_p,integ5,cases_int);
				weierstrass_int_J3(tp2,tobs,dd,del,-a_m,b4,
						index_p5,PI2_p,integ15,cases_int);
				// equations (21) (72) of Yang and Wang (2013).
				PI2_phi = ( PI2_p / ( one - pt->mu_tp12 )
					+ ( integ5[2] * c_add - integ15[2] * c_m )
					/ two ) * pt->lambda + c_phi * PI2_p;
			}
			p2_phi = PI2_phi + pp_phi;
		} else {
			if ( PI2_phi == zero ) PI2_phi = c_phi * PI2_p;
			p2_phi = PI2_phi + pp_phi;
		}
                
		if ( PI2_time == zero ) {
			cases_int = 4;
			weierstrass_int_J3(tp2,tobs,dd,del,h,b4,index_p5,PI2_p,integ,cases_int);
			// equations (20) (71) of Yang and Wang (2013).
			PI2_sig = ( PI2_p * pt->mu_tp12 + integ[2] * pt->mu_tp1 * b0
				/ two + integ[4] * b0 * b0 / sixteen )
				* pt->a2 + pt->r2 * PI2_p;
		// equations (22) of Yang and Wang (2013).
			PI2_time = PI2_sig + c_time * PI2_p;
		}
		p2_sig = PI2_sig + pp_sig;
		p2_t = PI2_time + pp_t;
	}



	if ( pt->mobseqmtp ) {
		if ( pt->muobs == pt->mu_tp1 ) {
			// equations (52) of Yang and Wang (2013).
			*phyt = -pp_phi + two * ( t1*p1_phi + t2*p2_phi );
			*timet = -pp_t + two * ( t1*p1_t + t2*p2_t );
			*sigmat = -pp_sig + two * ( t1*p1_sig + t2*p2_sig );
		} else {
			// equations (52) of Yang and Wang (2013).
			*phyt = pp_phi + two*( t1*p1_phi + t2*p2_phi );
			*timet = pp_t + two*( t1*p1_t + t2*p2_t );
			*sigmat = pp_sig + two*( t1*p1_sig + t2*p2_sig );
		}
	} else {
		if ( kt < zero ) {
			// equations (52) of Yang and Wang (2013).
			*phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi);	
			*timet=pp_t+two*(t1*p1_t+t2*p2_t);
			*sigmat=pp_sig+two*(t1*p1_sig+t2*p2_sig);
		}
		if ( kt > zero) {
			// equations (52) of Yang and Wang (2013).		
			*phyt = -pp_phi + two*( t1*p1_phi + t2*p2_phi );
			*timet = -pp_t + two*( t1*p1_t + t2*p2_t );
			*sigmat = -pp_sig + two*( t1*p1_sig + t2*p2_sig );
		}
	}
	if ( pt->mu_tp1 == one ) *phyt = *phyt + ( t1 + t2) * pi;
	// phyt = fmod(phyt, twopi);
	// if ( phyt < zero ) phyt = phyt + twopi;
	*tt1 = t1;
	*tt2 = t2;
}






/*
!*
!*
!*     PURPOSE:   To compute the constants of motion for the spherical motion of a photon. 
!*     INPUTS:    r_sm--------------the radius of the spherical motion.
!*                theta_max---------the maximum or minimum value the \theta coordinate of the 
!*                                  spherical motion, which also the turning point of the motion
!*                                  in theta coordinate.  
!*                signs-------------which determine the direction of the motion of a photon with respect 
!*                                  to the black hole. signs>0 is for co-rotating/prograde orbiting,
!*                                  signs< is for counter-rotating/retrograde orbiting.
!*                a_spin------------the black hole spin.
!*     OUTPUTS:   lamda,q-----------constants of motion.    
!*                theta_min---------the minimum value of the theta coordinate the orbiting.                         
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: 
!*
!*
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2023-01-04.
!*
!*
!*
!*/
void lambdaq_sphericalm( ptcl *pt, double r_sm, double *theta_min, double *rmin, double *rmax )
{
	double r_max, c_temp, r_min;
	double r_sm2, r_sm3;

	r_max = radiusofsphericalmotion( pt->a_spin, zero );
	c_temp = 90.0;
	r_min = radiusofsphericalmotion( pt->a_spin, c_temp );
	*rmin = r_min;
	*rmax = r_max;
	if ( pt->a_spin < zero ) {
		c_temp = r_max;
		r_max = r_min;
		r_min = c_temp;
	}

      
	if ( pt->a_spin != zero ) {
		if ( r_sm < r_min || r_sm > r_max || r_sm - r_min <= -1.e-6
		|| r_sm - r_max >= 1.e-6) {
			printf("lambdaq_sphericalm(): For spin =%f the \n\
			allowed range for radial coordinate of the spherical \n\
			motion of a photon is between rmin = %f and rmax = %f. \n\
			The radius you input is out of the range and \n\
			the code must be stopped.\n", pt->a_spin, r_min, r_max );
			exit(0);
		}
	} else {
		if ( r_sm != three ) {
			printf("For a = 0, the radius of the spherical motion of \n\
			a photon is 3 for all inclination anges of the orbit \n\
			with respect to the equatorial plane.\n");
			r_sm = three;
		}
	}
	r_sm2 = pow( r_sm, 2 );
	r_sm3 = pow( r_sm, 3 );

	double mutp2, sintp2;
	double a1, b1, c1;
	if ( pt->a_spin != zero ) {
		a1 = sq( pt->a2 ) * sq( r_sm - one );
		b1 = two * pt->a2 * r_sm * ( two * pt->a2 - three * r_sm + r_sm3 );
		c1 = r_sm3 * ( - four * pt->a2 + r_sm * sq( r_sm - three ) );
	    
		//double mu_starm;
		mutp2 = ( - b1 + sqrt( b1 * b1 - four * a1 * c1 ) ) / two / a1;
		sintp2 = one - mutp2;
		//mu_starm = ( - b1 - sqrt( b1 * b1 - four * a1 * c1 ) ) / two / a1;
		*theta_min = acos( sqrt( mutp2 ) ) / dtor;
	} else {
		printf(" For a = 0, we need you to input the minimum of the \
			theta coordinate:");
		scanf(" Please input theta_min: %lf", theta_min);
		if ( *theta_min != 90.0 ) {
			sintp2 = sq( sin( *theta_min * dtor ) );
			mutp2 = sq( cos( *theta_min * dtor ) );
		} else {
			sintp2 = one;
			mutp2 = zero;
		}
	}


	double Delta, AL, AE;
	//cotmax = costhemax / sinthemax;
	//Sigma = r_sm2 + pt->a2 * mutp2;
	Delta = r_sm2 - two * r_sm + pt->a2;

	AL = Delta * ( three * r_sm2*r_sm2 + pt->a2 * r_sm2
		* ( mutp2 + one ) - sq( pt->a2 ) * mutp2 );
	AE = Delta * ( r_sm2 - pt->a2 * mutp2 );

	if ( pt->a_spin != zero )
		pt->lambda = sign( one, pt->a_spin ) * sqrt( sintp2 * AL / AE );
	else
		pt->lambda = sqrt( sintp2 * AL / AE );

	pt->lam2 = sq( pt->lambda );
	pt->q = mutp2 * ( AL / AE - pt->a2 );
	//printf("AL AE = %f %f %f \n", AL, AE, Delta);
}




static double Denomenator_r( double r, double theta, double a_spin )
{
	double ant, sinthe, costhe, a2;

	a2 = a_spin * a_spin;
	if ( theta != 90.0 ) {
		sinthe = sin( theta * dtor );
		costhe = cos( theta * dtor );
	} else {
 		sinthe = one;
		costhe = zero;
	}
	ant =pow( ( - ( -three * r * r + ( r + one ) * a2
		* sq( costhe ) + two * a_spin * sinthe
		* sqrt( r * ( r * r - a2 * sq( costhe ) ) ) ) ), (one/three) ); 
	return ant;
}



static double radiusofsphericalmotion( double a, double c )
{
	double r1, r2, Dr;

	r1 = 30.0;
	do {
		r2 = Denomenator_r( r1, c, a);
		Dr = r2 - r1;
		r1 = r2;
	} while( fabs( Dr ) >= 1.e-10 );

	return r1;
}





/*
!*
!*     PURPOSE:   To compute the integral \int^y_x z^2/sqrt(1-(z/z0)^2)dz, and z0<1. 
!*     INPUTS:    x,y-----------the integral limits.    
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS:
!*
!*
!*
!*/
static double mveone_int( double x, double y, double z0 )
{
	double xt, yt, Iy, Ix, z02, z03;

	xt = x;
	yt = y;
	if (xt == yt) return zero;
	z02 = z0 * z0;
	z03 = z02 * z0;

	Iy = z03 * half * ( asin( yt / z0 ) - y / z02 * sqrt( z02 - y * y ) );
	Ix = z03 * half * ( asin( xt / z0 ) - x / z02 * sqrt( z02 - x * x ) );
	return ( Iy - Ix );
}




/*
!*
!*
!*
!*     PURPOSE:   This routine computs the four B_L coordiants r,\mu,\phi,t and affine 
!*                parameter \sigma of the spherical motion when black hole spin a is zero. 
!*     INPUTS:   p--------------the independent variable.
!*               kt-------------p_\theta, the \theta component of four momentum of the photon.
!*               kp-------------p_\phi,  the \phi component of four momentum of the photon.
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or the initial position of photon. 
!*               thetamax-------the maximum or minimum value of the \theta coordinate of the geodesic.
!*                              which also the \theta turning point of the spherical motion.       
!*     OUTPUTS:  phyt-----------the \phi coordinat of the photon.
!*               timet----------the t coordinats of the photon.
!*               sigmat---------the affine parameter \sigma.
!*               mucos----------the \theta coordinates of the photon, and mucos=cos(\theta).  
!*               t1,t2----------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively.      
!*     ROUTINES CALLED: metricg, circ_mb, weierstrass_int_J3, weierstrassP
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: 
!*
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2023-01-04.
!*
!*
!*
!*/
static void SPINZERO( ptcl *pt, double p, double kp, double kt, double theta_max, 
		double *phyt, double *timet, double *sigmat,
		double *mucos, int *tt1, int *tt2 )
{
	double pp_sigma, pp_time, pp_phi;
	double p1_sigma, p1_time, p1_phi;
	double p2_sigma, p2_time, p2_phi;
	static double AA, BB;
	static double PI1_phi, PI2_phi, PI1_time, 
		PI2_time, PI1_sigma, PI2_sigma;
	double p_mu = zero, Ptotal, PI1, PI2;
	double pp, p1, p2;
	int t1, t2;

	t1 = 0;
	t2 = 0;
	pt->mobseqmtp = false;



	if ( pt->q > zero ) {
		double mu;
		AA = sqrt( ( pt->q + pt->lam2 ) / pt->q );
		BB = sqrt( pt->q );

		if ( kt > zero )
			mu = sin( asin( pt->muobs * AA ) - p * BB * AA ) / AA;	
		else {
			if ( kt == zero )
				mu = cos( p * AA * BB ) * pt->muobs;
			else			      
				mu = sin( asin( pt->muobs * AA ) + p * AA * BB ) / AA;
		}
		*mucos = mu;

		if ( kt != zero ) {
			pt->mu_tp1 = sqrt( pt->q / ( pt->lam2 + pt->q ) );
			pt->mu_tp2 = - pt->mu_tp1;
		} else {
			pt->mu_tp1 = fabs( pt->muobs );
			pt->mu_tp2 = - pt->mu_tp1;
			pt->mobseqmtp = true;
		}
		if ( fabs( pt->muobs ) == one ) pt->mobseqmtp = true;

		if ( pt->mu_tp1 == zero ) {
			// photons are confined in the equatorial plane, 
			// so the integrations about !\theta are valished.
			*timet = zero;
			*sigmat = zero;
			*phyt = zero;
			*tt1 = 0;
			*tt2 = 0;
			return;
		}



		PI1 = (halfpi - asin( pt->muobs / pt->mu_tp1 ) ) * pt->mu_tp1 / BB;
		Ptotal = pi * pt->mu_tp1 / BB;
		PI2 = Ptotal - PI1;
		pp = ( asin( mu / pt->mu_tp1 ) - asin( pt->muobs / pt->mu_tp1 ) ) * pt->mu_tp1 / BB;
		p1 = PI1 - pp;
		p2 = Ptotal - p1;
		PI1_phi = zero;
		PI2_phi = zero;
		PI1_time = zero;
		PI2_time = zero;
		PI1_sigma = zero;
		PI2_sigma = zero;

		for ( int j = 0; j < 1e6; j++ ) {
			for ( int i = j; i <= j + 1; i++ ) {
				if ( pt->mobseqmtp ) {
					if( pt->muobs == pt->mu_tp1 ) {
						t1 = j;
						t2 = i;
						p_mu = -pp + two * ( t1 * p1 + t2 * p2 );
					} else {
						t1 = i;
						t2 = j;
						p_mu = pp + two * ( t1 * p1 + t2 * p2 );
					}
				} else {
					if ( kt < zero ) {
						t1 = i;
						t2 = j;
						p_mu = pp + two * ( t1 * p1 + t2 * p2 );
					}
					if ( kt > zero ) {
						t1 = j;
						t2 = i;
						p_mu = -pp + two * ( t1 * p1 + t2 * p2 );
					}
				}  
				if ( fabs( p - p_mu ) < 1.e-6 ) goto temp_one;
			}
		}
temp_one:
		// *************** p0 part *************** !
		Delta = pt->r2 + pt->a2 - two * pt->robs;
		c_phi = pt->a_spin * ( pt->robs * two - ( pt->lambda * pt->a_spin ) ) / Delta;
		c_time = ( two * pow( pt->robs, 3 ) + ( two * pt->a_spin
			* ( pt->a_spin - pt->lambda ) ) * pt->robs ) / Delta;


		pp_sigma = pt->a2 * mveone_int( pt->muobs, mu, one / AA ) / BB + pt->r2 * pp;
		pp_time = pp_sigma + c_time * pp;
		pp_phi = pt->lambda * Schwarz_integral( pt->muobs, mu, AA ) / BB + c_phi * pp;


		// *************** p1 part *************** !
		if ( t1 == 0 ) {
			p1_sigma = zero;
			p1_time = zero;
			p1_phi = zero;
		} else {
			if ( PI1_time == zero || PI1_sigma == zero) {
				PI1_sigma = pt->a2 * mveone_int( pt->muobs, 
					pt->mu_tp1, one / AA ) / BB + pt->r2 * PI1;
				PI1_time = PI1_sigma + c_time * PI1;
			}
			if ( PI1_phi == zero )
				PI1_phi = pt->lambda * Schwarz_integral( pt->muobs,
					pt->mu_tp1, AA ) / BB + c_phi * PI1;

			p1_sigma = PI1_sigma - pp_sigma;
			p1_time = PI1_time - pp_time;
			p1_phi = PI1_phi - pp_phi;
		}
 
		// *************** p2 part *************** !
		if ( t2 == 0 ) {
			p2_sigma = zero;
			p2_time = zero;
			p2_phi = zero;
		} else {
			if ( PI2_time == zero || PI2_sigma == zero ) {
				PI2_sigma = pt->a2 * mveone_int( pt->mu_tp2,
					pt->muobs, one / AA ) / BB + pt->r2 * PI2;
				PI2_time = PI2_sigma + c_time * PI2;
			}
			if ( PI2_phi == zero )
				PI2_phi = pt->lambda * Schwarz_integral( pt->mu_tp2,
					pt->muobs, AA ) / BB + c_phi * PI2;

			p2_sigma = PI2_sigma + pp_sigma;
			p2_time = PI2_time + pp_time;
			p2_phi = PI2_phi + pp_phi;
		}


		if ( pt->mobseqmtp ) {
			if ( pt->muobs == pt->mu_tp1 ) {
				// equations (52) of Yang and Wang (2013).	
				*sigmat = -pp_sigma+two*(t1*p1_sigma+t2*p2_sigma);
				*timet = -pp_time+two*(t1*p1_time+t2*p2_time);
				*phyt = -pp_phi+two*(t1*p1_phi+t2*p2_phi);
			} else {
				// equations (52) of Yang and Wang (2013).	
				*sigmat = pp_sigma+two*(t1*p1_sigma+t2*p2_sigma);
				*timet = pp_time+two*(t1*p1_time+t2*p2_time);
				*phyt = pp_phi+two*(t1*p1_phi+t2*p2_phi);
			} 
		} else {
		         if ( kt < zero ) {
				// equations (52) of Yang and Wang (2013).	
		              *sigmat = pp_sigma+two*(t1*p1_sigma+t2*p2_sigma);
		              *timet = pp_time+two*(t1*p1_time+t2*p2_time);
		              *phyt = pp_phi+two*(t1*p1_phi+t2*p2_phi);
			}
			if (kt > zero) {
				// equations (52) of Yang and Wang (2013).	
		              *sigmat = -pp_sigma+two*(t1*p1_sigma+t2*p2_sigma);
		              *timet = -pp_time+two*(t1*p1_time+t2*p2_time);
		              *phyt = -pp_phi+two*(t1*p1_phi+t2*p2_phi);
			}
		}

		*tt1 = t1;
		*tt2 = t2;
		if ( theta_max == zero || theta_max == 180.0 )
			*phyt = *phyt + ( t1 + t2) * pi;
	} else {
		printf(" phyt_schwatz(): q<0, which is a affending value, \
			the program should be stoped! and q = %f\n ", pt->q );
		exit(0);
		*mucos = pt->muobs;
		t1 = 0;
		t2 = 0;
		*phyt = zero;
		*timet = zero;
		*sigmat = zero;
		*tt1 = t1;
		*tt2 = t2;
	}
}




/*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-17.
!*
!*/
void YNOGK_sphmotion( ptcl *pt, double p, double theta_max, double kp, 
	double kt, double *phyt, double *timet, double *sigmat, 
	double *mucos, double *sint, int *t1, int *t2 )
{
	static unsigned int count_num = 1;
 
sphm_restarting:
	if ( count_num == 1 ) {

		Spherical_motion_Settings( pt, theta_max, kt );
		count_num++;

		kp_1 = kp;
		kt_1 = kt;
		lambda_1 = pt->lambda;
		q_1 = pt->q;
		muobs_1 = pt->muobs;
		sinobs_1 = pt->sinobs;
		a_spin_1 = pt->a_spin;
		robs_1 = pt->robs;
		thetamax_1 = theta_max;
 
		Get_Sphmot_Results( pt, p, theta_max, kp, kt, phyt,
				timet, sigmat, mucos, t1, t2 );

		*sint = sqrt( one - sq( *mucos ) ); 

    	} else {
		if ( kp_1 == kp &&
		kt_1 == kt &&
		lambda_1 == pt->lambda &&
		q_1 == pt->q &&
		muobs_1 == pt->muobs &&
		sinobs_1 == pt->sinobs &&
		a_spin_1 == pt->a_spin &&
		robs_1 == pt->robs &&
		thetamax_1 == theta_max ) {
 
			Get_Sphmot_Results( pt, p, theta_max, kp, kt, phyt,
					timet, sigmat, mucos, t1, t2 );
			*sint = sqrt( one - sq( *mucos ) );

		} else {
			count_num = 1;
			goto sphm_restarting;
		}

	}
}













