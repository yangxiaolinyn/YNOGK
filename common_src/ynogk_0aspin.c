/*
 * common_src/ynogk_0aspin.c
 *
 * This module aim on solve cubic and quartic polynomial equations.
 * One can use these subroutine root3 and root4 to find roots of cubic and 
 * quartic equations respectively. 
 *
 * Author       Yang Xiao-lin
 *
 * Date         2022-12-17
 *
 * City         Kun-ming Yunnan Provice, China
 *
 * 2022-12-11   Stating the writting of the module.
 */
 

#include "ynogk_0aspin.h"



/*****************************************************************************\
!*
!*  Private Variables for this module.
!*
\*****************************************************************************/


static double a_spin, lambda, lam2, q, muobs, sinobs, scal, f2, f3;
static int mobseqmtp;


static double mu_tp1, mu_tp2;

static double AA, BB;

static double pp, p1, p2;
static double PI1, PI2, Pt;
static double PI1_phi, PI2_phi;
/*****************************************************************************/




static double a_spin_1;
static double scal_1;
static double lambda_1;
static double q_1;
static double f3_1;
static double f2_1;
static double muobs_1;
static double sinobs_1;

/*
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! To determine the number of time_0 N_t1, N_t2 that the particle meets the
! two turn points mu_tp1, mu_tp2 respectively.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-17.
!*/
static void Get_tm1_tm2_New( double p, double PI1, double PI2, double Pt, 
		int mobseqmtp, double f2, 
		double f3, double muobs, double mu_tp1, 
		double mu_tp2, int *t1, int *t2 )
{
	double p1_temp, p2_temp;
	int N_temp;

	if ( mobseqmtp ) {
		p1_temp = zero;
		p2_temp = Pt;
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
					p2_temp = p2_temp + Pt;
					*t1 = *t2;
					*t2 = N_temp - *t1;
				}
			} while (true);
		} else if ( muobs == mu_tp2) {
			do {
				if ( p1_temp <= p && p <= p2_temp )
					break;
				else {
					N_temp = N_temp + 1;
					p1_temp = p2_temp;
					p2_temp = p2_temp + Pt;
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
		if ( f2 > zero ) {
			p2_temp = PI2;
			do {
				if ( p1_temp <= p && p <= p2_temp )
					break;
				else {
					N_temp = N_temp + 1;
					p1_temp = p2_temp;
					p2_temp = p2_temp + Pt;
					*t1 = *t2;
					*t2 = N_temp - *t1;
				}
			} while (true);
		}
		if ( f2 < zero ) {
			p2_temp = PI1;
			do {
				if( p1_temp <= p && p <= p2_temp )
					break;
				else {
					N_temp = N_temp + 1;
					p1_temp = p2_temp;
					p2_temp = p2_temp + Pt;
					*t2 = *t1;
					*t1 = N_temp - *t2;
				}
			} while (true);
		}
	}
}


/*
!*     PURPOSE:  Computes \int^x_y dt/(1-t^2)/sqrt(1-AA^2*t^2) and AA .gt. 1  
!*     INPUTS:   components of above integration.      
!*     OUTPUTS:  valve of integral.             
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-17.
!*
!*/
/*************************************************************************/
double Schwarz_integral( double y, double x, double AA ) 
/*************************************************************************/ 
{
	double yt, xt, schwatz_int, ppx, ppy, A2;

	xt = x;
	yt = y;
	if ( yt == xt ) {
		return (zero);
	}

	if ( fabs( AA ) != one ) {
		A2 = sq( AA );
		ppx = atan( sqrt( A2 - one ) * xt / sqrt( fabs( one - A2*xt*xt ) ) );
		ppy = atan( sqrt( A2 - one ) * yt / sqrt( fabs( one - A2*yt*yt ) ) );
		schwatz_int = ( ppx - ppy ) / sqrt( A2 - one );
	} else {
		if ( fabs(xt) == one )
			schwatz_int = infinity;
		else {
			if ( fabs(yt) == one )
				schwatz_int = -infinity;
			else {
				ppx = xt / sqrt( fabs( one - xt*xt ) );
				ppy = yt / sqrt( fabs( one - yt*yt ) );
				schwatz_int = ppx - ppy;
			}
		}
	}
	return (schwatz_int);
}



/*
!*     PURPOSE:  Computes \mu part of integrals in coordinates \phi, expressed by equation (72) 
!*               in Yang & Wang (2012) with zero spin of black hole.    
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f2-------------p_\theta, which is the \theta component of four momentum of a photon 
!*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
!*               f3-------------p_\phi, which is the \phi component of four momentum of a photon 
!*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  phyc_schwatz-----------value of integral \phi_\theta expressed by equation (71) in
!*                              Yang & Wang (2012).   
!*               mucos----------value of function \mu(p) with zero spin.
!*               t1,t2----------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively.            
!*     ROUTINES CALLED: schwatz_int
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: 
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-17.
!*
!*/
/*****************************************************************************/
int phit_Schwarzschild( ptcl *pt, double p, double *phit_Schwarz,
		double *mucos, double *sign_pth, int *t1, int *t2 )
/*****************************************************************************/
{
	double u;
	double pp_phi, p1_phi, p2_phi;


	lambda = pt->lambda;
	lam2 = pt->lam2;
	q = pt->q;
	sinobs = pt->sinobs;
	muobs = pt->muobs;
	scal = pt->scal;

	f2 = pt->f1234[2];
	f3 = pt->f1234[3];

	a_spin_1 = zero;
	lambda_1 = lambda;
	q_1 = q;
	muobs_1 = muobs;
	sinobs_1 = sinobs;
	scal_1 = scal;
	f3_1 = f3;
	f2_1 = f2;
	*t1 = 0;
	*t2 = 0;
	

	mobseqmtp = false;
	if ( q > zero ) {
		AA = sqrt( ( lam2 + q ) / q );
		BB = sqrt( q );

		if ( f2 < zero ) {
			u = asin( muobs * AA ) + p * BB * AA;
			// mucos=dsin(dasin(muobs*AA)+p*BB*AA)/AA;
			*mucos = sin( u ) / AA;
			u = fmod( u + halfpi, twopi );
			if ( zero <= u && u <= pi )
		              *sign_pth = - one;
			else
		              *sign_pth = one;

		} else {
			if (f2 == zero) {
				u = p * AA * BB;
				// mucos=dcos(p*AA*BB)*muobs;
				*mucos = cos(u) * muobs;
				u = fmod( u, twopi );
				if ( zero <= u && u <= pi )
					*sign_pth = sign( one, muobs );
				else
					*sign_pth = - sign( one, muobs );

			} else {
				u = asin( muobs * AA ) - p * BB * AA;
				//mucos = sin( asin( muobs * AA ) - p * AA * BB ) / AA;
				*mucos = sin( u ) / AA;
				u = fmod( u - halfpi, twopi );
				if ( -pi <= u && u <= zero )
					*sign_pth = one;
				else
					*sign_pth = -one;
     
			}
		}
		if ( f2 != zero ) {
			mu_tp1 = sqrt( q / ( lam2 + q ) );
			mu_tp2 = - mu_tp1;
		} else {
			mu_tp1 = fabs( muobs );
			mu_tp2 = - mu_tp1;
			mobseqmtp = true;
		}
		if ( fabs( muobs ) == one )
			mobseqmtp= true;

		if ( mu_tp1 == zero ) {
			// photons are confined in the equatorial plane, 
			// so the integrations about !\theta are valished.
			*sign_pth = zero;
			*phit_Schwarz = zero;
			return 0;
		}

		PI1 = ( halfpi - asin( muobs / mu_tp1 ) ) * mu_tp1 / BB;
		Pt = pi * mu_tp1 / BB;
		PI2 = Pt - PI1;
		pp = ( asin( *mucos / mu_tp1 ) - asin( muobs / mu_tp1 ) ) * mu_tp1 / BB;
		p1 = PI1 - pp;
		p2 = Pt - p1;
		PI1_phi = zero;
		PI2_phi = zero;

		/* To get the number of turnning points: t1 and t2. */
		Get_tm1_tm2_New( p, PI1, PI2, Pt, mobseqmtp, f2, 
		f3, muobs, mu_tp1, mu_tp2, t1, t2 );


		if ( lambda == zero ) {
			*phit_Schwarz = zero;
			return 0;
		}

		pp_phi = lambda * Schwarz_integral( muobs, *mucos, AA ) / BB;
		if ( *t1 == 0 )
			p1_phi = zero;
		else {
			if ( PI1_phi == zero )
				PI1_phi = lambda * 
				Schwarz_integral( muobs, mu_tp1, AA ) / BB;

			p1_phi = PI1_phi - pp_phi;
		}

		if ( *t2 == 0 )
			p2_phi = zero;
		else {
			if ( PI2_phi == zero )
				PI2_phi = lambda * Schwarz_integral( mu_tp2, muobs, AA) / BB;

			p2_phi = PI2_phi + pp_phi;
		}

		if ( mobseqmtp ) {
			if ( muobs == mu_tp1 )
				*phit_Schwarz = - pp_phi + two * ( *t1 * p1_phi + *t2 * p2_phi );
			else
				*phit_Schwarz = pp_phi + two * ( *t1 * p1_phi + *t2 * p2_phi );
       
		} else {
			if ( f2 < zero )
				*phit_Schwarz = pp_phi + two * ( *t1 * p1_phi + *t2 * p2_phi );

			if ( f2 > zero)
		              *phit_Schwarz = -pp_phi + two * ( *t1 * p1_phi + *t2 * p2_phi );
   
		}

	} else {
              /*write(unit=6,fmt=*)'phyt_schwatz(): q<0, which is an affending',&
              !                'value, the program should be',&  
              !                'stoped! and q = ',q
              !stop */
              *mucos = muobs;
              *t1 = 0;
              *t2 = 0;
              *phit_Schwarz = zero;
              *sign_pth = zero;
	}
	return 0;
}



/*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-17.
!*
!*/
/*****************************************************************************/
static void phit_Schwarzschild_Settings( ptcl *pt )
/*****************************************************************************/
{
	a_spin = pt->a_spin;
	lambda = pt->lambda;
	lam2 = pt->lam2;
	q = pt->q;
	sinobs = pt->sinobs;
	muobs = pt->muobs;
	scal = pt->scal;

	f2 = pt->f1234[2];
	f3 = pt->f1234[3];

	a_spin_1 = a_spin;
	lambda_1 = lambda;
	q_1 = q;
	muobs_1 = muobs;
	sinobs_1 = sinobs;
	scal_1 = scal;
	f3_1 = f3;
	f2_1 = f2;

	mobseqmtp = false;
	if ( q > zero ) {
		AA = sqrt( ( lam2 + q ) / q );
		BB = sqrt( q );
 
		if ( f2 != zero ) {
			mu_tp1 = sqrt( q / ( lam2 + q ) );
			mu_tp2 = - mu_tp1;
		} else {
			mu_tp1 = fabs( muobs );
			mu_tp2 = - mu_tp1;
			mobseqmtp = true;
		}
		if ( fabs( muobs ) == one )
			mobseqmtp= true;
 

		PI1 = ( halfpi - asin( muobs / mu_tp1 ) ) * mu_tp1 / BB;
		Pt = pi * mu_tp1 / BB;
		PI2 = Pt - PI1;

		PI1_phi = zero;
		PI2_phi = zero;

		/* To get the number of turnning points: t1 and t2. */
		//Get_tm1_tm2_New( p, PI1, PI2, Pt, mobseqmtp, f2, 
		//f3, muobs, mu_tp1, mu_tp2, t1, t2 );

	} else {
	}
}



/*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-17.
!*
!*/
/*****************************************************************************/
static int Get_phit_Schwarzschild( ptcl *pt, double p, double *phit_Schwarz,
		double *mucos, double *sign_pth, int *t1, int *t2 )
/*****************************************************************************/
{
	double u;
	double pp_phi, p1_phi, p2_phi;


	if ( q > zero ) {
		if ( f2 < zero ) {
			u = asin( muobs * AA ) + p * BB * AA;
			// mucos=dsin(dasin(muobs*AA)+p*BB*AA)/AA;
			*mucos = sin( u ) / AA;
			u = fmod( u + halfpi, twopi );
			if ( zero <= u && u <= pi )
		              *sign_pth = - one;
			else
		              *sign_pth = one;

		} else {
			if (f2 == zero) {
				u = p * AA * BB;
				// mucos=dcos(p*AA*BB)*muobs;
				*mucos = cos(u) * muobs;
				u = fmod( u, twopi );
				if ( zero <= u && u <= pi )
					*sign_pth = sign( one, muobs );
				else
					*sign_pth = - sign( one, muobs );

			} else {
				u = asin( muobs * AA ) - p * BB * AA;
				//mucos = sin( asin( muobs * AA ) - p * AA * BB ) / AA;
				*mucos = sin( u ) / AA;
				u = fmod( u - halfpi, twopi );
				if ( -pi <= u && u <= zero )
					*sign_pth = one;
				else
					*sign_pth = -one;
     
			}
		}

		if ( mu_tp1 == zero ) {
			// photons are confined in the equatorial plane, 
			// so the integrations about !\theta are valished.
			*sign_pth = zero;
			*phit_Schwarz = zero;
			return 0;
		}

		pp = ( asin( *mucos / mu_tp1 ) - asin( muobs / mu_tp1 ) ) * mu_tp1 / BB;
		p1 = PI1 - pp;
		p2 = Pt - p1;

		/* To get the number of turnning points: t1 and t2. */
		Get_tm1_tm2_New( p, PI1, PI2, Pt, mobseqmtp, f2, 
		f3, muobs, mu_tp1, mu_tp2, t1, t2 );


		if ( lambda == zero ) {
			*phit_Schwarz = zero;
			return 0;
		}

		pp_phi = lambda * Schwarz_integral( muobs, *mucos, AA ) / BB;
		if ( *t1 == 0 )
			p1_phi = zero;
		else {
			if ( PI1_phi == zero )
				PI1_phi = lambda * 
				Schwarz_integral( muobs, mu_tp1, AA ) / BB;

			p1_phi = PI1_phi - pp_phi;
		}

		if ( *t2 == 0 )
			p2_phi = zero;
		else {
			if ( PI2_phi == zero )
				PI2_phi = lambda * Schwarz_integral( mu_tp2, muobs, AA) / BB;

			p2_phi = PI2_phi + pp_phi;
		}

		if ( mobseqmtp ) {
			if ( muobs == mu_tp1 )
				*phit_Schwarz = - pp_phi + two * ( *t1 * p1_phi + *t2 * p2_phi );
			else
				*phit_Schwarz = pp_phi + two * ( *t1 * p1_phi + *t2 * p2_phi );
       
		} else {
			if ( f2 < zero )
				*phit_Schwarz = pp_phi + two * ( *t1 * p1_phi + *t2 * p2_phi );

			if ( f2 > zero)
		              *phit_Schwarz = -pp_phi + two * ( *t1 * p1_phi + *t2 * p2_phi );
   
		}

	} else {
              /*write(unit=6,fmt=*)'phyt_schwatz(): q<0, which is an affending',&
              !                'value, the program should be',&  
              !                'stoped! and q = ',q
              !stop */
              *mucos = muobs;
              *t1 = 0;
              *t2 = 0;
              *phit_Schwarz = zero;
              *sign_pth = zero;
	}
	return 0;
}





/*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-17.
!*
!*/
/*****************************************************************************/
void Get_phit_Integrals_Schwarzschild( ptcl *pt, double p, 
		double *phit_Schwarz, double *mucos, 
		double *sign_pth, int *t1, int *t2 )
/*****************************************************************************/
{
	static unsigned int count_num = 1;

restarting:
	if ( count_num == 1 ) {

		phit_Schwarzschild_Settings( pt );

		count_num += 1;

		Get_phit_Schwarzschild( pt, p, phit_Schwarz, mucos, 
			sign_pth, t1, t2 );

    	} else {
		if ( f3_1 == pt->f1234[3]
		&& f2_1 == pt->f1234[2]
		&& lambda_1 == lambda
		&& q_1 == q
		&& sinobs_1 == sinobs
		&& muobs_1 == muobs
		&& a_spin_1 == a_spin
		&& scal_1 == scal ) {

			Get_phit_Schwarzschild( pt, p, phit_Schwarz, mucos, 
				sign_pth, t1, t2 );
		} else {
			count_num = 1;
			goto restarting;
		}

	}
	
}


























