/*
 * common_src/root34.c
 *
 * This module aim on solve cubic and quartic polynomial equations.
 * One can use these subroutine root3 and root4 to find roots of cubic and 
 * quartic equations respectively. 
 *
 * Author       Yang Xiao-lin
 *
 * Date         2022-11-16
 *
 * City         Kun-ming Yunnan Provice, China
 *
 * 2022-12-11   Stating the writting of the module.
 */
 

#include "ynogkBL.h"




/*
!*     PURPOSE:  Computes four Boyer-Lindquist coordinates (r,\mu,\phi,t) and affine parameter 
!*               \sigma as functions of parameter p, i.e. functions r(p), \mu(p), \phi(p), t(p)
!*               and \sigma(p). Cf. discussions in Yang & Wang (2012).    
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f1234----------array of p_1, p_2, p_3, p_4, which are the components of four-
!*                              momentum of a photon measured under the LNRF frame. This array 
!*                              can be computed by subroutine lambdaq(...), see below.     
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or initialposition of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  radi-----------value of function r(p). 
!*               mu-------------value of function \mu(p). 
!*               phi------------value of function \phi(p).
!*               time-----------value of function t(p).
!*               sigma----------value of function \sigma(p).
!*               sign_pr--------sign of r component of 4-momentum of the photon.
!*               sign_pth-------sign of \theta component of 4-momentum of the photon.
!*               tm1,tm2--------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively. 
!*               tr1,tr2--------number of times of photon meets turning points r_tp1 and r_tp2
!*                              respectively.            
!*     ROUTINES CALLED: INTRPART, INTTPART.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS:
!*
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-11.
!*
!*/
void YNOGK( ptcl *p, double pm, double *radi, double *mu, double *phi, 
		double *time, double *sigma, double *sign_pr, double *sign_pth )
{
	double time_r, time_t, aff_r, phi_r, phi_t;
	double Rab;
	int rotate;
	int tm1, tm2;
  
	Get_Integrals_For_R_Part( p, pm, radi, sign_pr, &aff_r, &time_r, &phi_r );
	Get_Integrals_For_Theta_Part( p, pm, &phi_t, &time_t, mu, sign_pth, &tm1, &tm2 );  

 
// time coordinate value, equation (74) of Yang & Wang (2012).
	//printf("here11 = time_r = %f  time_t = %f \n", time_r, time_t);
        *time = time_r + time_t;
// affine parameter value, equation (74) of Yang & Wang (2012).

	//printf("here22 = aff_r = %f  time_t = %f \n", aff_r, time_t);
        *sigma = aff_r + time_t;
// phi coordinate value.
        rotate = false;
 

	if ( fabs( p->muobs ) != one ) {
		// equation (74) of Yang & Wang (2012).
		*phi = - ( phi_r + phi_t );
		if ( p->f1234[3] == zero )
			*phi = *phi + ( tm1 + tm2 ) * pi;

		*phi = fmod( *phi, twopi );
		if ( *phi < zero )
			*phi = *phi + twopi;

		//printf( "phi2 = %f 2pi = %f \n", phi_r, phi_t );
	} else {
		// equation (74) of Yang & Wang (2012).
		*phi = - ( phi_t + phi_r + ( tm1 + tm2 ) * pi );
		Rab = sqrt( sq( p->f1234[3] ) + sq( p->f1234[2] ) );
		if ( *phi != zero )
			rotate = true;
 
		if ( Rab != zero ) {
			// a muobs was multiplied to control the rotate direction
			if ( ( p->f1234[3] >= zero) && ( p->f1234[2] > zero) )
				*phi = p->muobs * *phi + asin( p->f1234[2] / Rab );

		        if( ( p->f1234[3] < zero ) && ( p->f1234[2] >= zero ) )
				*phi = p->muobs * *phi + pi - asin( p->f1234[2] / Rab );

		        if( ( p->f1234[3] <= zero) && ( p->f1234[2] < zero ) )
				*phi = p->muobs * *phi + pi - asin( p->f1234[2] / Rab );

		        if( ( p->f1234[3] > zero) && ( p->f1234[2] <= zero ) )
				*phi = p->muobs * *phi + twopi + asin( p->f1234[2] / Rab );

		} else
			*phi = zero;

		if ( rotate ) {
			*phi = fmod( *phi, twopi );
			if ( *phi < zero)
		           *phi += + twopi;
		}
	}
}














