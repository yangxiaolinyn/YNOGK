/*
 * common_src/ellfunctions.c
 *
 *              This module includes supporting functions and subroutines to compute 
!*              Weierstrass' and Jacobi's elliptical integrals and functions by Carlson's 
!*              integral method. Those codes mainly come from Press (2007) and geokerr.f of
!*              Dexter & Agol (2009).
 *
 * Author       Yang Xiao-lin
 *
 * Date         2022-11-16
 *
 * City         Kun-ming Yunnan Provice, China
 *
 * 2022-11-20   Starting the writting of the module.
 */
 

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "ConstantsUnits.h"
#include "Carlsons.h"
#include "ellfunctions.h"
  

/*
!*    PURPOSE:   Returns the Jacobian elliptic functions sn(u|k^2), cn(u|k^2), 
!*               and dn(u|k^2). Here uu=u, while emmc=1-k^2. 
!*    RETURN:    sn, cn, dn----Jacobian elliptic functions sn(u|k^2), 
!*               cn(u|k^2), dn(u|k^2). 
!*    AUTHOR:    Press et al. (2007) 
*/ 
int sncndn1(double uu, double emmc, double *sn, double *cn, double *dn)
{    
	double const CA = 3.e-8;
	int l;
	double a, b, c, d, emc, u, em[13], en[13];
	int bo;

	emc = emmc;
	u = uu;
	if ( emc != 0.0 ) {
		bo = (emc < 0.0);
		if ( bo ) {
			d = 1.0-emc; //t'=t*k, u'=k*u, k'^2=1./k^2,  
			emc = - emc / d;
			d = sqrt(d);
			u = d * u;
		} 
		a = 1.0;
		*dn = 1.0;
		for ( int i = 0; i < 13; i++ ) { 
			l = i;
			em[i] = a;
			emc = sqrt(emc);
			en[i] = emc;
			c = 0.50 * ( a + emc );
			if ( fabs( a - emc ) <= CA * a ) break;
			emc = a * emc;
			a = c;
		}

		u = c*u;
		*sn = sin(u);  
		*cn = cos(u);               
		if ( *sn == 0.0 ) {
			if (bo) {
				a = *dn;
				*dn = *cn;
				*cn = a;
				*sn = *sn / d;
			}
			return 0;
		}

		a = (*cn) / (*sn);
		c = a * c;
		for ( int ii = l; ii > 0; ii-- ) {
			b = em[ii];
			a = c * a;
			c = *dn * c;
			*dn = ( en[ii] + a ) / ( b + a );
			a = c / b;
		}
		a = 1.0 / sqrt( c*c + 1.0 );
		if ( *sn < 0.0 ) *sn = -a;
		else *sn = a;
		*cn = c * (*sn);
		return 0;
	} else {
		*cn = one / cosh(u);
		*dn = *cn;
		*sn = tanh(u);
		return 0;
	}
}

 
#define CA 0.0000000000000000001
/* Special Functions The accuracy is the square of CA.*/

/*
 *     Returns the Jacobian elliptic functions sn(u, kc ), cn(u, kc ), 
 *             and dn(u, kc ). Here uu = u, while emmc = kc2 .
 *     AUTHOR:    Press et al. (2007) 
 */
void sncndn(double uu, double emmc, double *sn, double *cn, double *dn)
{
	double a,b,c,d,emc,u;
	double em[14],en[14];
	int i,ii,l,bo;

	emc=emmc;
	u=uu;
	if (emc) {
		bo=(emc < 0.0);
		if (bo) {
			d = 1.0 - emc;
			emc /= -1.0 / d;
			u *= (d=sqrt(d));
		}
		a = 1.0;
		*dn = 1.0;
		for ( i=1; i<=13; i++ ) {
			l = i;
			em[i] = a;
			en[i] = (emc = sqrt(emc));
			c=0.5*(a+emc);
			//printf("tt1 =%40.35f %40.35f %40.35f %40.35f %40.35f \n", emc, a, fabs( a - emc ), CA*a, c);
			if ( fabs( a - emc ) <= CA*a ) break;
			emc *= a;
			a = c;
		}

		u *= c;
		*sn = sin(u);
		*cn = cos(u);

		if ( *sn != zero ) {
			a = (*cn) / (*sn);
			c *= a;
			for ( ii=l; ii>=1; ii-- ) {
				b = em[ii];
				a *= c;
				c *= (*dn);
				*dn = (en[ii] + a) / (b + a);
				a = c/b;
			}

			a = 1.0 / sqrt( c*c + 1.0 );
			*sn = (*sn >= 0.0 ? a : -a);
			*cn = c*(*sn);
		} else {
			if (bo) {
				a = (*dn);
				*dn = (*cn);
				*cn = a;
				*sn /= d;
			}
		}
	} else {
		*cn = 1.0 / cosh(u);
		*dn = (*cn);
		*sn = tanh(u);
	}
}



/*
!*    PURPOSE:   Returns the Jacobian elliptic functions sn(u|k^2), cn(u|k^2), 
!*            and dn(u|k^2). Here uu=u, while emmc=1-k^2. 
!*    RETURN:    sn, cn, dn----Jacobian elliptic functions sn(u|k^2), cn(u|k^2), dn(u|k^2). 
!*    AUTHOR:    Press et al. (2007)
!*
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-24.
!*/
void sncndn_old(double uu, double emmc, double *sn, double *cn, double *dn)
{
	double a,b,c,d,emc,u;
	double em[14],en[14];
	int i, ii, l, bo;

	emc = emmc;
	u = uu;
	if ( emc != 0.0 ) {
		//bo = ( emc < 0.0);
		if ( emc < 0.0 )
			bo = true;
		else
			bo = false;

		if (bo) {
			d = 1.0 - emc;   /* t'=t*k, u'=k*u, k'^2=1./k^2, */
			emc = - emc / d;
			d = sqrt( d );
			u = d * u;
		}
		a = 1.0;
		*dn = 1.0;
		for ( i = 1; i < 14; i++ ) {
			l = i;
			em[i] = a;
			emc = sqrt( emc );
			en[i] = emc;
			c = 0.50 * ( a + emc );
			//printf("tt1 =%40.35f %40.35f %40.35f %40.35f %40.35f \n", emc, a, fabs( a - emc ), CA*a, c);
			if ( fabs( a - emc ) <= CA*a )
				break;
			emc = a * emc;
			a = c;
		}
		u = c * u;
		*sn = sin(u);
		*cn = cos(u);
		//printf("c = %40.35f \n u = %40.35f \n sn = %40.35f \n cn = %40.35f \n", c, u, *sn, *cn);
		if ( *sn == 0.0 ) {
			if ( bo ) {
				a = *dn;
				*dn = *cn;
				*cn = a;
				*sn = (*sn) / d;
			}
		} else {
			a = (*cn) / (*sn);
			c = a * c;
			for ( ii = l; ii >= 1; ii-- ) {
			      b = em[ii];
			      a = c * a;
			      c = (*dn) * c;
			      *dn = ( en[ii] + a ) / ( b + a );
			      a = c / b;
			}
			a = 1.0 / sqrt( c*c + 1.0 );
			if ( *sn < 0.0 )
				*sn = - a;
			else
				*sn = a;

			*cn = c * (*sn);
		}
	} else {
		  *cn = 1.0 / cosh(u);
		  *dn = *cn;
		  *sn = tanh(u);
	}
}



/*
!*    PURPOSE:   to compute Weierstrass' elliptical function \wp(z;g_2,g_3) and all of 
!*               this function involved are real numbers.   
!*    INPUTS:    z----------the independent variable value.
!*               g_2, g_3---two parameters.
!*               r1(1:3)----an array which is the roots of equation W(t)=4t^3-g_2t-g_3=0.
!*               del--------number of real roots among r1(1:3). 
!*    RETURN:    weierstrassP----the value of function \wp(z;g_2,g_3).  
!*    ROUTINES CALLED:  sncndn
!*    AUTHOR:     Yang, Xiao-lin & Wang, Jian-cheng (2012)
!*    DATE WRITTEN:  1 Jan 2012 */

double weierstrassP(double z, double g2, double g3, double _Complex *r1, int del) 
{
	double ant;
	double z1, e1, e2, e3, sn, u, k2, alp, bet, sig2, lamb, cn, dn;
        
	if ( z == zero ) {
		ant = infinity;
	} else {
		z1 = fabs(z); 
		if ( del == 3 ) {
			e1 = creal(r1[1]);
			e2 = creal(r1[2]);
			e3 = creal(r1[3]);
			k2 = (e2-e3) / (e1-e3);
			u = z1 * sqrt( e1 - e3 );
			sncndn( u, 1.0 - k2, &sn, &cn, &dn );
			ant = e3 + ( e1 - e3 ) / sn / sn;
			//printf(" Wp= %40.35f tp = %40.35f  pp = %40.35f df = %40.35f \n", e3, e1 - e3, sn, ant);
			//printf(" Wp2= %40.35f \n %40.35f \n %40.35f\n", z1, sqrt( e1 - e3 ), e1 - e3);
		} else {
			//alp=-real(r1(1))/two
			alp = creal(r1[2]);
			bet = fabs( cimag(r1[2]) );
			//sig=(9.D0*alp**2+bet**2)**(one/four)
			sig2 = sqrt( 9.0* pow( alp, 2 ) + pow( bet, 2 ) );
			lamb = ( sig2 - three * alp ) / bet;
			//k2=0.5D0+1.5D0*alp/sig**2
			k2 = one / ( one + pow( lamb, 2 ) );
			u = two * sqrt( sig2 ) * z1;
			sncndn( u, 1.0 - k2, &sn, &cn, &dn );
			if (cn > zero) {
				ant = two * sig2 * ( one + cn ) /
					pow( sn, 2 ) - two * alp - sig2;
			} else {
				ant = two * sig2 / ( one - cn ) -
						two * alp - sig2;
			}     
		}
	}

	return ant;
}

 


/*
!*    PURPOSE:   to compute the semi period of Weierstrass' elliptical 
!*               function \wp(z;g_2,g_3) and all of 
!*               this function involved are real numbers.   
!*    INPUTS:    g_2, g_3---two parameters. 
!*               r1(1:3)----an array which is the roots of equation W(t)=4t^3-g_2t-g_3=0.
!*               del--------number of real roots among r1(1:3). 
!*    RETURN:    halfperiodwp----the semi period of function \wp(z;g_2,g_3).  
!*    ROUTINES CALLED:  rf
!*    AUTHOR:    Yang, Xiao-lin
!*    DATE WRITTEN:   2022-11-20
*/
 
double halfperiodwp(double g2, double g3, double _Complex *r1, int del)
{
	double halfperiodwps, e1, e2, e3, EK, alp, bet, sig2, k2;

	if ( del == 3 ) {
		e1 = creal(r1[1]);
		e2 = creal(r1[2]);
		e3 = creal(r1[3]);
		k2 = ( e2 - e3 ) / ( e1 - e3 );
		EK = rf( zero, one - k2, one );
		halfperiodwps = EK / sqrt( e1 - e3 );
	} else {
		alp = creal(r1[2]);
		bet = fabs( cimag(r1[2]) );
		sig2 = sqrt( 9.0 * alp*alp + bet*bet );
		k2 = one / two + 1.5 * alp / sig2;
		EK = rf( zero, one - k2, one );
		halfperiodwps = EK / sqrt( sig2 );
	}
	return halfperiodwps;
}


/*
!*     PURPOSE: calculate Legendre's first kind elliptic integral: 
!*              F(t,k2)=\int_0^t dt/sqrt{(1-t^2)*(1-k2*t^2)}.  
!*     ARGUMENTS: t, k2
!*     ROUTINES CALLED:  RF  
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
 */
double EllipticF(double t, double k2)
{
	double x1, y1, z1;  
      
	x1 = one - t*t;
	y1 = one - k2*t*t;
	z1 = one;
	/*Press et al. 2007 (6.12.19) */
	return t * rf(x1,y1,z1);
}

/*
!*     PURPOSE: calculate Legendre's second kind elliptic integrals: 
!*              E(t,k2)=\int_0^t sqrt{1-k2*t^2}/sqrt{(1-t^2)}dt.
!*     ARGUMENTS: t, k2
!*     ROUTINES CALLED:  RF, RD 
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS: 
 */
double EllipticE(double t, double k2)
{
	double x1, y1, z1;
      
	x1=1.0 - t*t;
	y1=1.0 - k2*t*t;
	z1=1.0;
	/* Press et al. 2007 (6.12.20) */
	return t * rf(x1, y1, z1) - 1.0 / 3.0 * k2* pow(t, 3) * rd(x1, y1, z1);
}


/*
!*     PURPOSE: calculate Legendre's third kind elliptic integrals: 
!*              PI(t,n,k2)=\int_0^t /(1+nt^2)/sqrt{(1-k2*t^2)(1-t^2)}dt. 
!*     ARGUMENTS: t, k2
!*     ROUTINES CALLED:  RF, RJ
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
 */
double EllipticPI(double t, double n, double k2)
{
	double x1, y1, z1, w1, t2;

	t2 = pow(t, 2);
	x1 = 1.0 - t2;
	y1 = 1.0 - k2 * t2;
	z1 = 1.0;
	w1 = 1.0 + n * t2;
	/* Press et al. 2007 (6.12.20) */
	return t * rf(x1,y1,z1) - 1.0 / 3.0*n*pow(t, 3)*rj(x1,y1,z1,w1);
}



 
