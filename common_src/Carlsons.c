/*
 * common_src/Carlsons.c
 *
 *              This module includes supporting functions and subroutines of
 *              Carlson's integral method. Those codes mainly come from 
 *              Press (2007) Numerical Recipes.
 *
 * Author       Yang Xiao-lin
 *
 * Date         2022-11-16
 *
 * City         Kun-ming Yunnan Provice, China
 *
 * 2022-11-20   Starting the writting of the module.
 */

#include "root34.h"

#include <complex.h>
#include <math.h>
#include "ConstantsUnits.h"
#include <stdlib.h>
#include <Carlsons.h>
  


double rf(double x, double y, double z)
/* Computes Carlson’s elliptic integral of the ﬁrst kind, 
 * RF (x, y, z). x, y, and z must be nonneg-ative, and at 
 * most one can be zero. TINY must be at least 5 times the 
 * machine underﬂow limit, BIG at most one ﬁfth the machine 
 * overﬂow limit.
 */
{
	double alamb, ave, delx, dely, delz, e2, e3, sqrtx, sqrty, sqrtz, xt, yt, zt;

	double static ERRTOL = 0.0001;
	//double static TINY = 1.5e-38;
	//double static BIG = 3.0e37;
	double static THIRD = (1.0/3.0);
	double static C1 = (1.0/24.0);
	double static C2 = 0.1;
	double static C3 = (3.0/44.0);
	double static C4 = (1.0/14.0);


	/*if ( fmin(fmin(x,y),z) < 0.0 || fmin(fmin(x + y,x + z),y + z) < TINY ||
		fmax(fmax(x,y),z) > BIG) {
		printf("invalid arguments in rf\n");
		printf("x = %f \n y = %f \n z = %f \n ", x, y, z);
	}*/

	if (x < zero) x = zero;
	if (y < zero) y = zero;
	if (z < zero) z = zero;

	xt = x;
	yt = y;
	zt = z;
	do {
		sqrtx = sqrt(xt);
		sqrty = sqrt(yt);
		sqrtz = sqrt(zt);
		alamb = sqrtx * ( sqrty + sqrtz ) + sqrty * sqrtz;
		xt = 0.25*(xt + alamb);
		yt = 0.25*(yt + alamb);
		zt = 0.25*(zt + alamb);
		ave = THIRD*(xt + yt + zt);
		delx = (ave-xt)/ave;
		dely = (ave-yt)/ave;
		delz = (ave-zt)/ave;
	} while (fmax(fmax(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
	e2 = delx*dely-delz*delz;
	e3 = delx*dely*delz;
	return (1.0 + (C1*e2-C2-C3*e3)*e2 + C4*e3)/sqrt(ave);
}




/*
 * Computes Carlson’s elliptic integral of the second kind,
 * RD (x, y, z). x and y must be non-negative, and at most 
 * one can be zero. z must be positive. TINY must be at 
 * least twice the negative 2/3 power of the machine 
 * overﬂow limit. BIG must be at most 0.1 × ERRTOL times
 * the negative 2/3 power of the machine underﬂow limit.
 */
double rd(double x, double y, double z)
{
	double alamb, ave, delx, dely, delz, ea, eb, ec, ed, 
		ee, fac, sqrtx, sqrty, sqrtz, sum, xt, yt, zt;

	double static ERRTOL = 0.00015;
	//double static TINY = 1.0-25;
	//double static BIG = 4.5e21;
	double static C1 = (3.0/14.0);
	double static C2 = (1.0/6.0);
	double static const C3 = (9.0/22.0);
	double static const C4 = (3.0/26.0);
	double static C5 = (0.25*C3);
	double static C6 = (1.5*C4);

	/*if (fmin(x,y) < 0.0 || fmin(x + y,z) < TINY || fmax(fmax(x,y),z) > BIG)
		printf("invalid arguments in rd\n x = %f \n y = %f \n", x, y);*/

	if (x < zero) x = zero;
	if (y < zero) y = zero;

	xt = x;
	yt = y;
	zt = z;
	sum = 0.0;
	fac = 1.0;
	do {
		sqrtx = sqrt(xt);
		sqrty = sqrt(yt);
		sqrtz = sqrt(zt);
		alamb = sqrtx*(sqrty + sqrtz) + sqrty*sqrtz;
		sum +=  fac/(sqrtz*(zt + alamb));
		fac = 0.25*fac;
		xt = 0.25*(xt + alamb);
		yt = 0.25*(yt + alamb);
		zt = 0.25*(zt + alamb);
		ave = 0.2*(xt + yt + 3.0*zt);
		delx = (ave-xt) / ave;
		dely = (ave-yt) / ave;
		delz = (ave-zt) / ave;
	} while (fmax(fmax(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
	ea = delx*dely;
	eb = delz*delz;
	ec = ea-eb;
	ed = ea-6.0*eb;
	ee = ed + ec + ec;
	return 3.0*sum + fac*(1.0 + ed*(-C1 + C5*ed-C6*delz*ee)
		+ delz*(C2*ee + delz*(-C3*ec + delz*C4*ea))) / (ave*sqrt(ave));
}


 


/*
 * Computes Carlson’s elliptic integral of the third kind, 
 * RJ (x, y, z, p). x, y, and z must be nonnegative, and at 
 * most one can be zero. p must be nonzero. If p < 0, the 
 * Cauchy principal value is returned. TINY must be at least 
 * twice the cube root of the machine underﬂow limit,
 * BIG at most one ﬁfth the cube root of the machine overﬂow limit. 
*/
double rj(double x, double y, double z, double p)
{
	double alamb, alpha, ans, ave, a = zero, b = zero, beta, delp, 
		delx, dely, delz, ea, eb, ec, ed, ee, fac, 
		pt, rcx = zero, rho, sqrtx, sqrty, sqrtz, 
		sum, tau, xt, yt, zt;

	double static ERRTOL = 0.00015;
	//double static TINY = 2.5e-13;
	//double static BIG  = 9.0e11;
	double static C1 = (3.0/14.0);
	double static const C2 = (1.0/3.0);
	double static const C3 = (3.0/22.0);
	double static const C4 = (3.0/26.0);
	double static C5 = (0.75*C3);
	double static C6 = (1.5*C4);
	double static C7 = (0.5*C2);
	double static C8 = (C3+C3);

	/*if (fmin(fmin(x,y),z) < 0.0 || fmin(fmin(x+y,x+z),fmin(y+z,fabs(p))) < TINY
		|| fmax(fmax(x,y),fmax(z,fabs(p))) > BIG)
		printf("invalid arguments in rj\n x = %f \n y = %f \n z = %f \n", x, y, z);*/

	if (x < zero) x = zero;
	if (y < zero) y = zero;
	if (z < zero) z = zero;

	sum=0.0;
	fac=1.0;
	if (p > 0.0) {
		xt=x;
		yt=y;
		zt=z;
		pt=p;
	} else {
		xt=fmin(fmin(x,y),z);
		zt=fmax(fmax(x,y),z);
		yt=x+y+z-xt-zt;
		a=1.0/(yt-p);
		b=a*(zt-yt)*(yt-xt);
		pt=yt+b;
		rho=xt*zt/yt;
		tau=p*pt/yt;
		rcx=rc(rho,tau);
	}

	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		alpha=pow(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz, 2);
		beta=pt*pow(pt+alamb, 2);
		sum += fac*rc(alpha,beta);
		fac=0.25*fac;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		pt=0.25*(pt+alamb);
		ave=0.2*(xt+yt+zt+pt+pt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
		delp=(ave-pt)/ave;
	} while (fmax(fmax(fabs(delx),fabs(dely)), 
		fmax(fabs(delz),fabs(delp))) > ERRTOL);

	ea=delx*(dely+delz)+dely*delz;
	eb=delx*dely*delz;
	ec=delp*delp;
	ed=ea-3.0*ec;
	ee=eb+2.0*delp*(ea-ec);
	ans=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))
			+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave));
	if (p <= 0.0) ans=a*(b*ans+3.0*(rcx-rf(xt,yt,zt)));
	return ans;
}




/*
 * Computes Carlson’s degenerate elliptic integral, 
 * RC (x, y). x must be nonnegative and y must be nonzero. 
 * If y < 0, the Cauchy principal value is returned. 
 * TINY must be at least 5 times the machine underﬂow limit, 
 * BIG at most one ﬁfth the machine maximum overﬂow limit.
 */
double rc(double x, double y)
{
	double alamb, ave, s, w, xt, yt;

	static double ERRTOL = 0.00012;
	static double const TINY = 1.69e-38;
	static double const SQRTNY = 1.3e-19;
	static double const BIG = 3.e37;
	static double const TNBG = TINY*BIG;
	static double COMP1 = 2.236/SQRTNY;
	static double COMP2 = TNBG*TNBG/25.0;
	//static double THIRD = 1.0/3.0;
	static double C1 = 0.3;
	static double C2 = 1.0/7.0;
	static double C3 = 0.375;
	static double C4 = 9.0/22.0;

	long int i = 0;



	if (x < 0.0 || y == 0.0 || (x+fabs(y)) < TINY || (x+fabs(y)) > BIG ||
		(y<-COMP1 && x > 0.0 && x < COMP2)){
		//printf("invalid arguments in rc \n x = %f \n y = %f \n", x, y);
		if ( y == 0.0 )
			return infinity;
	}

	if (y > 0.0) {
		xt=x;
		yt=y;
		w=1.0;
	} else {
		xt = x-y;
		yt = -y;
		w = sqrt(x)/sqrt(xt);
	}

	i = 0;
	do {
		alamb=2.0*sqrt(xt)*sqrt(yt)+yt;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		ave=(xt+yt+yt)/3.0;
		s=(yt-ave)/ave;
		i += 1;
		/*printf( "i = %ld \n", i);
		printf( "xt = %f \n  yt = %f \n ", xt, yt );
		printf( "w = %f \n  s = %f \n  ave = %20.15f \n alamb = %f \n ", w, s, ave, alamb );*/
	} while ( fabs(s) > ERRTOL );
	//if ( y == 0.0 ) //exit(0);
		//printf( "i = %ld  ave = %f  s = %f rc = %f \n", i, ave, s, w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave));
	return w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave);
}









