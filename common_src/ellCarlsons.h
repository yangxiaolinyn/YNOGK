/*
 * common_src/ellCarlsons.h
 *
 *              This module includes integrals and functions for Carlson's 
!*              integral method. Those codes mainly come from geokerr.f of
!*              Dexter & Agol (2009).
 *
 * Author       Yang Xiao-lin
 *
 * City         Kun-ming Yunnan Provice, China
 *
 * 2022-11-20   Starting the writting of the module.
 */


#ifndef ELLCARLSONS_H
#define ELLCARLSONS_H

#include "root34.h"
 
#include <complex.h>
#include <math.h>
#include "ConstantsUnits.h"
#include "stdlib.h"
#include "Carlsons.h"
#include "ellfunctions.h"


double sign( double y, double x );

/*int elldoublecomplexs( int *index_p5, double f1, double g1, double h1, 
		double f2, double g2, double h2, double a5, double b5,
		double y, double x, double rff_p, double *integ, int cases);*/

int ellcubiccomplexs(int *index_p4, double a1, double b1, double a4, 
		double b4, double f, double g, double h, double y, 
		double x, double rff_p, double *integ, int cases);

int ellcubicreals( int *index_p4, double a1, double b1, double a2, 
		double b2, double a3, double b3, double a4, double b4,
		double y, double x, double rff_p, double *integ, int cases);

int weierstrass_int_J3(double y, double x, double _Complex *bb, int del, 
			double a4, double b4, int *p4, 
			double rff_p, double *integ, int cases);


int carlson_doublecomplex5(double y, double x, double f1, double g1, 
			double h1, double f2, double g2, double h2,
			double a5, double b5,int *p5, double rff_p, 
			double *integ, int cases);


#endif

/*
!*     PURPOSE: Computes integrals: J_k(h)=\int^x_y (b4*t+a4)^(k/2)*(4*t^3-g_2*t-g_3)^(-1/2)dt.
!*              Where integer index k can be 0, -2, -4 and 2. (75) and (76) of Yang & Wang (2012).   
!*     INPUTS:  x,y -- limits of integral.
!*              bb(1:3) -- Roots of equation 4*t^3-g_2*t-g_3=0 solved by routine root3.
!*              del -- Number of real roots in bb(1:3).
!*              p4,rff_p,integ -- p4(1:4) is an array which specifies the value of index k of J_k(h). 
!*                 If p4(1)=0, then J_0 was computed and sent to integ(1).
!*                 If p4(1)!=0, then J_0 was replaced by parameter p, and rff_p=p.  
!*                 If p4(2)=-2, then J_{-2} was computed and sent to integ(2).                      
!*                 If p4(3)=2, then J_{2} was computed and sent to integ(3).
!*                 If p4(4)=-4, then J_{-4} was computed and sent to integ(4).
!*              cases -- If cases=1, then only J_0 was computed.
!*                       If cases=2, then only J_0 and J_{-2} are computed.
!*                       If cases=3, then only J_0, J_{-2} and J_{2} are computed.    
!*                       If cases=4, then J_0, J_{-2}, J_{2} and J_{-4} are computed.            
!*     OUTPUTS: integ -- is an array saved the results of J_k(h). 
!*     ROUTINES CALLED:  ellcubicreals, ellcubiccomplexs     
!*     ACCURACY:   Machine.
!*     REMARKS:
!*     AUTHOR:     Yang Xiao-lin
!*     DATE WRITTEN:  2022-11-20
!*
int weierstrass_int_J3(double y, double x, double _Complex *bb, int del, 
			double a4, double b4, int *p4, 
			double rff_p, double *integ, int cases)
 

!*     PURPOSE: Computes integrals: J_k(h)=\int^x_y (b5*r+a5)^(k/2)*
!*                                  [(h1*r^2+g1*r+f1)(h2*r^2+g2*r+f2)]^(-1/2)dr.
!*              Where integer index k can be 0, -2, -4, 2 and 4. 
!*                                  (77) of Yang & Wang (2012).   
!*     INPUTS:  x,y -- limits of integral.  
!*              p5,rff_p,integ -- p5(1:5) is an array which specifies 
!*                                the value of index k of J_k(h). 
!*                 If p5(1)=0, then J_0 was computed and sent to integ(1).
!*                 If p5(1)!=0, then J_0 was replaced by parameter p, and rff_p=p.  
!*                 If p5(2)=-2, then J_{-2} was computed and sent to integ(2).                      
!*                 If p5(3)=2, then J_{2} was computed and sent to integ(3).
!*                 If p5(4)=-4, then J_{-4} was computed and sent to integ(4).
!*                 If p5(4)=4, then J_{4} was computed and sent to integ(5).
!*              cases -- If cases=1, then only J_0 will be computed.
!*                       If cases=2, then only J_0 and J_{-2} will be computed.
!*                       If cases=3, then only J_0, J_{-2} and J_{2} will be computed.
!*                       If cases=4, then J_0, J_{-2}, J_{2} and J_{-4} will be computed. 
!*                       If cases=5, then J_0, J_{-2}, J_{2} and J_{4} will be computed.     
!*     OUTPUTS: integ -- is an array saved the results of J_k(h). 
!*     ROUTINES CALLED:  elldoublecomplexs
!*     ACCURACY:   Machine.
!*     AUTHOR:     Yang Xiao-lin
!*     DATE WRITTEN:  2022-11-21
!*
int carlson_doublecomplex5(double y, double x, double f1, double g1, 
			double h1, double f2, double g2, double h2,
			double a5, double b5,int *p5, double rff_p, 
			double *integ, int cases)    
 

!*     PURPOSE: Computes J_k(h)=\int_y^x dt (b4*t+a4)^(k/2)[(b1*t+a1)*(b2*t+a2)*(b3*t+a3)]^{-1/2}. 
!*              It is the case of equation W(t)=4*t^3-g_2*t-g_3=0 has three real roots. 
!*              Equation (61) of Yang & Wang (2012).   
!*     INPUTS:  Arguments for above integral. If index_p4(1)=0, then J_0 will be computed, 
!*              else J_0 will be replaced by parameter p, and rff_p=p. 
!*     OUTPUTS:  Value of integral J_k(h).
!*     ROUTINES CALLED: RF,RJ,RC,RD
!*     ACCURACY:   Machine.
!*     AUTHOR:     Dexter & Agol (2009)
!*     MODIFIED:   Yang Xiao-lin  2022-11-21
!*     DATE WRITTEN:  4 Mar 2009
!*     REVISIONS: 
!*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int ellcubicreals( int *index_p4, double a1, double b1, double a2, 
		double b2, double a3, double b3, double a4, double b4,
		double y, double x, double rff_p, double *integ, int cases)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 


!*     PURPOSE: Computes J_k(h)=\int_y^x dt (b4*t+a4)^(k/2)[(b1*t+a1)*(h*t^2+g*t+f)]^{-1/2}. 
!*              It is the case of equation W(t)=4*t^3-g_2*t-g_3=0 has one real root. 
!*              Equation (62) of Yang & Wang (2012).  
!*     INPUTS:  Arguments for above integral.
!*     OUTPUTS:  Value of integral. If index_p4(1)=0, then J_0 will be computed, 
!*               else J_0 will be replaced by parameter p, and rff_p=p. 
!*     ROUTINES CALLED: RF,RJ,RC,RD
!*     ACCURACY:   Machine.
!*     AUTHOR:     Dexter & Agol (2009)
!*     DATE WRITTEN:  4 Mar 2009
!*     MODIFIED:   Yang & Wang (2012)
!*     C version:  Yang Xiao-lin   2022-11-22
!*
int ellcubiccomplexs(int *index_p4, double a1, double b1, double a4, \
		double b4, double f, double g, double h, double y, \
		double x, double rff_p, double *integ, int cases)
 


!*     PURPOSE: Computes J_k(h)=\int_y^x dt (f_1+g_1t+h_1t^2)^{p_1/2} 
!*                       (f_2+g_2t+h_2t^2)^{p_2/2} (a_5+b_5t)^{p_5/2}. 
!*     INPUTS:  Arguments for above integral.
!*     OUTPUTS:  Value of integral. If index_p5(1)=0, then J_0 will be computed, 
!*               else J_0 will be replaced by parameter p, and rff_p=p. 
!*     ROUTINES CALLED: RF,RJ,RC,RD
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Dexter & Agol (2009)
!*     DATE WRITTEN:  4 Mar 2009
!*     MODIFIED:   Yang & Wang  (2012)
!*     C VERSION:   Yang Xiao-lin   2022-11-22
!*     REVISIONS: [-1,-1,-1,-1,0]=integ(1),[-1,-1,-1,-1,-2]=integ(2),
!*                [-1,-1,-1,-1,-4]=integ(4),[-1,-1,-1,-1,2]=integ(3),
!*                [-1,-1,-1,-1,4]=integ(5)
!*
int elldoublecomplexs( int *index_p5, double f1, double g1, double h1, 
		double f2, double g2, double h2, double a5, double b5,
		double y, double x, double rff_p, double *integ, int cases)
*/





