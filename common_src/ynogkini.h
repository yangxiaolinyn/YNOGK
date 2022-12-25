/*
 * common_src/particle.h
 *
 *     PURPOSE: This module aims on computing the four Boyer-Lindquist 
 *              coordinates (r, \mu = cos\theta, \phi, t) and the affine 
 *              parameter \sigam of the Kerr space-time.    
 *   
 *     AUTHOR:        Yang & Wang (2012)  
 *     DATE WRITTEN:  4 Jan 2012
 *
 *     C language Version:   Yang Xiao-lin  2022-11-22.
 *
 * Author       Yang Xiao-lin
 *
 * City         Kun-ming, Yunnan Provice, China.
 *
 * 2022-11-22   Starting the writting of the module.
 */


#ifndef YNOGKINI_H
#define YNOGKINI_H

#include "BLcoordinates_new.h"


/* functions for the int_r_part() */
//static void Get_t1_t2( double p, double f1234r, int robs_eq_rtp, double p_tp1_tp2, 
//		double robs, double r_tp1, double r_tp2, int *t1, int *t2 );
//static void Get_t1_t2_old( double p, double f1234r, int robs_eq_rtp, double pp, 
//		double p1, double p2, double robs, double r_tp1, int *t1, int *t2 );


/*void Set_Initializations_For_int_r_part( ptcl * p, double pm, 
		double *affr, double *timer, double *phir );*/
 
/*void Set_Initializations_For_Integrations_of_R_Part( ptcl * p );*/

/* void Get_R_Integral_Situations( ptcl *p, double pem );
void Integral_r_part( ptcl * p, double pm, double *affr, 
		double *timer, double *phir ); */



/* void Get_Results_of_Integrations_For_R_Part( ptcl * p, double pm, 
		double *r_coord, double *sign_pr, double *affr, 
		double *timer, double *phir ); */


void Get_Integrals_For_R_Part( ptcl *pt, double p, 
		double *r_coord, double *sign_pr, double *affr, 
		double *timer, double *phir );

#endif
 




