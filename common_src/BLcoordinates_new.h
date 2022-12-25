/*
 * common_src/BLcoordinates.h
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
 
#ifndef BLCOORDINATE_NEW_H
#define BLCOORDINATE_NEW_H
 
#include "rmu_tps.h"

 
typedef struct {
	double u, v, w, L1, L2, m2;
	double u2, v2, w2;
	double cc, cr, dr, pinf; 
	double PI0, PI0_total, PI0_inf_obs;
	double PI0_obs_hori, PI01, PI0_total_2;
	double sqrt_L1;
} out_data2;


//void mucos_set( ptcl *this );
double mucos( ptcl *this, double p );

//double radius_preparation( ptcl *this, double p );
double radius_settings( ptcl *p );
double radius( ptcl *this, double p );
double r2p( ptcl *p, double rend, int t1, int t2 );
//double radius( ptcl *p, double pem, int first );
 
double Pemdisk( ptcl * p, double mu, double rout, double rin, double *re );
double Pemdisk_all( ptcl * p, double mu, double rout, double rin, double *re );

void Get_data( out_data2 * tmp );
double p_total( ptcl *p );

#endif
 




