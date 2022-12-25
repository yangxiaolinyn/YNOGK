/*
 * common_src/rmu_tps.h
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
 * 2022-11-29   Starting the writting of the module.
 */


 


#ifndef RMU_TPS_H
#define RMU_TPS_H
 
#include "particle.h"


int mutp( ptcl* this );
int radiustp( ptcl * this );
 
double mu2p( ptcl * p, double mu, int t1, int t2 );
double mu2p_Schwarzschild( ptcl *p, double mu, int t1, int t2 );


#endif
 




