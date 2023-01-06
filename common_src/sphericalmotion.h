#ifndef SPHERICALMOTION_H
#define SPHERICALMOTION_H
 
#include "ynogkBL.h"


void lambdaq_sphericalm( ptcl *pt, double r_sm, double *theta_min, 
	double *rmin, double *rmax );
void YNOGK_sphmotion( ptcl *pt, double p, double theta_max, double kp, 
	double kt, double *phyt, double *timet, double *sigmat, 
	double *mucos, double *sint, int *t1, int *t2 );
 
#endif
