#ifndef YNOGK_0ASPIN_H
#define YNOGK_0ASPIN_H

#include "BLcoordinates_new.h"


double Schwarz_integral( double y, double x, double AA );
void Get_phit_Integrals_Schwarzschild( ptcl *pt, double p, 
		double *phit_Schwarz, double *mucos, 
		double *sign_pth, int *t1, int *t2 );

 
#endif
