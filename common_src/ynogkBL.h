#ifndef YNOGKBL_H
#define YNOGKBL_H

#include "ynogk_theta_part.h"
#include "ynogkini.h"

//void ynogk( ptcl *this );
void YNOGK( ptcl *p, double pm, double *radi, double *mu, double *phi, 
		double *time, double *sigma, double *sign_pr, double *sign_pth );
 
#endif
