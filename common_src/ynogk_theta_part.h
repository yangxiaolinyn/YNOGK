#ifndef YNOGK_THETA_PART_H
#define YNOGK_THETA_PART_H

#include "ynogk_0aspin.h"
 

/*int Integration_Theta_part( ptcl *tp, double p, double *phit, double *timet, 
		double *mucos, double *sign_pth, int *t1, int *t2 );*/

  
//static void Get_mu_t1_t2_New( double p, double f12342, int mobseqmtp, double p_mt1_mt2, 
//		double muobs, double mu_tp1, double mu_tp2, int *t1, int *t2 );


//int Integration_Theta_part_Settings( ptcl *pt );

/*int Get_Integrations_of_Theta_part( ptcl *pt, double p, double *phit, double *timet, 
		double *mucos, double *sign_pth, int *tm1, int *tm2 );*/


void Get_Integrals_For_Theta_Part( ptcl *pt, double p, double *phit, double *timet, 
		double *mucos, double *sign_pth, int *tm1, int *tm2 );

#endif
