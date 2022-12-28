#ifndef PEMFINDINGUNIT_H
#define PEMFINDINGUNIT_H
 
#include "ynogkBL.h"
 

void pemfinds( ptcl *p, int *cases, double rin, double rout, double muup, 
	double mudown, double phy1, double phy2, int caserange, double Fp, 
	double paras, double orir, double oricosth, int bisection, 
	double pemfind, int NN );
 
#endif
