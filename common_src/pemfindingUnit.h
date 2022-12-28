#ifndef PEMFINDINGUNIT_H
#define PEMFINDINGUNIT_H
 
#include "ynogkBL.h"
 

void pemfinds( ptcl *p, int *cases, double rins, double routs, double muups, 
	double mudowns, double phy1s, double phy2s, int caseranges, 
	double *parass, int bisections, double *pemfind, int NNs );
 
#endif
