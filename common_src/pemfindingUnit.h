#ifndef PEMFINDINGUNIT_H
#define PEMFINDINGUNIT_H
 
#include "ynogkBL.h"
 

void set_parameters( double rins, double routs, double muups, double mudowns,
	double phy1s, double phy2s, int bisections, int caseranges, int NNs );
void pemfinds( ptcl *p, double *pemfind );

 
#endif
