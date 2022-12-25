#ifndef ROOT34_H
#define ROOT34_H
 
#include <math.h>
#include "ConstantsUnits.h"
#include <stdlib.h>
#include <complex.h>
#include <stdio.h>

void sort( double a1, double a2, double a3, double *s);
int root3(double b, double c, double d, double _Complex *roots, unsigned int *del);
int root4(double b, double c, double d, double e, 
	double _Complex *r4, unsigned int *reals);

 
#endif
