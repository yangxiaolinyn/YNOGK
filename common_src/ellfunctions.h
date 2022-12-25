/*
 * common_src/ellfunctions.h
 *  
 * Author       Yang Xiao-lin
 *
 * Date         2022-11-16
 *
 * City         Kun-ming Yunnan Provice, China
 *
 * 2022-11-20   Starting the writting of the module.
 */


#ifndef ELLFUNCTIONS_H
#define ELLFUNCTIONS_H
 
 

void sncndn(double uu, double emmc, double *sn, double *cn, double *dn);
int sncndn1(double uu, double emmc, double *sn, double *cn, double *dn);
double weierstrassP(double z, double g2, double g3, double _Complex *r1, int del);
double halfperiodwp(double g2, double g3, double _Complex *r1, int del);
double EllipticF(double t, double k2);
double EllipticE(double t, double k2);
double EllipticPI(double t, double n, double k2);

 
#endif


