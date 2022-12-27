/*
 * common_src/particle.h
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
 * 2022-11-22   Starting the writting of the module.
 */


/*typedef struct
{
	double somiga;
	double expnu;
	double exppsi;
	double expmu1;
	double expmu2;

	double expnu2;
	double exppsi2;
	double expmu12;
	double expmu22;
} metrics; */



typedef struct
{
	double somiga;
	double expnu;
	double exppsi;
	double expmu1;
	double expmu2;

	double expnu2;
	double exppsi2;
	double expmu12;
	double expmu22;
} metrics;


typedef struct
{
	double somiga;
	double expnu;
	double exppsi;
	double expmu1;
	double expmu2;

	double expnu2;
	double exppsi2;
	double expmu12;
	double expmu22;
} ini_parats;


#ifndef PARTICLE_H
#define PARTICLE_H
       
#include "ellCarlsons.h"
 

typedef struct
{
	double velocity_ini[4];
	double robs;  
	double r2;   
	double sinobs;  
	double sin2;    
	double muobs;  
	double mu2;  
	double a_spin;
	double a2;   
	double scal;

	double alpha;
	double beta;
	double LNRF_Pthe;
	double LNRF_Pphi;

	double LNRF_P[5];
	double f1234[5];

	metrics mt;
} bvs;  /* bvs = basic variables */



typedef struct
{
	double lambda;
	double lam2;  /* lbd2 = lambda * lambda */
	double q;

	double alphac;
	double betac;

	double LNRF_P[5];
	double f1234[5];

	metrics mt;
} out_data;  /* bvs = basic variables */





typedef struct {
	double r;
	double mucos;
	double phi;
	double t;
	double sigma;
 
	double alpha;
	double beta;
	double velocity_ini[4];  /* the velocity of the observer. */

	double LNRF_P[5];
	double f1234[5];
	double LNRF_Pthe;
	double LNRF_Pphi;
	double lambda;
	double lam2;  /* lam2 = lambda * lambda */
	double q;
	double robs;  /* robs = radius of observer or = radius of initial r */
	double r2;   /* r2 = robs * robs */
	double sinobs; /* sinobs = sin(theta_ini) or sinobs = sin(theta_obs) */
	double sin2;   /* sin2 = sinobs * sinobs */
	double muobs; /* muobs = cos(theta_ini) or muobs = cos(theta_obs) */
	double mu2;   /* mu2 = muobs * muobs */
	double a_spin;
	double a2;  /* a2 = a_spin * a_spin */
	double rhorizon;
	double scal;


	double alphac;
	double betac;

	metrics mt;

	/* variables for the rmu_tps.c module */
	double mu_tp1;
	double mutp2;    /* mutp2 = mu_tp1^2. */
	double mutp3;    /* mutp3 = mu_tp1^3. */
	double mu_tp2;
	unsigned int mu_reals;
	int mobseqmtp;

	double r_tp1;
	double r_tp2;
	unsigned int r_reals;
	int robs_eq_rtp;
	int indrhorizon;
	int cases;       /* If r_tp2=infinity, then cases=1, else cases=2. */
	double _Complex R_roots[5];
	//double _Complex bb[5];


	/* variables for the function mucos(p) */
	double sign_pth;   /* the sign of p_\theta of the four momentum. */
	double _Complex gdd[4];
	unsigned int gdel;
	double half_period_wp;
	double period_wp;
	double fzero;
	double AA, BB;
	double g2, g3, b0, b1;


	/* variables for the function radius(p) */
	double _Complex rbb[5], rdd[4];
	double radius, sign_pr;
	double periodwp, half_periodwp, p_ini_r_tp2;
	double rg2, rg3, rb0, rb1, rb2, rb3;
	unsigned int rdel;

	//double PI0, PI01, PI0_total, PI0_inf_obs, PI0_obs_hori, PI0_total_2;
	int cases_int, index_p4[5];
	double tinf, thorizon, tinf1, tp2, t_inf, rff_p, tp;
	int first_time_call;

	unsigned int Situations;

	double r_p, mu_p, phi_p, time_p, sigma_p, sign_pr, sign_pth;
} ptcl;


/*extern double a_spin;
extern double a2;
extern double robs;
extern double r2;
extern double muobs;
extern double mu2;
extern double sinobs;
extern double sin2;
extern double scal; */



//struct ptcl;



// Memory allocator
ptcl* particle_new();


void particle_construct( ptcl *this, double a_spin, double rini, 
		double mucos, double sinobs, double scal, double * );

void Set_alphabeta( ptcl *this, double alpha, double beta );
void Set_alpha( ptcl *this, double alpha );
void Set_beta( ptcl *this, double beta );


void lambdaq( ptcl *th );
void show_lambdaq( ptcl *th );

void ini_direction2lamdaq( ptcl *this, double pr, double ptheta, double pphi );
void metricg( ptcl *this );
int metricgij( double robs, double muobs, double sinobs, double a_spin, 
	double *somiga, double *expnu, double *exppsi, double *expmu1, double *expmu2 );
void center_of_image( ptcl *this );

void transport_data_out( ptcl *this, out_data * od );


double rms( double a_spin );

#endif
 




