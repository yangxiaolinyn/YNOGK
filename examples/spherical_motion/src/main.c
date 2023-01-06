 
#include "sphericalmotion.h"

void disk( ptcl * p, double mudisk, double rdisk_out );

typedef struct
{
	double n3;
	double rin;
	double rout; 
} wpdisk_paras;


double n1, n2, n3;
double gama0, rin, rout;

/*
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2023-01-03.
!*
!*
!*
!*/ 
 


int main( int argc, char *argv[])
{
	double a_spin, rini, muobs, sinobs, scal, Vobs[3], theta_obs;
	double deltap;
	double p;
	double xcc, ycc, zcc, theta_max;
	int m, n;
	a_spin = 0.980;
	rini = 1.403;
	theta_obs = 90.0;
	muobs = cos( theta_obs * dtor );
	sinobs = sin( theta_obs * dtor );
	scal = 1.0;

	Vobs[0] = zero;
	Vobs[1] = zero;
	Vobs[2] = zero;
	
	ptcl* pt;
	pt = particle_new();
 
	particle_construct( pt, a_spin, rini, muobs, sinobs, scal, Vobs );
 
	m = 20;
	n = 20;
 
	deltap = 0.006;

	if ( pt->a_spin == zero )
		theta_max = 40.0;

	double ra, ra2;
	FILE *fp = fopen("./plot/sphmotion.txt", "w");

	double phyt, timet, sigmat, mut, sin_theta;
	double kp, kt;
	double rmin, rmax;
	int t1, t2;

	kp = one;
	kt = one;

	for ( int i = m; i <= m; i++ ) {
		for ( int j = n; j <= n; j++ ) {
			ra = rini;
			ra2 = ra * ra;
			lambdaq_sphericalm( pt, rini, &theta_max, &rmin, &rmax );
			printf(" theta_max = %f  %f %f \n", pt->lambda, pt->q, theta_max );
			for ( int k = 0; k <= 20000; k++ ) {
				p = k * deltap;
				YNOGK_sphmotion( pt, p, theta_max, kp, 
					kt, &phyt, &timet, &sigmat, &mut, &sin_theta, &t1, &t2 );
					//printf(" r1 = %d  %f%f  %f  %f  %f  %f \n", k, ra, phyt, 
					//	timet, sigmat, mut, sin_theta); 
 
					xcc = sqrt( ra2 + pt->a2 ) * sin_theta * cos(phyt);
					ycc = sqrt( ra2 + pt->a2 ) * sin_theta * sin(phyt);
					zcc = sqrt( ra2 + pt->a2 ) * mut;
					fprintf( fp, "%20.16f %20.16f %20.16f \n", xcc, ycc, zcc );
			}
		}
	}
	fclose( fp );

	free(pt);

	return 0;
}














