
//#include <mpi.h>
//#include <stdlib.h> 
//#include "ConstantsUnits.h"
 
//#include "disk.h"
#include "ynogkBL.h"

#define nn 500

 
//unsigned long long int icount = 0;
//unsigned long long mydutyphot, Total_Smple_Num;



int main( int argc, char *argv[])
{
	double a_spin, robs, muobs, sinobs, scal, Velocity_obs[3];
	double somiga, expnu, exppsi, expmu1, expmu2;
	a_spin = 0.93750;
	robs = rms( a_spin );
	muobs = cos( 90.0 * dtor );
	sinobs = sin( 90.0 * dtor );
	scal = 1.0;


	metricgij( robs, muobs, sinobs, a_spin, &somiga, &expnu, &exppsi, &expmu1, &expmu2 );

	Velocity_obs[0] = zero;
	Velocity_obs[1] = zero;
	Velocity_obs[2] = exppsi / expnu * ( one / ( pow(robs, (three/two)) + a_spin ) - somiga );
	
	
	ptcl* pt;
	pt = particle_new();

	//photon_ctor( &pt, a_spin, robs, muobs, sinobs, scal );
	particle_construct( pt, a_spin, robs, muobs, sinobs, scal, Velocity_obs );

	center_of_image( pt );
	printf( " beta c = %f \t alphac = %f \n", pt->betac, pt->alphac );

	int m = 20;
	int n = 20;
	int Num = 400;
	double delta_the, delta_phi, deltap, ptotal;
 
	delta_the = pi / (m * 1.0);
	delta_phi = twopi / (n * 1.0);
 

	FILE *fp = fopen("./pyplot/rays.txt", "w");

	double pem;
	double r_em, mu_em, phi_em, time_em, sigma_em, sign_pr, sign_pth;
	double theta, phi;
	double pphi, prad, pthe;

	for ( int i = 0; i < m + 1; i++ ) {
		theta = i * delta_the;
		for ( int j = 0; j < n + 1; j++ ) {
			phi = j * delta_phi;

			pphi = cos(theta);
			prad = sin(theta) * sin(phi);
			pthe = sin(theta) * cos(phi);

			ini_direction2lamdaq( pt, prad, pthe, pphi );
			//show_lambdaq( pt );
	
			ptotal = p_total( pt );
			deltap = ptotal / Num;
			for ( int k = 0; k <= Num; k++ ) { 
				pem = k * deltap;
				//printf("k = %d pem = %f \n", i, pem);
				YNOGK( pt, pem, &r_em, &mu_em, &phi_em, 
					&time_em, &sigma_em, &sign_pr, &sign_pth );
				//printf("**************************************\n");
 
				if ( r_em >= 1.10 * pt->rhorizon && r_em < 5000.0 ) {
					fprintf( fp, "%20.16e \t %20.16e \t %20.16e \n",
						sqrt( r_em*r_em + pt->a2 ) * sqrt( one - mu_em*mu_em ) * cos(phi_em),
						-sqrt( r_em*r_em + pt->a2 ) * sqrt( one - mu_em*mu_em ) * sin(phi_em),
						r_em * mu_em );
				} else {
		                        fprintf( fp, "%20.16f \t %20.16f \t %20.16f \n",
						sqrt( - 100.0), sqrt(- 100.0), sqrt(- 100.0) );
				}
			}
		}
	}
 	fclose(fp);
 
 

	free(pt);

	return 0;
}














