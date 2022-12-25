
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
	double a_spin, rini, muobs, sinobs, scal, Velocity_obs[3];
	a_spin = 0.950;
	rini = 4.e70;
	muobs = cos( 60.0 * dtor );
	sinobs = sin( 60.0 * dtor );
	scal = 1.0;

	Velocity_obs[0] = zero;
	Velocity_obs[1] = zero;
	Velocity_obs[2] = zero;
	
	
	ptcl* pt;
	pt = particle_new();

	//photon_ctor( &pt, a_spin, rini, muobs, sinobs, scal );
	particle_construct( pt, a_spin, rini, muobs, sinobs, scal, Velocity_obs );

	center_of_image( pt );
	printf( " beta c = %f \t alphac = %f \n", pt->betac, pt->alphac );
 
	int m = 70;
	int n = 70;
	double delta_the, delta_phi, deltap;
 
	delta_the = pi / (m * 1.0);
	delta_phi = twopi / (n * 1.0);

        deltap = 1.50 / 200.0;

	double pem;
	double r_em, mu_em, phi_em, time_em, sigma_em, sign_pr, sign_pth;
 
	double theta, phi;
	double pphi, prad, pthe;

	FILE *fp = fopen("./pyplot/rayxx.txt", "w");

	for ( int i = 0; i < m + 1; i++ ) {
		theta = i * delta_the;
		for ( int j = 0; j < n + 1; j++ ) {
			phi = j * delta_phi;

			pphi = cos(theta);
			prad = sin(theta) * sin(phi);
			pthe = sin(theta) * cos(phi);

			ini_direction2lamdaq( pt, prad, pthe, pphi );
			show_lambdaq( pt );

			//pem = Pemdisk( pt, mudisk, rdisk_out, rdisk_in , &r_em );
			//pem = Pemdisk_all( p, mudisk, rdisk_out, rms1, &r_em );

			if ( pem != -one && pem != -two ) {
 
				YNOGK( pt, pem, &r_em, &mu_em, &phi_em, 
					&time_em, &sigma_em, &sign_pr, &sign_pth );

				//printf(" x  = %20.16f \n y =  %20.16f \n", r_em, phi_em );
				fprintf( fp, "%20.16f \t %20.16f \n", r_em * cos(phi_em), r_em * sin(phi_em) );

			} else {
				fprintf( fp, "%20.16f \t %20.16f \n", sqrt(-1.0), sqrt(-1.0) );
			}
		}
	}
 	fclose(fp);
 
	free(pt);

	return 0;
}














