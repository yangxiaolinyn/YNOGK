
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

	int m = 40;
	int n = 400;
	double La, Lb;
	double dam, dan;
	double alpha, beta;

	La = 10.0;
	Lb = 10.0;
	dam = La * 2.0 / (m * 1.0);
	dan = Lb * 2.0 / (n * 1.0);

	FILE *fp = fopen("./pyplot/rayxx.txt", "w");

	double pem, mudisk, rdisk_out, rdisk_in;
	double r_em, mu_em, phi_em, time_em, sigma_em, sign_pr, sign_pth;

	rdisk_out = 15.0;
	rdisk_in = 1.40;
	mudisk = zero;

	for ( int i = 0; i < m + 1; i++ ) {
		beta = pt->betac - Lb + dam * i;
		Set_beta( pt, beta );
		for ( int j = 0; j < n + 1; j++ ) {
			alpha = pt->alphac - La + dan * j;
			Set_alpha( pt, alpha );
			
			lambdaq( pt );

			pem = Pemdisk( pt, mudisk, rdisk_out, rdisk_in , &r_em );
			//pem = Pemdisk_all( p, mudisk, rdisk_out, rms1, &r_em );
 
			/*printf(" pem = %f \t rem = %f  \n", pem, r_em);
			printf(" i = %d \t j = %d  \n", i, j);
			printf("*********************************************\n");*/

			//printf(" alphas = %f \t beta = %f  \n", alpha, beta );
			//printf(" betac = %f \t Lb = %f db = %f  i = %d \n", pt->betac, - Lb, db, i);


              		/*printf(" lam = %f   q = %f  a = %f \n ", pt->lambda, pt->q, pt->a_spin);
              		printf("f12341 = %f  f12342 = %f  f12343 = %f  f12344 = %f \n", 
				pt->f1234[1], pt->f1234[2], pt->f1234[3], pt->f1234[4]);*/

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

	FILE *fp1 = fopen("./pyplot/rayyy.txt", "w");

	for ( int i = 0; i < m + 1; i++ ) {
		alpha = pt->alphac - La + dam * i;
		Set_alpha( pt, alpha );
		for ( int j = 0; j < n + 1; j++ ) {
			beta = pt->betac - Lb + dan * j;
			Set_beta( pt, beta );
			lambdaq( pt );

			pem = Pemdisk( pt, mudisk, rdisk_out, rdisk_in , &r_em );
			//pem = Pemdisk_all( p, mudisk, rdisk_out, rms1, &r_em );
 
			/*printf(" pem = %f \t rem = %f  \n", pem, r_em);
			printf(" i = %d \t j = %d  \n", i, j);
			printf("*********************************************\n");*/

			//printf(" alphas = %f \t beta = %f  \n", alpha, beta );
			//printf(" betac = %f \t Lb = %f db = %f  i = %d \n", pt->betac, - Lb, db, i);


              		/*printf(" lam = %f   q = %f  a = %f \n ", pt->lambda, pt->q, pt->a_spin);
              		printf("f12341 = %f  f12342 = %f  f12343 = %f  f12344 = %f \n", 
				pt->f1234[1], pt->f1234[2], pt->f1234[3], pt->f1234[4]);*/

			if ( pem != -one && pem != -two ) {
 
				YNOGK( pt, pem, &r_em, &mu_em, &phi_em, 
					&time_em, &sigma_em, &sign_pr, &sign_pth );

				//printf(" x  = %20.16f \n y =  %20.16f \n", r_em, phi_em );
				fprintf( fp1, "%20.16f \t %20.16f \n", r_em * cos(phi_em), r_em * sin(phi_em) );

			} else {
				fprintf( fp1, "%20.16f \t %20.16f \n", sqrt(-1.0), sqrt(-1.0) );
			}
		}
	}
 	fclose(fp1);
 

	free(pt);

	return 0;
}














