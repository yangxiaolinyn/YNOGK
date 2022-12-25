
//#include <mpi.h>
#include "ynogkBL.h"
 

int main( int argc, char *argv[])
{
	double a_spin, rini, theta_ini, muobs, sinobs, scal, Vobs[3];
	double r_em, mu_em, phi_em, time_em, sigma_em, sign_pr, sign_pth;
	double deltax, deltay;
	int m;



	a_spin = 0.9980;
	rini = 1.e6;
	theta_ini = 90.0;
	muobs = cos( theta_ini * dtor );
	sinobs = sin( theta_ini * dtor );
	scal = 1.0;

	/******* Set the velocities of the observer ********/
	// (111) in Yang & Wang (2012).
	//Vobs[0] = this->robs/(robs**(three/two)+a_spin)*zero;
	//Vobs[1] = this->robs/(robs**(three/two)+a_spin)*cos(45.D0*dtors)*zero;
	//Vobs[2] = this->robs/(robs**(three/two)+a_spin)*sin(45.D0*dtors)*zero;
	Vobs[0] = zero;
	Vobs[1] = zero;
	Vobs[2] = zero;
	
	ptcl* pt;
	pt = particle_new();

	particle_construct( pt, a_spin, rini, muobs, sinobs, scal, Vobs );

	m = 800;

	center_of_image( pt );
	printf( " beta c = %f \t alphac = %f \n", pt->betac, pt->alphac );

	deltax = 12.0 / m;
	deltay = 12.0 / m;

	FILE *fp = fopen("./plot/shadow.txt", "w");

	int t1 = 0;
	int t2 = 0;
	double rend = zero;
	double p_end = zero;
	double beta, alpha;

	for ( int i = 0; i <= m; i++ ) {
		beta = pt->betac-i*deltay + 6.0;
		Set_beta( pt, beta );
		for ( int j = 0; j <= m; j++ ) {
			alpha = pt->alphac - j * deltax + 4.0;
			Set_alpha( pt, alpha );
			lambdaq( pt );

			radiustp( pt );
			//printf( "alpha = %f beta = %f \n", alpha, beta );
			//printf( "pt->r_tp1 = %f pt->r_tp2 = %f \n", pt->r_tp1, pt->r_tp2 );
   
			if ( pt->r_tp1 <= pt->rhorizon ) {
				t1 = 0;
				t2 = 0;
				rend = pt->rhorizon;
				p_end = r2p( pt, rend, t1, t2 );
				YNOGK( pt, p_end, &r_em, &mu_em, &phi_em, 
					&time_em, &sigma_em, &sign_pr, &sign_pth );

				//printf( "p_end1 = %f timeem = %f sigma = %f \n", p_end, time_em, sigma_em );
				fprintf( fp, "%20.16f \n", sigma_em );
			} else {
				t1 = 1;
				t2 = 0;
				rend = pt->robs;
				p_end = r2p( pt, rend, t1, t2 );
				YNOGK( pt, p_end, &r_em, &mu_em, &phi_em, 
					&time_em, &sigma_em, &sign_pr, &sign_pth );

				//printf( " p_end2 = %f timeem = %f sigma = %f \n", p_end, time_em, sigma_em );
				//printf( "pem2 = %f sigma = %f \n", p_end, sigma_em );
				fprintf( fp, "%20.16f \n", sigma_em*half );
			}
		}
	}
	fclose(fp);

	free(pt);

	return 0;
}














