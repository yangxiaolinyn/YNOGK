
//#include <mpi.h>
#include "ynogkBL.h"
#include "pemfindingUnit.h"

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
!*     C VERSION:  Yang Xiao-lin    2022-12-27.
!*
!*
!*
!*/
/****************************************************************/
double Fp( ptcl *pt, double p )
/****************************************************************/
{
	double beta, gamma0;
       
	//n3 = paras(1); //!0.95D0      
	//rin = paras(2); //!4.23D0
	//rout = paras(3); //!50.D0  
	//YNOGK( p, pem, &r_em, &mu_em, &phi_em, 
	//		&time_em, &sigma_em, &sign_pr, &sign_pth );
	YNOGKC( pt, p );

	// eq. (126) of Yang & Wang (2013).
	beta = n3 * sin( halfpi * ( pt->r_p - rin ) / ( rout - rin ) );
	// eq. (127) of Yang & Wang (2013).
	gamma0 = gama0 * dtors + n1 * pi * exp( n2 * 
			( rin - pt->r_p ) / ( rout - rin ) );
	//write(unit=6,fmt=*)pt->r_p,mua,phya,Fp
	// eq. (128) of Yang & Wang (2013).   
	return tan( beta ) * cos( pt->phi_p - gamma0 ) + 
		pt->mu_p / sqrt( one - pt->mu_p * pt->mu_p );
}


void warpeddisk( ptcl * pt, double mudisk, double rdisk_out )
{ 
	double beta, alpha;
	double deltax, deltay, pem, ut_em;
	double somiga_obs, expnu_obs, exppsi_obs, expmu1_obs, expmu2_obs;
	double somiga_em, expnu_em, exppsi_em, expmu1_em, expmu2_em;
	double g;
	int m, caserange, bisection;

	metricgij( pt->robs, pt->muobs, pt->sinobs, pt->a_spin,
		&somiga_obs, &expnu_obs, &exppsi_obs, &expmu1_obs, &expmu2_obs );

	// (111) in Yang & Wang (2012).
	//velocity(1)=-0.D0!expmu1_obs/expnu_obs*robs/(robs**(three/two)+a_spin)!*zero
	//velocity(2)=0.D0!expmu2_obs/expnu_obs/(robs**(three/two)+a_spin)*cos(45.D0*dtor)!*zero
	//velocity(3)=0.D0!exppsi_obs/expnu_obs*(one/(robs**(three/two)+a_spin)*sin(45.D0*dtor)-somiga_obs)!*zero   

	center_of_image( pt );
	printf( " beta c = %f \t alphac = %f \n", pt->betac, pt->alphac );

	m = 400;
	deltax = 110.0 / m;
	deltay = 110.0 / m;


	//double rin, rout, gama0, n1, n2, n3;
	n1 = 4.0;
	n2 = 4.0;
	n3 = 0.950; 
	rin = rms( pt->a_spin ); /* rms gives the radius of the 
				* ISCO (inner most stable circular orbits) */
	rout = 50.0;
	gama0 = -50.0;
	// caserange=4, means muup, mudown and phy1, phy2 are not provided.
	caserange = 4;
	// parameters to describe the curved surface of warped disk and sended to
	// function Fp as dummy variable.

	bisection = 0;

	FILE *fp = fopen("./plot/warpdiskg.txt", "w"); 

	set_parameters( rin, rout, one, -one, zero, twopi, bisection, caserange, 40 );
	

	double beta1, gamma0, V_theta, V_phi, theta_dot, phi_dot;
	double V, phy0, sin_phy0, cos_phy0, V_the_bar, V_phi_bar;
	double sq_Theta;
	//double dmu;

	for ( int i = 0; i <= m; i++ ) {
		beta = pt->betac-i*deltay + 55.0;
		Set_beta( pt, beta );
		printf("i = %d \t %f \n", i, beta);
		for ( int j = 0; j <= m; j++ ) {
			alpha = pt->alphac - j * deltax + 55.0;
			Set_alpha( pt, alpha );
			//printf("j = %d \t %f \n", j, alpha);
			lambdaq( pt );

			//pem = Pemdisk_all( p, mudisk, rdisk_out, rms1, &r_em );
			//pem = pemfind(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&                              
                        //rin,rout,muup,mudown,phy1,phy2,caserange,Fp,paras,bisection)
			//printf("pem1 = %f \n", pem);
			pemfinds( pt, &pem );
			//printf("pem2 = %30.25f \n", pem);

			if ( pem != -one && pem != -two ) {
				YNOGKC( pt, pem );

				beta1 = n3 * sin( halfpi * ( pt->r_p - rin ) / ( rout - rin ) );
				gamma0 = gama0 * dtor + n1 * pi * 
					exp( n2 * ( rin - pt->r_p ) / ( rout - rin ) );

				phy0 = atan( tan( pt->phi_p - gamma0 ) * cos( beta1 ) );
				sin_phy0 = sign( sin(phy0), sin( pt->phi_p - gamma0 ) );
				cos_phy0 = sign( cos(phy0), cos( pt->phi_p - gamma0 ) );
				V = pt->r_p / ( pow( pt->r_p, (three/two) ) + pt->a_spin );


				V_theta = V * ( - cos( pt->phi_p - gamma0 ) * cos( beta1 )
					* sin_phy0 * pt->mu_p + sin( pt->phi_p - gamma0 )
					* cos_phy0 * pt->mu_p - sin( beta1 ) * sin_phy0 * pt->sin_p );
 

				V_phi = V * ( sin( pt->phi_p - gamma0 ) * cos( beta1 )
					* sin_phy0 + cos( pt->phi_p - gamma0 ) * cos_phy0 );

				theta_dot = V_theta / pt->r_p;
				phi_dot = V_phi / pt->r_p / pt->sin_p;
 
				metricgij( pt->r_p, pt->mu_p, pt->sin_p, pt->a_spin,
					&somiga_em, &expnu_em, &exppsi_em, &expmu1_em, &expmu2_em );
      
				V_the_bar = theta_dot * expmu2_em / expnu_em;
				V_phi_bar = exppsi_em / expnu_em * ( phi_dot - somiga_em );
				ut_em = one / expnu_em / sqrt( one - sq(V_the_bar) - sq(V_phi_bar) );

				sq_Theta = sqrt( pt->q + pt->a2 * sq( pt->mu_p )
					- pt->lam2 * sq( pt->mu_p / pt->sin_p ) );
				// Eq. (110) of Yang & Wang (2012). 

				//dmu = mucos( pt, pem + 1.e-6) - mucos( pt, pem );

				/*g = one / expnu_obs * ( one - pt->lambda * somiga_obs )
					/ pt->f1234[4] / ut_em / ( one + 
					* sign( sq_Theta, - pt->sign_pth_p ) * V_theta / pt->r_p
					- pt->lambda * V_phi / pt->sin_p / pt->r_p );*/
 
				g = one / expnu_obs * ( one - pt->lambda * somiga_obs )
					/ pt->f1234[4] / ut_em / ( one + 
					sign( sq_Theta, pt->sign_pth_p ) * V_theta / pt->r_p
					- pt->lambda * V_phi / pt->sin_p / pt->r_p );
 
				fprintf( fp, "%20.16f \n", g );

			} else {

				if ( pem == -one )
					fprintf( fp, "%20.16f \n", zero );
				else
					fprintf( fp, "%20.16f \n", 0.10 );

			}
		}
	}
	fclose(fp);
} 
 


int main( int argc, char *argv[])
{
	double a_spin, rini, muobs, sinobs, scal, Vobs[3], theta_obs;
	a_spin = 0.9980;
	rini = 4.e10;
	theta_obs = 0.0;
	muobs = cos( theta_obs * dtor );
	sinobs = sin( theta_obs * dtor );
	scal = 1.0;

	Vobs[0] = zero;
	Vobs[1] = zero;
	Vobs[2] = zero;
	
	ptcl* pt;
	pt = particle_new();

	//photon_ctor( &pt, a_spin, rini, muobs, sinobs, scal );
	particle_construct( pt, a_spin, rini, muobs, sinobs, scal, Vobs );

	warpeddisk( pt, zero, 30.0 );

	free(pt);

	return 0;
}














