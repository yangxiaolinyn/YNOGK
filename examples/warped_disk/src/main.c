
//#include <mpi.h>
#include "ynogkBL.h"

void disk( ptcl * p, double mudisk, double rdisk_out );

typedef struct
{
	double n3;
	double rin;
	double rout; 
} wpdisk_paras;


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
double Fp( ptcl *pt, double p, double n3, double rin, 
	double paras4, double paras5, double paras6 )
/****************************************************************/
{
	double beta, gamma0;
       
	//n3 = paras(1); //!0.95D0      
	//rin = paras(2); //!4.23D0
	//rout = paras(3); //!50.D0  
	//YNOGK( p, pem, &r_em, &mu_em, &phi_em, 
	//		&time_em, &sigma_em, &sign_pr, &sign_pth );
	YNOGKC( ptcl *p, double p );

	// eq. (126) of Yang & Wang (2013).
	beta = n3 * sin( halfpi * ( pt->r_p - rin ) / ( rout - rin ) );
	// eq. (127) of Yang & Wang (2013).
	gamma0 = paras4 * dtors + paras5 * pi * exp( paras6 * 
			( rin - pt->r_p ) / ( rout - rin ) );
	//write(unit=6,fmt=*)pt->r_p,mua,phya,Fp
	// eq. (128) of Yang & Wang (2013).   
	return tan( beta ) * cos( pt->phi_p - gamma0 ) + 
		pt->mu_p / sqrt( one - pt->mu_p * pt->mu_p );
}


void disk( ptcl * p, double mudisk, double rdisk_out )
{ 
	double beta, alpha;
	double deltax, deltay, rms1, pem, bomiga, ut_em;
	double somiga_em, expnu_em, exppsi_em, expmu1_em, expmu2_em;
	double r_em; //, mu_em, phi_em, time_em, sigma_em, sign_pr, sign_pth;
	double g;
	int m;

	metricg( p );
	rms1 = rms( p->a_spin );  /* rms gives the radius of the 
				   * ISCO (inner most stable circular orbits) */
	m = 1200;
	// (111) in Yang & Wang (2012).
	//this->velocity_ini[1] = this->robs/(robs**(three/two)+a_spin)*zero;
	//this->velocity_ini[2] = this->robs/(robs**(three/two)+a_spin)*cos(45.D0*dtors)*zero;
	//this->velocity_ini[3] = this->robs/(robs**(three/two)+a_spin)*sin(45.D0*dtors)*zero;

	center_of_image( p );
	printf( " beta c = %f \t alphac = %f \n", p->betac, p->alphac );

	deltax = 80.0 / m;
	deltay = 30.0 / m;
	//fopen(unit=15,file='tdiskg.txt',status="replace");
	FILE *fp = fopen("./plot/diskg5.txt", "w"); 


	//for ( int i = 330; i <= 330; i++ ) {
	for ( int i = 0; i <= m; i++ ) {
		beta = p->betac-i*deltay + 15.0;
		Set_beta( p, beta );
		//for ( int j = 448; j <= 449; j++ ) {
		for ( int j = 0; j <= m; j++ ) {
			alpha = p->alphac - j * deltax + 40.0;
			Set_alpha( p, alpha );
			lambdaq( p );
			// call pemdisk just gives the direct image.      
			//pem = Pemdisk( p, mudisk, rdisk_out, rms1, &r_em );
			// while call pemdisk_all will gives the direct and high-order images.
			pem = Pemdisk_all( p, mudisk, rdisk_out, rms1, &r_em );

			if ( pem != -one && pem != -two ) {
				r_em = radius( p, pem );

				//YNOGK( p, pem, &r_em, &mu_em, &phi_em, 
				//	&time_em, &sigma_em, &sign_pr, &sign_pth ); 
				 
				if ( r_em >= rms1 ) {
					// Keplerian velocity of the particle.    
					bomiga = one / ( p->a_spin + pow(r_em, (three/two)) );
					metricgij(r_em, zero, one, p->a_spin, &somiga_em, 
						&expnu_em, &exppsi_em, &expmu1_em, &expmu2_em);
					// time component of four momentum of particle.
					ut_em = one / expnu_em / sqrt( one - 
						sq(exppsi_em / expnu_em * (bomiga-somiga_em) ) );
					// redshift formula of (113) of Yang and Wang (2012).
					g = ( one - p->mt.somiga * p->lambda ) / p->mt.expnu / 
						p->f1234[4] / ( one - bomiga * p->lambda) / ut_em;
				  
					//printf("g = %f \n", g);
					//printf("ut_em = %f \n", ut_em);
					fprintf( fp, "%20.16f \n", g );
				} else {
					fprintf( fp, "%20.16f \n", zero );
				}
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
	double a_spin, rini, muobs, sinobs, scal, Vobs[3];
	a_spin = 0.9980;
	rini = 4.e16;
	muobs = cos( 86.0 * dtor );
	sinobs = sin( 86.0 * dtor );
	scal = 1.0;

	Vobs[0] = zero;
	Vobs[1] = zero;
	Vobs[2] = zero;
	
	ptcl* pt;
	pt = particle_new();

	//photon_ctor( &pt, a_spin, rini, muobs, sinobs, scal );
	particle_construct( pt, a_spin, rini, muobs, sinobs, scal, Vobs );

	disk( pt, zero, 30.0 );

	free(pt);

	return 0;
}














