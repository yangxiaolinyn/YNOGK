
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
double Fp( ptcl *pt, double p, double n3, double rin, double rout,
	double paras4, double paras5, double paras6 )
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
	gamma0 = paras4 * dtors + paras5 * pi * exp( paras6 * 
			( rin - pt->r_p ) / ( rout - rin ) );
	//write(unit=6,fmt=*)pt->r_p,mua,phya,Fp
	// eq. (128) of Yang & Wang (2013).   
	return tan( beta ) * cos( pt->phi_p - gamma0 ) + 
		pt->mu_p / sqrt( one - pt->mu_p * pt->mu_p );
}


void warpeddisk( ptcl * p, double mudisk, double rdisk_out )
{ 
	double beta, alpha;
	double deltax, deltay, pem, bomiga, ut_em;
	double somiga_em, expnu_em, exppsi_em, expmu1_em, expmu2_em;
	double r_em; //, mu_em, phi_em, time_em, sigma_em, sign_pr, sign_pth;
	double g;
	int m, caserange, bisection;

	metricg( p );

	// (111) in Yang & Wang (2012).
	//velocity(1)=-0.D0!expmu1_obs/expnu_obs*robs/(robs**(three/two)+a_spin)!*zero
	//velocity(2)=0.D0!expmu2_obs/expnu_obs/(robs**(three/two)+a_spin)*cos(45.D0*dtor)!*zero
	//velocity(3)=0.D0!exppsi_obs/expnu_obs*(one/(robs**(three/two)+a_spin)*sin(45.D0*dtor)-somiga_obs)!*zero   

	center_of_image( p );
	printf( " beta c = %f \t alphac = %f \n", p->betac, p->alphac );

	m = 400;
	deltax = 110.0 / m;
	deltay = 110.0 / m;


	double rin, rout, gama0, n1, n2, n3;
	n1 = 4.0;
	n2 = 4.0;
	n3 = 0.950; 
	rin = rms( p->a_spin ); /* rms gives the radius of the 
				* ISCO (inner most stable circular orbits) */
	rout = 50.0;
	// caserange=4, means muup, mudown and phy1, phy2 are not provided.
	caserange = 4;
	// parameters to describe the curved surface of warped disk and sended to
	// function Fp as dummy variable.

	bisection = false;

	FILE *fp = fopen("./plot/warpdiskg.txt", "w"); 


	//for ( int i = 330; i <= 330; i++ ) {
	for ( int i = 0; i <= m; i++ ) {
		beta = p->betac-i*deltay + 55.0;
		Set_beta( p, beta );
		//for ( int j = 448; j <= 449; j++ ) {
		for ( int j = 0; j <= m; j++ ) {
			alpha = p->alphac - j * deltax + 55.0;
			Set_alpha( p, alpha );
			lambdaq( p );

			pem = Pemdisk_all( p, mudisk, rdisk_out, rms1, &r_em );
			//pem = pemfind(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&                              
                        //rin,rout,muup,mudown,phy1,phy2,caserange,Fp,paras,bisection)   

			if ( pem != -one && pem != -two ) {
				YNOGKC( pt, pem );
				sinp=sqrt(one-mua**two)  
				//write(unit=6,fmt=*)re,mua,phya
				beta=n3*sin(PI/two*(re-rin)/(rout-rin))
				gamma0=paras(4)*dtor+paras(5)*PI*exp(paras(6)*(rin-re)/(rout-rin)) 

				phy0=atan(tan(phya-gamma0)*cos(beta))
				sin_phy0=sign(sin(phy0),sin(phya-gamma0))
				cos_phy0=sign(cos(phy0),cos(phya-gamma0))
				V=re/(re**(three/two)+a_spin)
				V_theta=V*(-cos(phya-gamma0)*cos(beta)*sin_phy0*mua+sin(phya-gamma0)*cos_phy0*mua&
							-sin(beta)*sin_phy0*sinp)
				V_phi=V*(sin(phya-gamma0)*cos(beta)*sin_phy0+cos(phya-gamma0)*cos_phy0) 

				theta_dot=V_theta/re
				phi_dot=V_phi/re/sinp

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

	warpeddisk( pt, zero, 30.0 );

	free(pt);

	return 0;
}














