
#include "disk.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*\
*  Functions' declarations defined in my code.
\*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

 
void disk( ptcl * p, double mudisk, double rdisk_out )
{ 
	double beta, alpha, betac, alphac;
	double deltax, deltay, rms1, pem, bomiga, ut_em;
	double somiga_em, expnu_em, exppsi_em, expmu1_em, expmu2_em;
	double r_em, mu_em, phi_em, time_em, sigma_em, sign_pr, sign_pth;
	double g;
	int m;

	metricg( p );
	rms1 = 1.30; //rms(a_spin)
	m = 800;
	// (111) in Yang & Wang (2012).
	//this->velocity_ini[1] = this->robs/(robs**(three/two)+a_spin)*zero;
	//this->velocity_ini[2] = this->robs/(robs**(three/two)+a_spin)*cos(45.D0*dtors)*zero;
	//this->velocity_ini[3] = this->robs/(robs**(three/two)+a_spin)*sin(45.D0*dtors)*zero;
	// (104) and (105) in Yang & Wang (2012).
	//center_of_image( robs,sinobs,muobs,a_spin,scal,velocity,alphac,betac );
	betac = zero;
	alphac = zero;
	//write(*,*) alphac, betac

	deltax = 80.0 / m;
	deltay = 30.0 / m;
	//fopen(unit=15,file='tdiskg.txt',status="replace");
	FILE *fp = fopen("./examples/ynogkc/plot/diskg4.txt", "w"); 


	//for ( int i = 330; i <= 330; i++ ) {
	for ( int i = 0; i <= m; i++ ) {
		beta = betac-i*deltay + 15.0;
		Set_beta( p, beta );
		//for ( int j = 448; j <= 449; j++ ) {
		for ( int j = 0; j <= m; j++ ) {
			alpha = alphac - j * deltax + 40.0;
			Set_alpha( p, alpha );
			lambdaq( p );
			// call pemdisk just gives the direct image.      
			//pem = Pemdisk( p, mudisk, rdisk_out, rms1, &r_em );
			pem = Pemdisk_all( p, mudisk, rdisk_out, rms1, &r_em );
			//pem = Pemdisk(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,mudisk,rdisk_out,rms1);
			// while call pemdisk_all will gives the direct and high-order images.            
			//pem = Pemdisk_all(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,mudisk,rdisk_out,rms1) 
			//printf(" j = %d \t \n lam = %f \n  q = %f \n", j, p->lambda, p->q );
			//printf("j = %d \t \n  pem = %f \n \%f \n %f \n %f \n", j, 
				//p->f1234[1], p->f1234[2], p->f1234[3], p->f1234[4]);
			//printf(" j = %d \t \n pem = %f \t %f \n", j, pem, alpha );

			//if (i == 330) {
				//fprintf( fp, "%20.16f \n", zero );

			//}

			if ( pem != -one && pem != -two ) {
				//printf("~~~~~~~~~~~~~~~2222222~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
				//radius_settings( p );
				r_em = radius( p, pem );

				//YNOGK( p, pem, &r_em, &mu_em, &phi_em, 
				//	&time_em, &sigma_em, &sign_pr, &sign_pth );
				//printf("~~~~~~~~~~~~~~~2222222~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
				//re = radius(pem,f1234(1),lambda,q,a_spin,robs,scal);
				//YNOGK(pem, f1234, lambda, q, sinobs, muobs, a_spin, robs, scal, 
				//	r_em, mu_em, phi_em, time_em, sigma_em);
				//printf("pem = %f \t r_em = %f \n", pem, r_em);

				//x_em = dsqrt( a_spin**2 + r_em**2 ) * dcos(phi_em);
				//y_em = dsqrt( a_spin**2 + r_em**2 ) * dsin(phi_em);
				 
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
 
 

 
