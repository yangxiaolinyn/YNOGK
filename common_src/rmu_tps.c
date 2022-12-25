/*
 * common_src/rmu_tps.c
 *
 *     PURPOSE: th module aims on computing the four Boyer-Lindquist 
 *              coordinates (r,\theta,\phi,t) and the affine 
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
 * 2022-11-29   Starting the writting of the module.
 */
 

#include "rmu_tps.h"
 

/*
!*     PURPOSE: Returns the coordinates of turning points \mu_tp1 and \mu_tp2 of poloidal motion, judges
!*                whether the initial poloidal angle \theta_{obs} is one of turning points, if 
!*                it is true then mobseqmtp=.TRUE..  
!*     INPUTS:   f12342---------p_2, the \theta component of four momentum of the photon measured 
!*                              under the LNRF, see equation (84) in Yang & Wang (2012).
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*     OUTPUTS:  mu_tp1, mu_tp2----the turning points, between which the poloidal motion of 
!*                                 the photon was confined, and mu_tp2 <= mu_tp1. 
!*               reals------number of real roots of equation \Theta_\mu(\mu)=0.
!*               mobseqmtp---If mobseqmtp=.TRUE., then muobs equals to be one of the turning points. 
!*     ROUTINES CALLED: NONE
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS:
!*     C VERSION:  Yang Xiao-lin    2022-11-30.
!*/
int mutp( ptcl* this ) 
{
	this->mobseqmtp = false;
	if ( this->a_spin == zero ) {
		if( this->f1234[2] != zero ) {
			this->mu_tp1 = sqrt( this->q / 
				( this->lam2 + this->q ) );
			this->mu_tp2 = - this->mu_tp1;
		} else {
			this->mu_tp1 = fabs( this->muobs );
			this->mu_tp2 = - this->mu_tp1;
			this->mobseqmtp = true;
		}
		this->mu_reals = 2;    
	} else {
		if ( this->lambda != zero ) {
			double delta, sqrt_delta, tmp;
			delta = sq( this->a2 - this->lam2 - 
				this->q ) + four * this->a2 * this->q;

			sqrt_delta = sqrt(delta);

			this->mu_tp1 = sqrt( fabs( ( sqrt_delta - 
				( this->lam2 + this->q - this->a2) ) / two ) ) 
				/ fabs( this->a_spin );

			tmp = sqrt_delta + ( this->lam2 + this->q - this->a2 );
			if ( tmp <= zero ) {
				this->mu_tp2 = sqrt( - tmp / two ) / fabs( this->a_spin );
				if ( this->f1234[2] == zero ) {
					if ( fabs( this->muobs - this->mu_tp1 ) <= 1.e-4 )
						this->mu_tp1 = fabs( this->muobs );
					else       
						this->mu_tp2 = fabs( this->muobs );
					this->mobseqmtp = true;
				}
				this->mu_reals = 4;
			} else {
				if ( this->f1234[2] != zero ) {
					this->mu_tp2 = - this->mu_tp1;
				} else {
					this->mu_tp1 = fabs( this->muobs );
					this->mu_tp2 = - this->mu_tp1;
					this->mobseqmtp = true;
				}
				this->mu_reals = 2;
			}
		} else {
			if ( fabs( this->muobs ) != one ) {
				if ( this->q <= zero ) {
					if ( this->f1234[2] != zero ) {
				        	this->mu_tp2 = sqrt(- this->q )
						/ fabs( this->a_spin );
					} else {
						/* ! a = B = zero. */
						this->mu_tp2 = fabs( this->muobs );
						this->mobseqmtp = true;
					}
					this->mu_reals = 4;
				} else {
					this->mu_tp2 = - one;
					this->mu_reals = 2;
				}
				this->mu_tp1 = one;
			} else {
				this->mu_tp1 = one;
				if ( this->q <= zero && ( sq( this->f1234[2] ) + 
				sq( this->f1234[3] ) ) != zero ) {
					this->mu_tp2 = sqrt( - this->q )
						/ fabs( this->a_spin );
					this->mu_reals = 4;
				} else {
					this->mu_tp2 = - one;
					this->mu_reals = 2;
				}
			}
		}
	}

	if ( fabs( this->muobs ) == one ) this->mobseqmtp = true;
	if ( ( this->muobs < zero ) && ( this->mu_reals == 4 ) ) {
		double  mutemp;
		mutemp = this->mu_tp1;
		this->mu_tp1 = - this->mu_tp2;
		this->mu_tp2 = - mutemp;
	}

	this->mutp2 = pow( this->mu_tp1, 2 );
	this->mutp3 = pow( this->mu_tp1, 3 );
	return 0;
}  



/*
!*     PURPOSE: Returns the coordinates of turning points r_tp1 and r_tp2 of radial motion, judges
!*                whether the initial radius robs is one of turning points, if 
!*                it is true then robs_eq_rtp=.TRUE.. And if r_tp1 less or equal r_horizon,
!*                then indrhorizon=.TRUE. Where r_horizon is the radius of the event horizon.  
!*     INPUTS:   f12341---------p_r, the r component of four momentum of the photon measured 
!*                              under the LNRF, see equation (83) in Yang & Wang (2012).
!*               a_spin---------spin of black hole, on interval (-1,1).
!*               robs-----------radial coordinate of the observer. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*     OUTPUTS:  r_tp1, r_tp2----the turning points, between which the radial motion of 
!*                                 the photon was confined, and r_tp2 >= r_tp1.
!*               bb(1:4)----roots of equation R(r)=0.                
!*               reals------number of real roots of equation R(r)=0.
!*               robs_eq_rtp---If robs_eq_rtp=.TRUE., then robs equal to be one of turning points. 
!*               cases-------If r_tp2=infinity, then cases=1, else cases=2.
!*               indrhorizon----if r_tp1 less or equals r_horizon, indrhorizon=.TRUE.. 
!*     ROUTINES CALLED: root4
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: 
!*
!*     C VERSION:  Yang Xiao-lin    2022-11-30.
!*/
int radiustp( ptcl * this )
{
	double f12341, b1, c1, d1, e1;

	f12341 = this->f1234[1];
 
	this->robs_eq_rtp = false;
	this->indrhorizon = false;

	b1 = zero;
	c1 = this->a2 - this->lam2 - this->q;
	d1 = two * ( this->q + sq( this->lambda - this->a_spin ) );
	e1 = -this->q * this->a2;
	root4( b1, c1, d1, e1, this->rbb, &(this->r_reals) );

	/*printf("a = %f \n lambda = %f \n q = %f \n", this->a_spin, this->lambda, this->q);

	printf("111 = %d  robs = %f \n", this->r_reals, this->robs);
	printf("%20.16f + %20.16fi \n", creal(this->rbb[1]), cimag(this->rbb[1]));
	printf("%20.16f + %20.16fi \n", creal(this->rbb[2]), cimag(this->rbb[2]));
	printf("%20.16f + %20.16fi \n", creal(this->rbb[3]), cimag(this->rbb[3]));
	printf("%20.16f + %20.16fi \n", creal(this->rbb[4]), cimag(this->rbb[4]));
	printf("b1 = %f \n c1 = %f \n d1 = %f \n e1 = %f \n", b1, c1, d1, e1);*/

	if ( this->r_reals == 4 ) {
		if ( f12341 == zero ) {
			if ( fabs( this->robs - creal( this->rbb[4] ) ) <= 1.e-4 ) {
				this->r_tp1 = this->robs;
				this->r_tp2 = infinity;
				this->cases = 1;
			}
			if ( fabs( this->robs - creal( this->rbb[2]) ) <=  1.e-4) {
				this->r_tp1 = this->robs;  
				this->r_tp2 = creal( this->rbb[3] );
				this->cases = 2;
			}
			if ( fabs( this->robs - creal( this->rbb[3] ) ) <= 1.e-4 ) {
				this->r_tp1 = creal( this->rbb[2] );
				this->r_tp2 = this->robs;
				this->cases = 2;
			}
			if ( fabs( this->robs - creal( this->rbb[1] ) ) <= 1.e-4) {
				this->r_tp1 = - infinity;
				this->r_tp2 = this->robs;
				this->cases = 3;
				printf("radiustp(): wrong! 4 roots, cases = 3");
				exit(0);  
			}
			this->robs_eq_rtp = true;
		} else {
			if ( this->robs >= creal( this->rbb[4] ) ) {
				this->r_tp1 = creal( this->rbb[4] );
				this->r_tp2 = infinity;
				this->cases = 1;
				/*printf( "tttgg = %f \n %f \n %f \n ", this->robs, 
					creal( this->rbb[4] ), infinity );
				printf("%20.16f + %20.16fi \n", creal( this->rbb[1] ), cimag( this->rbb[1] ));
				printf("%20.16f + %20.16fi \n", creal( this->rbb[2] ), cimag( this->rbb[2] ));
				printf("%20.16f + %20.16fi \n", creal( this->rbb[3] ), cimag( this->rbb[3] ));
				printf("%20.16f + %20.16fi \n", creal( this->rbb[4] ), cimag( this->rbb[4] ));
				printf( "r_tp1 = %f \n %f \n ", this->r_tp1, 
					this->r_tp2 ); */

			} else if ( fabs( this->robs - creal( this->rbb[4] ) ) <= 1.e-10 ) {
				this->r_tp1 = creal( this->rbb[4] );
				this->r_tp2 = infinity;
				this->cases = 1;
				this->robs_eq_rtp = true;
				printf("radiustp(): p_r_LNRF === ', f12341, robs, rbb");
				exit(0);
			} else {
				if( this->robs >= creal( this->rbb[2] ) && 
				this->robs <= creal( this->rbb[3] ) ) {
					this->r_tp1 = creal( this->rbb[2] );
					this->r_tp2 = creal( this->rbb[3] );
					this->cases = 2;
				} else {
					if ( creal( this->rbb[1] ) > this->rhorizon 
					&& this->robs <= creal( this->rbb[1] ) ) {
						printf("radiustp(): wrong! 4 roots,cases = 3'");
						exit(0);
						this->r_tp2 = creal( this->rbb[1] );
						this->r_tp1 = - infinity;
					} else {
						printf("radiustp(): wrong! 4 roots',robs,rbb");
						exit(0);
					}
				}
			}
		}
	} else if ( this->r_reals == 2 ) {
		double r1[3];
		int j;
		j = 1;
		for ( int i = 1; i < 5; i++ ) {
			if ( cimag( this->rbb[i] ) == zero ) {
				r1[j] = creal( this->rbb[i] );
				j = j + 1;
			}
		}

		if ( f12341 == zero ) {
			if ( fabs( this->robs - r1[2] ) <= 1.e-4 ) {
				this->r_tp1 = this->robs;
				this->r_tp2 = infinity;
				this->cases = 1;
			}
			if ( fabs( this->robs - r1[1] ) <= 1.e-4 ) {
				this->r_tp1 = - infinity;
				this->r_tp2 = this->robs;
				printf(" radiustp(): wrong! 2 roots, cases = 3'");
				exit(0);
			}
			this->robs_eq_rtp = true;
		} else if ( fabs( this->robs - r1[2] ) <= 1.e-10 ) {
			this->r_tp1 = this->robs;
			this->r_tp2 = infinity;
			this->cases = 1;
		} else {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			if ( this->robs >= r1[2] ) {
				this->r_tp1 = r1[2];     
				this->r_tp2 = infinity;
				this->cases = 1;
			} else {
				if ( r1[1] >= this->rhorizon && this->robs <= r1[1] ) {
					printf("radiustp(): wrong! 2 roots, cases = 3");
					exit(0);
				}
			}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		}
	} else if ( this->r_reals == 0 ) {
		this->r_tp1 = zero;
		this->r_tp2 = infinity;
		this->cases = 1;
	} 
 
	//printf( "eee = %f \n %f \n %f \n r_reals = %d \n  = %f \n %d \n", this->rhorizon, 
	//	this->r_tp1, this->r_tp2, this->r_reals, this->f1234[1], this->cases );
	if ( this->rhorizon > this->r_tp1 && this->rhorizon < this->r_tp2 )
		this->indrhorizon = true;

	return 0;
}


 
/*
!*     PURPOSE:  Computes the value of parameter p from \mu coordinate. In other words, to compute 
!*               the \mu part of integral of equation (24), using formula (54) in Yang & Wang (2012).
!*               (54) is: p=-sign(p_\theta)*p_0+2*t1*p_2+2*t2*p_2. where p_\theta is initial \theta
!*               component of 4 momentum of photon.
!*     INPUTS:   f12342---------p_\theta, which is the \theta component of four momentum of a photon 
!*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
!*               f12343---------p_\phi, which is the \phi component of four momentum of a photon 
!*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               t1,t2----------Number of photon meets the turning points \mu_tp1 and \mu_tp2
!*                              respectively.
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.     
!*               mu-------------\mu coordinate of photon.     
!*     OUTPUTS:  value of \mu part of integral of (24). 
!*     ROUTINES CALLED: mu2p_schwartz, mutp, root3, weierstrass_int_J3
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: 
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-08.
!*/
double mu2p( ptcl * p, double mu, int t1, int t2 )
{
	double mu2p = zero, tposition, tp2, 
		b0, b1, b2, b3, g2, g3, tinf, 
		p1, p2, pp, a4, b4, integ4[5], rff_p = zero;

	int index_p4[5], cases;
	unsigned int del;
	double _Complex dd[4];

	if ( p->f1234[3] == zero && p->f1234[2] == zero && fabs(p->muobs) == one) {
		return zero;
	}
	if ( p->a_spin == zero ) {
		return (mu2p_Schwarzschild( p, mu, t1, t2 ));
	}

	a4 = zero;
	b4 = one;
	p->mobseqmtp = false;
	mutp( p );
	// equatorial plane motion.
	if ( p->mu_tp1 == zero ) {
		return zero;
	}

	// equations (26)-(29) in Yang & Wang (2012).

	b0 = - four * p->a2 * p->mutp3 + two * p->mu_tp1 * (p->a2 - p->lam2 - p->q);
	b1 = - two * p->a2 * p->mutp2 + one / three * (p->a2 - p->lam2 - p->q);
	b2 = - four / three * p->a2 * p->mu_tp1;
	b3 = - p->a2;
 
	g2 = three / four * ( sq(b1) - b0*b2 );
	g3 = one / 16.0 * ( three * b0 * b1 * b2 - two * sq3(b1) - sq(b0) * b3 );

	//equation (30) in Yang & Wang (2012).
	if ( fabs( mu - p->mu_tp1 ) != zero )
		tposition = b0 / ( four * ( mu - p->mu_tp1 ) ) + b1 / four;
	else
		tposition = infinity;

	if ( p->muobs != p->mu_tp1 ) 
		tinf = b0 / four / ( p->muobs - p->mu_tp1 ) + b1 / four;
	else
		tinf = infinity;  

	root3( zero, -g2/four, -g3/four, dd, &del );
	index_p4[1] = 0;
	cases = 1;

	if ( mu > p->mu_tp1 || mu < p->mu_tp2 ) {
		return (-one);
	}
     
	// equation (30) in Yang & Wang (2012).
	tp2 = b0 / four / ( p->mu_tp2 - p->mu_tp1 ) + b1 / four;
	if ( t1 == 0 ) 
		p1 = zero;
	else {
		// equation (53) in Yang & Wang (2012).
		weierstrass_int_J3( tposition, infinity, dd, del, a4, b4,
				index_p4, rff_p, integ4, cases );
		p1 = integ4[1];
	}
	if ( t2 == 0 )
		p2 = zero;
	else {
		weierstrass_int_J3( tp2, tposition, dd, del, a4, b4, 
				index_p4, rff_p, integ4, cases);
		p2 = integ4[1];
	}
	weierstrass_int_J3( tinf, tposition, dd, del, a4, b4, 
			index_p4, rff_p, integ4, cases);
	pp = integ4[1];

	//equation (54) in Yang & Wang (2012).
	if ( p->mobseqmtp ) {
		if ( p->muobs == p->mu_tp1 ) 
			mu2p = -pp + two*(t1*p1 + t2*p2);
		else
			mu2p = pp + two*(t1*p1 + t2*p2);  
	} else {
		if( p->f1234[2] < zero ) 
			mu2p = pp + two*(t1*p1 + t2*p2);
		if( p->f1234[2] > zero ) 
			mu2p = -pp + two*(t1*p1 + t2*p2);
	}
	return (mu2p);
}



/*
!*     PURPOSE:  Computes the value of parameter p from \mu coordinate. In other words, to compute 
!*               the \mu part of integral of equation (24), using formula (54) in Yang & Wang (2012).
!*               (54) is: p=-sign(p_\theta)*p_0+2*t1*p_2+2*t2*p_2. where p_\theta is initial \theta
!*               component of 4 momentum of photon.
!*               And black hole spin is zero. 
!*     INPUTS:   f12342---------p_\theta, which is the \theta component of four momentum of a photon 
!*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
!*               f12343---------p_\phi, which is the \phi component of four momentum of a photon 
!*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               t1,t2----------Number of photon meets the turning points \mu_tp1 and \mu_tp2
!*                              respectively.
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.     
!*               mu-------------\mu coordinate of photon.     
!*     OUTPUTS:  mu2p-----------value of \mu part of integral of (24). 
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: 
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-10.
!*
!*/
double mu2p_Schwarzschild( ptcl *p, double mu, int t1, int t2 )
{
	double  pp, p1, p2, BB, ant;

	if ( p->f1234[3] == zero && p->f1234[2] == zero ) {
		return (- two);
		/* this Theta=q(1-mu^2),so if B=0,then q=0.
		so Theta_mu=0 for ever.But we do not need to
		consider it,for q=0,so the next part means 
		that it will return zero value. */
	}

	p->mobseqmtp = false;
	if ( p->q > zero ) {
		BB = sqrt(p->q);
		if ( p->f1234[2] != zero ) {
			p->mu_tp1 = sqrt( p->q / ( p->lam2 + p->q ) );
			p->mu_tp2 = - p->mu_tp1;
		} else {
			p->mu_tp1 = p->muobs;
			p->mu_tp2 = - p->mu_tp1;
			p->mobseqmtp = true;
		}
		if ( fabs( p->muobs ) == one )
			p->mobseqmtp = true;

		pp = ( asin( mu / p->mu_tp1 ) - asin( p->muobs / p->mu_tp1 ) )
			* p->mu_tp1 / BB;
		if( t1 == 0 )
			p1 = zero;
		else      
			p1 = ( halfpi - asin( mu / p->mu_tp1 ) ) * p->mu_tp1 / BB;

		if ( t2 == 0 )
			p2 = zero;
		else      
			p2 = ( asin( mu / p->mu_tp1 ) + halfpi ) * p->mu_tp1 / BB;

		if ( p->mobseqmtp ) {
			if ( p->muobs == p->mu_tp1 )
				ant = -pp + two * (t1*p1 + t2*p2);
			else
				ant =  pp + two * (t1*p1 + t2*p2);           
		} else
			ant = sign( -p->f1234[2] )*pp + two * (t1*p1 + t2*p2);

	} else {
		ant = zero;
	} 
	return (ant);
}




















