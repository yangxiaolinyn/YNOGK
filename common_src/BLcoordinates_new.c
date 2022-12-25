/*
 * common_src/BLcoordinates.c
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
 * 2022-11-22   Starting the writting of the module.
 */

 

#include "BLcoordinates_new.h"
 
 
static double b0, b1, b2, b3, g2, g3;
static double u, v, w, L1, L2, m2;
static double u2, v2, w2;
static double sn, cn, dn;
static double cc, cr, dr, pinf; 
static double PI0, PI0_total, PI0_inf_obs;
static double PI0_obs_hori, PI01, PI0_total_2;
static double sqrt_L1;
//static double t_inf;
 

static double f12342, muobs, a_spin, q, lam2, a2;
static double b0, b1, b2, b3, g2, g3, tinf;
static double a4, b4, integ4[5], rff_p = zero;

static int index_p4[5];




static double a_spin_1;
static double scal_1;
static double lambda_1;
static double q_1;
static double robs_1;
static double f12341_1;
static double f12342_1;
static double f12343_1;
static double muobs_1;
static double sinobs_1;



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!*     PURPOSE:  Computes function \mu(p) defined by equation (32) in Yang & Wang (2012). That is
!*               \mu(p)=b0/(4*\wp(p+PI0;g_2,g_3)-b1)+\mu_tp1. \wp(p+PI0;g_2,g_3) is the Weierstrass'
!*               elliptic function.  
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f12342---------p_2, \theta component of four momentum of a photon measured under a LNRF.
!*               f12343---------p_3, \phi component of four momentum of a photon measured under a LNRF..
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  mucos----------\mu coordinate of photon corresponding to a given p. 
!*     ROUTINES CALLED: weierstrass_int_J3, mutp, root3.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS:
!*     C VERSION:  Yang Xiao-lin    2022-12-01.
!*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!*
!*
!*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!*/  
static void mucos_set( ptcl *this )
{
	f12342 = this->f1234[2];
	lam2 = this->lam2;
	q = this->q;
	muobs = this->muobs;
	a_spin = this->a_spin;
	a2 = this->a2;
              
	/* spin is zero. */
	if ( a_spin == zero ) {
		if ( q > zero ) {
			this->AA = sqrt( ( lam2 + q ) / q );
			this->BB = sqrt( q );
		}
	} else {
		int cases;

		mutp( this );

		a4 = zero;
		b4 = one;
		/* equations (26)-(29) in Yang & Wang (2012). */
		b0 = - four * a2 * this->mutp3 + two * this->mu_tp1 * (a2 - lam2 - q);
		b1 = - two * a2 * this->mutp2 + one / three * (a2 - lam2 - q);
		b2 = - four / three * a2 * this->mu_tp1;
		b3 = - a2;

		this->b0 = b0;
		this->b1 = b1;

		/* equation (31) in Yang & Wang (2012). */
		g2 = three / four * ( sq(b1) - b0*b2 );
		g3 = one / 16.0 * ( three * b0 * b1 * b2 - two * sq3(b1) - sq(b0) * b3 );
		this->g2 = g2;
		this->g3 = g3;

		root3( zero, -g2/four, -g3/four, this->gdd, &this->gdel );
		index_p4[1] = 0;
		cases = 1;

		tinf = b0 / ( four*( this->mu_tp2 - this->mu_tp1 ) ) + b1 / four;
		weierstrass_int_J3( tinf, infinity, this->gdd, this->gdel, a4, b4,
					index_p4, rff_p, integ4, cases );

		this->half_period_wp = integ4[1];
		//pp0 = halfperiodwp(g2, g3, dd, del) 
		this->period_wp = two * this->half_period_wp;

		/* equation (33) in Yang & Wang (2012). */
		if ( muobs != this->mu_tp1 ) {
			tinf = b0 / ( four*( muobs - this->mu_tp1 ) ) + b1 / four;
			weierstrass_int_J3( tinf, infinity, this->gdd, this->gdel, a4, b4, 
					index_p4, rff_p, integ4, cases);
			this->fzero = integ4[1];

		} else
                  	this->fzero=zero;

		if ( f12342 < zero )
			this->fzero = - this->fzero; 
	}
}


/*     C VERSION:  Yang Xiao-lin    2022-12-01.   */
double mucos_tmp( ptcl *this, double p )
{
	double mucos, f12343, f12342, muobs, a_spin, q; 
	double u;
 


	f12342 = this->f1234[2];
	f12343 = this->f1234[3];
	q = this->q;
	muobs = this->muobs;
	a_spin = this->a_spin; 

	if ( (f12343 == zero) && (f12342 == zero) && (fabs(muobs) == one) ) {
		/* this is because that mu==1 for ever,this because */
		/* that Theta_mu=-a^2(1-mu^2). So, mu must = +1 or -1 */
		/* for ever, and q=-a^2, X=lambda/sin(theta)=0 .*/
		/* So Theta_mu=q+a^2mu^2-X^2mu^4=-a^2(1-mu^2) */
		this->sign_pth = zero;   
		return muobs;
	}              
	/* spin is zero. */
	if ( a_spin == zero ) {
		if ( q > zero ) {
			if ( f12342 < zero ) {
				u = asin( muobs * this->AA ) + p * this->BB * this->AA;
				/* mucos=dsin(asin(muobs*this->AA) + 
				p*this->BB*this->AA)/this->AA */
				mucos = sin( u ) / this->AA;
				u = fmod( u + halfpi, twopi );
				if ( zero <= u && u <= pi )
					this->sign_pth = - one;
				else
					this->sign_pth = one;
			} else {
				if ( f12342 == zero ) {
					u = p * this->AA * this->BB;
					mucos = cos( u ) * muobs;
					u = fmod(u, twopi);
					if ( zero <= u && u <= pi )
						this->sign_pth = sign( muobs );
					else
						this->sign_pth = - sign( muobs );
				} else {
					u = asin( muobs * this->AA ) - p * this->BB * this->AA;
					/* mucos=dsin(dasin(muobs*this->AA) - 
					p*this->AA*this->BB)/this->AA */
					mucos = sin( u ) / this->AA;
					u = fmod(u - halfpi, twopi);
					if ( -pi <= u && u <= zero )
						this->sign_pth = one;
					else
						this->sign_pth = -one;     
				}
			}
		} else {
			this->sign_pth = zero;
			mucos = muobs;
		}
	} else {
		/* Equatorial plane motion. */
		if ( muobs == zero && q == zero ) {
			this->sign_pth = zero;
			mucos = zero;
			return (mucos);
		} 
 
		/* equation (32) in Yang & Wang (2012). */
		u = p + this->fzero;
		mucos = this->mu_tp1 + this->b0 / ( four * 
			weierstrassP( u, this->g2, this->g3, this->gdd, 
					this->gdel ) - this->b1 );

		/*printf("t1 = %f \n  t2 = %f \n, b0 =  %f \n, b1 =  %f \n ", 
			this->mu_tp1, this->mu_tp2,  this->b0, this->b1 );

		printf("Weis = ", weierstrassP( u, this->g2, this->g3, this->gdd, 
					this->gdel ));*/

		if ( u <= zero )
			this->sign_pth = -one;
		else {
			u = fmod(u, this->period_wp);
			if ( u <= this->half_period_wp )
				this->sign_pth = one;
			else
				this->sign_pth = - one;
		}

		if( mucos  <  this->mu_tp2) {
			//printf("mucos = %f \n  mu_tp2 = %f \n", mucos, this->mu_tp2);
			exit(0);
		}
	} 
	return (mucos);
}



/*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-19.
!*
!*/
double mucos( ptcl *pt, double p )
{

	static unsigned int count_num = 1;

mu_restarting:
	if ( count_num == 1 ) {

		mucos_set( pt );
		count_num += 1;

		a_spin_1 = pt->a_spin;
		scal_1 = pt->scal;
		lambda_1 = pt->lambda;
		q_1 = pt->q;
		f12341_1 = pt->f1234[1];
		f12342_1 = pt->f1234[2];
		f12343_1 = pt->f1234[3];
		muobs_1 = pt->muobs;
		sinobs_1 = pt->sinobs;
		robs_1 = pt->robs;

		return (mucos_tmp( pt, p ));

    	} else {
		if ( f12343_1 == pt->f1234[3]
		&& f12342_1 == pt->f1234[2]
		&& lambda_1 == pt->lambda
		&& q_1 == pt->q
		&& sinobs_1 == pt->sinobs
		&& muobs_1 == pt->muobs
		&& a_spin_1 == pt->a_spin
		&& scal_1 == pt->scal
		&& robs_1 == pt->robs ) {
 
			return (mucos_tmp( pt, p ));

		} else {
			count_num = 1;
			goto mu_restarting;
		}

	}
}






/*
!*     PURPOSE:  Computes function r(p) defined by equation (41) and (49) in Yang & Wang (2012). That is
!*               r(p)=b0/(4*\wp(p+PIr;g_2,g_3)-b1)+r_tp1. \wp(p+PIr;g_2,g_3) is the Weierstrass'
!*               elliptic function; Or r=r_+, r=r_-. 
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f1234r---------p_1, r components of four momentum of a photon measured under a LNRF. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.  
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  radius---------radial coordinate of photon corresponding to a given p.
!*               sign_pr--------the sign of r component of 4-momentum of the photon. 
!*     ROUTINES CALLED: weierstrass_int_J3, mutp, root3.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS:
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-03.
!*
!*/
double radius_preparation( ptcl *this, double p )
{
	double f1234r, lambda, q, a_spin, robs;
	double up, a2, lam2;
	double integ4[5], integ04[5], integ14[5];
	double cc, cr, dr, pinf, a4, b4, 
             t_inf, rff_p = zero;

	//double b0, b1, b2, b3, g2, g3;
	double tinf, thorizon, tinf1, tp2; 
	//double u, v, w, L1, L2, m2;


	int cases_int, index_p4[5];

	/* save :: f1234r_1,lambda_1,q_1,a_spin_1,robs_1,scal_1,r_tp1,r_tp2,reals,&
                robs_eq_rtp,indrhorizon,cases,bb,dd,del,cc,tinf,tp2,&
                thorizon,tinf1, t_inf,pinf,a4,b4
		, half_periodwp, periodwp, p_ini_r_tp2 */


	f1234r = this->f1234[1];
	lambda = this->lambda;
	lam2 = this->lam2;
	q = this->q;
	a_spin = this->a_spin;
	a2 = this->a2;
	robs = this->robs;

	a4 = zero;
	b4 = one;
	cc = a2 - lam2 - q;
	this->robs_eq_rtp = false;
	this->indrhorizon = false;

	radiustp( this );

	if (this->r_reals != 0 ) {
		/* equations (35)-(38) in Yang & Wang (2012). */
		b0 = four * sq3(this->r_tp1) + two * (a2 - lam2 - q ) * this->r_tp1 + 
			two * ( q + sq( lambda - a_spin ) );
		b1 = two * sq(this->r_tp1) + one / three * ( a2 - lam2 - q );
		b2 = four / three * this->r_tp1;
		b3 = one;
		g2 = three / four * ( sq(b1) - b0 * b2 );
		g3 = one / 16.0 * ( 3.0 * b0 * b1 * b2 - 2.0 * sq3(b1) - sq(b0) * b3 );

		this->rb0 = b0;
		this->rb1 = b1;
		this->rg2 = g2;
		this->rg3 = g3;
		/* equation (39) in Yang & Wang (2012). */
		if ( robs - this->r_tp1 != zero )
			tinf = b0 / four / ( robs - this->r_tp1 ) + b1 / four;
		else
			tinf = infinity;

		if ( this->rhorizon - this->r_tp1 != zero ) 
			thorizon = b1 / four + b0 / four / ( this->rhorizon - this->r_tp1 );
		else
			thorizon = infinity;

		tp2 = b0 / four / ( this->r_tp2 - this->r_tp1 ) + b1 / four;
		tinf1 = b1 / four;

		root3( zero, -g2/four, -g3/four, this->rdd, &this->rdel );

		index_p4[1] = 0;
		cases_int = 1;

		/* equation (42) in Yang & Wang (2012). */
		weierstrass_int_J3( tinf, infinity, this->rdd, this->rdel, zero, one, 
				index_p4, rff_p, integ04, cases_int );
		PI0 = integ04[1];

		if ( this->cases == 1 ) {  /* in this case, r_tp2 = infinity */
			if ( !this->indrhorizon ) {
				if ( f1234r  <  zero ) {
					weierstrass_int_J3( tinf1, infinity, 
						this->rdd, this->rdel, zero, one, 
						index_p4, rff_p, integ14, cases_int );

					PI0_total = PI0 + integ14[1];
					if ( p < PI0_total ) {
						/* equation (41) in Yang & Wang (2012).*/
						this->radius = this->r_tp1 + b0 / ( four * 
							weierstrassP( p - PI0, g2, g3, 
							this->rdd, this->rdel )- b1 );
						if ( p < PI0 )
							this->sign_pr = -one;
						else
							this->sign_pr = one;
					} else {
						this->radius = infinity;  /* Goto infinity, far away. */
						this->sign_pr = one;
					}
				} else {
					weierstrass_int_J3( tinf1, tinf, this->rdd, this->rdel, zero, one,
						index_p4, rff_p, integ04, cases_int );
					PI0_inf_obs = integ04[1];
					if ( p < PI0_inf_obs ) 
						/* equation (41) in Yang & Wang (2012). */
						this->radius=this->r_tp1+b0/(four * 
						weierstrassP(p + PI0, g2, g3, this->rdd, this->rdel)-b1);
				        else
						this->radius=infinity;  /*Goto infinity, far away.*/
				        this->sign_pr = one;
				}
			} else {
				if ( f1234r < zero ) {
					weierstrass_int_J3(tinf,thorizon, this->rdd, this->rdel, zero, one,
				                           index_p4,rff_p,integ04, cases_int);
				        PI0_obs_hori = integ04[1];
				        if ( p < PI0_obs_hori )
						/* equation (41) in Yang & Wang (2012). */
						this->radius = this->r_tp1+b0/(four*
						weierstrassP(p - PI0,g2,g3, this->rdd, this->rdel)-b1);
				        else
						this->radius = this->rhorizon;  /*Fall into black hole. */
				        this->sign_pr = -one;
				} else {
					weierstrass_int_J3(tinf1,tinf, this->rdd, this->rdel, zero, one, 
						index_p4,rff_p,integ04,cases_int);
					PI0_inf_obs = integ04[1];
				        if ( p < PI0_inf_obs )
					/* equation (41) in Yang & Wang (2012). */
						this->radius=this->r_tp1 + b0 /
						(four*weierstrassP(p + PI0,g2,g3, this->rdd, this->rdel)-b1);
				        else
						this->radius=infinity; /*Goto infinity, far away.*/
				        this->sign_pr = one;
				}
			}
		} else if ( this->cases == 2 ) {   //  case(2)
			if ( !this->indrhorizon ) {
				if ( f1234r < zero ) 
					PI01 = - PI0;
				else
					PI01 = PI0;
				/* equation (41) in Yang & Wang (2012). */
				this->radius = this->r_tp1 + b0 / ( four* weierstrassP( p + PI01, 
					g2, g3, this->rdd, this->rdel ) - b1 );
				/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				// this part of the code aims on obtaining the sign of the r component of  
				// 4-momentum of the photon.                                               
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
				weierstrass_int_J3(tp2,infinity, this->rdd, this->rdel,a4,b4,          /**/
							index_p4,rff_p,integ14,cases_int);             /**/
				this->half_periodwp = integ14[1];                                            /**/
				this->periodwp = two * this->half_periodwp;                                        /**/
				if (f1234r > zero) {                                                   /**/
					up = fmod(p + PI01, this->periodwp);                                    /**/
					if ( up <= this->half_periodwp )                                       /**/
						this->sign_pr = one;                                                /**/
					else                                                             /**/
						this->sign_pr = - one;                                              /**/
				} else if (f1234r < zero) {                                              /**/
					up = fmod(p + PI01 + this->half_periodwp, this->periodwp);              /**/
					if (up < this->half_periodwp)                                          /**/
						this->sign_pr = -one;                                               /**/
					else                                                             /**/
						this->sign_pr = one;                                               /**/
				} else {                                                                  /**/
					if (robs == this->r_tp1) {                                              /**/
						if ( fmod(p, this->periodwp) <= this->half_periodwp )                     /**/
							this->sign_pr = one;                                            /**/
						else                                                         /**/
							this->sign_pr = -one;                                    /**/
					} else {                                                             /**/ 
						if ( fmod(p + this->half_periodwp, this->periodwp)                       /**/
				                                   <= this->half_periodwp )                    /**/
							this->sign_pr = one;                           /**/
						else                                                         /**/
							this->sign_pr = -one;                            /**/ 
				        }                                                            /**/
				}                                                                /**/
			/*============================================================================
			!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			} else {
				if( f1234r <= zero ) {
					weierstrass_int_J3( tinf, thorizon, this->rdd, this->rdel, zero, one,
				                           index_p4, rff_p, integ14, cases_int);
					PI0_obs_hori = integ14[1]; 
					if ( p < PI0_obs_hori )
					/* equation (41) in Yang & Wang (2012). */
						this->radius = this->r_tp1 + b0 / ( four * 
						weierstrassP( p - PI0, g2, g3, this->rdd, this->rdel) - b1 );
					else
						this->radius = this->rhorizon; //Fall into black hole
				        this->sign_pr = -one;
				} else {
					weierstrass_int_J3( tp2, thorizon, this->rdd, this->rdel, a4, b4, index_p4, 
						rff_p,integ14,cases_int);
					weierstrass_int_J3( tp2, tinf, this->rdd, this->rdel, a4, b4, 
						index_p4, rff_p, integ4, cases_int);
					this->p_ini_r_tp2 = integ4[1];
					PI0_total_2 = integ14[1] + integ4[1];
					if ( p < PI0_total_2 )
						// equation (41) in Yang & Wang (2012). 
						this->radius = this->r_tp1 + b0 / ( four * 
						weierstrassP( p + PI0, g2, g3, this->rdd, this->rdel ) - b1 );
				        else
						this->radius = this->rhorizon; //Fall into black hole. 
					//!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					if (p <= this->p_ini_r_tp2)                  /**/
						this->sign_pr = one;                       /**/
				        else                                    /**/
						this->sign_pr = -one;                      /**/                                /**/
					//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
				}
			}
		} // end casese.
		if ( a_spin == zero ) {
			if ( cc == zero ) {
				if ( f1234r < zero ) {
					if ( p < one / this->rhorizon - one / this->robs )
						this->radius = this->robs / ( this->robs * p + one );
					else
						this->radius = this->rhorizon;
					this->sign_pr = -one;
				} else {
					if ( p < one/robs )
						this->radius = this->robs / ( one - this->robs * p );
					else
						this->radius = infinity;
				        this->sign_pr = one;
				}
			}
			if ( cc == -27.0 ) {
				if (f1234r < zero) {
					cr = - three*fabs((sqrt(robs*(robs+6.0))+(three+two*robs)/sqrt3) 
							/(three-robs))*exp(three*sqrt3*p)-sqrt3;
					dr = - fabs((sqrt(robs*(robs+6.0))+(three+two*robs)/sqrt3)/ 
							(three-robs))*exp(three*sqrt3*p)+two/sqrt3;
					if ( p != zero )
						this->radius = (three+cr*dr+sqrt(9.0+6.0*cr*dr+ cr*cr))/(dr*dr-one);
					else
						this->radius = this->robs; //infinity;

				        this->sign_pr = - one;
				} else {
				        cr = -three*fabs((sqrt(robs*(robs+6.0))+(three+two*robs)/sqrt3)/ 
				                      (three-robs))*exp(-three*sqrt3*p)-sqrt3;
				        dr = -fabs((sqrt(robs*(robs+6.0))+(three+two*robs)/sqrt3)/ 
				                      (three-robs))*exp(-three*sqrt3*p)+two/sqrt3;
				        PI0 = log(fabs((sqrt(robs*(robs+6.0))+(three+two*robs)/sqrt3)/ 
				                      (robs-three)))/three/sqrt3-log(one+two/sqrt3)/three/sqrt3;
				        if ( p < PI0 )
				            this->radius = (three+cr*dr+sqrt(9.0+6.0*cr*dr+cr*cr))/(dr*dr-one);
				        else
				            this->radius = infinity;
				        this->sign_pr = one;
				}
			}
		}
	} else {
            u=creal(this->rbb[4]);
            w=fabs( cimag(this->rbb[4]));
            v=fabs( cimag(this->rbb[2]));
            if (u != zero) {
// equation (45) in Yang & Wang (2012). 
                L1=(four*u*u+w*w+v*v+sqrt(sq(four*u*u+w*w+v*v)-four*w*w*v*v))/(two*w*w);
                L2=(four*u*u+w*w+v*v-sqrt(sq(four*u*u+w*w+v*v)-four*w*w*v*v))/(two*w*w);
// equation (46) in Yang & Wang (2012). 
                thorizon=sqrt((L1-one)/(L1-L2))*(this->rhorizon-u*(L1+one)/
                                  (L1-one))/sqrt(sq(this->rhorizon-u)+w*w);
// equation (48) in Yang & Wang (2012). 
                m2=(L1-L2)/L1;
                tinf=sqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/
                          (L1-one))/sqrt(sq(robs-u) + w*w);
                t_inf=sqrt((L1-one)/(L1-L2));
// equation (50) in Yang & Wang (2012). 
                pinf=EllipticF(tinf,m2)/w/sqrt(L1);
                sncndn(p*w*sqrt(L1) + sign( f1234r )*pinf*w*sqrt(L1),one-m2, &sn, &cn, &dn);
                if(f1234r < zero) {
                    PI0=pinf-EllipticF(thorizon,m2)/(w*sqrt(L1));
                    if (p < PI0) 
// equation (49) in Yang & Wang (2012), and p_r <0, r=r_{+}
                        this->radius=u+(-two*u+w*(L1-L2)*sn*fabs(cn))/((L1-L2)*sn*sn-(L1-one));
                    else
                        this->radius=this->rhorizon; 
                    this->sign_pr = - one;
                } else {
                    PI0 = EllipticF(t_inf,m2)/(w*sqrt(L1))-pinf;
                    if ( p < PI0 ) 
// equation (49) in Yang & Wang (2012), and p_r >0, r=r_{-}
                        this->radius=u+(-two*u-w*(L1-L2)*sn*fabs(cn))/((L1-L2)*sn*sn-(L1-one));
                    else
                        this->radius=infinity; 
                    this->sign_pr = one;
		}
            } else {
		        if (f1234r < zero) {
		            if ( p < ( tan(robs/w)- tan(this->rhorizon/w))/w ) 
		                this->radius=w*tan( tan(robs/w)-p*w);
		            else
		                this->radius=this->rhorizon; 
		            this->sign_pr = - one;
		        } else {
		            if(p < (pi/two- tan(robs/w))/w) 
		                this->radius=w*tan( tan(robs/w)+p*w);
		            else
		                this->radius=infinity;
		            this->sign_pr = one;
			}
		}
	}

	return this->radius;
}





/*
!*     PURPOSE:  Computes function r(p) defined by equation (41) and (49) in Yang & Wang (2012). That is
!*               r(p)=b0/(4*\wp(p+PIr;g_2,g_3)-b1)+r_tp1. \wp(p+PIr;g_2,g_3) is the Weierstrass'
!*               elliptic function; Or r=r_+, r=r_-. 
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f1234r---------p_1, r components of four momentum of a photon measured under a LNRF. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.  
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  radius---------radial coordinate of photon corresponding to a given p.
!*               sign_pr--------the sign of r component of 4-momentum of the photon. 
!*     ROUTINES CALLED: weierstrass_int_J3, mutp, root3.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS:
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-03.
!*
!*/
double radius_settings( ptcl *p )
{
	double integ4[5], integ04[5], integ14[5];

	radiustp( p );

	if (p->r_reals != 0 ) {
		/* equations (35)-(38) in Yang & Wang (2012). */
		p->rb0 = four * sq3(p->r_tp1) + two * (p->a2 - p->lam2 - p->q ) * p->r_tp1 + 
			two * ( p->q + sq( p->lambda - p->a_spin ) );
		p->rb1 = two * sq(p->r_tp1) + one / three * ( p->a2 - p->lam2 - p->q );
		p->rb2 = four / three * p->r_tp1;
		p->rb3 = one;
		p->rg2 = three / four * ( sq(p->rb1) - p->rb0 * p->rb2 );
		p->rg3 = one / 16.0 * ( 3.0 * p->rb0 * p->rb1 * p->rb2 - 2.0 * sq3(p->rb1) - sq(p->rb0) * p->rb3 );

		/* equation (39) in Yang & Wang (2012). */
		if ( p->robs - p->r_tp1 != zero )
			p->tinf = p->rb0 / four / ( p->robs - p->r_tp1 ) + p->rb1 / four;
		else
			p->tinf = infinity;

		if ( p->rhorizon - p->r_tp1 != zero ) 
			p->thorizon = p->rb1 / four + p->rb0 / four / ( p->rhorizon - p->r_tp1 );
		else
			p->thorizon = infinity;

		p->tp2 = p->rb0 / four / ( p->r_tp2 - p->r_tp1 ) + p->rb1 / four;
		p->tinf1 = p->rb1 / four;

		root3( zero, - (p->rg2) / four, - (p->rg3) / four, p->rdd, &p->rdel );

		p->index_p4[1] = 0;
		p->cases_int = 1;

		/* equation (42) in Yang & Wang (2012). */
		weierstrass_int_J3( p->tinf, infinity, p->rdd, p->rdel, zero, one, 
				p->index_p4, p->rff_p, integ04, p->cases_int );
		PI0 = integ04[1];
		//p->PI0 = PI0;
		//printf(" ss = %f \t %f \t %f \t \n  ", p->rg2, p->rg3, PI0 );

		if ( p->cases == 1 ) {
			if ( !p->indrhorizon ) {
				if ( p->f1234[1]  <  zero ) {
					weierstrass_int_J3( p->tinf1, infinity, 
						p->rdd, p->rdel, zero, one, 
						p->index_p4, p->rff_p, integ14, p->cases_int );

					PI0_total = PI0 + integ14[1];
				} else {
					weierstrass_int_J3( p->tinf1, p->tinf, p->rdd, p->rdel, zero, one,
						p->index_p4, p->rff_p, integ04, p->cases_int );

					PI0_inf_obs = integ04[1];
				}
			} else {
				//printf(" pr = %f \t \n  ", p->f1234[1] );
				if ( p->f1234[1] < zero ) {
					weierstrass_int_J3( p->tinf, p->thorizon, p->rdd, p->rdel, zero, one,
				                           p->index_p4, p->rff_p, integ04, p->cases_int);

				        PI0_obs_hori = integ04[1];
				} else {
					weierstrass_int_J3( p->tinf1, p->tinf, p->rdd, p->rdel, zero, one, 
						p->index_p4, p->rff_p, integ04, p->cases_int);

					PI0_inf_obs = integ04[1];
				}
			}
		} else if ( p->cases == 2 ) { // case(2)
			if ( !p->indrhorizon ) {
				if ( p->f1234[1] < zero ) 
					PI01 = - PI0;
				else
					PI01 = PI0;
				/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				// p part of the code aims on obtaining the sign of the r component of  
				// 4-momentum of the photon.                                               
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
				weierstrass_int_J3( p->tp2, infinity, p->rdd, p->rdel, zero, one,        /**/
						p->index_p4, p->rff_p, integ14, p->cases_int);           /**/
				p->half_periodwp = integ14[1];                                           /**/
				p->periodwp = two * p->half_periodwp;                                    /**/
			} else {
				if( p->f1234[1] <= zero ) {
					weierstrass_int_J3( p->tinf, p->thorizon, p->rdd, p->rdel, zero, one,
				                           p->index_p4, p->rff_p, integ14, p->cases_int);

					PI0_obs_hori = integ14[1];
				} else {
					weierstrass_int_J3( p->tp2, p->thorizon, p->rdd, p->rdel, zero, one, p->index_p4, 
						p->rff_p, integ14, p->cases_int);
					weierstrass_int_J3( p->tp2, p->tinf, p->rdd, p->rdel, zero, one, 
						p->index_p4, p->rff_p, integ4, p->cases_int);
					p->p_ini_r_tp2 = integ4[1];
					PI0_total_2 = integ14[1] + integ4[1];
				}
			}
		} // end casese.
		if ( p->a_spin == zero ) {
			cc = p->a2 - p->lam2 - p->q;
			if ( cc == zero ) { 
			} else if ( cc == - 27.0 ) {
				if ( p->f1234[1] < zero ) {
					cr = - three * fabs( (sqrt( p->robs*( p->robs + 6.0 ) )
						+ ( three + two * p->robs ) / sqrt3 ) / ( three - p->robs ) );
					dr = - fabs( ( sqrt( p->robs * ( p->robs+6.0 ) )
						+ ( three + two * p->robs ) / sqrt3 ) / ( three - p->robs ) );

				} else {
					cr = - three * fabs( ( sqrt( p->robs * ( p->robs + 6.0 ) )
						+ ( three + two * p->robs ) / sqrt3 ) / ( three - p->robs ) );
					dr = -fabs( ( sqrt( p->robs * ( p->robs + 6.0 ) )
						+ ( three + two * p->robs ) / sqrt3 ) / ( three - p->robs ) );
					PI0 = log( fabs( ( sqrt( p->robs * (p->robs + 6.0 ) )
						+ ( three + two * p->robs ) / sqrt3 ) / ( p->robs - three ) ) )
						/ three / sqrt3 - log( one + two / sqrt3 ) / three / sqrt3;
				}
			}
		}
	} else {
		u = creal( p->rbb[4] );
		u2 = u * u;
		w = fabs( cimag( p->rbb[4] ) );
		w2 = w * w;
		v = fabs( cimag( p->rbb[2] ) );
		v2 = v * v;
		if ( u != zero ) {
			// equation (45) in Yang & Wang (2012). 
			L1 = ( four * u2 + w2 + v2 + sqrt( sq( four*u2 + w2 + v2 )
				- four * w2 * v2 ) ) / ( two * w2 );
			L2 = ( four * u2 + w2 + v2 - sqrt( sq( four*u2 + w2 + v2 )
				- four * w2 * v2 ) ) / ( two * w2 );
			// equation (46) in Yang & Wang (2012). 
			p->thorizon = sqrt( ( L1-one ) / ( L1 - L2 ) )
				* ( p->rhorizon - u * ( L1 + one ) / ( L1 - one ) )
				/ sqrt( sq( p->rhorizon - u ) + w2 );
			// equation (48) in Yang & Wang (2012). 
			m2 = ( L1-L2 ) / L1;
			p->tinf = sqrt( ( L1 - one ) / ( L1 - L2 ) )
				* ( p->robs - u * ( L1 + one ) / ( L1 - one ) )
				/ sqrt( sq( p->robs-u ) + w2);
			p->t_inf = sqrt( ( L1 - one ) / ( L1 - L2 ) );
			// equation (50) in Yang & Wang (2012).
			sqrt_L1 = sqrt(L1);
			pinf = EllipticF( p->tinf, m2 ) / w / sqrt_L1;
			if( p->f1234[1] < zero ) {
				PI0 = pinf - EllipticF( p->thorizon, m2 ) / ( w * sqrt_L1 );
			} else {
				PI0 = EllipticF( p->t_inf,m2) / ( w * sqrt_L1 ) - pinf; 
			}
		} else {
			PI01 = atan( p->robs / w );
			if ( p->f1234[1] < zero )
				PI0 = ( PI01 - atan( p->rhorizon / w ) ) / w;
			else
				PI0 = ( pi / two - atan( p->robs / w ) ) / w;
		}
	}

	return p->radius;
}




/*
!*     PURPOSE:  Computes function r(p) defined by equation (41) and (49) in Yang & Wang (2012). That is
!*               r(p)=b0/(4*\wp(p+PIr;g_2,g_3)-b1)+r_tp1. \wp(p+PIr;g_2,g_3) is the Weierstrass'
!*               elliptic function; Or r=r_+, r=r_-. 
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f1234r---------p_1, r components of four momentum of a photon measured under a LNRF. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.  
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  radius---------radial coordinate of photon corresponding to a given p.
!*               sign_pr--------the sign of r component of 4-momentum of the photon. 
!*     ROUTINES CALLED: weierstrass_int_J3, mutp, root3.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS:
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-03.
!*
!*/
static double radius_tmp( ptcl *p, double pem )
{
	if (p->r_reals != 0 ) { 
		/* equation (39) in Yang & Wang (2012). */ 
		if ( p->cases == 1 ) {
			if ( !p->indrhorizon ) {
				if ( p->f1234[1]  <  zero ) {
					if ( pem < PI0_total ) {
						/* equation (41) in Yang & Wang (2012).*/
						p->radius = p->r_tp1 + p->rb0 / ( four * 
							weierstrassP( pem - PI0, p->rg2, p->rg3, 
							p->rdd, p->rdel )- p->rb1 );
						if ( pem < PI0 )
							p->sign_pr = -one;
						else
							p->sign_pr = one;
					} else {
						p->radius = infinity;  /* Goto infinity, far away. */
						p->sign_pr = one;
					}
				} else {
					if ( pem < PI0_inf_obs ) 
						/* equation (41) in Yang & Wang (2012). */
						p->radius=p->r_tp1 + p->rb0/(four * 
							weierstrassP( pem + PI0, p->rg2, 
							p->rg3, p->rdd, p->rdel) - p->rb1 );
				        else
						p->radius=infinity;  /*Goto infinity, far away.*/
				        p->sign_pr = one;
				}
			} else {
				//printf(" pr = %f \t \n  ", p->f1234[1] );
				if ( p->f1234[1] < zero ) {
				        if ( pem < PI0_obs_hori ) {
						/* equation (41) in Yang & Wang (2012). */
						p->radius = p->r_tp1 + p->rb0 / ( four *
						weierstrassP( pem - PI0, p->rg2, p->rg3, 
							p->rdd, p->rdel) - p->rb1 );
				        } else
						p->radius = p->rhorizon;  /*Fall into black hole. */
				        p->sign_pr = -one;
				} else {
				        if ( pem < PI0_inf_obs )
					/* equation (41) in Yang & Wang (2012). */
						p->radius=p->r_tp1 + p->rb0 /
						(four*weierstrassP( pem + PI0, p->rg2, 
							p->rg3, p->rdd, p->rdel) - p->rb1 );
				        else
						p->radius=infinity; /*Goto infinity, far away.*/
				        p->sign_pr = one;
				}
			}
		} else if ( p->cases == 2 ) {   //  case(2)
			if ( !p->indrhorizon ) {
				/* equation (41) in Yang & Wang (2012). */
				p->radius = p->r_tp1 + p->rb0 / ( four* weierstrassP( pem + PI01, 
					p->rg2, p->rg3, p->rdd, p->rdel ) - p->rb1 );
				/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				// p part of the code aims on obtaining the sign of the r component of  
				// 4-momentum of the photon.                                               
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
				double up;                                                                 //
				if ( p->f1234[1] > zero ) {                                                //
					up = fmod(pem + PI01, p->periodwp);                                //
					if ( up <= p->half_periodwp )                                      //
						p->sign_pr = one;                                          //
					else                                                               //
						p->sign_pr = - one;                                        //
				} else if ( p->f1234[1] < zero ) {                                         //
					up = fmod(pem + PI01 + p->half_periodwp, p->periodwp);             //
					if ( up < p->half_periodwp )                                       //
						p->sign_pr = -one;                                         //
					else                                                               //
						p->sign_pr = one;                                          //
				} else {                                                                   //
					if ( p->robs == p->r_tp1 ) {                                       //
						if ( fmod(pem, p->periodwp) <= p->half_periodwp )          //
							p->sign_pr = one;                                  //
						else                                                       //
							p->sign_pr = -one;                                 //
					} else {                                                           //
						if ( fmod(pem + p->half_periodwp, p->periodwp)             //
				                                   <= p->half_periodwp )                   //
							p->sign_pr = one;                                  //
						else                                                       //
							p->sign_pr = -one;                                 //
				        }                                                                  //
				}                                                                          //
				/*=========================================================================*/
			} else {
				if( p->f1234[1] <= zero ) {
					if ( pem < PI0_obs_hori )
					/* equation (41) in Yang & Wang (2012). */
						p->radius = p->r_tp1 + p->rb0 / ( four * 
							weierstrassP( pem - PI0, p->rg2, 
							p->rg3, p->rdd, p->rdel) - p->rb1 );
					else
						p->radius = p->rhorizon; //Fall into black hole
				        p->sign_pr = -one;
				} else {
					if ( pem < PI0_total_2 )
						// equation (41) in Yang & Wang (2012). 
						p->radius = p->r_tp1 + p->rb0 / ( four * 
							weierstrassP( pem + PI0, p->rg2, 
							p->rg3, p->rdd, p->rdel ) - p->rb1 );
				        else
						p->radius = p->rhorizon; //Fall into black hole. 
					//!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					if (pem <= p->p_ini_r_tp2)                  /**/
						p->sign_pr = one;                   /**/
				        else                                        /**/
						p->sign_pr = -one;                  /**/
					//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
				}
			}
		} // end casese.
		if ( p->a_spin == zero ) {
			if ( cc == zero ) {
				if ( p->f1234[1] < zero ) {
					if ( pem < one / p->rhorizon - one / p->robs )
						p->radius = p->robs / ( p->robs * pem + one );
					else
						p->radius = p->rhorizon;
					p->sign_pr = -one;
				} else {
					if ( pem < one / p->robs )
						p->radius = p->robs / ( one - p->robs * pem );
					else
						p->radius = infinity;
				        p->sign_pr = one;
				}
			} else if ( cc == -27.0 ) {
				double tmp1, tmp2;
				if ( p->f1234[1] < zero ) {
					tmp1 = cr * exp( three * sqrt3 * pem ) - sqrt3;
					tmp2 = dr * exp( three * sqrt3 * pem ) + two / sqrt3;
					if ( pem != zero )
						p->radius = ( three + tmp1 * tmp2 + sqrt( 9.0 + 6.0 * tmp1 * tmp2
							+ tmp1 * tmp1 ) ) / ( tmp2 * tmp2 - one );
					else
						p->radius = p->robs; //infinity;

				        p->sign_pr = - one;
				} else {
					tmp1 = cr * exp( - three * sqrt3 * pem ) - sqrt3;
					tmp2 = dr * exp( - three * sqrt3 * pem ) + two / sqrt3;

					if ( pem < PI0 )
						p->radius = ( three + tmp1 * tmp2 + sqrt( 9.0 + 6.0 * tmp1 * tmp2
							+ tmp1 * tmp1 ) ) / ( tmp2 * tmp2 - one );
					else
						p->radius = infinity;
					p->sign_pr = one;
				}
			}
		}
	} else {
		if ( u != zero ) {
			sncndn( pem * w * sqrt_L1 + sign( p->f1234[1] ) * 
				pinf * w * sqrt_L1, one - m2, &sn, &cn, &dn);
			if( p->f1234[1] < zero ) {
				if ( pem < PI0 ) 
					// equation (49) in Yang & Wang (2012), and p_r <0, r=r_{+}
					p->radius = u + ( - two * u + w * (L1-L2) * sn * fabs(cn) )
						/ ( ( L1 - L2 ) * sn * sn - ( L1 - one ) );
				else
					p->radius = p->rhorizon;

				p->sign_pr = - one;
			} else {
				if ( pem < PI0 ) 
					// equation (49) in Yang & Wang (2012), and p_r >0, r=r_{-}
					p->radius = u + ( - two * u - w * ( L1 - L2 ) * sn * fabs(cn) )
						/ ( ( L1 - L2 ) * sn * sn - ( L1 - one ) );
				else
					p->radius = infinity;

				p->sign_pr = one;
			}
		} else {
			if ( p->f1234[1] < zero ) {
				if ( pem < PI0 ) 
					p->radius = w * tan( PI01 - pem * w );
				else
					p->radius = p->rhorizon;

				p->sign_pr = - one;
			} else {
				if( pem < PI0 ) 
					p->radius = w * tan( PI01 + pem * w );
				else
					p->radius = infinity;

				p->sign_pr = one;
			}
		}
	}

	return p->radius;
}




/*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-19.
!*
!*/
double radius( ptcl *pt, double p )
{

	static unsigned int count_num = 1;

restarting:
	if ( count_num == 1 ) {

		radius_settings( pt );
		count_num += 1;

		a_spin_1 = pt->a_spin;
		scal_1 = pt->scal;
		lambda_1 = pt->lambda;
		q_1 = pt->q;
		f12341_1 = pt->f1234[1];
		f12342_1 = pt->f1234[2];
		f12343_1 = pt->f1234[3];
		muobs_1 = pt->muobs;
		sinobs_1 = pt->sinobs;
		robs_1 = pt->robs;

		return (radius_tmp( pt, p ));

    	} else {
		if ( f12343_1 == pt->f1234[3]
		&& f12342_1 == pt->f1234[2]
		&& lambda_1 == pt->lambda
		&& q_1 == pt->q
		&& sinobs_1 == pt->sinobs
		&& muobs_1 == pt->muobs
		&& a_spin_1 == pt->a_spin
		&& scal_1 == pt->scal
		&& robs_1 == pt->robs ) {
 
			return (radius_tmp( pt, p ));

		} else {
			count_num = 1;
			goto restarting;
		}

	}
}



/*
!*     PURPOSE:  Computes the value of parameter p from radial coordinate. In other words, to compute 
!*               the r part of integral of equation (24), using formula (58) in Yang & Wang (2012).
!*               (58) is: p=-sign(p_r)*p_0+2*t1*p_2+2*t2*p_2. where p_r is initial radial
!*               component of 4 momentum of photon. 
!*     INPUTS:   f1234r---------f_1, which was defined by equation (106) in Yang & Wang (2012). 
!*               a_spin---------spin of black hole, on interval (-1,1).
!*               robs-----------radial coordinate of the observer. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
!*               scal-----------a dimentionless parameter to control the size of the images. 
!*                              Which is usually be set to 1.D0.   
!*               t1,t2----------Number of photon meets the turning points r_tp1 and r_tp2
!*                              respectively in radial motion.
!*     OUTPUTS:  r2p------------value of r part of integral (24) in Yang & Wang (2012).
!*     ROUTINES CALLED: radiustp, root3, weierstrass_int_j3, EllipticF.
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
double r2p( ptcl *p, double rend, int t1, int t2 )
{
	double integ4[5], integ04[5], integ14[5], pp, p1, p2;
	double ant = zero;

	radiustp( p );

	if (p->r_reals != 0 ) {
		if ( rend < p->r_tp1 || rend > p->r_tp2 )
			return (-one);

		/* equations (35)-(38) in Yang & Wang (2012). */
		p->rb0 = four * sq3(p->r_tp1) + two * (p->a2 - p->lam2 - p->q ) * p->r_tp1 + 
			two * ( p->q + sq( p->lambda - p->a_spin ) );
		p->rb1 = two * sq(p->r_tp1) + one / three * ( p->a2 - p->lam2 - p->q );
		p->rb2 = four / three * p->r_tp1;
		p->rb3 = one;
		p->rg2 = three / four * ( sq(p->rb1) - p->rb0 * p->rb2 );
		p->rg3 = one / 16.0 * ( 3.0 * p->rb0 * p->rb1 * p->rb2 - 2.0 * sq3(p->rb1) - sq(p->rb0) * p->rb3 );

		/* equation (39) in Yang & Wang (2012). */
		if ( p->robs - p->r_tp1 != zero )
			p->tinf = p->rb0 / four / ( p->robs - p->r_tp1 ) + p->rb1 / four;
		else
			p->tinf = infinity;

		if ( p->rhorizon - p->r_tp1 != zero ) 
			p->thorizon = p->rb1 / four + p->rb0 / four / ( p->rhorizon - p->r_tp1 );
		else
			p->thorizon = infinity;

		if ( rend - p->r_tp1 != zero ) 
			p->tp = p->rb1 / four + p->rb0 / four / ( rend - p->r_tp1 );
		else
			p->tp = infinity;

		p->tp2 = p->rb0 / four / ( p->r_tp2 - p->r_tp1 ) + p->rb1 / four;
		p->tinf1 = p->rb1 / four;

		root3( zero, - (p->rg2) / four, - (p->rg3) / four, p->rdd, &p->rdel );

		p->index_p4[1] = 0;
		p->cases_int = 1;

		/* equation (42) in Yang & Wang (2012). */
		weierstrass_int_J3( p->tinf, infinity, p->rdd, p->rdel, zero, one, 
				p->index_p4, p->rff_p, integ04, p->cases_int );
		PI0 = integ04[1]; 
		//printf(" ss = %f \t %f \t %f \t \n  ", p->rg2, p->rg3, PI0 );

		if ( p->cases == 1 ) {
			if ( !p->indrhorizon ) {
				if ( p->f1234[1]  <  zero ) {
					weierstrass_int_J3( p->tinf, p->tp, 
						p->rdd, p->rdel, zero, one, 
						p->index_p4, p->rff_p, integ14, 
						p->cases_int );

					pp = integ14[1];
					if ( t1 != 0 ) {
						weierstrass_int_J3( p->tp, infinity,
							p->rdd, p->rdel, zero, one,
							p->index_p4, p->rff_p, 
							integ14, p->cases_int );
						p1 = integ14[1];
					} else
						p1 = zero;

					ant = pp + two * p1 * t1;
					//printf("ant pp = %40.35f \n p1 = %40.35f \n  t1 = %d \n", pp, p1, t1);
					//printf("ant pp = %40.35f \n p1 = %40.35f \n ", p->tinf, p->tp);
				} else {
					weierstrass_int_J3( p->tinf, p->tp, 
						p->rdd, p->rdel, zero, one,
						p->index_p4, p->rff_p, integ04, 
						p->cases_int );

					ant = - integ04[1];
				}
			} else {
				//printf(" pr = %f \t \n  ", p->f1234[1] );
				if ( p->f1234[1] < zero ) {
					if ( rend <= p->rhorizon ) {
						weierstrass_int_J3( p->tinf, p->thorizon, 
							p->rdd, p->rdel, zero, one,
							p->index_p4, p->rff_p, integ04, 
							p->cases_int );
						ant = integ04[1];
					} else {
						weierstrass_int_J3( p->tinf, p->tp, 
							p->rdd, p->rdel, zero, one,
							p->index_p4, p->rff_p, integ04, 
							p->cases_int );
						ant = integ04[1];
					}
				} else {
					if ( rend < infinity ) {
						weierstrass_int_J3( p->tinf, p->tp, 
							p->rdd, p->rdel, zero, one,
							p->index_p4, p->rff_p, integ04, 
							p->cases_int );
						ant = - integ04[1];
					} else {
						weierstrass_int_J3( p->tinf1, p->tinf, 
							p->rdd, p->rdel, zero, one, 
							p->index_p4, p->rff_p, integ04, 
							p->cases_int );

						ant = integ04[1];
					}
				}
			}
		} else if ( p->cases == 2 ) { // case(2)
			if ( !p->indrhorizon ) {
				weierstrass_int_J3( p->tinf, p->tp, p->rdd, p->rdel, zero, one,
						p->index_p4, p->rff_p, integ4, p->cases_int);
				pp = integ4[1];
				if ( t1 == 0 )
					p1 = zero;
				else {
					weierstrass_int_J3( p->tp, infinity, p->rdd, 
						p->rdel, zero, one, p->index_p4, 
						p->rff_p, integ4, p->cases_int );
					p1 = integ4[1];
				}

				if ( t2 == 0 )
					p2 = zero;
				else {
					weierstrass_int_J3( p->tp2, p->tp, p->rdd, 
						p->rdel, zero, one, p->index_p4, 
						p->rff_p, integ4, p->cases_int );
					p2 = integ4[1];
				}

				if ( p->f1234[1] != zero )
					ant = sign( - p->f1234[1] ) * pp + 
						two * ( t1 * p1 + t2 * p2 );
				else {
					if ( p->robs == p->r_tp1 )
						ant = - pp + two * ( t1 * p1 + t2 * p2 );
					else
						ant =   pp + two * ( t1 * p1 + t2 * p2 );
				}
			} else {
				if( p->f1234[1] <= zero ) {
					if ( rend <= p->rhorizon ) {
						weierstrass_int_J3( p->tinf, p->thorizon, p->rdd, 
							p->rdel, zero, one, p->index_p4, 
							p->rff_p, integ4, p->cases_int );
						pp = integ4[1];
					} else {
						weierstrass_int_J3( p->tinf, p->tp, p->rdd, 
							p->rdel, zero, one, p->index_p4, 
							p->rff_p, integ4, p->cases_int );
						pp = integ4[1];
					}
				} else {
					weierstrass_int_J3( p->tinf, p->tp, p->rdd, 
						p->rdel, zero, one, p->index_p4, 
						p->rff_p, integ4, p->cases_int);
					pp = integ4[1];

					if ( t2 == 0 )
						p2 = zero;
					else {
						weierstrass_int_J3( p->tp2, p->tp, p->rdd, 
							p->rdel, zero, one, p->index_p4, 
							p->rff_p, integ4, p->cases_int );
						p2 = integ4[1];
					}
					ant = -pp + two * t2 * p2;
				}
			}
		} // end casese.
		if ( p->a_spin == zero ) {
			cc = p->a2 - p->lam2 - p->q;
			if ( cc == zero ) {
				if ( p->f1234[1] < zero ) {
					if ( rend <= p->rhorizon )
						ant = one / p->rhorizon - one / p->robs;
					else
						ant = one / rend - one / p->robs;
				} else {
					if ( rend <= infinity )
						ant = one / p->robs - one / rend;
					else
						ant = one / p->robs;
				}
			} else if ( cc == - 27.0 ) {
				if ( p->f1234[1] < zero ) {
					if ( rend > p->rhorizon )
						ant = - log( fabs( ( sqrt( p->robs
							* ( p->robs + 6.0 ) ) + ( three + two * p->robs )
							/ (sqrt3) ) / ( p->robs - three ) ) )
							/ ( three * sqrt3 ) + log( fabs( ( 
							sqrt( rend * ( rend + 6.0 ) ) + ( three + two * rend )
							/ ( sqrt3 ) ) / ( rend - three ) ) ) / ( three * sqrt3 );
					else
						ant = - log( fabs( ( sqrt( p->robs
							* ( p->robs + 6.0 ) ) + ( three + two * p->robs )
							/ (sqrt3) ) / ( p->robs - three ) ) )
							/ ( three * sqrt3 ) + log( fabs( ( 
							sqrt( p->rhorizon * ( p->rhorizon + 6.0 ) )
							+ ( three + two * p->rhorizon )
							/ ( sqrt3 ) ) / ( p->rhorizon - three ) ) )
							/ ( three * sqrt3 );

				} else {
					if ( rend < infinity )
						ant = - log( fabs( ( sqrt( rend * ( rend + 6.0 ) )
							+ ( three + two * rend ) / ( sqrt3 ) )
							/ ( rend - three ) ) ) / ( three * sqrt3 )
							+ log( fabs( ( sqrt( p->robs * ( p->robs + 6.0 ) )
							+ ( three + two * p->robs ) / ( sqrt3 ) )
							/ ( p->robs - three ) ) ) / ( three * sqrt3 );
					else
						ant = - log( one + two / sqrt3 ) / three / sqrt3
							+ log( fabs( ( sqrt( rend * ( rend + 6.0 ) )
							+ ( three + two * rend ) / ( sqrt3 ) )
							/ ( rend - three ) ) ) / ( three * sqrt3 );
				}
			}
		}
	} else {
		u = creal( p->rbb[4] );
		u2 = u * u;
		w = fabs( cimag( p->rbb[4] ) );
		w2 = w * w;
		v = fabs( cimag( p->rbb[2] ) );
		v2 = v * v;
		if ( u != zero ) {
			// equation (45) in Yang & Wang (2012). 
			L1 = ( four * u2 + w2 + v2 + sqrt( sq( four*u2 + w2 + v2 )
				- four * w2 * v2 ) ) / ( two * w2 );
			L2 = ( four * u2 + w2 + v2 - sqrt( sq( four*u2 + w2 + v2 )
				- four * w2 * v2 ) ) / ( two * w2 );
			// equation (46) in Yang & Wang (2012). 
			p->thorizon = sqrt( ( L1-one ) / ( L1 - L2 ) )
				* ( p->rhorizon - u * ( L1 + one ) / ( L1 - one ) )
				/ sqrt( sq( p->rhorizon - u ) + w2 );
			// equation (48) in Yang & Wang (2012). 
			m2 = ( L1-L2 ) / L1;
			p->tinf = sqrt( ( L1 - one ) / ( L1 - L2 ) )
				* ( p->robs - u * ( L1 + one ) / ( L1 - one ) )
				/ sqrt( sq( p->robs-u ) + w2);
			p->t_inf = sqrt( ( L1 - one ) / ( L1 - L2 ) );
			// equation (50) in Yang & Wang (2012).
			sqrt_L1 = sqrt(L1);
			pinf = EllipticF( p->tinf, m2 ) / w / sqrt_L1;
			if( p->f1234[1] < zero ) {
				if ( rend <= p->rhorizon )
					ant = pinf - EllipticF( p->thorizon, m2 )
						/ ( w * sqrt( L1 ) );
				else
					ant = pinf - EllipticF( p->tp, m2)
						/ ( w * sqrt( L1 ) );
			} else {
				if ( rend < infinity )
					ant = EllipticF( p->tp, m2 ) / ( w * sqrt( L1 ) ) - pinf;
				else
					ant = EllipticF( p->t_inf, m2 ) / ( w * sqrt( L1 ) ) - pinf;
			}
		} else {
			if ( p->f1234[1] < zero )
				if ( rend <= p->rhorizon )
					ant = ( atan( p->robs / w ) - atan( p->rhorizon / w ) ) / w;
				else
					ant = ( atan( p->robs / w ) - atan( rend / w ) ) / w;
			else
				if ( rend < infinity )
					ant = ( atan( rend / w ) - atan( p->robs / w ) ) / w;
				else
					ant = ( halfpi - atan( p->robs / w ) ) / w;
		}
	}

	return ant;
}


/*
!*     PURPOSE:  Solves equation \mu(p)=mu, i.e. to search the value p_{em} of
!*               parameter p, corresponding to the intersection point of geodesic with with 
!*               disk, i.e. a surface has a constants inclination angle with respect to 
!*               equatorial plane of black hole.  
!* 
!*     INPUTS:   f1234(1:4)-----array of p_r, p_theta, p_phi, p_t, which are defined by equation 
!*                              (82)-(85) in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon.
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               mu-------------mu=\cos(\pi/2-\theta_disk), where \theta_disk is the inclination 
!*                              angle of disk surface with respect to the equatorial plane of 
!*                              black hole.    
!*               rin, rout------inner and outer radius of disk.        
!*     OUTPUTS:  pemdisk--------value of root of equation \mu(p)= mu for p.  
!*                              pemdisk=-1.D0, if the photon goto infinity.
!*                              pemdisk=-2.D0, if the photon fall into event horizon.       
!*     REMARKS:                 This routine just search the intersection points of geodesic with 
!*                              up surface of disk. Following routine Pemdisk_all will searches 
!*                              intersection points of geodesic with up and down surface of disk.      
!*     ROUTINES CALLED: mutp, mu2p.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: 
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-08.
!*/
double Pemdisk( ptcl * p, double mu, double rout, double rin, double *re )
{
	double pm;

	int t1, t2;

	mutp( p );
 
	if ( p->mu_reals == 2 ) {
		if ( p->mobseqmtp ) {
		        t1 = 0;
		        t2 = 0;
		} else {
			if ( p->muobs > zero ) {
				if ( p->f1234[2] < zero ) {
				        t1 = 1;
				        t2 = 0;
				}
				if ( p->f1234[2] > zero ) {
				        t1 = 0;
				        t2 = 0;
				}
			} else {
				if ( p->muobs == zero ) {
					if ( p->f1234[2] < zero ) {
				            t1 = 1;
				            t2 = 0;
					}
				        if ( p->f1234[2] > zero ) {
				            t1 = 0;
				            t2 = 1;
					}                       
				} else {
				        if ( p->f1234[2] < zero ) {
				            t1 = 0;
				            t2 = 0;
					}
					if ( p->f1234[2] > zero ) {
				            t1 = 0;
				            t2 = 1;
					}
				}
			}
		}
		//printf("mu = %f  \t \n %d \t %d \n",mu, t1, t2);
		pm = mu2p( p, mu, t1, t2 );
		//printf("~~~~~~~~~~~~~~~1111111111~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
		*re = radius( p, pm );
		//printf("~~~~~~~~~~~~~~~1111111111~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
		//printf("mu = %20.15f  \t  %d \t %d \t %20.15f \n", pm, t1, t2, *re);
		if (*re <= rout && *re >= rin) {
			return (pm);
		} else if ( *re > rout ) {
			return (-one);
		} else if ( *re < rin ) {
			return (-two);
		}
	} else {
		//printf("mu = %20.15f  \t  %d \t %d \t %20.15f \n", pm, t1, t2, *re);
		return (- one);
	}
	return (pm);
	//printf("mu = %f  \t \n %d \t %d \n",mu, t1, t2);
}


/*
!*     PURPOSE:  Solves equation \mu(p)=mu, where \mu(p)=\mu(p), i.e. to search the value p_{em} of
!*               parameter p, corresponding to the intersection point of geodesic with 
!*               disk, i.e. a surface has a constants inclination angle with respect to 
!*               equatorial plane of black hole.
!* 
!*     INPUTS:   f1234(1:4)-----array of f_1, f_2, f_3, f_0, which are defined by equation 
!*                              (106)-(109) in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               mu-------------mu=\cos(\pi/2-\theta_disk), where \theta_disk is the inclination 
!*                              angle of disk surface with respect to the equatorial plane of 
!*                              black hole.  
!*               rin, rout------inner and outer radius of disk.        
!*     OUTPUTS:  pemdisk_all--------value of root of equation \mu(p)= mu for p.  
!*                              pemdisk=-1.D0, if the photon goto infinity.
!*                              pemdisk=-2.D0, if the photon fall into event horizon.        
!*     REMARKS:                 This routine will searches intersection points of 
!*                              geodesic with double surfaces of disk.      
!*     ROUTINES CALLED: mutp, mu2p, radius.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: 
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-10.
!*/
double Pemdisk_all( ptcl * p, double mu, double rout, double rin, double *re )
{
	double pm;
	int t1, t2;

	mutp( p );

	t1 = 0;
	t2 = 0;
	if ( p->mu_reals == 2 ) {
		for ( int j = 0; j < 11; j++ ) {
			for ( int i = j; i <= j + 1; i++ ) {
				if ( p->mobseqmtp ) {
					if ( p->muobs == p->mu_tp1 ) {
						t1 = j;
						t2 = i;
					} else {
						t1 = i;
						t2 = j;
					}
				} else {
					if ( p->muobs > zero ) {
						if ( p->f1234[2] < zero ) {
							t1 = i;
							t2 = j;
						}
						if ( p->f1234[2] > zero ) {
							t1 = j;
							t2 = i;
						}
					else
						if ( p->muobs == zero ) {
							if ( p->f1234[2] < zero ) {
								t1 = i;
								t2 = j;
							}
							if ( p->f1234[2] > zero ) {
								t1 = j;
								t2 = i;
							}                      
						} else {
							if ( p->f1234[2] < zero ) {
								t1 = i;
								t2 = j;
							}
							if ( p->f1234[2] > zero ) {
								t1 = j;
								t2 = i;
							}
						}       
					}
				}
				pm = mu2p( p, mu, t1, t2 );
				//printf("~~~~~~~~~~~~~~tt1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
				//printf( "mu = %d \t pm = %f \n", p->mu_reals, pm );
				if ( pm  <=  zero ) continue;
				radius_settings( p );
				*re = radius( p, pm );
				//printf("~~~~~~~~~~~~~~~1111111111~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
				//printf("mu = %20.15f  \t  %d \t %d \t %20.15f \n", pm, t1, t2, *re);
				if ( *re <= rout && *re >= rin )
					return pm;
				else {
					if ( *re >= infinity ) {
						return (-one);
				        } else {
						if( *re  <=  p->rhorizon ) {
							return (-two);
						}
					}
				}
			}
		}
	} else {
		//printf("~~~~~~~~~~~~~~~22222~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
		//printf("mu = %20.15f  \t  %d \t %d \t %20.15f \n", pm, t1, t2, *re);
		return (-two);
	}
        return pm;
}



void Get_data( out_data2 * tmp )
{
	(*tmp).u = u;
	(*tmp).v = v;
	(*tmp).w = w;
	(*tmp).u2 = u2;
	(*tmp).v2 = v2;
	(*tmp).w2 = w2;

	(*tmp).L1 = L1;
	(*tmp).L2 = L2;
	(*tmp).m2 = m2;


	(*tmp).cc = cc;
	(*tmp).cr = cr;
	(*tmp).dr = dr;
	(*tmp).pinf = pinf;

	(*tmp).PI0 = PI0;
	(*tmp).PI0_total = PI0_total;
	(*tmp).PI0_inf_obs = PI0_inf_obs;

	(*tmp).PI0_obs_hori = PI0_obs_hori;
	(*tmp).PI01 = PI01;
	(*tmp).PI0_total_2 = PI0_total_2;

	(*tmp).sqrt_L1 = sqrt_L1;
}



/*
!*     REVISIONS:
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-03.
!*
!*/
double p_total( ptcl *p )
{
	radius_settings( p );

	double ptotal = zero;
	if (p->r_reals != 0 ) {
		if ( p->cases == 1 ) {
			if ( !p->indrhorizon ) {
				if ( p->f1234[1]  <  zero ) {
					ptotal = PI0_total;
				} else {
					ptotal = PI0_inf_obs;
				}
			} else {
				//printf(" pr = %f \t \n  ", p->f1234[1] );
				if ( p->f1234[1] < zero ) {
				        ptotal = PI0_obs_hori;
				} else {
					ptotal = PI0_inf_obs;
				}
			}
		} else if ( p->cases == 2 ) { // case(2)
			if ( !p->indrhorizon ) {
				ptotal = p->half_periodwp * 4.0;
			} else {
				if( p->f1234[1] <= zero ) {
					ptotal = PI0_obs_hori;
				} else {
					ptotal = PI0_total_2;
				}
			}
		} // end casese.
	} else {
		if ( u != zero ) {
			if( p->f1234[1] < zero ) {
				ptotal = PI0;
			} else {
				ptotal = PI0;
			}
		} else {
			if ( p->f1234[1] < zero )
				ptotal = ( atan( p->robs / w ) - atan( p->rhorizon / w ) ) / w;
			else
				ptotal = ( halfpi - atan( p->robs / w ) ) / w;
		}
	}
	return ptotal;
}


























