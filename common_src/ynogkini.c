/*
 * common_src/particle.c
 *
 *     PURPOSE: th module aims on computing the four Boyer-Lindquist 
 *              coordinates (r,\theta,\phi,t) and the affine 
 *              parameter \sigam of the Kerr space-time.    
 *   
 *     AUTHOR:        Yang & Wang (2013)  
 *     DATE WRITTEN:  4 Jan 2013
 *
 *     C language Version:   Yang Xiao-lin  2022-11-22.
 *
 * Author       Yang Xiao-lin
 *
 * City         Kun-ming, Yunnan Provice, China.
 *
 * 2022-11-22   Starting the writting of the module.
 */


#include "ynogkini.h"

  
//********************************************
// variables for int_r_part functions.
//********************************************
static double PI1_phi;
static double PI2_phi;
static double PI1_time;
static double PI2_time;
static double PI1_aff;
static double PI2_aff;
static double Ap, Am;
static double B_add, B_m;
static double r_add, r_m;
static double r_add2, r_m2;
static double rhorizon;

static double a4, b4;
static double a5, b5;
static double t_inf;
static double f1, g1, h1;
static double f2, g2, h2;



//static double cc;
int robs_eq_rtp;
int indrhorizon;


double r_tp1, r_tp1sq;
double r_tp2;
double robs;
double a_spin;
double lambda;
double muobs, sinobs;
unsigned int r_reals;

int r_cases;       /* If r_tp2=infinity, then cases=1, else cases=2. */
double _Complex R_roots[5];
double _Complex bb[5], dd[4];
//********************************************


 
//********************************************
static double b0, b0sq, b1, b2, b3, g2, g3;
static double u, v, w, L1, L2, m2;
static double u2, v2, w2;
static double sn, cn, dn;
static double pinf; 
static double PI0, PI0_total, PI0_inf_obs, PI0_obs_inf;
static double PI0_obs_hori, PI01, PI0_total_2, PI0_obs_tp2;
static double sqrt_L1;
//static double t_inf;


static double h;
static double E_add, E_m;
static double D_add, D_m;
static double hp, hm;
static double tobs;
static double thorizon;
static double tp2;
static double tinf;
static double tp;

static double wp, wm, wbarp, wbarm;
static double rff_p, integ04[5], integ4[5], integ14[5];
static double integ5[6], integ15[6];


static int cases_int, index_p4[5], index_p5[6], del;
int t1, t2;

static double r_coord, sign_pr;
static double pp, p1, p2, p_tp1_tp2;
static double PI1_p, PI2_p;
static double half_periodwp, periodwp;


static double x, y;
//********************************************
 
static void Integral_r_part( ptcl * p, double pm, double *affr, 
				double *timer, double *phir );


/*
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The Get_t1_t2() Subroutine was firstly written at 2017--11--16 am 11:23 in Fortran,
! and transformed into C version by Yang Xiao-lin at 2022--12--15 am 10:00.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-15.
!*
!*/
static void Get_t1_t2( double p, double f1234r, int robs_eq_rtp, double p_tp1_tp2, 
		double robs, double r_tp1, double r_tp2, int *t1, int *t2 )
{
	double p1_temp, p2_temp;
	int N_temp;
	if ( robs_eq_rtp ) {
		p1_temp = zero;
		p2_temp = p_tp1_tp2;
		*t1 = 0;
		*t2 = 0;
		N_temp = 0;
		if ( robs == r_tp1 ) {
			do {
				if ( p1_temp <= p && p <= p2_temp )
					break;
				else {
                                      N_temp = N_temp + 1;
                                      p1_temp = p2_temp;
                                      p2_temp = p2_temp + p_tp1_tp2;
                                      *t1 = *t2;
                                      *t2 = N_temp - *t1;
				}
			} while (true);
		} else if (robs == r_tp2) {
			do {
				if ( p1_temp <= p && p <= p2_temp )
                                      break;
				else {
                                      N_temp = N_temp + 1;
                                      p1_temp = p2_temp;
                                      p2_temp = p2_temp + p_tp1_tp2;
                                      *t2 = *t1;
                                      *t1 = N_temp-*t2;
				}
			} while(true);
		}
	} else {
		p1_temp = zero;
		*t1 = 0;
		*t2 = 0;
		N_temp = 0;
		if ( f1234r > zero ) {
			p2_temp = PI2_p;
			do {
				if (p1_temp <= p && p <= p2_temp)
					break;
				else {
					N_temp = N_temp + 1;
					p1_temp = p2_temp;
					p2_temp = p2_temp + p_tp1_tp2;
					*t1 = *t2;
					*t2 = N_temp - *t1;
				}
			} while(true);
		}
		if ( f1234r < zero ) {
			p2_temp = PI1_p;
			do {
				if ( p1_temp <= p && p <= p2_temp )
					break;
				else {
					N_temp = N_temp + 1;
					p1_temp = p2_temp;
					p2_temp = p2_temp + p_tp1_tp2;
					*t2 = *t1;
					*t1 = N_temp - *t2;
				}
			}  while(true);
		}
	}
}


/*
!*     C VERSION:  Yang Xiao-lin    2022-12-15.
!*/
static void Get_t1_t2_old( double p, double f1234r, int robs_eq_rtp, double pp, 
		double p1, double p2, double robs, double r_tp1, int *t1, int *t2 )
{
	double p_temp = zero;
	//int tt1, tt2;
	for ( int j = 0; j < 101; j++ ) {
		for ( int i = j; i <= j + 1; i++ ) {
			if ( robs_eq_rtp ) {
				if ( robs == r_tp1 ) {
					*t1 = j;
					*t2 = i;
					p_temp = - pp + two
					* ( (*t1) * p1 + (*t2) * p2 );
				} else {
					*t1 = i;
					*t2 = j;
					p_temp = pp + two
					* ( (*t1) * p1 + (*t2) * p2 );
				}
			} else {
				if ( f1234r > zero ) {
					*t1 = j;
					*t2 = i;
					p_temp = - pp + two
					* ( (*t1) * p1 + (*t2) * p2 );
				}
				if ( f1234r < zero ) {
					*t1 = i;
					*t2 = j;
					p_temp = pp + two*( (*t1)*p1 + (*t2) * p2 );
				}
			}
			if ( fabs( p - p_temp ) < 1e-4 )
				goto here;
		}
	}
here:
	//tt1 = *t1;
	//tt2 = *t2;
}

/*
|*
|*
|*
|*
!*     C VERSION:  Yang Xiao-lin    2022-12-13.
|*/
void Set_Initializations_For_int_r_part( ptcl * p, double pm, 
		double *affr, double *timer, double *phir )
{
	out_data2 tmp;
	Get_data( &tmp );

	robs = p->robs;
	lambda = p->lambda;
	a_spin = p->a_spin;
	muobs = p->muobs;
	sinobs = p->sinobs;

	rhorizon = p->rhorizon;
	// equation (64) in Yang & Wang (2013).
	r_add = rhorizon;
	r_m = one - sqrt( one - p->a2 );
	r_add2 = sq( r_add );
	r_m2 = sq( r_m );   
	// equation (64) in Yang & Wang (2013).  
	B_add = ( two * r_add - p->a_spin * p->lambda ) / ( r_add - r_m );
	B_m = ( two * r_m - p->a_spin * p->lambda) / ( r_add - r_m );
	// equation (64) in Yang & Wang (2013).
	Ap = ( r_add * ( four - p->a_spin * p->lambda )
		- two * p->a2 ) / sqrt( one - p->a2 );
	Am = ( r_m * ( four - p->a_spin * p->lambda ) 
		- two * p->a2 ) / sqrt( one - p->a2 );
	b4 = one;
	a4 = zero;
	//cc = p->a2 - p->lam2 - p->q;
	
	radiustp( p );
	r_tp1 = p->r_tp1;
	r_tp1sq = sq( r_tp1 );
	r_tp2 = p->r_tp2;
	r_reals = p->r_reals;
	robs_eq_rtp = p->robs_eq_rtp;
	indrhorizon = p->indrhorizon;
	r_cases = p->cases;

	for ( int i = 1; i < 5; i++ )
		bb[i] = p->rbb[i];

	PI1_phi = zero;
	PI2_phi = zero;
	PI1_time = zero;
	PI2_time = zero;
	PI1_aff = zero;
	PI2_aff = zero;

	if ( p->r_reals != 0 ) {
		/* equations (35)-(38) in Yang & Wang (2013). */
		b0 = p->rb0;
		b0sq = sq( b0 );
		b1 = p->rb1;
		b2 = p->rb2;
		b3 = p->rb3;
		g2 = p->rg2;
		g3 = p->rg3;

		tobs = p->tinf;

		thorizon = p->thorizon;

		tp2 = p->tp2;
		tinf = p->tinf1;
		h = - b1 / four;
		// equation (64), (66) and (70) in Yang & Wang (2013).       
		E_add = b0 / ( four * ( r_add - r_tp1 ) ) + b1 / four;
		E_m = b0 / ( four * ( r_m - r_tp1 ) ) + b1 / four;
		D_add = b0 / ( four * sq( r_tp1 - r_add ) );
		D_m = b0 / ( four * sq( r_tp1 - r_m ) );

		wp = one / ( r_tp1 - r_add );
		wm = one / ( r_tp1 - r_m );
		wbarp = b0 / four / sq( r_tp1 - r_add );
		wbarm = b0 / four / sq( r_tp1 - r_m );
		hp = b0 / four / (r_add - r_tp1 ) + b1 / four;
		hm = b0 / four / (r_m - r_tp1 ) + b1 / four;

		for ( int i = 1; i < 4; i++ )
			dd[i] = p->rdd[i];
		del = p->rdel;
		index_p4[1] = 0;
		cases_int = 1;

		PI0 = tmp.PI0;

		if ( r_cases == 1 ) {
			if ( p->f1234[1] >= zero ) {
				PI0_obs_inf = tmp.PI0_inf_obs;
				if ( pm < PI0_obs_inf ) {
					//equation (41) in Yang & Wang (2013).    
					tp = weierstrassP( pm + PI0, g2, g3, dd, del);
					r_coord = r_tp1 + b0 / ( four * tp - b1);
					pp = - pm;
				} else {
					tp = tinf;  // Goto infinity, far away. 
					r_coord = infinity;
					pp = - PI0_obs_inf;
				}
				t1 = 0;
				t2 = 0;
				sign_pr = one;
			} else { 
				if ( !indrhorizon ) {
					PI0_total = tmp.PI0_total;
					t2 = 0;
					if ( pm <= PI0 ) {
						t1 = 0;
						pp = pm;
						// equation (41) in Yang & Wang (2013).
						tp = weierstrassP( pm - PI0, g2, g3, dd, del );
						r_coord = r_tp1 + b0 / ( four * tp - b1 );
						sign_pr = - one;
					} else {
						t1 = 1;
						PI1_p = PI0;
						if ( pm < PI0_total ) {
							// equation (41) in Yang & Wang (2013).     
							tp = weierstrassP(pm - PI0, g2, g3, dd, del );
							r_coord = r_tp1 + b0 / ( four * tp - b1 );
							pp = two * PI0 - pm;
							p1 = fabs( pm - PI0 );
						} else {
							tp = tinf;   //Goto infinity, far away.
							r_coord = infinity;
							pp = - PI0_total + two * PI0;
							p1 = PI0_total - PI0;
						}
						sign_pr = one;
					}
				} else {
					// f1234r<0, photon will fall into black 
					// hole unless something encountered.
					PI0_obs_hori = tmp.PI0_obs_hori;
					if ( pm < PI0_obs_hori ) {
						// equation (41) in Yang & Wang (2013).   
		                            tp = weierstrassP(pm - PI0, g2, g3, dd, del );
		                            r_coord = r_tp1 + b0 / ( four * tp - b1 );
		                            pp = pm;
					} else {
		                            tp = thorizon; //Fall into black hole.
		                            r_coord = rhorizon;
		                            pp=PI0_obs_hori;
					}
					t1 = 0;
					t2 = 0;
					sign_pr = - one;
				}
			}
		} else if ( r_cases == 2 ) {
			if ( !indrhorizon ) {
				if( p->f1234[1] < zero )
					PI01 = - PI0;
				else
					PI01 = PI0;
				// equation (41) in Yang & Wang (2013).
				tp = weierstrassP( pm + PI01, g2, g3, dd, del);
				r_coord = r_tp1 + b0 / ( four * tp - b1 );

				half_periodwp = p->half_periodwp;
				periodwp = p->periodwp;

				weierstrass_int_J3( tobs, tp, dd, del, zero, one,
						index_p4, rff_p, integ4, cases_int);
				pp = integ4[1];

				// equation (57) in Yang & Wang (2013).
				p_tp1_tp2 = half_periodwp;
				PI2_p = p_tp1_tp2 - PI0;
				PI1_p = PI0;
				p1 = PI0 - pp;
				p2 = p_tp1_tp2 - p1;

				// equation (58) in Yang & Wang (2013).
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				// To determine t1 and t2, which are the times 
				// of the photon meets the two turnning 
				// points: r_tp1 and r_tp2 respectively.
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				Get_t1_t2( pm, p->f1234[1], robs_eq_rtp, p_tp1_tp2, 
						p->robs, r_tp1, r_tp2, &t1, &t2 );
				Get_t1_t2_old( pm, p->f1234[1], robs_eq_rtp, pp, 
						p1, p2, p->robs, r_tp1, &t1, &t2 );
			} else {  // photon has probability to fall into black hole.
				if ( p->f1234[1] <= zero ) {
					PI0_obs_hori = tmp.PI0_obs_hori;
					if ( pm < PI0_obs_hori ) {
						// equation (41) in Yang & Wang (2013).     
						tp = weierstrassP( pm - PI0, g2, g3, dd, del);
						r_coord = r_tp1 + b0 / ( four * tp - b1 );
						pp = pm;
					} else {
						tp = thorizon;   // Fall into black hole.
						r_coord = rhorizon;
						pp = PI0_obs_hori;
					}
		                        t1 = 0;
		                        t2 = 0;
		                        sign_pr = -one;
				} else {
					// p_r>0, photon will meet the r_tp2 &
					// turning point and turn around then goto vevnt horizon.     
					index_p4[1] = 0;
					cases_int = 1;
					weierstrass_int_J3( tp2, tobs, dd, del, a4, b4,
							index_p4, rff_p, integ04, cases_int );
					weierstrass_int_J3( tp2, thorizon, dd, del, a4, b4,
							index_p4, rff_p, integ14, cases_int );
		                        PI0_obs_tp2 = integ04[1];
		                        PI2_p = PI0_obs_tp2;
		                        PI0_total = integ14[1] + PI0_obs_tp2;

					if (pm <= PI0_obs_tp2) {
						t1 = 0;
						t2 = 0;
						pp = -pm;
						// equation (41) in Yang & Wang (2013).
						tp=weierstrassP( pm + PI0, g2, g3, dd, del );
						r_coord = r_tp1 + b0 / ( four * tp - b1 );
						sign_pr = one;
					} else {
						t1 = 0;
						t2 = 1;
						if (pm < PI0_total) {
							// equation (41) in Yang & Wang (2013).      
							tp = weierstrassP( pm + PI0, g2, g3, dd, del );
							r_coord = r_tp1 + b0 / ( four * tp - b1 );
							pp = pm - two * PI0_obs_tp2;
							p2 = pm - PI0_obs_tp2;
						} else {
							tp = thorizon;    // Fall into black hole. 
							r_coord = rhorizon;
							pp = PI0_total - two * PI0_obs_tp2;
							p2 = PI0_total - PI0_obs_tp2;
						}
						sign_pr = -one;
					}
				}
			}

		}
		Integral_r_part( p, pm, affr, timer, phir );
	} else {
		u = tmp.u;
		v = tmp.v;
		w = tmp.w;
		u2 = tmp.u2;
		v2 = tmp.v2;
		w2 = tmp.w2;
		double t_inf;
		double f1, g1, h1;
		double f2, g2, h2;
		double a5, b5;

		if ( u != zero ) {
			// equation (45) in Yang & Wang (2012).       
                        L1 = tmp.L1;
                        L2 = tmp.L2;
			// equation (46) in Yang & Wang (2012).       
                        thorizon = p->thorizon;
			// equation (48) in Yang & Wang (2012).   
                        m2 = tmp.m2;
                        tinf = p->tinf;
                        t_inf = p->t_inf;
			// equation (50) in Yang & Wang (2012).       
                        pinf = tmp.pinf;
			sqrt_L1 = tmp.sqrt_L1;
			//sncndn(p*w*dsqrt(L1)+dsign(one,f1234r)*&
                         //        pinf*w*dsqrt(L1),one-m2,sn,cn,dn)
			sncndn( pm * w * sqrt_L1 + sign( p->f1234[1] ) * 
				pinf * w * sqrt_L1, one - m2, &sn, &cn, &dn);
                        f1 = u2+w2;
                        g1 = -two*u;
                        h1 = one;
                        f2 = u2+v2;
                        g2 = -g1;
                        h2 = one;
                        a5 = zero;
                        b5 = one;
                        index_p5[1] = -1;
                        index_p5[2] = -2;
                        index_p5[3] = 2;
                        index_p5[4] = -4;
                        index_p5[5] = 4;

			if ( p->f1234[1] < zero ) {
				PI0 = pinf - EllipticF( thorizon, m2 ) / ( w * sqrt_L1 );
				if ( pm < PI0 ) {
					// equation (49) in Yang & Wang (2012).       
					y = u + ( -two * u + w * ( L1 - L2 ) * sn * fabs( cn ) )
						/ ( ( L1 - L2 ) * sq(sn) - ( L1 - one ) );
					r_coord = y;
					pp = pm;
				} else {
					y = rhorizon;
					r_coord = y;
					pp = PI0;
				}
				x = robs;
				sign_pr = - one;
			} else {
				PI0 = EllipticF( t_inf, m2 ) / ( w * sqrt_L1 ) - pinf;
				if ( pm < PI0 ) {
					// equation (49) in Yang & Wang (2012).       
		                        x = u + ( - two * u - w * ( L1 - L2 ) * sn * fabs( cn ) )
						/ ( ( L1 - L2 ) * sq(sn) - ( L1 - one ) );
		                        r_coord = x;
		                        pp = pm;
				} else {
		                        x = infinity;
		                        r_coord = x;
		                        pp = PI0;
				}
				y = robs;
				sign_pr = one;
			}
			//**********affine parameter part integration**********
			cases_int = 5;
			carlson_doublecomplex5( y, x, f1, g1, h1, f2, g2, 
				h2, a5, b5, index_p5, fabs(pp), integ5, cases_int );
			// equation (67) in Yang & Wang (2012).   
			*affr = integ5[5];
			// equation (68) in Yang & Wang (2012).   
			*timer = two * integ5[3] + four * pp + *affr;
			cases_int = 2;
			carlson_doublecomplex5( y, x, f1, g1, h1, f2, g2,
				h2, -r_add, b5, index_p5, fabs(pp), integ5, cases_int );
			carlson_doublecomplex5( y, x, f1, g1, h1, f2, g2,
				h2, -r_m, b5, index_p5, fabs(pp), integ15, cases_int);
			//******************** phi part ***********************
			// equation (68) in Yang & Wang (2012).   
			*timer = *timer + Ap * integ5[2] - Am * integ15[2];
			// equation (69) in Yang & Wang (2012).   
			*phir = a_spin * ( B_add * integ5[2] - B_m * integ15[2] );
			if ( muobs == zero && p->f1234[2] == zero ) {
			// equation (18) in Yang & Wang (2012).   
				*phir = *phir + pp * lambda;
			}
		} else {
			double rdius, pp_time;
			if ( p->f1234[1] < zero ) {
				PI0 = ( atan( robs / w ) - atan( rhorizon / w ) ) / w;
				if( pm < PI0 ) {
					rdius = w * tan( atan( robs / w ) - pm * w );
					r_coord = rdius;
				} else {
					rdius = rhorizon;
					r_coord = rdius;
				}
				//~~~~~~~~~~~~ timer part ~~~~~~~~~~~~~~~~~~~~~
				y = rdius;
				x = robs;
				sign_pr = - one;
			} else {
				PI0 = ( halfpi - atan( robs / w ) ) / w;
				if ( pm < PI0 ) {
					rdius = w * tan( atan( robs / w ) + pm * w );
					r_coord = rdius;
				} else {
					rdius = infinity;
					r_coord = rdius;
				}
				//~~~~~~~~~~~~ timer part ~~~~~~~~~~~~~~~~~~~~~
				y = robs;
				x = rdius;
				sign_pr = one;
			}
			pp_time = (x-y) + atan( x / w ) * ( -w + four / w - r_add * Ap / w
				/ ( w2 + r_add2 ) + r_m * Am / w / ( w2 + r_m2 ) ) -
				atan( y / w ) * ( -w + four / w - r_add * Ap / w
				/(w2 + r_add2 ) + r_m * Am / w / ( w2 + r_m2 ) );

			pp_time = pp_time + log( x*x + w2 ) * ( one - Ap / two
				/ ( w2 + r_add2 ) + Am / two / ( w2 + r_m2 ) )
				- ( log( y*y + w2 ) * ( one - Ap / two / ( w2 + r_add2 )
				+ Am / two / ( w2 + r_m2 ) ) );
			*timer = pp_time + Ap * log( fabs( x -r_add ) ) / ( w2 + r_add2 )
				- Am * log( fabs( x - r_m ) ) / ( w2 + r_m2 )
				- ( Ap * log( fabs( y - r_add ) ) / ( w2 + r_add2 )
				- Am * log( fabs( y - r_m ) ) / ( w2 + r_m2 ) );
			// affine parameter part **************************************
			*affr = ( x - y ) - w * atan( x / w ) + w * atan( y / w );
			// phy part ***************************************************
			if ( a_spin != zero ) {
				*phir = ( -B_add * r_add / w / ( r_add2 + w2 ) + B_m * r_m / w / ( r_m2 + w2 ) )
						* ( atan( x / w ) - atan( y / w ) ) +
				log( fabs( x - r_add ) / sqrt( x*x + w2 ) ) * B_add / ( r_add * two + w2 ) -
				log( fabs( y - r_add ) / sqrt( y*y + w2 ) ) * B_add / ( r_add * two + w2 ) -
				log( fabs( x - r_m ) / sqrt( x*x + w2 ) ) * B_m / ( r_m * two + w2 ) +
				log( fabs( y - r_m ) / sqrt( y*y + w2 ) ) * B_m / ( r_m * two + w2 );
				*phir = *phir * a_spin;
			} else
                            	*phir = zero;
      
			if ( muobs == zero && p->f1234[2] == zero )
				*phir = *phir + lambda * ( atan( x / w ) - atan( y / w ) ) / w;

		}
	}
}




/*
|*
|*
|*
|*
!*     C VERSION:  Yang Xiao-lin    2022-12-13.
|*/
static void Set_Initializations_For_Integrations_of_R_Part( ptcl * p )
{
	out_data2 tmp;
	Get_data( &tmp );

	robs = p->robs;
	lambda = p->lambda;
	a_spin = p->a_spin;
	muobs = p->muobs;
	sinobs = p->sinobs;

	rhorizon = p->rhorizon;
	// equation (64) in Yang & Wang (2013).
	r_add = rhorizon;
	r_m = one - sqrt( one - p->a2 );
	r_add2 = sq( r_add );
	r_m2 = sq( r_m );   
	// equation (64) in Yang & Wang (2013).  
	B_add = ( two * r_add - p->a_spin * p->lambda ) / ( r_add - r_m );
	B_m = ( two * r_m - p->a_spin * p->lambda) / ( r_add - r_m );
	// equation (64) in Yang & Wang (2013).
	Ap = ( r_add * ( four - p->a_spin * p->lambda )
		- two * p->a2 ) / sqrt( one - p->a2 );
	Am = ( r_m * ( four - p->a_spin * p->lambda ) 
		- two * p->a2 ) / sqrt( one - p->a2 );
	b4 = one;
	a4 = zero;
	//cc = p->a2 - p->lam2 - p->q;
	
	radiustp( p );
	r_tp1 = p->r_tp1;
	r_tp1sq = sq( r_tp1 );
	r_tp2 = p->r_tp2;
	r_reals = p->r_reals;
	robs_eq_rtp = p->robs_eq_rtp;
	indrhorizon = p->indrhorizon;
	r_cases = p->cases;

	for ( int i = 1; i < 5; i++ )
		bb[i] = p->rbb[i];

	PI1_phi = zero;
	PI2_phi = zero;
	PI1_time = zero;
	PI2_time = zero;
	PI1_aff = zero;
	PI2_aff = zero;

	if ( p->r_reals != 0 ) {
		/* equations (35)-(38) in Yang & Wang (2013). */
		b0 = p->rb0;
		b0sq = sq( b0 );
		b1 = p->rb1;
		b2 = p->rb2;
		b3 = p->rb3;
		g2 = p->rg2;
		g3 = p->rg3;

		tobs = p->tinf;

		thorizon = p->thorizon;

		tp2 = p->tp2;
		tinf = p->tinf1;
		h = - b1 / four;
		// equation (64), (66) and (70) in Yang & Wang (2013).       
		E_add = b0 / ( four * ( r_add - r_tp1 ) ) + b1 / four;
		E_m = b0 / ( four * ( r_m - r_tp1 ) ) + b1 / four;
		D_add = b0 / ( four * sq( r_tp1 - r_add ) );
		D_m = b0 / ( four * sq( r_tp1 - r_m ) );

		wp = one / ( r_tp1 - r_add );
		wm = one / ( r_tp1 - r_m );
		wbarp = b0 / four / sq( r_tp1 - r_add );
		wbarm = b0 / four / sq( r_tp1 - r_m );
		hp = b0 / four / (r_add - r_tp1 ) + b1 / four;
		hm = b0 / four / (r_m - r_tp1 ) + b1 / four;

		for ( int i = 1; i < 4; i++ )
			dd[i] = p->rdd[i];
		del = p->rdel;
		index_p4[1] = 0;
		cases_int = 1;

		PI0 = tmp.PI0;

		if ( r_cases == 1 ) {
			if ( p->f1234[1] >= zero ) {
				PI0_obs_inf = tmp.PI0_inf_obs;
			} else { 
				if ( !indrhorizon ) {
					PI0_total = tmp.PI0_total;
				} else {
					// f1234r<0, photon will fall into black 
					// hole unless something encountered.
					PI0_obs_hori = tmp.PI0_obs_hori;
				}
			}
		} else if ( r_cases == 2 ) {
			if ( !indrhorizon ) {
				if( p->f1234[1] < zero )
					PI01 = - PI0;
				else
					PI01 = PI0;
				// equation (41) in Yang & Wang (2013).
				half_periodwp = p->half_periodwp;
				periodwp = p->periodwp;

				weierstrass_int_J3( tobs, tp, dd, del, zero, one,
						index_p4, rff_p, integ4, cases_int);
				pp = integ4[1];

				// equation (57) in Yang & Wang (2013).
				p_tp1_tp2 = half_periodwp;
				PI2_p = p_tp1_tp2 - PI0;
				PI1_p = PI0;

				// equation (58) in Yang & Wang (2013).
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				// To determine t1 and t2, which are the times 
				// of the photon meets the two turnning 
				// points: r_tp1 and r_tp2 respectively.
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			} else {  // photon has probability to fall into black hole.
				if ( p->f1234[1] <= zero ) {
					PI0_obs_hori = tmp.PI0_obs_hori;
				} else {
					// p_r>0, photon will meet the r_tp2 &
					// turning point and turn around then goto vevnt horizon.     
					index_p4[1] = 0;
					cases_int = 1;
					weierstrass_int_J3( tp2, tobs, dd, del, a4, b4,
							index_p4, rff_p, integ04, cases_int );
					weierstrass_int_J3( tp2, thorizon, dd, del, a4, b4,
							index_p4, rff_p, integ14, cases_int );
		                        PI0_obs_tp2 = integ04[1];
		                        PI2_p = PI0_obs_tp2;
		                        PI0_total = integ14[1] + PI0_obs_tp2;
				}
			}

		}
		//Integral_r_part( p, pm, affr, timer, phir );
	} else {
		u = tmp.u;
		v = tmp.v;
		w = tmp.w;
		u2 = tmp.u2;
		v2 = tmp.v2;
		w2 = tmp.w2;

		if ( u != zero ) {
			// equation (45) in Yang & Wang (2012).       
                        L1 = tmp.L1;
                        L2 = tmp.L2;
			// equation (46) in Yang & Wang (2012).       
                        thorizon = p->thorizon;
			// equation (48) in Yang & Wang (2012).   
                        m2 = tmp.m2;
                        tinf = p->tinf;
                        t_inf = p->t_inf;
			// equation (50) in Yang & Wang (2012).       
                        pinf = tmp.pinf;
			sqrt_L1 = tmp.sqrt_L1;
			//sncndn(p*w*dsqrt(L1)+dsign(one,f1234r)*&
                        //	pinf*w*dsqrt(L1),one-m2,sn,cn,dn)
			//sncndn( pm * w * sqrt_L1 + sign( p->f1234[1] ) * 
			//	pinf * w * sqrt_L1, one - m2, &sn, &cn, &dn);
                        f1 = u2+w2;
                        g1 = -two*u;
                        h1 = one;
                        f2 = u2+v2;
                        g2 = -g1;
                        h2 = one;
                        a5 = zero;
                        b5 = one;
                        index_p5[1] = -1;
                        index_p5[2] = -2;
                        index_p5[3] = 2;
                        index_p5[4] = -4;
                        index_p5[5] = 4;

			if ( p->f1234[1] < zero ) {
				PI0 = pinf - EllipticF( thorizon, m2 ) / ( w * sqrt_L1 );
			} else {
				PI0 = EllipticF( t_inf, m2 ) / ( w * sqrt_L1 ) - pinf;
			}

		} else {
			if ( p->f1234[1] < zero ) {
				PI0 = ( atan( robs / w ) - atan( rhorizon / w ) ) / w;
			} else {
				PI0 = ( halfpi - atan( robs / w ) ) / w;
			}
		}
	}
}



/*
|*
|*
|*
|*
!*     C VERSION:  Yang Xiao-lin    2022-12-13.
|*/
static void Get_Results_of_Integrations_For_R_Part( ptcl * p, double pm, 
		double *r_coord, double *sign_pr, double *affr, 
		double *timer, double *phir )
{
	if ( p->r_reals != 0 ) {
		if ( r_cases == 1 ) {
			if ( p->f1234[1] >= zero ) {
				if ( pm < PI0_obs_inf ) {
					//equation (41) in Yang & Wang (2013).    
					tp = weierstrassP( pm + PI0, g2, g3, dd, del);
					*r_coord = r_tp1 + b0 / ( four * tp - b1);
					pp = - pm;
				} else {
					tp = tinf;  // Goto infinity, far away. 
					*r_coord = infinity;
					pp = - PI0_obs_inf;
				}
				t1 = 0;
				t2 = 0;
				*sign_pr = one;
			} else { 
				if ( !indrhorizon ) {
					t2 = 0;
					if ( pm <= PI0 ) {
						t1 = 0;
						pp = pm;
						// equation (41) in Yang & Wang (2013).
						tp = weierstrassP( pm - PI0, g2, g3, dd, del );
						*r_coord = r_tp1 + b0 / ( four * tp - b1 );
						*sign_pr = - one;
					} else {
						t1 = 1;
						PI1_p = PI0;
						if ( pm < PI0_total ) {
							// equation (41) in Yang & Wang (2013).     
							tp = weierstrassP(pm - PI0, g2, g3, dd, del );
							*r_coord = r_tp1 + b0 / ( four * tp - b1 );
							pp = two * PI0 - pm;
							p1 = fabs( pm - PI0 );
						} else {
							tp = tinf;   //Goto infinity, far away.
							*r_coord = infinity;
							pp = - PI0_total + two * PI0;
							p1 = PI0_total - PI0;
						}
						*sign_pr = one;
					}
					//printf("pm = %40.35f \n PI0 = %40.35f \n  tp = %40.35f \n", pm, PI0, tp);
				} else {
					// f1234r<0, photon will fall into black 
					// hole unless something encountered.

					if ( pm < PI0_obs_hori ) {
						// equation (41) in Yang & Wang (2013).   
		                            tp = weierstrassP(pm - PI0, g2, g3, dd, del );
		                            *r_coord = r_tp1 + b0 / ( four * tp - b1 );
		                            pp = pm;
					} else {
		                            tp = thorizon; //Fall into black hole.
		                            *r_coord = rhorizon;
		                            pp=PI0_obs_hori;
					}
					t1 = 0;
					t2 = 0;
					*sign_pr = - one;
				}
			}
		} else if ( r_cases == 2 ) {
			if ( !indrhorizon ) {
				// equation (41) in Yang & Wang (2013).
				tp = weierstrassP( pm + PI01, g2, g3, dd, del);
				*r_coord = r_tp1 + b0 / ( four * tp - b1 );

				weierstrass_int_J3( tobs, tp, dd, del, zero, one,
						index_p4, rff_p, integ4, cases_int);
				pp = integ4[1];
				// equation (57) in Yang & Wang (2013).
				p1 = PI0 - pp;
				p2 = p_tp1_tp2 - p1;

				// equation (58) in Yang & Wang (2013).
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				// To determine t1 and t2, which are the times 
				// of the photon meets the two turnning 
				// points: r_tp1 and r_tp2 respectively.
				//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				Get_t1_t2( pm, p->f1234[1], robs_eq_rtp, p_tp1_tp2, 
						p->robs, r_tp1, r_tp2, &t1, &t2 );
				/* Get_t1_t2_old( pm, p->f1234[1], robs_eq_rtp, pp, 
						p1, p2, p->robs, r_tp1, &t1, &t2 ); */
			} else {  // photon has probability to fall into black hole.
				if ( p->f1234[1] <= zero ) {
					if ( pm < PI0_obs_hori ) {
						// equation (41) in Yang & Wang (2013).     
						tp = weierstrassP( pm - PI0, g2, g3, dd, del);
						*r_coord = r_tp1 + b0 / ( four * tp - b1 );
						pp = pm;
					} else {
						tp = thorizon;   // Fall into black hole.
						*r_coord = rhorizon;
						pp = PI0_obs_hori;
					}
		                        t1 = 0;
		                        t2 = 0;
		                        *sign_pr = -one;
				} else {
					// p_r>0, photon will meet the r_tp2 &
					// turning point and turn around then goto vevnt horizon.
					if (pm <= PI0_obs_tp2) {
						t1 = 0;
						t2 = 0;
						pp = - pm;
						// equation (41) in Yang & Wang (2013).
						tp=weierstrassP( pm + PI0, g2, g3, dd, del );
						*r_coord = r_tp1 + b0 / ( four * tp - b1 );
						*sign_pr = one;
					} else {
						t1 = 0;
						t2 = 1;
						if (pm < PI0_total) {
							// equation (41) in Yang & Wang (2013).      
							tp = weierstrassP( pm + PI0, g2, g3, dd, del );
							*r_coord = r_tp1 + b0 / ( four * tp - b1 );
							pp = pm - two * PI0_obs_tp2;
							p2 = pm - PI0_obs_tp2;
						} else {
							tp = thorizon;    // Fall into black hole. 
							*r_coord = rhorizon;
							pp = PI0_total - two * PI0_obs_tp2;
							p2 = PI0_total - PI0_obs_tp2;
						}
						*sign_pr = -one;
					}
				}
			}

		}
		Integral_r_part( p, pm, affr, timer, phir );
	} else {

		//printf(" u = %f \n", u);
		if ( u != zero ) {

			sncndn( pm * w * sqrt_L1 + sign( p->f1234[1] ) * 
				pinf * w * sqrt_L1, one - m2, &sn, &cn, &dn);

			if ( p->f1234[1] < zero ) {
				if ( pm < PI0 ) {
					// equation (49) in Yang & Wang (2012).       
					y = u + ( -two * u + w * ( L1 - L2 ) * sn * fabs( cn ) )
						/ ( ( L1 - L2 ) * sq(sn) - ( L1 - one ) );
					*r_coord = y;
					pp = pm;
				} else {
					y = rhorizon;
					*r_coord = y;
					pp = PI0;
				}
				x = robs;
				*sign_pr = - one;
			} else {
				if ( pm < PI0 ) {
					// equation (49) in Yang & Wang (2012).       
		                        x = u + ( - two * u - w * ( L1 - L2 ) * sn * fabs( cn ) )
						/ ( ( L1 - L2 ) * sq(sn) - ( L1 - one ) );
		                        *r_coord = x;
		                        pp = pm;
				} else {
		                        x = infinity;
		                        *r_coord = x;
		                        pp = PI0;
				}
				y = robs;
				*sign_pr = one;
			}
			//**********affine parameter part integration**********
			cases_int = 5;
			carlson_doublecomplex5( y, x, f1, g1, h1, f2, g2, 
				h2, a5, b5, index_p5, fabs(pp), integ5, cases_int );
			// equation (67) in Yang & Wang (2012).   
			*affr = integ5[5];
			//printf(" x = %f y = %f int 5 = %f \n", x, y, integ5[5]);
			// equation (68) in Yang & Wang (2012).   
			*timer = two * integ5[3] + four * pp + *affr;
			//printf(" time1 = %f y = %f int 5 = %f \n", integ5[3], pp, integ5[5]);
			cases_int = 2;
			carlson_doublecomplex5( y, x, f1, g1, h1, f2, g2,
				h2, -r_add, b5, index_p5, fabs(pp), integ5, cases_int );
			carlson_doublecomplex5( y, x, f1, g1, h1, f2, g2,
				h2, -r_m, b5, index_p5, fabs(pp), integ15, cases_int);
			//******************** phi part ***********************
			// equation (68) in Yang & Wang (2012).   
			*timer = *timer + Ap * integ5[2] - Am * integ15[2];
			//printf(" time2 %f \n y = %f \n int 5 = %f \n = %f \n", Ap, integ5[2], Am, integ15[2]);
			//printf(" time3 = %f  = %f \n", -r_add, -r_m );
			// equation (69) in Yang & Wang (2012).   
			*phir = a_spin * ( B_add * integ5[2] - B_m * integ15[2] );
			if ( muobs == zero && p->f1234[2] == zero ) {
			// equation (18) in Yang & Wang (2012).   
				*phir = *phir + pp * lambda;
			}
		} else {
			double pp_time;
			if ( p->f1234[1] < zero ) {
				PI0 = ( atan( robs / w ) - atan( rhorizon / w ) ) / w;
				if( pm < PI0 ) {
					y = w * tan( atan( robs / w ) - pm * w );
					*r_coord = y;
				} else {
					y = rhorizon;
					*r_coord = y;
				}
				//~~~~~~~~~~~~ timer part ~~~~~~~~~~~~~~~~~~~~~
				x = robs;
				*sign_pr = - one;
			} else {
				PI0 = ( halfpi - atan( robs / w ) ) / w;
				if ( pm < PI0 ) {
					x = w * tan( atan( robs / w ) + pm * w );
					*r_coord = x;
				} else {
					x = infinity;
					*r_coord = x;
				}
				//~~~~~~~~~~~~ timer part ~~~~~~~~~~~~~~~~~~~~~
				y = robs;
				*sign_pr = one;
			}
			pp_time = (x-y) + atan( x / w ) * ( -w + four / w - r_add * Ap / w
				/ ( w2 + r_add2 ) + r_m * Am / w / ( w2 + r_m2 ) ) -
				atan( y / w ) * ( -w + four / w - r_add * Ap / w
				/(w2 + r_add2 ) + r_m * Am / w / ( w2 + r_m2 ) );

			pp_time = pp_time + log( x*x + w2 ) * ( one - Ap / two
				/ ( w2 + r_add2 ) + Am / two / ( w2 + r_m2 ) )
				- ( log( y*y + w2 ) * ( one - Ap / two / ( w2 + r_add2 )
				+ Am / two / ( w2 + r_m2 ) ) );
			*timer = pp_time + Ap * log( fabs( x -r_add ) ) / ( w2 + r_add2 )
				- Am * log( fabs( x - r_m ) ) / ( w2 + r_m2 )
				- ( Ap * log( fabs( y - r_add ) ) / ( w2 + r_add2 )
				- Am * log( fabs( y - r_m ) ) / ( w2 + r_m2 ) );
			// affine parameter part **************************************
			*affr = ( x - y ) - w * atan( x / w ) + w * atan( y / w );
			// phy part ***************************************************
			if ( a_spin != zero ) {
				*phir = ( -B_add * r_add / w / ( r_add2 + w2 ) + B_m * r_m / w / ( r_m2 + w2 ) )
						* ( atan( x / w ) - atan( y / w ) ) +
				log( fabs( x - r_add ) / sqrt( x*x + w2 ) ) * B_add / ( r_add * two + w2 ) -
				log( fabs( y - r_add ) / sqrt( y*y + w2 ) ) * B_add / ( r_add * two + w2 ) -
				log( fabs( x - r_m ) / sqrt( x*x + w2 ) ) * B_m / ( r_m * two + w2 ) +
				log( fabs( y - r_m ) / sqrt( y*y + w2 ) ) * B_m / ( r_m * two + w2 );
				*phir = *phir * a_spin;
			} else
                            	*phir = zero;
      
			if ( muobs == zero && p->f1234[2] == zero )
				*phir = *phir + lambda * ( atan( x / w ) - atan( y / w ) ) / w;

		}
	}
}

/*
!*
!*
!*
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-15.
!*
!*
!*/
static void Integral_r_part( ptcl * p, double pm, double *affr, double *timer, double *phir )
{
	double pp_aff, pp_phi, pp_time, time_temp;
	double p1_aff, p1_phi, p1_time;
	double p2_aff, p2_phi, p2_time;

	index_p4[1] = -1;
	index_p4[2] = -2;
	index_p4[3] = 0;
	index_p4[4] = -4;
	//**************** pp part ********************************************
	cases_int = 4;
	/*printf( "11111111111111111111111 tobs = %f  tp = %f  h = %f  b4 = %f \n", tobs, tp, h, b4 );
	printf( "11111111111111111111111 g2 = %f  g3 = %f pp = %f infs = %d  r_tp1 = %f \n", g2, g3, pp, indrhorizon, r_tp1 );
	printf("1111 = %d  fr = %f \n", r_cases, p->f1234[1]);
	printf(" pm = %f PI0 = %f PI0_tatol = %f \n", pm, PI0, PI0_total);*/
	weierstrass_int_J3( tobs, tp, dd, del, h, b4,
			index_p4, fabs(pp), integ4, cases_int );

	// equation (62) in Yang & Wang (2013).
	pp_aff = integ4[4] * b0sq / sixteen + integ4[2] * b0 * r_tp1 / two + pp * r_tp1sq;
	/*printf("pp_aff 1111 = %f  %f  %f  %f  %f  %f  %f \n", integ4[4], b0sq, integ4[2], b0, r_tp1, pp, r_tp1sq );
	printf(" tobs = %40.35f tp = %40.35f  pp = %40.35f df = %40.35f \n", tobs, tp, fabs(pp), tobs - tp);
	printf(" dd1  = %f  + %f * I  \n", creal( dd[1] ), cimag( dd[1] ) );
	printf(" dd2  = %f  + %f * I  \n", creal( dd[2] ), cimag( dd[2] ) );
	printf(" dd3  = %f  + %f * I  \n", creal( dd[3] ), cimag( dd[3] ) );
	exit(0);*/
	// equation (63) in Yang & Wang (2013).
	pp_time = integ4[2] * b0 / two + pp * ( two * r_tp1 + four + Ap * wp ) + pp_aff;
	time_temp = pp * ( - Am * wm );
	
	//printf("pp_time1 = %f  %f  %f  %f \n", integ4[2], b0, pp * ( two * r_tp1 + four + Ap * wp ), pp_aff );

	//printf( "phi3 = %f 2pi = %f affr = %d  2pi = %f affr = %f \n", tobs, tp, del, h, b4 );
	//printf( "t1 = %f \t t2 =  %d \n", pp, cases_int );
	//printf( " I1 = %f I2 = %f I3 = %f I4 = %f \n", integ4[1], integ4[2], integ4[3], integ4[4] );
        
	cases_int = 2;
	b4 = 1.0;
	weierstrass_int_J3( tobs, tp, dd, del, -E_add,
			b4, index_p4, fabs(pp), integ4, cases_int );
	// equation (63) in Yang & Wang (2013).
	pp_time = pp_time - Ap * wbarp * integ4[2];
	//printf("pp_time23 = %f  %f  %f  %f\n", Ap, wbarp, integ4[2], -E_add );
	//exit(0);
	if ( p->a_spin != zero ) {
		weierstrass_int_J3( tobs, tp, dd, del, -E_m,
			b4, index_p4, fabs(pp), integ14, cases_int );
		// equation (63) in Yang & Wang (2013).
		pp_time = pp_time + Am * wbarm * integ14[2] + time_temp;
		// equation (65) in Yang & Wang (2013).               
		pp_phi = pp * p->a_spin * ( B_add / ( r_tp1 - r_add ) - B_m / ( r_tp1 - r_m ) )
			- p->a_spin * B_add * D_add * integ4[2] + p->a_spin * B_m * D_m * integ14[2];
	} else
		pp_phi=zero;

              
	if ( p->muobs == zero && p->f1234[2] == zero )
		// equation (18) in Yang & Wang (2013).
		pp_phi = pp_phi + pp * p->lambda;

	// ******************** p1 part****************************************
	if (t1 == 0) {
		p1_phi = zero;
		p1_time = zero;
		p1_aff = zero;
	} else {
		if ( PI1_aff == zero && PI1_time == zero ) {
			cases_int = 4;
			weierstrass_int_J3( tobs, infinity, dd, del, 
				h, b4, index_p4, PI0, integ4, cases_int );
			// equation (62) in Yang & Wang (2013).
			PI1_aff = integ4[4] * b0sq / sixteen + integ4[2] *
				b0 * r_tp1 / two + PI0 * r_tp1sq;
 
			//printf("PI1_aff = %f  %f  %f  %f  %f  %f  %f \n", integ4[4], 
			//		b0sq, integ4[2], b0, r_tp1, PI0, r_tp1sq );
			// equation (63) in Yang & Wang (2013).     
			PI1_time = integ4[2] * b0 / two + PI0 * 
				( two * r_tp1 + four + Ap * wp ) + PI1_aff;
			time_temp = PI0 * ( - Am * wm );
        
			cases_int = 2;
			weierstrass_int_J3( tobs, infinity, dd, del, -E_add,
				b4, index_p4, PI0, integ4, cases_int );
			// equation (63) in Yang & Wang (2013).
			PI1_time = PI1_time - Ap * wbarp * integ4[2];
			if ( p->a_spin != zero) {
				weierstrass_int_J3( tobs, infinity, dd, del, -E_m,
					b4, index_p4, PI0, integ14, cases_int );
				// equation (63) in Yang & Wang (2013).       
				PI1_time = PI1_time + Am * wbarm * integ14[2] + time_temp;
				// equation (65) in Yang & Wang (2013).
				PI1_phi = PI0 * p->a_spin * ( B_add / ( r_tp1 - r_add ) - B_m / ( r_tp1 - r_m ) )
					- p->a_spin * B_add * D_add * integ4[2] + p->a_spin * B_m * D_m * integ14[2];
			} else
                                PI1_phi = zero;

			if ( p->muobs == zero && p->f1234[2] == zero ) {
				// equation (18) in Yang & Wang (2013).       
				PI1_phi += PI0 * lambda;
			}
		}
		// equation (55) in Yang & Wang (2013).       
		p1_aff = PI1_aff - pp_aff;
		p1_time = PI1_time - pp_time;
		p1_phi = PI1_phi - pp_phi;
	}

	// ******************** p2 part****************************************
	if ( t2 == zero ) {
                        p2_phi = zero;
                        p2_time = zero;
                        p2_aff = zero;
	} else {
		if ( PI2_aff == zero && PI2_time == zero ) {
			cases_int = 4;
			weierstrass_int_J3( tp2, tobs, dd, del, h, b4,
					index_p4, PI2_p, integ4, cases_int );
			// equation (62) in Yang & Wang (2013).       
			PI2_aff = integ4[4] * b0*b0 / sixteen + integ4[2] *
					b0 * r_tp1 / two + PI2_p * r_tp1sq;
			// equation (63) in Yang & Wang (2013).       
			PI2_time = integ4[2] * b0 / two + PI2_p * 
				( two * r_tp1 + four + Ap * wp ) + PI2_aff;
			time_temp = PI2_p * ( - Am * wm );

			cases_int = 2;
			weierstrass_int_J3( tp2, tobs, dd, del, -E_add,
					b4, index_p4, PI2_p, integ4, cases_int );
			// equation (63) in Yang & Wang (2013).
			PI2_time = PI2_time - Ap * wbarp * integ4[2];
			if ( p->a_spin != zero) {
				weierstrass_int_J3( tp2, tobs, dd, del, -E_m, 
					b4, index_p4, PI2_p, integ14, cases_int );
				// equation (63) in Yang & Wang (2013).       
				PI2_time = PI2_time + Am * wbarm * integ14[2] + time_temp;
				// equation (65) in Yang & Wang (2013).
				PI2_phi = PI2_p * p->a_spin * ( B_add / ( r_tp1 - r_add ) - B_m / ( r_tp1 - r_m ) )
					- p->a_spin * B_add * D_add * integ4[2] + p->a_spin * B_m * D_m * integ14[2];
			} else {
				PI2_phi = zero;
			}
			if ( p->muobs == zero && p->f1234[2] == zero ) {
				// equation (18) in Yang & Wang (2013).       
				PI2_phi = PI2_phi + PI2_p * p->lambda;
			}
		}
		// equation (55) in Yang & Wang (2013).       
		p2_aff = PI2_aff + pp_aff;
		p2_time = PI2_time + pp_time;
		p2_phi = PI2_phi + pp_phi;
	}

	//***************** phi, aff,time part ********************************
	// equation (56) in Yang & Wang (2013).       
	if ( p->f1234[1] != zero ) {
		*phir = sign( -p->f1234[1] ) * pp_phi + two * (t1 * p1_phi + t2 * p2_phi );
		*timer = fabs( sign( -p->f1234[1] ) * pp_time + two * (t1 * p1_time + t2 * p2_time ) );
		*affr = sign( -p->f1234[1] ) * pp_aff + two * (t1 * p1_aff + t2 * p2_aff );
		//printf("pp_time = %f  %f  %f \n t1 = %d \n t2 = %d \n", pp_time, p1_time, p2_time, t1, t2 );
		//printf("pp_aff = %f  %f  %f  \n", pp_aff, p1_aff, p2_aff );
	} else {
		if ( robs == r_tp1 ) {
			*phir = -pp_phi + two * ( t1 * p1_phi + t2 * p2_phi );
			*timer = -pp_time + two * ( t1 * p1_time + t2 * p2_time );
			*affr = -pp_aff + two * ( t1 * p1_aff + t2 * p2_aff );
		} else {
			*phir = pp_phi + two * ( t1 * p1_phi + t2 * p2_phi );
			*timer = pp_time + two * (t1 * p1_time + t2 * p2_time );
			*affr = pp_aff + two * ( t1 * p1_aff + t2 * p2_aff );
		}    
	}
}



/*
!*     C VERSION:  Yang Xiao-lin    2022-12-10.
!*/
void Get_R_Integral_Situations( ptcl *p, double pem )
{
	if ( p->r_reals != 0 ) { 
		if ( p->cases == 1 ) {
			if ( !p->indrhorizon && ( p->f1234[1] < zero ) ) {
				if ( pem < PI0_total ) {
					p->Situations = 11;
					if ( pem < PI0 )
						p->sign_pr = -one;
					else
						p->sign_pr = one;
				} else {
					p->Situations = 12;
					p->sign_pr = one;
				}
			} else if ( !p->indrhorizon && ( p->f1234[1] > zero ) ) {
				p->sign_pr = one;
				if ( pem < PI0_inf_obs ) 
					/* equation (41) in Yang & Wang (2013). */
					p->Situations = 13;
				else
					p->Situations = 14;  /*Goto infinity, far away.*/
			} else if ( p->indrhorizon && ( p->f1234[1] < zero ) ) {
				p->sign_pr = - one;
				if ( pem < PI0_obs_hori ) {
					/* equation (41) in Yang & Wang (2013). */
					p->Situations = 15;
				} else
					p->Situations = 16;  /*Fall into black hole. */
			} else if ( p->indrhorizon && ( p->f1234[1] > zero ) ) {
				p->sign_pr = one;
				if ( pem < PI0_inf_obs )
					/* equation (41) in Yang & Wang (2013). */
					p->Situations = 17;
				else
					p->Situations = 18; /*Goto infinity, far away.*/
			}
		} else if ( p->cases == 2 ) {   //  case(2)
			double up;
			if ( !p->indrhorizon && p->f1234[1] > zero ) {

				up = fmod(pem + PI01, p->periodwp);
				if ( up <= p->half_periodwp ) {
					p->sign_pr = one;
					p->Situations = 21;
				} else {
					p->sign_pr = - one;
					p->Situations = 22;
				}

			} else if ( !p->indrhorizon && p->f1234[1] < zero ) {

				up = fmod(pem + PI01 + p->half_periodwp, p->periodwp);
				if ( up < p->half_periodwp ) {
					p->sign_pr = -one;
					p->Situations = 23;
				} else {
					p->sign_pr = one;
					p->Situations = 24;
				}

			} else if ( !p->indrhorizon && p->robs == p->r_tp1 ) {

				if ( fmod(pem, p->periodwp) <= p->half_periodwp ) {
					p->sign_pr = one;
					p->Situations = 25;
				} else {
					p->sign_pr = -one;
					p->Situations = 26;
				}

			} else if ( !p->indrhorizon && p->robs == p->r_tp2 ) {

				if ( fmod(pem + p->half_periodwp, p->periodwp)
							<= p->half_periodwp ) {
					p->sign_pr = one;
					p->Situations = 27;
				} else {
					p->sign_pr = -one;
					p->Situations = 28;
				}

			} else if ( p->indrhorizon && p->f1234[1] <= zero ) {

				p->sign_pr = - one;
				if ( pem < PI0_obs_hori )
					p->Situations = 29;
				else
					p->Situations = 30;

			} else if ( p->indrhorizon && p->f1234[1] > zero ) {

				if ( pem < PI0_total_2 )
					p->Situations = 31;
				else
					p->Situations = 32;

				if (pem <= p->p_ini_r_tp2) {
					p->sign_pr = one;
					p->Situations = 33;
				} else {
					p->sign_pr = -one;
					p->Situations = 34;
				}

			}
		} // end casese.
	} else {
		if ( u != zero && p->f1234[1] < zero ) {

			p->sign_pr = - one;
			if ( pem < PI0 )
				p->Situations = 41;
			else
				p->Situations = 42;

		} else if ( u != zero && p->f1234[1] > zero ) {

			p->sign_pr = one;
			if ( pem < PI0 )
				p->Situations = 43;
			else
				p->Situations = 44;

		} else if ( u == zero && p->f1234[1] < zero ) {

			p->sign_pr = - one;
			if ( pem < PI0 ) 
				p->Situations = 45;
			else
				p->Situations = 46;

		} else if ( u == zero && p->f1234[1] > zero ) {

			p->sign_pr = one;
			if( pem < PI0 ) 
				p->Situations = 47;
			else
				p->Situations = 48;

		}
	}
}





static double a_spin_1;
static double scal_1;
static double lambda_1;
static double q_1;
static double f12341_1;
static double f12342_1;
static double muobs_1;
static double sinobs_1;


/*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-17.
!*
!*/
void Get_Integrals_For_R_Part( ptcl *pt, double p, 
		double *r_coord, double *sign_pr, double *affr, 
		double *timer, double *phir )
{
	static unsigned int count_num = 1;
 
rint_restarting:
	if ( count_num == 1 ) {

		radius_settings( pt );
		Set_Initializations_For_Integrations_of_R_Part( pt );
		count_num += 1;

		a_spin_1 = pt->a_spin;
		scal_1 = pt->scal;
		lambda_1 = pt->lambda;
		q_1 = pt->q;
		f12341_1 = pt->f1234[1];
		f12342_1 = pt->f1234[2];
		muobs_1 = pt->muobs;
		sinobs_1 = pt->sinobs;
 
		Get_Results_of_Integrations_For_R_Part( pt, p, 
			r_coord, sign_pr, affr, timer, phir );

		//printf("here332 %d \n ", count_num);
		//printf( "phir = %f timer = %f affr = %f \n", *phir, *timer, *affr );

    	} else {
		if ( f12341_1 == pt->f1234[1]
		&& f12342_1 == pt->f1234[2]
		&& lambda_1 == pt->lambda
		&& q_1 == pt->q
		&& sinobs_1 == pt->sinobs
		&& muobs_1 == pt->muobs
		&& a_spin_1 == pt->a_spin
		&& scal_1 == pt->scal ) {
 
			Get_Results_of_Integrations_For_R_Part( pt, p, 
				r_coord, sign_pr, affr, timer, phir );
		} else {
			count_num = 1;
			goto rint_restarting;
		}

	}
}




 

