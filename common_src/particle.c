/*
 * common_src/particle.c
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


#include "particle.h"
 




// Memory allocator
ptcl* particle_new() {
	return (ptcl *)malloc( sizeof( ptcl ) );
}


void particle_construct( ptcl *this, double a_spin, double rini, 
		double mucos, double sinobs, double scal, double *Vobs )
{
	this->a_spin = a_spin;
	this->a2 = sq( this->a_spin );
	this->rhorizon = one + sqrt( one - this->a2 );

	this->robs = rini;
	this->r2 = sq( this->robs );
	this->muobs = mucos;
	this->mu2 = sq( this->muobs );
	this->sinobs = sinobs;
	this->sin2 = sq( this->sinobs );
	this->scal = scal;

	this->velocity_ini[1] = Vobs[0];
	this->velocity_ini[2] = Vobs[1];
	this->velocity_ini[3] = Vobs[2];
}


void Set_alphabeta( ptcl *this, double alpha, double beta )
{
	this->alpha = alpha;
	this->beta = beta;
}

void Set_alpha( ptcl *this, double alpha )
{
	this->alpha = alpha;
}

void Set_beta( ptcl *this, double beta )
{
	this->beta = beta;
}







/*
!*     PURPOSE:  Computes constants of motion from impact parameters alpha and beta by using 
!*               formulae (86) and (87) in Yang & Wang (2012).    
!*     INPUTS:   alpha,beta-----Impact parameters.
!*               robs-----------radial coordinate of observer or the initial position of photon.
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               velocity(1:3)--Array of physical velocities of the observer or emitter with respect to
!*                              LNRF.        
!*     OUTPUTS:  f1234(1:4)-----array of p_r, p_theta, p_phi, p_t, which are the components of 
!*                              four momentum of a photon measured under the LNRF frame, and 
!*                              defined by equations (82)-(85) in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.            
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  5 Jan 2012.
!*     REVISIONS: 
!*     C VERSION:    Yang Xiao-lin   2022-11-23. Kumming, Yunnan Observatory.
!*/
void lambdaq( ptcl *th )
{
	double A1, Delta, Sigma, At,Bt, Vr, Vt, Vp, 
		prt, ptt, ppt, /* expnu2, eppsi2, epmu12,epmu22, */
		bigA, gama, gama2, somiga;
 

	if ( fabs(th->beta) < 1.e-7  )th->beta = zero;
	if ( fabs(th->alpha) < 1.e-7 )th->alpha = zero;
	/* equations (94), (95) in Yang & Wang (2012). */
	At = th->alpha / th->scal / th->robs;
	Bt = th->beta / th->scal / th->robs;        
	Vr = th->velocity_ini[1];
	Vt = th->velocity_ini[2];
	Vp = th->velocity_ini[3];
	/* equation (90) in Yang & Wang (2012). */ 
	gama = one / sqrt( one - ( Vr*Vr + Vt*Vt + Vp*Vp ) );
	gama2 = pow( gama, 2);

	/* equations (97), (98), (99) in Yang & Wang (2012). */
	prt = - one / sqrt( one + At*At + Bt*Bt );
	ptt = Bt*prt;
	ppt = At*prt;
	/* equations (89), (90) and (91) in Yang & Wang (2012). */
	th->f1234[1] = ( gama * Vr - prt * ( one + gama2 * Vr * Vr / 
			( one + gama ) ) - ptt * gama2 * Vr * Vt / 
			( one + gama ) - ppt * gama2 * Vr * Vp / 
			( one + gama ) ) * th->robs * th->scal;
	th->f1234[2] = ( gama * Vt - prt * gama2 * Vt * Vr / ( one + gama) - 
			ptt * ( one + gama2 * Vt * Vt / ( one + gama ) ) - 
			ppt * gama2 * Vt * Vp / ( one + gama ) ) * 
			th->robs * th->scal;
	th->f1234[3] = ( gama * Vp - prt * gama2 * Vp * Vr / (one+gama) - 
			ptt * gama2 * Vp * Vt / ( one + gama ) - 
			ppt * ( one + gama2 * Vp * Vp / ( one + gama ) ) ) * 
			th->robs * th->scal;
	th->f1234[4] = gama * ( one - prt * Vr - ptt * Vt - ppt * Vp );

	/* Keep r component p_r of four momentum to be negative, 
		so the photon will go to the central black hole. */
	th->f1234[1] = - th->f1234[1];  /* making the r component of 
					velocity to be negative. */
	for (int i = 1; i < 5; i++) {
		if ( fabs( th->f1234[i] ) < 1.e-6 ) th->f1234[i] = zero;
	}

	/* equations (1), (2) in Yang & Wang (2012). */


	Delta = th->r2 - two*th->robs + th->a2;
	Sigma = th->r2 + pow((th->a_spin*th->muobs), 2);
	bigA = pow((th->r2 + th->a2), 2) - pow((th->a_spin*th->sinobs), 2) * Delta;
	somiga = two * th->a_spin * th->robs / bigA;

	//expnu2 = Sigma * Delta / bigA;
	//eppsi2 = pow(th->sinobs, 2) * bigA / Sigma;
	//epmu12 = Sigma / Delta;
	//epmu22 = Sigma;

	/* equations (86) and (87) in Yang & Wang (2012). */
	A1 = th->f1234[3] / ( sqrt(Delta) * Sigma / bigA * th->f1234[4] * 
		th->robs * th->scal + th->f1234[3] * somiga * th->sinobs );
	th->lambda = A1 * th->sinobs;
	th->lam2 = sq( th->lambda );
	th->q = (A1 * A1 - th->a2 ) * th->muobs * th->muobs + 
		pow( (th->f1234[2] / th->f1234[4] / th->robs / th->scal * 
			( one-th->lambda * somiga ) ), 2) * bigA / Delta;
}


void show_lambdaq( ptcl *th )
{
	printf("lambda = %20.16f \n q = %20.16f     \n", th->lambda, th->q);
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
}



/*
!*     PURPOSE:  Computes constants of motion from components of initial 4 momentum 
!*               of photon measured by emitter in its local rest frame, by using 
!*               formulae (86) and (87) in Yang & Wang (2012).    
!*     INPUTS:   pr,ptheta,pphi-----components of initial 4 momentum of photon measured by 
!*               emitter in its local rest frame.
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or the initial position of photon.
!*               velocity(1:3)--Array of physical velocities of the observer or emitter with respect to
!*                              LNRF.        
!*     OUTPUTS:  f1234(1:4)-----array of p_r, p_theta, p_phi, p_t, which are the components of 
!*                              four momentum of a photon measured under the LNRF frame, and 
!*                              defined by equations (82)-(85) in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.            
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  5 Jan 2012.
!*     REVISIONS:
!*     C VERSION:    Yang Xiao-lin   2022-11-23. Kumming, Yunnan Observatory.
!*/
void ini_direction2lamdaq( ptcl *this, double pr, double ptheta, double pphi )
{
	double Vr, Vt, Vp;
	double gama, gama_tp, gama_p, Delta, Sigma, bigA, somiga, A1;
 
	Vr = this->velocity_ini[1];
	Vt = this->velocity_ini[2];
	Vp = this->velocity_ini[3];

	/* equation (92) in Yang & Wang (2012). */
	gama = one / sqrt( one - ( Vr*Vr + Vt*Vt + Vp*Vp ) );
	gama_tp = one / sqrt( one - ( Vt*Vt + Vp*Vp ) );
	gama_p = one / sqrt( one - Vp*Vp );
	/* equation (106)-(109) in Yang & Wang (2012). */
	this->f1234[1] = - ( - gama * Vr + gama / gama_tp * pr );
	this->f1234[2] = - ( - gama * Vt + gama * gama_tp * Vr * Vt * pr + 
				gama_tp / gama_p * ptheta );
	this->f1234[3] = - ( - gama * Vp + gama * gama_tp * Vr * Vp * pr + 
			gama_tp * gama_p * Vt * Vp * ptheta + gama_p * pphi );
	this->f1234[4] = ( gama - gama * gama_tp * Vr * pr - 
		gama_tp * gama_p * Vt * ptheta - gama_p * Vp * pphi );

	for (int i = 1; i < 5; i++)
		if ( fabs( this->f1234[i] ) < 1.e-7)
			this->f1234[i] = zero;

	/* equations (1), (2) in Yang & Wang (2012). */
        Delta = this->r2 - two * this->robs + this->a2;
        Sigma = this->r2 + this->a2 * this->mu2;
        bigA = pow( ( this->r2 + this->a2 ), 2) - this->a2 * this->sin2 *Delta;
        somiga = two * this->a_spin * this->robs / bigA;
        //eppsi2 = this->sin2 * bigA / Sigma;
        //epmu12 = Sigma / Delta;
        //epmu22 = Sigma;
	/* equations (86) and (87) in Yang & Wang (2012). */
	//printf("he1 = %f  he1 = %f  he1 = %f  he1 = %f  \n", Delta, this->r2, this->robs, this->a2);
        A1 = this->f1234[3] / ( sqrt(Delta) * Sigma / bigA * this->f1234[4] + 
			this->f1234[3] * somiga * this->sinobs );
        this->lambda = A1 * this->sinobs;
        this->lam2 = sq(this->lambda);
        this->q = ( A1 * A1 - this->a2 ) * this->mu2 + 
		pow( ( this->f1234[2] / this->f1234[4] * 
                   ( one - this->lambda * somiga ) ), 2) * bigA / Delta;
}



/*
!*     PURPOSE:  Computes Kerr metric, exp^\nu, exp^\psi, exp^mu1, exp^\mu2, and omiga at position:
!*               r_obs, \theta_{obs}.     
!*     INPUTS:   robs-----------radial coordinate of observer or the initial position of photon. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).         
!*     OUTPUTS:  somiga,expnu,exppsi,expmu1,expmu2------------Kerr metrics under Boyer-Lindquist coordinates.
!*     ROUTINES CALLED: root3, weierstrass_int_J3, radiustp, weierstrassP, EllipticF, carlson_doublecomplex5 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  5 Jan 2012.
!*     REVISIONS: 
!*     C VERSION:  Yang Xiao-lin    2022-11-23.
!*/
void metricg( ptcl *this )
{
	double sigma, Delta, bigA;      

	/* equations (1) and (2) in Yang & Wang (2012). */
	Delta = this->r2 - two * this->robs + this->a2;
	sigma = this->r2 + this->a2 * this->mu2;
	bigA = pow( ( this->r2 + this->a2 ), 2 ) - this->a2 * this->sin2 * Delta;
	this->mt.somiga = two * this->a_spin * this->robs / bigA;

	this->mt.expnu2 = sigma * Delta / bigA;
	this->mt.exppsi2 = this->sin2 * bigA / sigma;
	this->mt.expmu12 = sigma / Delta;
	this->mt.expmu22 = sigma;

	this->mt.expnu = sqrt( this->mt.expnu2 );
	this->mt.exppsi = sqrt( this->mt.exppsi2 );
	this->mt.expmu1 = sqrt( this->mt.expmu12 );
	this->mt.expmu2 = sqrt( this->mt.expmu22 );
}



int metricgij( double robs, double muobs, double sinobs, double a_spin, 
	double *somiga, double *expnu, double *exppsi, double *expmu1, double *expmu2 )
{
	double Delta, bigA, sigma, robs2, a_spin2;       

	robs2 = robs * robs;
	a_spin2 = a_spin * a_spin;

        Delta = robs2 - two * robs + a_spin2;
        sigma = robs2 + sq(a_spin * muobs);
        bigA = sq( robs2 + a_spin2 ) - sq( a_spin * sinobs ) * Delta;
        *somiga = two * a_spin * robs / bigA;
        *expnu = sqrt( sigma * Delta / bigA );
        *exppsi = sinobs * sqrt( bigA / sigma );
        *expmu1 = sqrt( sigma / Delta );
        *expmu2 = sqrt( sigma );
        return 0;	
}


/*
!*     PURPOSE:  Solves equations f_3(alphac,betac)=0, f_2(alphac,betac)=0, of (100)  
!*               and (101) in Yang & Wang (2012). alphac, betac are the coordinates of 
!*               center point of images on the screen of observer.    
!*     INPUTS:   robs-----------radial coordinate of observer or the initial position of photon.
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.  
!*               velocity(1:3)--Array of physical velocity of observer or emitter with respect to
!*                              LNRF.        
!*     OUTPUTS:  alphac,betac---coordinates of center point of images on the screen of observer.            
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  9 Jan 2012.
!*     REVISIONS: 
!*     C VERSION:  Yang Xiao-lin    2022-11-23.
!*/
void center_of_image( ptcl *this )
{
	double Vr, Vt, Vp, gama, a1, b1, c1, alphap, alpham, betap, betam;
	double Vt2, Vp2, gama2;

        Vr = this->velocity_ini[1];
        Vt = this->velocity_ini[2];
        Vp = this->velocity_ini[3];
	Vt2 = Vt * Vt;
	Vp2 = Vp * Vp;
	/* equation (90) in Yang & Wang (2012). */
	gama = one / sqrt( one - ( Vr*Vr + Vt2 + Vp2 ) );
	gama2 = gama * gama;

	if ( Vt != zero ) {
		if ( Vp != zero ) {
			a1 = pow( ( one + gama * gama * ( Vt2 + Vp2 ) / 
				( one + gama ) ), 2) - gama2 * 
				( Vp2 + Vt2 );  
			b1 = two * gama2 * Vt * Vr * ( one + gama + gama2 * 
				( Vt * Vt + Vp2 ) ) / pow(( one + gama ), 2);
			c1 = pow( ( gama2 * Vt * Vr / ( one + gama ) ), 2) - gama2 * Vt2;
			betap = ( - b1 + sqrt( b1*b1 - four*a1*c1 ) ) / two / a1;
			betam = ( - b1 - sqrt( b1*b1 - four*a1*c1 ) ) / two / a1;                         
			if ( betap*Vp < zero )
				this->betac = betap;
			else
				this->betac = betam;
			this->alphac = Vp / Vt * this->betac;
		} else {
			this->alphac = zero;
			a1 = pow( (one + gama2 * Vt2 / (one + gama) ), 2) - gama2 * Vt2;
			b1 = two * gama2 * Vt * Vr * ( one + gama + gama2 * Vt2 ) / 
				pow( (one + gama ), 2);
			c1 = pow( ( gama2 * Vt * Vr / ( one + gama ) ), 2) - gama2 * Vt2;
			betap = (-b1 + sqrt( b1*b1 - four*a1*c1) ) / two / a1;
			betam = (-b1 - sqrt( b1*b1 - four*a1*c1) ) / two / a1;        
			if ( betap * Vt < zero )
				this->betac = betap;
			else
				this->betac = betam;
		}  
	} else {
		this->betac = zero;
		if ( Vp != zero ) {
			a1 = pow( (one + gama2 * Vp2 / ( one + gama ) ), 2) - 
				gama2 * Vp2;
			b1 = two * gama2 * Vp * Vr * ( one + gama + 
				gama2 * Vp2 ) / pow( ( one + gama ), 2);
			c1 = pow( ( gama2 * Vp * Vr / ( one + gama ) ), 2) - 
				gama2 * Vp2;
			alphap = ( -b1 + sqrt( b1*b1 - four*a1*c1 ) ) / two / a1;
			alpham = ( -b1 - sqrt( b1*b1 - four*a1*c1 ) ) / two / a1;
			if ( alphap * Vp < zero )
				this->alphac = alphap;
		        else
				this->alphac = alpham;
		} else
		        this->alphac=zero;
	}
        this->alphac = this->alphac * this->robs * this->scal;
        this->betac = this->betac * this->robs * this->scal;
}




void transport_data_out( ptcl *this, out_data * od )
{
	od->lambda = this->lambda; 
	od->lam2 = this->lam2; 
	od->q = this->q;
 
	od->alphac = this->alphac; 
	od->betac = this->betac;
 
	for ( int i = 1; i < 5; i++)
		od->f1234[i] = this->f1234[i];

	od->mt.expnu2 = this->mt.expnu2;
	od->mt.exppsi2 = this->mt.exppsi2;
	od->mt.expmu12 = this->mt.expmu12;
	od->mt.expmu22 = this->mt.expmu22;

	od->mt.expnu = this->mt.expnu;
	od->mt.exppsi = this->mt.exppsi;
	od->mt.expmu1 = this->mt.expmu1;
	od->mt.expmu2 = this->mt.expmu2;
}  



/*
|*
|*
|*
|*
|*/ 



/*
!*     PURPOSE: Computes inner most stable circular orbit r_{ms}. 
!*     INPUTS:   a_spin ---- Spin of black hole, on interval [-1,1].
!*     OUTPUTS:  radius of inner most stable circular orbit: r_{ms}
!*     ROUTINES CALLED: root4
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: 
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-21.
!*
!*/
double rms( double a_spin )
{
	double b, c, d, e;
	double _Complex rt[5];
	unsigned int reals;
	if ( a_spin == 0.0 ) {
		return 6.0;
	}

	b = 0.0;
	c = -6.0;
	d = 8.0 * a_spin;
	e = -3.0 * pow( a_spin, 2 );
	// Bardeen et al. (1972) 
	root4( b, c, d, e, rt, &reals);
	for ( int i = 4; i >= 1; i-- ) {
		if (cimag( rt[i] ) == 0.0 )
			return pow( creal( rt[i] ), 2 );
	}
	return zero;
}









 

