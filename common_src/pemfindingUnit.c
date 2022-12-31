/*
 * common_src/pemfindingUnit.c
 *
!********************************************************************************************
!      module pemfinding
!*******************************************************************************
!*     PURPOSE: This module aims on solving more general equation f(p)=0. For 
!*              detail definition of f(p), cf. Yang & Wang (2013).     
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2013).  
!*     DATE WRITTEN:  15 Jan 2013. 
!*******************************************************************************    
 *
 * Author       Yang Xiao-lin
 *
 * Date         2022-12-27
 *
 * City         Kun-ming Yunnan Provice, China
 *
 * 2022-11-16   Finished the fisrt writting of the module.
 */

#include "pemfindingUnit.h"


static double rin, rout;
static double muup, mudown;
static double phy1, phy2;

//static double oricosth;
static int bisection, caserange;
static int NN;


static void pemfindcase2( ptcl *p, double *pemfind );
static void pemfindcase3( ptcl *p, double *pemfind );
static void pemfindcase4( ptcl *p, double *pemfind );
static void pemfindcase5( ptcl *p, double *pemfind );
static double Sectionp( ptcl *pt, double p1, double p2 );
static double rootfind( ptcl *p, double p1, double p2 );
static double Bisectionp( ptcl *pt, double p1, double p2 );
static double NewRapson( ptcl *pt, double p1, double p2 );

extern double Fp( ptcl *pt, double p );



void set_parameters( double rins, double routs, double muups, double mudowns,
	double phy1s, double phy2s, int bisections, int caseranges, int NNs )
{
	rin = rins;
	rout = routs;
	muup = muups;
	mudown = mudowns;
	phy1 = phy1s;
	phy2 = phy2s;
	bisection = bisections;
	caserange = caseranges;
	NN = NNs;	
}



/*
!*
!*
!*     PURPOSE:  Searches for minimum root pem of equations f(p)=0.  s   
!* 
!*     INPUTS:   f1234(1:4)-----array of p_r, p_theta, p_\phi, p_0. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               rin,rout-------Inner and outer radius of emission region or emission objects.
!*               muup,mudown----Boundary \mu coordinates of emission region or objects.
!*                              Or the maximum and minimum of \mu coordinates of emission region or objects.
!*                              These parameter can be specified void value by setting
!*                              values of parameter caserange. 
!*               phy1,phy2------Boundary \phi coordinates of emission region or objects.
!*                              Or the maximum and minimum of \phi coordinates of emission region or objects.
!*                              These parameter can be specified void value by setting
!*                              values of parameter caserange. 
!*               caserange------Tell routine whether muup, mudown and phy1, phy2 are provided.
!*                              caserange = 1, muup, mudown and phy1, phy2 are provided.
!*                              caserange = 2, muup, mudown are provided, but phy1, phy2 not.
!*                              caserange = 3, muup, mudown are not provided, phy1, phy2 are provided.
!*                              caserange = 4, muup, mudown and phy1, phy2 not are provided. 
!*                              Provided means corresponding parameters have been set specific value.
!*                              Not provided means corresponding parameters have not been set specific value,
!*                              but these parameter should also be given as dummy parameter.
!*               Fp-------------name of function f(p). This routine to compute f(p) should be prvided by
!*                              user, and the dummy variable of Fp should have following form:
!*                              Fp(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras).
!*               paras(1:10)----Array of parameters to descirbe function f(p).  
!*               bisection------Logical variable, if TURE, then use bisection method to search the 
!*                              root of equation f(p)=0, else use Newton-Raphson method. 
!*               NN-------------In the subroutine pemfind there is an important parameter: NN, which is the number
!*                              of sections of the interval (p1 , p2 ) or (p3 , p4 ) has been divided when 
!*                              searching the roots. One needs to set a proper value for NN, since if 
!*                              NN is too small, the roots exit on the interval (p1 , p2 ) or (p3 , p4 ) 
!*                              maybe omitted, if NN is too big, the time consumed by the code will be
!*                              large.
!*     OUTPUTS:  pemfind--------value of root of equation f(p)=0 for p.  
!*                              pemdisk=-1.D0, if the photon goto infinity.
!*                              pemdisk=-2.D0, if the photon fall into event horizon.         
!*     REMARKS:  This routine will search root between interval (p1, p2). We will chose NN points on 
!*               this interval, and check one by one to wheter f(p) changing its sign, if so, the root
!*               must be on interval (p_{i-1}, p_{i}). Then we use Bisection or Newton-Raphson method 
!*               to find the roots. One should set NN propriately to guarantee no root missing and
!*               also quickly find the root, thus saving the running time of the code.
!*     ROUTINES CALLED: radiustp, r2p, rootfind, Sectionp.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS:
!*
!*
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-27.
!*
!*
!*
!*
!*/
void pemfinds( ptcl *p, double *pemfind )
{
	//paras = parass;

	NN = 20;
	int cases;
	radiustp( p );

	if ( rin > p->rhorizon ) {
		if( p->r_tp1 >= rout )
			cases = 1;
		else {
			if ( p->r_tp1 > rin )
				cases = 2;
			else {
				if ( p->r_tp1 > p->rhorizon )
					cases = 3;
				else
					cases = 4;
			}
		}
	} else {
		cases=5;
	}

	//printf("pem1 = %d \n", cases);
	switch ( cases ) {
		case 1:
			*pemfind = - one;
			break;
		case 2:
			do {
				pemfindcase2( p, pemfind );
				if ( *pemfind > zero || *pemfind == -one || *pemfind == -two ) break;
				NN = NN * 2;
			} while (1);
			break;
		case 3:
			do {
				pemfindcase3( p, pemfind );
				if ( *pemfind > zero || *pemfind == -one || *pemfind == -two ) break;
				NN = NN * 2;
			} while (1);
			break;
		case 4:
			do {
				//printf(" pfind1 = %f   NN = %d \n", *pemfind, NN );
				pemfindcase4( p, pemfind );
				if ( NN > 20000 )exit(0);
				//printf(" pfind2 = %f   NN = %d \n", *pemfind, NN );
				if ( *pemfind > zero || *pemfind == -one || *pemfind == -two ) break;
				NN = NN * 2;
			} while (1);
			break;
		case 5:
			do {
				pemfindcase5( p, pemfind );
				if ( *pemfind > zero || *pemfind == -one || *pemfind == -two ) break;
				NN = NN * 2;
			} while (1);
			break;
	}
} 




static void pemfindcase2( ptcl *p, double *pemfind )
{
	double mu1, mu2;
	double pr, p1, p2, r1;
	int t1, t2;
	int tr1, tr2;
	mutp( p );
	if ( p->muobs > muup ) {
		if ( p->f1234[2] > zero || ( p->mobseqmtp && p->muobs == p->mu_tp1 ) ) {
			t1 = 0;
			t2 = 0;
			mu1 = muup;
			mu2 = - muup;
			p1 = mu2p( p, mu1, t1, t2 );
			p2 = mu2p( p, mu2, t1, t2 );

			tr1 = 0;
			tr2 = 0;
			pr = r2p( p, rout, tr1, tr2 );       

			if ( pr <= p1 )
				*pemfind = rootfind( p, p1, p2 );
			else if ( p1 < pr && pr < p2 )
				*pemfind = rootfind( p, pr, p2 );
			else if ( pr >= p2 ) {
				*pemfind = - one;
			}

			if ( *pemfind == -one ) {
				t2 = 1;
				mu1 = - muup;
				mu2 = muup;
				p1 = mu2p( p, mu1, t1, t2 );
				p2 = mu2p( p, mu2, t1, t2 );
				r1 = radius( p, p1 );

				if ( rout <= r1 )
					*pemfind = - one;
				else if ( r1 < rout ) {
					*pemfind = rootfind( p, p1, p2 );
				}
			}
		} else {
			t1 = 1;
			t2 = 0;
			mu1 = muup;
			mu2 = - muup;
			p1 = mu2p( p, mu1, t1, t2 );
			p2 = mu2p( p, mu2, t1, t2 );
			r1 = radius( p, p1 );
			if ( r1 < rout )
				*pemfind = rootfind( p, p1, p2 );
			else
				*pemfind = - one;

			if ( *pemfind == -one ) {
				t2 = 1;
				mu1 = - muup;
				mu2 = muup;
				p1 = mu2p( p, mu1, t1, t2 );
				p2 = mu2p( p, mu2, t1, t2 );
				r1 = radius( p, p1 );
				if ( r1 < rout )
					*pemfind = rootfind( p, p1, p2 );
				else
					*pemfind = - one;
			}
		}
	} else {
		if ( p->f1234[2] > zero || ( p->mobseqmtp && p->muobs == p->mu_tp1 ) ) {
			t1 = 0;
			t2 = 0;
			//mu1 = muup
			mu2 = - fmin( muup, p->mu_tp1 );
			tr1 = 0;
			tr2 = 0;
			//p1 = mu2p(f1234(3),f1234(2),lambda,q,mu1,sinobs,muobs,a_spin,t1,t2,scal)
			p1 = r2p( p, rout, tr1, tr2 );
			p2 = mu2p( p, mu2, t1, t2 );
    
			if (p2 > p1) {
				*pemfind = rootfind( p, p1, p2 );
				if ( *pemfind == - one ) {
					t2 = 1;
					if ( muup > p->mu_tp1 )
						p1 = p2;
					else {
						mu1 = - muup;
						p1 = mu2p( p, mu1, t1, t2 );
					}
					mu2 = fmin( muup, p->mu_tp1 );
					p2 = mu2p( p, mu2, t1, t2 );
					*pemfind = rootfind( p, p1, p2 );
				}
			} else {
				t2 = 1;
				mu2 = fmin( muup, p->mu_tp1 );
				if ( muup < p->mu_tp1 ) {
					mu1 = - muup;
					p1 = mu2p( p, mu1, t1, t2 );
				}
				p2 = mu2p( p, mu2, t1, t2 );   
				// write(*,*)'here p1 p2 333 sdf= ',p1, p2
				*pemfind = rootfind( p, p1, p2 );
			}
		} else {
			t1 = 0;
			t2 = 0;
			mu1 = muup;
			mu2 = fmin(muup, p->mu_tp1);
			tr1 = 0;
			tr2 = 0;
			//p1 = mu2p(f1234(3),f1234(2),lambda,q,mu1,sinobs,muobs,a_spin,t1,t2,scal)
			p1 = r2p( p, rout, tr1, tr2 );
			p2 = mu2p( p, mu2, t1, t2 );

			if (p2 > p1) {
				*pemfind = rootfind( p, p1, p2 );
				if ( *pemfind == - one ) {
					t1 = 1;
					if ( muup > p->mu_tp1 )
						p1 = p2;
					else {
						mu1 = muup;
						p1 = mu2p( p, mu1, t1, t2 );
					}
					mu2 = - fmin( muup, p->mu_tp1 );
					p2 = mu2p( p, mu2, t1, t2 );
    
					*pemfind = rootfind( p, p1, p2 );
				}
			} else {
				t1 = 1;
				if ( muup < p->mu_tp1 ) {
					mu1 = muup;
					p1 = mu2p( p, mu1, t1, t2 ); 
				}
				mu2 = - fmin( muup, p->mu_tp1 );
				p2 = mu2p( p, mu2, t1, t2 );
      
				*pemfind = rootfind( p, p1, p2 );
			}
		}
	}
}



static void pemfindcase3( ptcl *p, double *pemfind )
{
	double p1, p2;
	int tr1, tr2;

	tr1 = 0;
	tr2 = 0;
	p1 = r2p( p, rout, tr1, tr2 );
	p2 = r2p( p, rin, tr1, tr2 );  
	if ( caserange == 4 )
		*pemfind = rootfind( p, p1, p2 );
	else
		*pemfind = Sectionp( p, p1, p2 );


	if ( *pemfind == - one ) {
		tr1 = 1;
		p1 = r2p( p, rin, tr1, tr2 );
		p2 = r2p( p, rout, tr1, tr2 );

		if ( caserange == 4 )
			*pemfind = rootfind( p, p1, p2 );      
		else      
			*pemfind = Sectionp( p, p1, p2 );
			// write(*,*)'ff=',p1,p2,pemfind
	}
}



static void pemfindcase4( ptcl *p, double *pemfind )
{
	double p1, p2;
	int tr1, tr2;
        
	tr1 = 0;
	tr2 = 0;
	p1 = r2p( p, rout, tr1, tr2 );
	p2 = r2p( p, rin, tr1, tr2 );

	if ( caserange == 4 ) {
		//printf("pemfindcase4 = %f  %f \n", p1, p2);
		*pemfind = rootfind( p, p1, p2 ); 
		//printf("pemfindcase42 = %f \n", *pemfind);      
	} else        
		*pemfind = Sectionp( p, p1, p2 );

	//the photon will fall into the black hole.
	if ( *pemfind == - one )
		*pemfind = - two;

}



static void pemfindcase5( ptcl *p, double *pemfind )
{
	double p1, p2;
	int tr1, tr2;
        
	tr1 = 0;
	tr2 = 0;
	p1 = r2p( p, rout, tr1, tr2 );
              
	if ( p->r_tp1 > p->rhorizon ) {
		tr1 = 1;
		p2 = r2p( p, rout, tr1, tr2 );
	} else
		p2 = r2p( p, p->rhorizon, tr1, tr2 );

	*pemfind = rootfind( p, p1, p2 );
}



/*
!*
!*
!*     PURPOSE:  To judge whether a geodesic intersects with the surface of object or emission region
!*               which are described by function f(p). The space range of surface of object or
!*               emission region is (rin, rout) in radius, (mudown,muup) in poloidal and 
!*               (phy1, phy2) in azimuthal. If geodesic has no intersections on those intervals
!*               then it will no intersects with the surface and emission region, a special       
!*               value -1.D0 will be returned.
!* 
!*     INPUTS:   f1234(1:4)-----array of p_r, p_theta, p_phi, p_t which are the four momentum of the photon. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               t1,t2----------Number of photon meets turning points r_tp1 and r_tp2 respectively.
!*               p1,p2----------roots may hidden inside interval (p1,p2).
!*               muup,mudown,phy1,phy2,
!*               caserange,Fp,paras,bisection--------see instructions in pemfind.  
!*
!*     OUTPUTS:  Sectionp-------value of root of equation f(p)=0 for p.  
!*                              If Sectionp=-1.D0, No roots were found no interval (p1, p2).     
!*     ROUTINES CALLED: phi, mucos, rootfind. 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: 
!*
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-28.
!*
!*/
double Sectionp( ptcl *pt, double p1s, double p2s )
{
	double p1, p2, phya1, phya2, rtn;
	double deltax = 5.e-5;
	double mu1, mu2;

	p1 = p1s - deltax;
	p2 = p2s + deltax;
	phya1 = zero;
	phya2 = zero;
	if ( caserange == 1 || caserange == 3 ) {
		YNOGKC( pt, p1 );
		phya1 = pt->phi_p;
		YNOGKC( pt, p2 );
		phya2 = pt->phi_p;
	}
	if ( caserange == 1 || caserange == 2 ) {
		mu1 = mucos( pt, p1 );
		mu2 = mucos( pt, p2 );
		if ( (mu1 - muup)*(mu1 - mudown) < zero || (mu2-muup)*(mu2-mudown) < zero || 
		(muup - mu1)*(muup - mu2) < zero || (mudown - mu1)*(mudown - mu2) < zero ) {

			if ( caserange == 1 || caserange == 3 ) {
				if ( ( phya1 - phy1 ) * ( phya1 - phy2 ) < zero
				  || ( phya2 - phy1 ) * ( phya2 - phy2 ) < zero
				  || ( phy1 - phya1 ) * ( phy1 - phya2 ) < zero
				  || ( phy2 - phya1 ) * ( phy2 - phya2 ) < zero )
					/* the geodesic intersecte the zone defined 
					by r1,r2,mu1,mu2,phy1,phy2,so it has the 
					possibility to hit the surface of the object. */
					rtn = rootfind( pt, p1, p2 ); //(1,1)
					// write(*,*)'ff=',p1,p2,Sectionp
		                else 
					/* the (phy1,phy2) and (phya1,phya2) has 
					no public point,so we don't consider it.*/
					rtn = - one;  // (1,1)

			} else
		                rtn = rootfind( pt, p1, p2 );  //!(1,0)

		} else {
			/* the internal of (mu1,mu2) and (muup,mudown) does not
			overfold each other,so the geodesic will not hit the
			surface of the object at the internal (p_rout,p_rout2),
			which also means the geodesic will never hit the object.
			so nothing needed to be done. */
			rtn = - one;  //    
		}
	} else {
		if ( caserange == 1 || caserange == 3 ) {
			if ( (phya1-phy1)*(phya1-phy2) < zero
			  || (phya2-phy1)*(phya2-phy2) < zero
			  || (phy1-phya1)*(phy1-phya2) < zero
			  || (phy2-phya1)*(phy2-phya2) < zero )
				rtn = rootfind( pt, p1, p2 );  // (0,1)
			else
				/* the (phy1,phy2) and (phya1,phya2) has no 
				public points,so we don't consider it. */
		                rtn = - one;  //(0,1)

		} else
			rtn = rootfind( pt, p1, p2 );        // (0,0)        

	}
        return rtn;
}



/*
!*
!*
!*     PURPOSE:  To search roots on interval (p1, p2). If no roots were found 
!*               a special value -1.D0 will return. 
!* 
!*     INPUTS:   f1234(1:4)-----array of p_r, p_theta, p_phi, p_t which are the four momentum of the photon. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               t1,t2----------Number of photon meets turning points r_tp1 and r_tp2 respectively.
!*               p1,p2----------roots may hidden inside interval (p1,p2). 
!*               Fp,paras,bisection--------see instructions in pemfind.  
!*
!*     OUTPUTS:  rootfind-------value of root of equation f(p)=0 for p.  
!*                              If rootfind=-1.D0, No roots were found no interval (p1, p2).     
!*     ROUTINES CALLED: phi, mucos, rootfind. 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS:
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-28.
!*
!*/
static double rootfind( ptcl *pt, double p1, double p2 )
{
	double rtfd, deltap, p, f_p, sp1, sp2;;
	//double const dp = 1.e-5;
	int NNf, k;
 
        //p1 = p1 + dp
	deltap = ( p2 - p1 ) / NN;
	NNf = NN;
	p = p1;
	f_p = Fp( pt, p );

	if ( f_p == 0.0 ) {
                rtfd = p;
                return (rtfd);
	}
	if ( f_p < zero ) {
                k = 0;
                do {
			if (f_p == zero)
                            return p;

			if (f_p > zero || k > NNf)
				break;
                        k += 1;
                        p = p1 + deltap * k;
                        f_p = Fp( pt, p );
		} while (true);
	} else {
		k = 0;
		do {
			if (f_p == zero) {
                            return p;
			}
			if ( f_p < zero || k > NNf )
				break;
			k = k + 1;
			p = deltap * k + p1;
			f_p = Fp( pt, p );
			//printf("rootfind2 = %d %d %f %f %f %f \n", k, NNf, f_p, Fp( pt, p1 ), Fp( pt, p ), Fp( pt, p2 ) );
		} while(true);
	}

	//printf("ttsf = %d  %f  %f \n", k, p, f_p);

	if (k <= NNf) {
		sp1 = p - deltap;
		sp2 = p;
		// Using bisection or Newton Raphson method to find roots on interval (sp1, sp2).   
		if ( bisection )
			rtfd = Bisectionp( pt, sp1, sp2 );
		else {
			//printf("rootfind2 = %f  %f \n", sp1, sp2);
			rtfd = NewRapson( pt, sp1, sp2 );
		}

	} else {
		// On interval (p1, p2) no roots are found!
		rtfd = - 1.0;
	}
        return rtfd;
}




/*
!*     PURPOSE:  To search roots on interval (p1, p2) by bisenction method for 
!*               equation Fp(p)=0.
!*     INPUTS:   Parameters descriptions ---------see instructions in pemfind.  
!*
!*     OUTPUTS:  Bisectionp-------value of root of equation f(p)=0 for p.        
!*     ROUTINES CALLED: phi, mucos, rootfind. 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: ****************************************** 
!*
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-28.
!*
!*
!*/
static double Bisectionp( ptcl *pt, double p1, double p2 )
{
	double pc, f1, f2, fc;
	int counter;
        
	pc = ( p1 + p2 ) / two;
	f1 = Fp( pt, p1 );
	f2 = Fp( pt, p2 );
	fc = Fp( pt, pc );

        counter = 0;
	do {
		if ( f1 * fc > zero ) {
			p1 = pc;
			f1 = fc;
			pc = ( p1 + p2 ) / two;
			fc = Fp( pt, pc );
		}
		if ( f2 * fc > zero ) {
			p2 = pc;
			f2 = fc;
			pc = ( p1 + p2 ) / two;
			fc = Fp( pt, pc );
		}
		if ( fc == zero ) break;
		if ( counter > 200 ) break;
		//printf("fc = %f %f %f %f %f %f  \n ", f1, fc, f2, p1, pc, p2);
		//printf("fc = %d \n", counter);
		counter += 1;
	} while( fabs( p2 - p1 ) / p1 > 1.e-5 );
 
        return pc;
}




/*
!*     PURPOSE:  To search roots on interval (p1, p2) by Newton Raphson method for 
!*               equation Fp(p)=0.
!*     INPUTS:   Parameters descriptions ---------see instructions in pemfind.  
!*
!*     OUTPUTS:  NewRapson-------value of root of equation f(p)=0 for p.        
!*     ROUTINES CALLED: phi, mucos, rootfind. 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*
!*
!*
!*     C VERSION:  Yang Xiao-lin    2022-12-28.
!*
!*/
static double NewRapson( ptcl *pt, double p1, double p2 )
{
	double ptn, dp, f_p, f_pj, dfp, h, temp;
	double pj;
	int k;
	double const kmax = 30, pacc = 1.e-2, EPS = 1.e-5;

        ptn = ( p2 + p1 ) / two;
	for ( k = 1; k <= kmax; k++ ) {
		//printf("tt1 = %d  %f \n", k, ptn);
		temp = ptn;
		h = EPS * fabs( temp );
		if ( h == zero ) h = EPS;
		pj = temp + h;
		h = pj - temp;
		//printf("tt2 = %f  %f  %f \n", pj, h, temp);
		if ( temp < zero || pj < zero )
		        return -fabs(temp);

		f_p = Fp( pt, temp );

		f_pj = Fp( pt, pj );


		dfp = ( f_pj - f_p ) / h;
		dp = f_p / dfp;
		ptn = ptn - dp;
		//If((ptn-p1)*(p2-ptn).lt.zero)then
		//write(unit=6,fmt=*)'ptn jumps out of brakets [p1, p2]!'
		//printf("fc2 = %f %f \n ", dp, ptn);
		if ( fabs( dp ) < pacc )
		        return ptn;
	}
	return ptn;
}












