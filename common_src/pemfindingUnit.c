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

 
static void pemfindcase2( ptcl *p );


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
void pemfinds( ptcl *p, int *cases, double rin, double rout, double muup, 
	double mudown, double phy1, double phy2, int caserange, double Fp, 
	double paras, double orir, double oricosth, int bisection, 
	double pemfind, int NN )
{
	radiustp( p );

	if ( rin > p->rhorizon ) {
		if( p->r_tp1 >= rout )
			*cases = 1;
		else {
			if ( p->r_tp1 > rin )
				*cases = 2;
			else {
				if ( p->r_tp1 > p->rhorizon )
					*cases = 3;
				else
					*cases = 4;
			}
		}
	} else {
		*cases=5;
	}

	switch ( *cases ) {
		case 1:
			pemfind = - one;
			break;
		case 2:
			break;
		case 3:
			break;
		case 4:
			break;
		case 5:
			break;
	}
} 




static void pemfindcase2( ptcl *p, double muup )
{
	double mu1, mu2;
	double pr, p1, p2;
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
				pemfind = rootfind(f1234,lambda,q,sinobs,muobs,a_spin,
					robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection);
			else if ( p1 < pr && pr < p2 )
				pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,
					robs,scal,tr1,tr2,pr,p2,NN,Fp,paras,orir,oricosth,bisection);
			else if (pr >= p2) {
				pemfind = - one;
			}

			if ( pemfind == -one ) {
				t2 = 1;
				mu1 = - muup;
				mu2 = muup;
				p1 = mu2p( p, mu1, t1, t2 );
				p2 = mu2p( p, mu2, t1, t2 );
				r1 = radius( p, p1 );

				if (rout <= r1)
					pemfind = - one;
				else if ( r1 < rout ) {
					pemfind = rootfind(f1234,lambda,q,sinobs,muobs,a_spin,
						robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection);
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
				pemfind = rootfind(f1234,lambda,q,sinobs,muobs,a_spin,
					robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection);
			else
				pemfind = - one;

			if ( pemfind == -one ) {
				t2=1
				mu1 = - muup
				mu2 = muup  
				p1 = mu2p(f1234(3),f1234(2),lambda,q,mu1,sinobs,muobs,a_spin,t1,t2,scal) 
				p2 = mu2p(f1234(3),f1234(2),lambda,q,mu2,sinobs,muobs,a_spin,t1,t2,scal) 
				r1 = radius(p1,f1234(1),lambda,q,a_spin,robs,scal) 
				if (r1 < rout) then          
					pemfind = rootfind(f1234,lambda,q,sinobs,muobs,a_spin,
						robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection) 
				else
					pemfind = - one;
			}
		}
	else
		      if (f1234(2)>zero .or. (mobseqmtp .and. muobs == mu_tp1)) then 
		        t1 = 0
		        t2 = 0
		        !mu1 = muup
		        mu2 = - min(muup, mu_tp1)  
		        tr1 = 0
		        tr2 = 0
		        !p1 = mu2p(f1234(3),f1234(2),lambda,q,mu1,sinobs,muobs,a_spin,t1,t2,scal) 
		        p1 = r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)        
		        p2 = mu2p(f1234(3),f1234(2),lambda,q,mu2,sinobs,muobs,a_spin,t1,t2,scal)  
		                 !write(*,*)'here p1 p2 1111= ',p1, p2      
		        if (p2 > p1) then           
		            pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
		                       robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection)  
		            If(pemfind.eq.-one)then
		                t2=1 
		                if (muup > mu_tp1) then
		                    p1 = p2
		                else
		                    mu1 = -muup
		                    p1 = mu2p(f1234(3),f1234(2),lambda,q,mu1,sinobs,muobs,a_spin,t1,t2,scal) 
		                endif
		                mu2 = min(muup, mu_tp1)  
		                p2 = mu2p(f1234(3),f1234(2),lambda,q,mu2,sinobs,muobs,a_spin,t1,t2,scal)    
		                 !write(*,*)'here p1 p2 2222= ',p1, p2      
		                pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
		                            robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection)                     
		            endif 
		        else 
		            t2=1 
		            mu2 = min(muup, mu_tp1)  
		            if (muup < mu_tp1) then
		                mu1 = -muup
		                p1 = mu2p(f1234(3),f1234(2),lambda,q,mu1,sinobs,muobs,a_spin,t1,t2,scal)  
		            endif
		            p2 = mu2p(f1234(3),f1234(2),lambda,q,mu2,sinobs,muobs,a_spin,t1,t2,scal)    
		           !write(*,*)'here p1 p2 333 sdf= ',p1, p2 
		            pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
		                        robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection)   
		        endif 
		      else 
		        t1 = 0
		        t2 = 0
		        mu1 = muup
		        mu2 = min(muup, mu_tp1)  
		        tr1 = 0
		        tr2 = 0
		        !p1 = mu2p(f1234(3),f1234(2),lambda,q,mu1,sinobs,muobs,a_spin,t1,t2,scal) 
		        p1 = r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)        
		        p2 = mu2p(f1234(3),f1234(2),lambda,q,mu2,sinobs,muobs,a_spin,t1,t2,scal)
		        !write(*,*)'here p1 p2 1111= ',p1, p2  
		        if (p2 > p1) then  
		            pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
		                      robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection)  
		            If(pemfind.eq.-one)then
		                t1=1 
		                if (muup > mu_tp1) then
		                    p1 = p2
		                else
		                    mu1 = muup
		                    p1 = mu2p(f1234(3),f1234(2),lambda,q,mu1,sinobs,muobs,a_spin,t1,t2,scal) 
		                endif
		                mu2 = -min(muup, mu_tp1)  
		                p2 = mu2p(f1234(3),f1234(2),lambda,q,mu2,sinobs,muobs,a_spin,t1,t2,scal)    
		                !write(*,*)'here p1 p2 2222= ',p1, p2      
		                pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
		                            robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection)                     
		            endif 
		        else
		            t1=1
		            if (muup < mu_tp1) then
		                mu1 = muup 
		                p1 = mu2p(f1234(3),f1234(2),lambda,q,mu1,sinobs,muobs,a_spin,t1,t2,scal)  
		            endif
		            mu2 = - min(muup, mu_tp1)  
		            p2 = mu2p(f1234(3),f1234(2),lambda,q,mu2,sinobs,muobs,a_spin,t1,t2,scal)  
		            !write(*,*)'here p1 p2 333 sdf= ',p1, p2           
		            pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
		                        robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection)  
		        endif 
		      end if
	}
}

























