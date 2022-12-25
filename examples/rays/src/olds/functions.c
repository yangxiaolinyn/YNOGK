#include "1RandUtils.h"
#include "ConstantsUnits.h"
#include <stdio.h>
#include <math.h>

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*\
*  Functions' declarations defined in my code.
\*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

 
double F_cos_mphi( long m, double phi, double am[] )
{
	double temp_v; 

	temp_v = zero;
	for(long i = 0; i <= m; i++)
		temp_v += am[i] * cos(i * phi);
	return temp_v;
} 
 
double F_sin_mphi( long m, double phi, double am[] )
{
	double temp_v; 

	temp_v = zero;
	for(long i = 0; i <= m; i++)
		temp_v += am[i] * sin(i * phi);
	return temp_v;
}   

 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
double cosm_phi_sampling( const int m, const double Phi1, const double rnum )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
{ 
	double temp_phi1, m2, xi; 
  
	if (m == 0) 
		return Phi1;
	else {   
		xi = rnum ;  
		if (  xi * two <= one + cos( m*Phi1 ) ) {
			temp_phi1 = Phi1; 
		} else {
			m2 = floor( Phi1 / ( pi/m ) );
			temp_phi1 = (two * m2 + one) * pi / m - Phi1; 
		}  
		return temp_phi1;
	} 
} 
 

 
