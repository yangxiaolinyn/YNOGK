/*
 * common_src/root34.c
 *
 * This module aim on solve cubic and quartic polynomial equations.
 * One can use these subroutine root3 and root4 to find roots of cubic and 
 * quartic equations respectively. 
 *
 * Author       Yang Xiao-lin
 *
 * Date         2022-11-16
 *
 * City         Kun-ming Yunnan Provice, China
 *
 * 2022-11-16   Finished the fisrt writting of the module.
 */

#include "root34.h"


/*
!* PURPOSE:  This subroutine aim on sorting a1, a2, a3 by decreasing way.
!* INPUTS:    a1,a2,a3----they are the number list required to bo sorted. 
!* OUTPUTS:   s1,s2,s3----sorted number list with decreasing way. 
!*      
!* AUTHOR:        Yang, Xiao-lin
!* DATE WRITTEN:  2022-11-16
!*/
void sort( double a1, double a2, double a3, double *s)
{
	double temp;
	int i,j;
   
	s[0] = a1;
	s[1] = a2;
	s[2] = a3;
   
	for (i = 0; i < 3; i++) {
		for (j = i + 1; j < 3; j++) {
			if ( s[i] < s[j] ) {
				temp = s[i];
				s[i] = s[j];
				s[j] = temp;
			}
		}
	}
}


 
/*
!* PURPOSE:   This subroutine aim on solving cubic equations: x^3+b*x^2+c*x+d=0.
!* INPUTS:    b, c, d-----they are the coefficients of the equation. 
!* OUTPUTS:   r1,r2,r3----roots of the equation, with complex number form.
!*            del---------the number of real roots among r1, r2, r3. 
!* ROUTINES CALLED:  sort    
!* This code comes from internet.
 */
int root3(double b, double c, double d, double _Complex *roots, unsigned int *del)
{
	double p, q, DD;
	//double _Complex roots[3] = {};

	/*Step 1: Calculate p and q */
	p  = c - pow(b, 2) / 3.0;
	q  = (2.0 * pow(b, 3) - 9.0 * b * c + 27.0 * d) / 27.0;

	/*Step 2: Calculate DD (discriminant) */
	DD = pow(p, 3) / 27.0 + pow(q, 2) / 4.0;

	double phi, y1, y2 = 0.0, y3, y2r, y2i, temp1, temp2, u, v;
	/*Step 3: Branch to different algorithms based on DD*/
	if ( DD < 0.0 ) {
		/* Step 3b:
		 * 3 real unequal roots -- 
		use the trigonometric formulation */
		phi = acos( -q / two / sqrt( fabs( pow(p, 3) ) / 27.0 ) );
		temp1 = two * sqrt( fabs(p) / 3.0 );
		y1 =  temp1 * cos( phi / 3.0 );
		y2 = -temp1 * cos( (phi + pi) / 3.0 );
		y3 = -temp1 * cos( (phi - pi) / 3.0 );
	} else {
		/* Step 3a:
		 * 1 real root & 2 conjugate complex roots 
		 * OR 3 real roots (some are equal)*/
		double sqrtDD;
		sqrtDD = sqrt( DD );
		temp1 = - q / two + sqrtDD;
		temp2 = - q / two - sqrtDD;
		u = pow( fabs(temp1), 1.0 / 3.0 ); //fabs(temp1)**(1.0/3.0);
		v = pow( fabs(temp2), 1.0 / 3.0 ); //fabs(temp2)**(1.0/3.0);
		if( temp1 < 0.0 ) u = -u;
		if( temp2 < 0.0 ) v = -v;
		y1  = u + v;
		y2r = - (u + v) / two;
		y2i =   (u - v) * sqrt3 / two;
	};

	/* Step 4: Final transformation */
	temp1 = b / 3.0;
	y1 -= temp1;
	y2 -= temp1;
	y3 -= temp1;
	y2r -= temp1;

	double arr[3];
	/* Assign answers */
	if ( DD < 0.0 ) {
		sort(y1, y2, y3, arr);
		roots[1] = arr[0];
		roots[2] = arr[1];
		roots[3] = arr[2];
		*del = 3;
	} else if ( DD == 0.0) {
		sort(y1, y2r, y2r, arr);
		roots[1] = arr[0];
		roots[2] = arr[1];
		roots[3] = arr[2];
		*del = 3;
	} else {
		roots[1] = y1;
		roots[2] = y2r + y2i*I;
		roots[3] = y2r - y2i*I;
		*del = 1;
	}

	return 0;
}



/*
!* PURPOSE:   This subroutine aims on solving quartic equations: 
!*            x^4+b*x^3+c*x^2+d*x+e=0.
!* INPUTS:    b, c, d, e-----they are the coefficients of equation. 
!* OUTPUTS:   r1,r2,r3,r4----roots of the equation, with complex number form.
!*            reals------------the number of real roots among r1,r2,r3,r4.  
!* ROUTINES CALLED:  root3   
!* AUTHOR:        Yang, Xiao-lin 
!* DATE WRITTEN:  2022-11-16
*/ 
int root4(double b, double c, double d, double e, 
	double _Complex *r4, unsigned int *reals)
{
	double q, r, s;
	double _Complex s3[3], s3t[4], temp[4], temp1;
	int i, j;
	unsigned del;
	double a1, b1, c1;

	q = c - 3.0 * pow(b, 2) / 8.0;
	r = d - b * c / 2.0 + pow(b, 3) / 8.0;
	s = e - b * d / 4.0 + pow(b, 2) * c / 16.0 - 3.0 * pow(b, 4) / 256.0;

	a1 = 2.0*q;
	b1 = pow(q, 2) - 4.0 * s;
	c1 = -pow(r, 2);

	root3( a1, b1, c1, s3t, &del );
 
	if ( del == 3 ) { 
		if ( creal(s3t[3]) >= 0.0 ) {
			*reals = 4;
			s3[0] = csqrt( s3t[1] );
			s3[1] = csqrt( s3t[2] );
			s3[2] = csqrt( s3t[3] );
		} else { 
			*reals = 0;  
			s3[0] = csqrt( s3t[1] );
			s3[1] = csqrt( s3t[2] );
			s3[2] = csqrt( s3t[3] );
		}
      } else {
		if ( creal(s3t[1]) >= 0.0 ) {
			*reals = 2;
			s3[0] = sqrt( creal(s3t[1]) ) + zero * I;
			s3[1] = csqrt( s3t[2] );
			s3[2] = creal( s3[1] ) - cimag( s3[1] ) * I;
		} else {
			*reals = 0;
			s3[0] = csqrt( s3t[1] );
			s3[1] = csqrt( s3t[2] );
			s3[2] = csqrt( s3t[3] );
		}
	}
  

	if ( creal( - s3[0]*s3[1]*s3[2] ) * r < 0.0 ) s3[0] = - s3[0];

	temp[0] = (s3[0] + s3[1] + s3[2]) / 2.0 - b / 4.0;
	temp[1] = (s3[0] - s3[1] - s3[2]) / 2.0 - b / 4.0;
	temp[2] = (s3[1] - s3[0] - s3[2]) / 2.0 - b / 4.0;
	temp[3] = (s3[2] - s3[1] - s3[0]) / 2.0 - b / 4.0;

	/* make the order of the results correct */
	for ( i = 0; i < 4; i++ ) {
		for ( j = i + 1; j < 4; j++ ) {
			if( creal( temp[i] ) > creal( temp[j] ) ) {
				temp1   = temp[i];
				temp[i] = temp[j];
				temp[j] = temp1;
			}             
		}                
	}
	r4[1] = temp[0];
	r4[2] = temp[1];
	r4[3] = temp[2];
	r4[4] = temp[3];

	return 0;
}
















