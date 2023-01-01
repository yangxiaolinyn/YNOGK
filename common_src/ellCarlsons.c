/*
 * common_src/ellCarlsons.c
 *
 *              This module includes integrals and functions for Carlson's 
!*              integral method. Those codes mainly come from geokerr.f of
!*              Dexter & Agol (2009).
 *
 * Author       Yang Xiao-lin
 *
 * City         Kun-ming Yunnan Provice, China
 *
 * 2022-11-20   Starting the writting of the module.
 */


#include <ellCarlsons.h> 



static int elldoublecomplexs( int *index_p5, double f1, double g1, double h1, 
		double f2, double g2, double h2, double a5, double b5,
		double y, double x, double rff_p, double *integ, int cases);



double sign( double y, double x )
{ 
	if ( x > 0.0 )
		return fabs(y);
	else if ( x < 0.0 )
		return (-fabs(y));
	else
		return 0.0;
}

/*
!*     PURPOSE: Computes integrals: J_k(h)=\int^x_y (b4*t+a4)^(k/2)*(4*t^3-g_2*t-g_3)^(-1/2)dt.
!*              Where integer index k can be 0, -2, -4 and 2. (75) and (76) of Yang & Wang (2012).   
!*     INPUTS:  x,y -- limits of integral.
!*              bb(1:3) -- Roots of equation 4*t^3-g_2*t-g_3=0 solved by routine root3.
!*              del -- Number of real roots in bb(1:3).
!*              p4,rff_p,integ -- p4(1:4) is an array which specifies the value of index k of J_k(h). 
!*                 If p4(1)=0, then J_0 was computed and sent to integ(1).
!*                 If p4(1)!=0, then J_0 was replaced by parameter p, and rff_p=p.  
!*                 If p4(2)=-2, then J_{-2} was computed and sent to integ(2).                      
!*                 If p4(3)=2, then J_{2} was computed and sent to integ(3).
!*                 If p4(4)=-4, then J_{-4} was computed and sent to integ(4).
!*              cases -- If cases=1, then only J_0 was computed.
!*                       If cases=2, then only J_0 and J_{-2} are computed.
!*                       If cases=3, then only J_0, J_{-2} and J_{2} are computed.    
!*                       If cases=4, then J_0, J_{-2}, J_{2} and J_{-4} are computed.            
!*     OUTPUTS: integ -- is an array saved the results of J_k(h). 
!*     ROUTINES CALLED:  ellcubicreals, ellcubiccomplexs     
!*     ACCURACY:   Machine.
!*     REMARKS:
!*     AUTHOR:     Yang Xiao-lin
!*     DATE WRITTEN:  2022-11-20
!*/
int weierstrass_int_J3(double y, double x, double _Complex *bb, int del, 
			double a4, double b4, int *p4, 
			double rff_p, double *integ, int cases)
{
	double yt, xt, a1, b1, a2, b2, a3, b3,
		tempt, f, g, a44, b44, sign_h;
	int i;
	//complex*16 bb[3] integ[4]
	int inverse;

	xt=x;
	yt=y;
	a44=a4;
	b44=b4;

	inverse=0;
	if ( fabs( xt - yt ) <= 1.e-9 ) {
		for ( i = 1; i <= 4; i++ )
			integ[i] = 0.0;
		return 0;
	}
	if ( yt > xt ) {
		tempt=xt;
		xt=yt;
		yt=tempt;
		inverse=1;
	}

	sign_h = sign( one, b44*xt+a44 );
	//printf(" yt =%f xt= %f \n ", yt, xt);
	//printf("weierstrass_int_J3 = %d case = %d \n", del, cases);
	if ( del == 3) {
		/* equation (75) of Yang & Wang (2013). */
		a44 = sign_h*a44;
		b44 = sign_h*b44;
		a1 = -creal(bb[1]);
		a2 = -creal(bb[2]);
		a3 = -creal(bb[3]);
		b1 = 1.0;
		b2 = 1.0;
		b3 = 1.0;
		ellcubicreals( p4, a1, b1, a2, b2, a3, b3, a44, b44, 
				yt, xt, rff_p*two, integ, cases );

		if (inverse) {
			integ[1] = - integ[1];
			integ[2] = - integ[2];
			integ[3] = - integ[3];
			integ[4] = - integ[4];
		}
		for (int i = 1; i <= 4; i++ ) {
			integ[i] = integ[i] / two;
			integ[i] = integ[i] * pow(sign_h, -p4[i]/2);
		}
		//write(*,*)'weierstrass_int_J3',integ
	} else {
		//equation (76) of Yang & Wang (2012).
		a44 = sign_h*a44;
		b44 = sign_h*b44;
		a1 = -creal(bb[1]);
		b1 = one;
		f = pow(creal(bb[2]), 2) + pow(cimag(bb[2]), 2);
		g = -two * creal(bb[2]);
 
		ellcubiccomplexs( p4, a1, b1, a44, b44, f, g, one, 
				yt, xt, rff_p*two, integ, cases);
		//printf("f = %f g = %f yt = %f xt = %f  g3 = %f g3 = %d \n", f, g, yt, xt, rff_p*two, cases);
		if (inverse) {
			integ[1] = - integ[1];
			integ[2] = - integ[2];
			integ[3] = - integ[3];
			integ[4] = - integ[4];
		}
		for (i = 1; i <= 4; i++) { 
			integ[i] = integ[i] / two;
			if ( sign_h != zero )
				integ[i] = integ[i] * pow(sign_h, -p4[i]/2);
		}
	}
	//printf("I1 = %f I2 = %f I3 = %f I4 = %f \n", integ[1], integ[2], integ[3], integ[4]);
	return 0;
}


/*
!*     PURPOSE: Computes integrals: J_k(h)=\int^x_y (b5*r+a5)^(k/2)*
!*                                  [(h1*r^2+g1*r+f1)(h2*r^2+g2*r+f2)]^(-1/2)dr.
!*              Where integer index k can be 0, -2, -4, 2 and 4. 
!*                                  (77) of Yang & Wang (2012).   
!*     INPUTS:  x,y -- limits of integral.  
!*              p5,rff_p,integ -- p5(1:5) is an array which specifies 
!*                                the value of index k of J_k(h). 
!*                 If p5(1)=0, then J_0 was computed and sent to integ(1).
!*                 If p5(1)!=0, then J_0 was replaced by parameter p, and rff_p=p.  
!*                 If p5(2)=-2, then J_{-2} was computed and sent to integ(2).                      
!*                 If p5(3)=2, then J_{2} was computed and sent to integ(3).
!*                 If p5(4)=-4, then J_{-4} was computed and sent to integ(4).
!*                 If p5(4)=4, then J_{4} was computed and sent to integ(5).
!*              cases -- If cases=1, then only J_0 will be computed.
!*                       If cases=2, then only J_0 and J_{-2} will be computed.
!*                       If cases=3, then only J_0, J_{-2} and J_{2} will be computed.
!*                       If cases=4, then J_0, J_{-2}, J_{2} and J_{-4} will be computed. 
!*                       If cases=5, then J_0, J_{-2}, J_{2} and J_{4} will be computed.     
!*     OUTPUTS: integ -- is an array saved the results of J_k(h). 
!*     ROUTINES CALLED:  elldoublecomplexs
!*     ACCURACY:   Machine.
!*     AUTHOR:     Yang Xiao-lin
!*     DATE WRITTEN:  2022-11-21
!*/
int carlson_doublecomplex5(double y, double x, double f1, double g1, 
			double h1, double f2, double g2, double h2,
			double a5, double b5, int *p5, double rff_p, 
			double *integ, int cases)    
{
	double xt, yt, tempt, a55, b55, sign_h;
	int i;
	int inverse;

	xt=x;
	yt=y;
	a55=a5;
	b55=b5;
	inverse = 0; 

	if ( fabs(xt-yt) <= 1.e-9 ) {
		for ( i = 1; i <= 5; i++ )
			integ[i] = 0.0;
		return 0;
	}

	if ( yt > xt ) {
		tempt=xt;
		xt=yt;
		yt=tempt;
		inverse=1;
	}

	sign_h = sign( one, b55 * xt + a55 );
	a55=sign_h*a55;
	b55=sign_h*b55;
	/* equation (77) of Yang & Wang (2012). */
	elldoublecomplexs( p5, f1, g1, h1, f2, g2, h2, a55, b55,
				yt, xt, rff_p, integ, cases);
	if ( inverse ) {
		integ[1]=-integ[1];
		integ[2]=-integ[2];
		integ[3]=-integ[3];
		integ[4]=-integ[4];
		integ[5]=-integ[5];
	}

	for ( i = 1; i <= 5; i++ )
            integ[i] = integ[i] * pow(sign_h, -p5[i]/2);
	return 0;
}

/*
!*     PURPOSE: Computes J_k(h)=\int_y^x dt (b4*t+a4)^(k/2)[(b1*t+a1)*(b2*t+a2)*(b3*t+a3)]^{-1/2}. 
!*              It is the case of equation W(t)=4*t^3-g_2*t-g_3=0 has three real roots. 
!*              Equation (61) of Yang & Wang (2012).   
!*     INPUTS:  Arguments for above integral. If index_p4(1)=0, then J_0 will be computed, 
!*              else J_0 will be replaced by parameter p, and rff_p=p. 
!*     OUTPUTS:  Value of integral J_k(h).
!*     ROUTINES CALLED: RF,RJ,RC,RD
!*     ACCURACY:   Machine.
!*     AUTHOR:     Dexter & Agol (2009)
!*     MODIFIED:   Yang Xiao-lin  2022-11-21
!*     DATE WRITTEN:  4 Mar 2009
!*     REVISIONS: 
!*/ 
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int ellcubicreals( int *index_p4, double a1, double b1, double a2, 
		double b2, double a3, double b3, double a4, double b4,
		double y, double x, double rff_p, double *integ, int cases)
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
{
	double d12, d13, d14, d24, d34, X1, X2, X3, X4, 
		Y1, Y2, Y3, Y4, U1c, U32, U22, W22, U12, Q22, P22,
		I1c, I3c, r12, r13, r24i, r34i, I2c, J2c, K2c, r14,
		r24, r34, rff;
 
	/* c (2.1) Carlson (1989) */
	d12=a1*b2-a2*b1;
	d13=a1*b3-a3*b1;
	d14=a1*b4-a4*b1;
	d24=a2*b4-a4*b2;
	d34=a3*b4-a4*b3;
	r14=a1/b1-a4/b4;
	r24=a2/b2-a4/b4;
	r34=a3/b3-a4/b4;
           
	//c (2.2) Carlson (1989)
	X1 = sqrt(fabs(a1+b1*x));
	X2 = sqrt(fabs(a2+b2*x));
	X3 = sqrt(fabs(a3+b3*x));
	X4 = sqrt(fabs(a4+b4*x));
	Y1 = sqrt(fabs(a1+b1*y));
	Y2 = sqrt(fabs(a2+b2*y));
	Y3 = sqrt(fabs(a3+b3*y));
	Y4 = sqrt(fabs(a4+b4*y));

	// (2.3) Carlson (1989)
	if ( x < infinity ) {
		U1c = (X1*Y2*Y3+Y1*X2*X3)/(x-y);
		U12 = pow( U1c, 2 );
		U22 = pow( (X2*Y1*Y3+Y2*X1*X3)/(x-y), 2);
		U32 = pow( (X3*Y1*Y2+Y3*X1*X2)/(x-y), 2);
	} else {
		U1c = sqrt( fabs( b2 * b3 ) ) * Y1;
		U12 = pow( U1c, 2 );
		U22 = b1 * b3 * pow(Y2, 2);
		U32 = b2 * b1 * pow(Y3, 2);
	}
	//c (2.4) Carlson (1989)
	W22 = U12;
	W22 = U12 - b4 * d12 * d13 / d14;
	//(2.5) Carlson (1989)
	if ( x < infinity ) {
		Q22 = pow((X4*Y4/X1/Y1), 2) * W22;
		//printf( "X41 = %f \n  Y4 = %f \n  X1 = %f \n Y1 = %f \n ", X4, Y4, X1, Y1 );
		//printf( "W2= %f \n  ", W22 );
	} else {
		Q22 = b4 / b1 * pow(Y4/Y1, 2) * W22;
		//printf( "X42 = %f \n  Y4 = %f \n  X1 = %f \n Y1 = %f \n ", b4, b1, Y4, Y1 );
		//printf( "W2= %f \n  ", W22 );
	}

	P22 = Q22 + b4 * d24 * d34 / d14;
       
	/* Now, compute the three integrals we need [-1,-1,-1],
	 * [-1,-1,-1,-2], and [-1,-1,-1,-4]: we need to calculate 
	 * the [-1,-1,-1,2] integral,we add it in this part. */
	if ( index_p4[1] == 0 ) {
		//c (2.21) Carlson (1989)
		rff = rf(U32,U22,U12);
		integ[1] = two * rff;
		//printf( "ellcubic = %40.35f \n %40.35f \n %40.35f \n %40.35f \n  ", rff, U32, U22, U12 );
		if ( cases == 1 ) return 0;
	} else {
		rff = rff_p / two;
	}
	//c (2.12) Carlson (1989)
	I1c = rff_p; //two*rff
	if ( index_p4[3] == 2 ) {
		//c (2.13) Carlson (1989)
		I2c = two / three * d12 * d13 * rd(U32,U22,U12) + 
			two * X1 * Y1 / U1c;
		//(2.39) Carlson (1989)
		integ[3] = (b4*I2c-d14*I1c)/b1;
		if ( cases == 3 ) return 0;
	}
	if ( X1 * Y1 != zero ) {
		// (2.14) Carlson (1989)
		I3c = two * rc(P22,Q22) - two * b1 * d12 * d13 / three / 
			d14 * rj(U32,U22,U12,W22);
	} else {
		// One can read the paragraph between 
		// (2.19) and (2.20) of Carlson (1989).
		I3c = - two * b1 * d12 * d13 / three / d14 * 
			rj(U32, U22, U12, W22);
	}

	if ( index_p4[2] == -2 ) {
		// (2.49) Carlson (1989)
		integ[2]=(b4*I3c-b1*I1c)/d14;
		//printf("int2 = %f  %f  %f  %f  %f \n", b4, I3c, b1, I1c, d14);
		if ( cases == 2 ) return 0;
	}

	if ( index_p4[4] == -4 ) {
		/* (2.1)  Carlson (1989) */
		r12=a1/b1-a2/b2;
		r13=a1/b1-a3/b3;
		r24i=b2*b4/(a2*b4-a4*b2);
		r34i=b3*b4/(a3*b4-a4*b3);
		if ( x < infinity ) {
			/* (2.17) Carlson (1989) */
			J2c = two / three * d12 * d13 * rd(U32,U22,U12) +
				two * d13 * X2 * Y2 / X3 / Y3 / U1c;
			/* (2.59) & (2.6) Carlson (1989) */
			K2c = b3 * J2c - two * d34 * ( X1 * X2 / X3 / pow(X4, 2) -
				Y1 * Y2 / Y3 / pow(Y4, 2) );
		} else {
			J2c = two / three * d12 * d13 * rd(U32,U22,U12) + 
				two * d13 * Y2 / b3 / Y3 / Y1;
			K2c = b3 * J2c + two * d34 * Y1 * Y2 / Y3 / pow(Y4, 2);
		}
		//c (2.62) Carlson (1989)
		integ[4] = - I3c * half / d14 * ( one / r14 + one / r24 + 
				one / r34) + half * b4 / d14 / d24 / d34 * K2c + 
				pow( (b1 / d14), 2) * ( one - half * 
				r12 * r13 * r24i * r34i ) * I1c;
		/*printf("int4= %f  %f  %f \n", - I3c * half / d14 * ( one / r14 + one / r24 + 
				one / r34), half * b4 / d14 / d24 / d34 * K2c,
				pow( (b1 / d14), 2) * ( one - half * 
				r12 * r13 * r24i * r34i ) * I1c);*/
	}  
	return 0;
}



/*
!*     PURPOSE: Computes J_k(h)=\int_y^x dt (b4*t+a4)^(k/2)[(b1*t+a1)*(h*t^2+g*t+f)]^{-1/2}. 
!*              It is the case of equation W(t)=4*t^3-g_2*t-g_3=0 has one real root. 
!*              Equation (62) of Yang & Wang (2012).  
!*     INPUTS:  Arguments for above integral.
!*     OUTPUTS:  Value of integral. If index_p4(1)=0, then J_0 will be computed, 
!*               else J_0 will be replaced by parameter p, and rff_p=p. 
!*     ROUTINES CALLED: RF,RJ,RC,RD
!*     ACCURACY:   Machine.
!*     AUTHOR:     Dexter & Agol (2009)
!*     DATE WRITTEN:  4 Mar 2009
!*     MODIFIED:   Yang & Wang (2012)
!*     C version:  Yang Xiao-lin   2022-11-22
!*/
int ellcubiccomplexs(int *index_p4, double a1, double b1, double a4, \
		double b4, double f, double g, double h, double y, \
		double x, double rff_p, double *integ, int cases)
{
	double X1, X4, Y1, Y4, d14, beta1, beta4, a11, c44, 
		a142, xi, eta, M2, Lp2, Lm2, I1c, U, U2, Wp2, W2, 
		Q2, P2, rho, I3c, I2c, r24Xr34, r12Xr13, N2c, K2c,
		rff, rdd, r14, Ay1;
 
	X1 = sqrt( fabs(a1+b1*x) );
	X4 = sqrt( fabs(a4+b4*x) );
	Y1 = sqrt( fabs(a1+b1*y) );
	Y4 = sqrt( fabs(a4+b4*y) );

	r14 = a1 / b1 - a4 / b4;
	d14 = a1 * b4 - a4 * b1;
	/* (2.2) Carlson (1991) */
	beta1 = g * b1 - two * h * a1;
	beta4 = g * b4 - two * h * a4;
	/*c (2.3) Carlson (1991) */
	a11 = sqrt( two*f*b1*b1 - two*g*a1*b1 + two*h*a1*a1 );
	c44 = sqrt( two*f*b4*b4 - two*g*a4*b4 + two*h*a4*a4 );
	a142 = two*f*b1*b4 - g*( a1*b4 + a4*b1 ) + two*h*a1*a4;
	/* c (2.4) Carlson (1991) */
	xi = sqrt( f + g*x + h*x*x );
	eta = sqrt( f + g*y + h*y*y );
	/* write(*,*)'ellcubiccomplexs= ',f,g,x,f+g*x+h*x*x,eta */
	/* (3.1) Carlson (1991): */
	if ( x < infinity )
		M2 = pow( (X1 + Y1) * sqrt( pow( (xi + eta), 2 ) - 
			h * pow( (x - y), 2) )/(x-y), 2);
	else 
		M2 = b1 * ( two * sqrt(h) * eta + g + two * h * y );

	/* (3.2) Carlson (1991): */
	Lp2 = M2 - beta1 + sqrt( two * h ) * a11;
	Lm2 = M2 - beta1 - sqrt( two * h ) * a11;

	if ( index_p4[1] == 0 ) {
		/* (1.2) Carlson (1991) */
		rff = rf(M2,Lm2,Lp2);
		integ[1] = four * rff;
		if ( cases == 1 ) return 0;
	} else
		rff = rff_p / four;

	/* (3.8)  Carlson (1991) */
	I1c = rff_p; //four*rff;
	/* (3.3) Carlson (1991) */
	if ( x < infinity )
		U = ( X1 * eta + Y1 * xi ) / ( x - y );
	else {
		U = sqrt( h ) * Y1;
		/* if(Y1 == zero) U = 1.D-9; */ 
	}

	/*c (3.5) Carlson (1991) */
	rho = sqrt( two * h ) * a11 - beta1;
	rdd = rd(M2, Lm2, Lp2);
	if ( index_p4[3] == 2 ) {
		/* (3.9) Carlson (1991) */
		I2c = a11 * sqrt( two / h ) / three * ( four * rho * rdd - 
			six * rff + three / U ) + two * X1 * Y1 / U;
		/* (2.39) Carlson (1989) */ 
		integ[3] = ( b4 * I2c - d14 * I1c ) / b1;
		if ( cases == 3 ) return 0;
	}

	U2 = U * U;
	Wp2 = M2 - b1 * ( a142 + a11 * c44 ) / d14;
	W2 = U2 - pow(a11, 2) * b4 / two / d14;
	//printf( "W2  U2 = %f \n  a11 = %f \n  b4 = %f \n d14 = %f \n ", U2, a11, b4, d14 );
	/* (3.4) Carlson (1991) */
	if (x < infinity ) {
		Q2 = pow( ( X4 * Y4 / X1 / Y1 ), 2) * W2;
		//printf( "X4 = %f \n  Y4 = %f \n  X1 = %f \n Y1 = %f \n ", X4, Y4, X1, Y1 );
		//printf( "W2= %f \n  ", W2 );
	} else {
		Q2 = ( b4 / b1 ) * pow( ( Y4 / Y1 ), 2) * W2;
		//printf( "X4 = %f \n  Y4 = %f \n  X1 = %f \n Y1 = %f \n ", b4, b1, Y4, Y1 );
		//printf( "W2= %f \n  ", W2 );
	}

	P2 = Q2 + pow(c44, 2) * b4 / two / d14;
	/* (3.9) 1991 */
	if ( X1*Y1 != 0.0 ) {
		I3c = ( two * a11 / three / c44 ) * ( ( - four * b1 / d14 ) * 
			( a142 + a11 * c44 ) * rj( M2, Lm2, Lp2, Wp2 ) 
			 - six * rff + three * rc(U2, W2) ) + two * rc(P2, Q2);
		//printf( "tmp1 = %f \n  tmp2 = %f \n tmp3 = %f \n ", rj( M2, Lm2, Lp2, Wp2 ), rc(U2, W2), rc(P2, Q2) );
		//printf( "U2 = %f \n  W2 = %f \n  P2 = %f \n Q2 = %f \n ", U2, W2, P2, Q2 );
	} else {
		I3c = ( two * a11 / three / c44 ) * ( ( - four * b1 / d14 ) * 
			( a142 + a11 * c44 ) * rj( M2, Lm2, Lp2, Wp2 )
                        - six * rff + three * rc(U2, W2) );
		//printf( "tmp11 = %f \n  tmp22 = %f \n ", rj( M2, Lm2, Lp2, Wp2 ), rc(U2, W2) );
	}


	if ( index_p4[2] == -2 ) {
		/* (2.49) Carlson (1989) */
		if ( Y4 != zero && X4 != zero )
			integ[2] = ( b4 * I3c - b1 * I1c ) / d14;
		else {
			integ[2] = infinity;
			//printf( "tmp1 = %f \n ", integ[2]);
		}

		//printf( "tmp1 = %f \n  tmp2 = %f \n  tmp3 = %f \n  tmp4 = %f \n  tmp5 = %f \n  \n", b4,  I3c, b1, I1c, d14 );

		if ( cases == 2 ) return 0;
	}

	if ( index_p4[4] == -4 ) {
		/* (2.19) Carlson (1991) */
		r24Xr34 = half * pow(c44, 2) / h / pow(b4, 2);
		r12Xr13 = half * pow(a11, 2) / h / pow(b1, 2);
		/* (3.11) Carlson (1991) */                 
		N2c = two / three * sqrt( two * h ) / a11 * 
			( four * rho * rdd - six * rff ); 
			/* +three/U)!+two/X1/Y1/U */
		if ( Y1 == zero ) {
			Ay1 = - two * b1 * xi / X1 + sqrt( two * h ) * a11 / U;
			if( x >= infinity )Ay1 = zero;
			/* two*d14*eta/pow(Y4, 2)/U+dsqrt(two*h)*a11/U */
		} else {
			Ay1 = ( two * d14 * eta / pow(Y4, 2) + 
				pow(a11, 2) / X1 / U ) / Y1 + 
				sqrt( two * h ) * a11 / U;
		}

		/* (2.5) & (3.12) Carlson (1991) */
		K2c = half * pow(a11, 2) * N2c - two * d14 * 
				( xi / X1 / pow(X4, 2) ) + Ay1; 
		/* (2.62) Carlson (1989) */
		if ( Y4 != zero && X4 != zero )
			integ[4] = -I3c * half / d14 * ( one / r14 + 
				two * b4 * beta4 / pow(c44, 2) ) + 
				half / d14 / ( h * b4 * r24Xr34 ) * K2c + 
				pow(( b1 / d14), 2) * ( one - half * 
				r12Xr13 / r24Xr34 ) * I1c;
		else
			integ[4] = infinity;
	}
	return 0;
}


/*
!*     PURPOSE: Computes J_k(h)=\int_y^x dt (f_1+g_1t+h_1t^2)^{p_1/2} 
!*                       (f_2+g_2t+h_2t^2)^{p_2/2} (a_5+b_5t)^{p_5/2}. 
!*     INPUTS:  Arguments for above integral.
!*     OUTPUTS:  Value of integral. If index_p5(1)=0, then J_0 will be computed, 
!*               else J_0 will be replaced by parameter p, and rff_p=p. 
!*     ROUTINES CALLED: RF,RJ,RC,RD
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Dexter & Agol (2009)
!*     DATE WRITTEN:  4 Mar 2009
!*     MODIFIED:   Yang & Wang  (2012)
!*     C VERSION:   Yang Xiao-lin   2022-11-22
!*     REVISIONS: [-1,-1,-1,-1,0]=integ(1),[-1,-1,-1,-1,-2]=integ(2),
!*                [-1,-1,-1,-1,-4]=integ(4),[-1,-1,-1,-1,2]=integ(3),
!*                [-1,-1,-1,-1,4]=integ(5)
!*/
static int elldoublecomplexs( int *index_p5, double f1, double g1, double h1, 
		double f2, double g2, double h2, double a5, double b5,
		double y, double x, double rff_p, double *integ, int cases)
{
	double xi1, xi2, eta1, eta2;
	double theta1, theta2, zeta1, zeta2, M, M2;
	double delta122, delta112, delta222, delta, deltap, Lp2, Lm2;
	double deltam, rff, U, U2, alpha15, beta15, alpha25, 
		beta25, Lambda, Omega2, psi, xi5, eta5;
	double gamma1, gamma2, Am111m1, XX;  //A1111m4
	double S, mu, T, V2, b2, a2, H, A1111m2, xi1p, B, G;
	double Sigma;
	double H0, S2, T2, eta1p, psi2, A1111;

	/* (2.1) Carlson (1992) */
	if ( x < infinity ) {
		xi1 = sqrt( f1 + g1*x + h1*x*x );
		xi2 = sqrt( f2 + g2*x + h2*x*x );
	} else {
		xi1 = x * sqrt( h1 );
		xi2 = x * sqrt( h2 );
	}
           
	eta1 = sqrt( f1 + g1 * y + h1 * y * y );
	eta2 = sqrt( f2 + g2 * y + h2 * y * y );
	/* (2.4) Carlson (1992) */
	if ( x < infinity ) {
		theta1=two*f1+g1*(x+y)+two*h1*x*y;
		theta2=two*f2+g2*(x+y)+two*h2*x*y;
	} else {
		theta1=(g1+two*h1*y)*x;
		theta2=(g2+two*h2*y)*x;
	}
	/* (2.5) Carlson (1992) */
	zeta1 = sqrt( two * xi1 * eta1 + theta1 );
	zeta2 = sqrt( two * xi2 * eta2 + theta2 );
	/* (2.6) Carlson (1992) */
	if ( x < infinity ) {
		M = zeta1 * zeta2 / ( x - y );
		M2 = M * M;
	} else
		M2 = ( two * sqrt( h1 ) * eta1 + g1 + two * h1 * y ) * 
			(two * sqrt( h2 ) * eta2 + g2 + two * h2 * y );


	/* (2.7) Carlson (1992) */
	delta122 = two*f1*h2 + two*f2*h1 - g1*g2;
	delta112 = four*f1*h1 - g1*g1;
	delta222 = four*f2*h2 - g2*g2;
	delta = sqrt( delta122 * delta122 - delta112 * delta222 );
	/* (2.8) Carlson (1992) */
	deltap = delta122 + delta;
	deltam = delta122 - delta;
	Lp2 = M2 + deltap;
	Lm2 = M2 + deltam;
 
	if ( index_p5[1] == 0 ) {
		rff = rf( M2, Lm2, Lp2 );
		/* (2.36) Carlson (1992) */
		integ[1] = four * rff;
		if( cases == 1 ) return 0;
	} else
		rff = rff_p / four;

	/* (2.6) Carlson (1992) */
	if ( x < infinity ) {
		U = ( xi1 * eta2 + eta1 * xi2 ) / ( x - y );
		U2 = U * U;
	} else {
             U = sqrt( h1 ) * eta2 + sqrt( h2 ) * eta1;
             U2 = U * U;
	}
        
	/* (2.11) Carlson (1992) */
	alpha15=two*f1*b5-g1*a5;
	alpha25=two*f2*b5-g2*a5;
	beta15=g1*b5-two*h1*a5;
	beta25=g2*b5-two*h2*a5;
	/* (2.12) Carlson (1992) */
        gamma1=half*(alpha15*b5-beta15*a5);
        gamma2=half*(alpha25*b5-beta25*a5);
	/* (2.13) Carlson (1992) */
        Lambda=delta112*gamma2/gamma1;
        Omega2=M2+Lambda;
        psi=half*(alpha15*beta25-alpha25*beta15);
        psi2=psi*psi;
	/* (2.15) Carlson (1992) */
	xi5=a5+b5*x;
	eta5=a5+b5*y;
	/* (2.16) Carlson (1992) */
	if ( x < infinity ) {
		Am111m1=one/xi1*xi2-one/eta1*eta2;
		//A1111m4=xi1*xi2/pow(xi5, 2)-eta1*eta2/pow(eta5, 2);
		A1111m2=xi1*xi2/xi5-eta1*eta2/eta5;
		XX = (xi5*eta5*theta1*half*Am111m1-xi1*xi2*pow(eta5, 2)+
                                  eta1*eta2*pow(xi5, 2))/pow((x-y), 2);
		mu=gamma1*xi5*eta5/xi1/eta1;
	} else {
		Am111m1 = sqrt(h2/h1)-one/eta1*eta2;
		//A1111m4 = sqrt(h1*h2)/b5/b5 - eta1*eta2/pow(eta5, 2);
		A1111m2 = sqrt(h1*h2)*x/b5-eta1*eta2/eta5;
		XX=b5*eta5*h1*y*Am111m1-pow(eta5, 2) * sqrt(h1*h2) + eta1*eta2*b5*b5;
		mu=gamma1*b5*eta5 / sqrt(h1)/eta1;
	}

	/* (2.17) Carlson (1992)
	 * (2.18) Carlson (1992) */
	S=half*(M2+delta122)-U2;
	S2=S*S;
	/* (2.19) Carlson (1992) */
        T=mu*S+two*gamma1*gamma2;
        T2=T*T;
        V2=mu*mu*(S2+Lambda*U2);
	/* (2.20) Carlson (1992) */
        b2=pow(Omega2, 2) * (S2/U2+Lambda);
        a2=b2+ pow(Lambda, 2)*psi2/gamma1/gamma2;
	/* (2.22) Carlson (1992) */
        H=delta112*psi*(rj(M2,Lm2,Lp2,Omega2)/three+ 
             half*rc(a2,b2))/pow(gamma1, 2)-XX*rc(T2,V2);
	//printf("H1 = %f \n %f \n %f \n %f \n %f \n %f \n %f \n ", delta112, 
	//	psi, rj(M2,Lm2,Lp2,Omega2), rc(a2,b2), sq(gamma1), XX, rc(T2,V2));
	//printf("H2 = %f \n %f \n %f \n", rc(T2,V2), T2, V2);
	if ( index_p5[3] == 2 || index_p5[5] == 4 ) {
		/* (2.23)--(2.29) Carlson (1992) */
		psi=g1*h2-g2*h1;
		Lambda=delta112*h2/h1;
		Omega2=M2+Lambda;
		Am111m1=one/xi1*xi2-one/eta1*eta2;
		A1111=xi1*xi2-eta1*eta2;
		XX=(theta1*half*Am111m1-A1111)/pow((x-y), 2);
		b2=pow(Omega2, 2) * (S2/U2+Lambda);
		a2=b2 + pow(Lambda, 2) * pow(psi, 2) /h1/h2;
		mu=h1/xi1/eta1;
		T=mu*S+two*h1*h2;
		T2=T*T;
		V2=mu * mu * (S2+Lambda*U2);
		H0=delta112*psi*(rj(M2,Lm2,Lp2,Omega2)/three+ 
			half*rc(a2,b2))/pow(h1, 2)-XX*rc(T2,V2);
		if ( index_p5[3] == 2 ) {
			/* (2.42) Carlson (1992) */
			integ[3] = two*b5*H0-two*beta15*rff/h1;
			if ( cases == 3 ) return 0;
		}

		if ( index_p5[5] == 4 ) {
			/* (2.2) Carlson (1992) */
			if (x < infinity)
				xi1p = half*(g1 + two*h1*x)/xi1;
			else
				xi1p = sqrt(h1);

			eta1p=half*(g1+two*h1*y)/eta1;
			/* (2.3) Carlson (1992) */
			B=xi1p*xi2-eta1p*eta2;
			/* (2.9) Carlson (1992) */
			if ( x < infinity ) {
				G=two/three*delta*deltap*rd(M2,Lm2,Lp2)+half*delta/U+ 
					(delta122*theta1-delta112*theta2)/four/xi1/eta1/U;
			} else {
				G=two/three*delta*deltap*rd(M2,Lm2,Lp2)+half*delta/U+ 
				(delta122*(g1+two*h1*y)-delta112*(g2+two*h2*y))/four/ sqrt(h1)/eta1/U;
			}
			/* (2.10) Carlson (1992) */  
			Sigma=G-deltap*rff+B;
			/* (2.44) Carlson (1992) */
			integ[5]=-b5*(beta15/h1+beta25/h2)*H0 + b5*b5* 
					Sigma/h1/h2 + pow(beta15, 2)*rff/pow(h1, 2);
			if ( cases == 5 ) return 0;
		}            
	}


	if ( index_p5[2] == -2 ) {
		/* (2.39) Carlson (1992) */
		integ[2] = -two*(b5*H+beta15*rff/gamma1);
		//printf("fffs = %f \n = %f \n = %f \n = %f \n = %f \n = %d \n = %f \n ", 
		//		b5, H, beta15, rff, gamma1, cases, integ[2]);
		if ( cases == 2 ) return 0;
	}

	if ( index_p5[4] == -4) {
		/* (2.2) Carlson (1992) */
		if ( x < infinity )
			xi1p=half*(g1+two*h1*x)/xi1;
                else
			xi1p= sqrt(h1);

                eta1p=half*(g1+two*h1*y)/eta1;
		/* (2.3) Carlson (1992) */
                B=xi1p*xi2-eta1p*eta2;
		/* (2.9) Carlson (1992) */
		if ( x < infinity )
			G=two/three*delta*deltap*rd(M2,Lm2,Lp2)+half*delta/U+
				(delta122*theta1-delta112*theta2)/four/xi1/eta1/U;
                else
			G=two/three*delta*deltap*rd(M2,Lm2,Lp2)+half*delta/U+
				(delta122*(g1+two*h1*y)-delta112*(g2+two*h2*y))/four / sqrt(h1)/eta1/U;

		/* (2.10) Carlson (1992) */
		Sigma=G-deltap*rff+B;
		/* (2.41) Carlson (1992) */
		integ[4]=b5*(beta15/gamma1+beta25/gamma2)*H + pow( beta15, 2) * rff / pow(gamma1, 2) +
                                   pow( b5, 2) * (Sigma-b5*A1111m2)/gamma1/gamma2;
	}
	return 0;
}






