/*  module RandUtils
*    //use MpiUtils
*    implicit none
*/
/*
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "ConstantsUnit.h"

#include<sys/time.h> */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "ConstantsUnits.h"

#include<sys/time.h>

#include "1RandUtils.h"


static double U[98], C, CD, CM, S, T;
static long I97, J97;

 
/*static long int rand_inst; 
static long Rand_Feedback;*/

//static double U[98], C, CD, CM, S, T;
//static long I97, J97;
 
//int rand_inst = 0;
//int const krand = sizeof 1.0;
//int Rand_Feedback = 1;

/*    interface RandRotation
    module procedure RandRotationS, RandRotationD
    end interface

    contains
*/
 
//void RMARIN( const long IJ, const long KL );

    //subroutine init_random_seed()
    //integer, allocatable :: seed(:)
    //integer :: i, n, un, istat, dt(8), pid, t(2), s
    //integer(8) :: count, tms
    //
    //call random_seed(size = n)
    //allocate(seed(n))
    //// First try if the OS provides a random number generator
    //// Fallback to XOR:ing the current time and pid. The PID is
    //// useful in case one launches multiple instances of the same
    //// program in parallel.
    //call system_clock(count)
    //if (count /= 0) then
    //    t = transfer(count, t)
    //else
    //    call date_and_time(values=dt)
    //    tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
    //    + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
    //    + dt(3) * 24 * 60 * 60 * 60 * 1000 &
    //    + dt(5) * 60 * 60 * 1000 &
    //    + dt(6) * 60 * 1000 + dt(7) * 1000 &
    //    + dt(8)
    //    t = transfer(tms, t)
    //end if
    //s = ieor(t(1), t(2))
    //pid = getpid() + 1099279 // Add a prime
    //s = ieor(s, pid)
    //if (n >= 3) then
    //    seed(1) = t(1) + 36269
    //    seed(2) = t(2) + 72551
    //    seed(3) = pid
    //    if (n > 3) then
    //        seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
    //    end if
    //else
    //    seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
    //end if
    //call random_seed(put=seed)
    //end subroutine init_random_seed


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	void RMARIN( const long IJ, const long KL )
	{ 
		if( (IJ < 0) || (IJ > 31328) || (KL < 0) || (KL > 30081) ) { 
			printf("The first random number seed must \
                               have a value between 0 and 31328.\n");
			printf("The second seed must have a value \
				between 0 and 30081. \n");
			printf("T = %ld \t %ld\n", IJ, KL);
			exit(100);
		};

		long I, J, K, L, M;
		I = ( IJ/177) % 177 + 2;
		J = IJ % 177 + 2;
		K = (KL/169) % 178 + 1;
		L = KL % 169;


		for(int II = 1; II < 98; II++) {
			S = 0.0;
			T = 0.5;
			for(int JJ = 1; JJ < 25; JJ++) {  
				M = ( ( (I*J) % 179 )*K ) % 179;
				I = J;
				J = K;
				K = M;
				L = (53*L+1) % 169;
				if ( ( (L*M) % 64 ) > 31 ) 
					S += T;
				T = 0.5 * T;
			}; 
			U[II] = S; 
		};
		C = 362436.0 / 16777216.0;
		CD = 7654321.0 / 16777216.0;
		CM = 16777213.0 /16777216.0;
		I97 = 97;
		J97 = 33;
		//U[0] = zero; 
		//for( int i=1; i < 98; i++)
		//   cout << U[i] << endl;
	}



int initRandom( int i, int i2 ) 
{ 
	int kl, ij;
	long int rand_inst = 1;  

	if ( i != -1 ) {
		if ( i2 != 9373 ) {
			kl=i2;
			if (i2 > 30081) {
				printf("initRandom: second seed too large");
				exit(0);
			}
		} else
			kl = 9373;
		ij = i;
	} else {
		struct timeval tv;

		gettimeofday(&tv, NULL);
 
		//long int tt;
		//tt = time(NULL);

		//tt = tv.tv_sec;
		//diff = tv.tv_usec;
		ij = (tv.tv_usec*tv.tv_usec + rand_inst*100) % 31328;
    
		kl = ( (tv.tv_usec * 1000) ) % 30081; 
		//printf("the seed numbers are: tt=%ld\n", tt);
		printf("tv.tv_usec=%ld\n  SEC: = %ld \n", tv.tv_usec, CLOCKS_PER_SEC);
	}

	printf("the seed numbers are: kl = %d \n ij=%d\n  \n", kl, ij);
	RMARIN( ij, kl ); 
	return 0;

    //end subroutine initRandom
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    double RANMAR()
    {
    // This is the random number generator proposed by George Marsaglia in
    // Florida State University Report: FSU-SCRI-87-50
    // It was slightly modified by F. James to produce an array of pseudorandom
    // numbers.
        //double U[98], C, CD, CM, S, T;
        //long I97, J97;
	//long I97, J97;
        double UNI;
 
        //     INTEGER IVEC
        UNI = U[I97] - U[J97];
        if( UNI < 0.0 )
            UNI += 1.0;
        U[I97] = UNI;
        I97 = I97 - 1;
        if(I97 == 0)
            I97 = 97;
        J97 = J97 - 1;
        if(J97 == 0)
            J97 = 97;
        C = C - CD;
        if( C < 0.0 )
            C = C + CM;
        UNI = UNI - C;
        if( UNI < 0.0 )
            UNI += 1.0; // bug?
        if( UNI == 0.0 || UNI == 1.0 )
            UNI = 0.2459428724394560; 
        return UNI;
    };













 
