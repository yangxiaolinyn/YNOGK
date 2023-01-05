/********************************************************************
!*    This module defines many constants often uesd in our code.
!*    One can use these constants through a command "use constants" in their
!*    own subroutines or functions. 
\********************************************************************/  

#ifndef CONSTANTS_H
#define CONSTANTS_H

#define max(a, b) ((a)>(b)?(a):(b))
#define sq(a) pow((a), 2)
#define sq3(a) pow((a), 3)

#define true 1
#define false 0


#define zero     0.0
#define one      1.0
#define two      2.0
#define three    3.0
#define four     4.0
#define five     5.0
#define six      6.0
#define seven    7.0
#define eight    8.0
#define nine     9.0
#define ten      10.0
#define sixteen  16.0
#define sqrt2    1.4142135623730951
#define sqrt3    1.7320508075688772
//unsigned const int krand1 sizeof 1.0
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define pi 3.14159265358979323848280
#define twopi 6.283185307179586476920
#define halfpi 1.570796326794896619231321691630
/*sqrtpi 1.772453850905516027298167483341145182797549456122387128213...*/
#define sqrtpi 1.772453850905516027290
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define half 0.50
#define half2 0.250
#define sigma_T 6.652448e-25
#define barn 1.e-24
#define mec2 0.5110
#define log10_mec2 log10( 0.5110 )  //In unit of MeV.
#define Boltzman_Costant_k 1.380662e-16
#define electron_mass 9.109534e-28
#define electron_Charge 4.803242e-10
#define ln_e 2.7182818284590452353600
#define infinity 1.e100 
#define Boltzman_k_ergK 1.380662e-16
#define Boltzman_k_evK 1.380662e-16 / 1.6021892e-12
#define dtors (3.14159265358979323848280/180.0)
#define dtor (3.14159265358979323848280/180.0)
#define mh 1.6726231e-24
#define mp 1.6726485e-24
#define hbar 1.0545887e-27    // h / 2 pi
#define planck_h 6.626178e-27   // erg * s
#define h_ev 4.13566743e-15   // eV·s
#define pho_v 2.99792458e10 
#define erg_of_one_ev 1.6021892e-12
#define Sigma_SB 5.67032e-5
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#define Gnewton 6.6720e-8
#define Msun 1.9891e33
#define Cv 2.99792458e10
#define rg_SUN 6.6720e-8 * 1.9891e33 / Cv / Cv
/*****************************************************************************/
/*****************************************************************************/
 
#endif
