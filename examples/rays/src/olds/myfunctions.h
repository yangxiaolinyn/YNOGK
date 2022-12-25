#include "1RandUtils.h"
#include <stdio.h>
#include <math.h>

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Functions' declarations defined in my code.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void print_time(int hour, int minute);
void print_time_globl(void);
int increment(int x);
int hour = 5, minute = 233;
//int total_minute = hour * 60 + minute;
//double total_hour = hour + minute / 60.0;
int x = 10;
//double pi = acos(-1.0);


void fool(void);
void print_parity(int x);


int is_leap_year(int year);
long long int factorial(long long int n);


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Functions' definitions in my code.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void print_time(int hour, int minute)
{
	printf("hour=%d \nminute=%d \n", hour, minute);
	int total_minute = hour * 60 + minute;
	printf("minutes=%d \n", total_minute);
	if (total_minute < 100) {
		printf("Total minutes are1:%d\n", total_minute);
	} else {
		printf("Total minutes are2:%d\n", total_minute);
		//int a = 1, b = 2, c = 3; 
	}
}



void print_time_globl(void)
{
	printf("hour=%d \nminute=%d \n", hour, minute);
	printf("minutes=%d \n", hour * 60 + minute);
}


int increment(int x)
{
	x = x + 1;
	return x;
}


void fool(void)
{
	int i = 899;
	printf("i = %d \n", i);
	i = 777;
}


void print_parity(int x)
{
	if (x % 2 == 0) {
		printf("The integer is even %d \n", x);
	} else
		printf("The integer is odd %d \n", x);
}


int is_leap_year(int year)
{
	if (year % 400 == 0) {
		printf("The year %d is leap year:\n", year);
		return 1;
	} else {
		if (year % 4 == 0 && year % 100 != 0) {
			printf("The year %d is leap year:\n", year);
			return 1;
		} else {
			printf("The year %d is not leap year:\n", year);
			return 0;
		}
	}
}


long long int factorial(long long int n)
{
	if (n == 0)
		return 1;
	else {
		if (n > 0)
			return n * factorial(n - 1);
		else
			return -1;
	}
}

long double ldfactorial(long double x)
{
	if (x == 0.0)
		return 1.0;
	else {
		if (x > 0)
			return x * factorial(x - 1.0);
		else
			return -1;
	}
}


long long int is_prime(long long int n)
{
	long long int i;
	for (i = 2; i <= sqrt(n); i++) {
		if (n % i == 0)
			return 0;
	}

	//printf(" i = %lld \n n=%lld\n", i, n);
	return 1;
}


struct complex_struct {
	double x, y;
};


void show_complex(struct complex_struct z)
{
	printf("complex number z = ( %f, %f )\n", z.x, z.y);
}


/**********************************\
*
*
*
\**********************************/
