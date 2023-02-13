// 
// complex.h
// 
// pheno
// Feb 09, 2023
// Feb 11, 2023
// 

#pragma once

#include <math.h>

#ifndef PI
#define PI (3.1415926535897932384626433832795028841971693993751058209749445923)
#endif

typedef struct {
	double real;
	double imag;
} complex;

typedef struct {
	double radial;
	double angular;
} polar;

polar	c2p(complex z);
double  cabs(complex z);
complex cacos(complex z);
complex cacosh(complex z);
complex cadd(complex x, complex y);
double  carg(complex z);
complex casin(complex z);
complex casinh(complex z);
complex catan(complex z);
complex catanh(complex z);
complex ccos(complex z);
complex ccosh(complex z);
complex cdiv(complex x, complex y);
complex cexp(complex z);
complex clog(complex z);
complex clog10(complex z);
complex clogn(complex z, double n);
complex cmul(complex x, complex y);
complex conj(complex z);
complex cpow(complex x, complex y);
complex csin(complex z);
complex csinh(complex z);
complex csub(complex x, complex y);
complex csqrt(complex z);
void	cswap(complex* a, complex* b);
complex ctan(complex z);
complex ctanh(complex z);
double  norm(complex z);
complex p2c(polar p);
polar	padd(polar a, polar b);
polar	pdiv(polar a, polar b);
polar	pmul(polar a, polar b);
polar	psub(polar a, polar b);
