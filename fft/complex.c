// 
// complex.c
// 
// pheno
// Feb 09, 2023
// Feb 11, 2023
// 

#include "complex.h"

polar	c2p(complex z) {
	polar p = { cabs(z), carg(z) };
	return p;
}

double  cabs(complex z) {
	return sqrt(norm(z));
}

complex cacos(complex z) {
	complex I = { 0., 1. };
	complex one = { 1., 0. };
	complex zsquare = cmul(z, z);
	return cmul(cdiv(one, I), clog(cadd(z, cmul(I, csqrt(csub(one, zsquare))))));
}
complex cacosh(complex z) {
	complex I = { 0., 1. };
	complex one = { 1., 0. };
	complex zsquare = cmul(z, z);
	complex minusone = { -1., 0. };
	return cmul(((z.real < 0) ? minusone : one), clog(cadd(z, csqrt(csub(zsquare, one)))));
}

complex cadd(complex x, complex y) {
	complex r = { x.real + y.real, x.imag + y.imag };
	return r;
}

double  carg(complex z) {
	if (z.real > 0) return atan(z.imag / z.real);
	else if (z.real == 0) return(z.imag < 0) ? -PI / 2 : PI / 2;
	else return (z.imag < 0) ? atan(z.imag / z.real) - PI : PI + atan(z.imag / z.real);
}

complex casin(complex z) {
	complex I = { 0., 1. };
	complex one = { 1., 0. };
	complex zsquare = cmul(z, z);
	return cmul(cdiv(one, I), clog(cadd(cmul(I, z), csqrt(csub(one, zsquare)))));
}

complex casinh(complex z) {
	complex I = { 0., 1. };
	complex one = { 1., 0. };
	complex zsquare = cmul(z, z);
	return clog(cadd(z, csqrt(cadd(one, zsquare))));
}

complex catan(complex z) {
	complex I = { 0., 1. };
	complex one = { 1., 0. };
	complex doublei = { 0., 2. };
	return cmul(cdiv(one, doublei), clog(cdiv(csub(I, z), cadd(I, z))));
}

complex catanh(complex z) {
	complex one = { 1., 0. };
	complex half = { .5, 0. };
	return cmul(half, clog(cdiv(cadd(one, z), csub(one, z))));
}

complex ccos(complex z) {
	complex r = { cos(z.real) * cosh(z.imag), -sin(z.real) * sinh(z.imag) };
	return r;
}

complex ccosh(complex z) {
	complex r = { cosh(z.real) * cos(z.imag), sinh(z.real) * sin(z.imag) };
	return r;
}

complex cdiv(complex x, complex y) {
	return p2c(pdiv(c2p(x), c2p(y)));
}

complex cexp(complex z) {
	complex r = { exp(z.real) * cos(z.imag), exp(z.real) * sin(z.imag) };
	return r;
}

complex clog(complex z) {
	polar p = c2p(z);
	complex r = { log(p.radial), p.angular };
	return r;
}

complex clog10(complex z) {
	complex t = clog(z);
	complex r = { t.real / log(10), t.imag / log(10) };
	return r;
}

complex clogn(complex z, double n) {
	complex t = clog(z);
	complex r = { t.real / log(n), t.imag / log(n) };
	return r;
}

complex cmul(complex x, complex y) {
	return p2c(pmul(c2p(x), c2p(y)));
}

complex conj(complex z) {
	complex r = { z.real, -z.imag };
	return r;
}

complex cpow(complex x, complex y) {
	return cexp(cmul(y, clog(x)));
}

complex csin(complex z) {
	complex r = { sin(z.real) * cosh(z.imag), cos(z.real) * sinh(z.imag) };
	return r;
}

complex csinh(complex z) {
	complex r = { sinh(z.real) * cos(z.imag), cosh(z.real) * sin(z.imag) };
	return r;
}

complex csub(complex x, complex y) {
	complex r = { x.real - y.real, x.imag - y.imag };
	return r;
}

complex csqrt(complex z) {
	if (z.imag != 0) {
		complex r = {
			sqrt((cabs(z) + z.real) / 2),
			z.imag / fabs(z.imag) * sqrt((cabs(z) - z.real) / 2)
		};
		return r;
	}
	else {
		complex r = { sqrt(z.real), 0. };
		return r;
	}
}

void cswap(complex* a, complex* b) {
	complex tmp = { 0., 0. };
	tmp = *a;
	*a = *b;
	*b = tmp;
	return;
}

complex ctan(complex z) {
	double v = cos(2 * z.real) + cosh(2 * z.imag);
	complex t = { sin(2 * z.real), sinh(2 * z.imag) };
	complex r = { t.real / v, t.imag / v };
	return r;
}

complex ctanh(complex z) {
	double v = cosh(2 * z.real) + cos(2 * z.imag);
	complex t = { sinh(2 * z.real), sin(2 * z.imag) };
	complex r = { t.real / v, t.imag / v };
	return r;
}

double  norm(complex z) {
	return z.real * z.real + z.imag * z.imag;
}

complex p2c(polar p) {
	complex r = { p.radial * cos(p.angular), p.radial * sin(p.angular) };
	return r;
}

polar	padd(polar a, polar b) {
	return c2p(cadd(p2c(a), p2c(b)));
}

polar	pdiv(polar a, polar b) {
	polar r = { a.radial / b.radial, a.angular - b.angular };
	return r;
}

polar	pmul(polar a, polar b) {
	polar r = { a.radial * b.radial, a.angular + b.angular };
	return r;
}

polar	psub(polar a, polar b) {
	return c2p(csub(p2c(a), p2c(b)));
}
