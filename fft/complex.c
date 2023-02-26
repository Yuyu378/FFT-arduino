// 
// _Dcomplex.c
// 
// 

#include "complex.h"

_Dpolar	c2p(_Dcomplex z) {
	_Dpolar p = { cabs(z), carg(z) };
	return p;
}

double  cabs(_Dcomplex z) {
	return sqrt(norm(z));
}

_Dcomplex cacos(_Dcomplex z) {
	_Dcomplex one = { 1., 0. };
	_Dcomplex zsquare = cmul(z, z);
	return cmul(cdiv(one, J), clog(cadd(z, cmul(J, csqrt(csub(one, zsquare))))));
}
_Dcomplex cacosh(_Dcomplex z) {
	_Dcomplex one = { 1., 0. };
	_Dcomplex zsquare = cmul(z, z);
	_Dcomplex minusone = { -1., 0. };
	return cmul(((z.real < 0) ? minusone : one), clog(cadd(z, csqrt(csub(zsquare, one)))));
}

_Dcomplex cadd(_Dcomplex x, _Dcomplex y) {
	_Dcomplex r = { x.real + y.real, x.imag + y.imag };
	return r;
}

double  carg(_Dcomplex z) {
	if (z.real > 0) return atan(z.imag / z.real);
	else if (z.real == 0) return(z.imag < 0) ? -PI / 2 : PI / 2;
	else return (z.imag < 0) ? atan(z.imag / z.real) - PI : PI + atan(z.imag / z.real);
}

_Dcomplex casin(_Dcomplex z) {
	_Dcomplex one = { 1., 0. };
	_Dcomplex zsquare = cmul(z, z);
	return cmul(cdiv(one, J), clog(cadd(cmul(J, z), csqrt(csub(one, zsquare)))));
}

_Dcomplex casinh(_Dcomplex z) {
	_Dcomplex one = { 1., 0. };
	_Dcomplex zsquare = cmul(z, z);
	return clog(cadd(z, csqrt(cadd(one, zsquare))));
}

_Dcomplex catan(_Dcomplex z) {
	_Dcomplex one = { 1., 0. };
	_Dcomplex doublei = { 0., 2. };
	return cmul(cdiv(one, doublei), clog(cdiv(csub(J, z), cadd(J, z))));
}

_Dcomplex catanh(_Dcomplex z) {
	_Dcomplex one = { 1., 0. };
	_Dcomplex half = { .5, 0. };
	return cmul(half, clog(cdiv(cadd(one, z), csub(one, z))));
}

_Dcomplex* ccopy(_Dcomplex* z, unsigned N) {
    _Dcomplex* r = 0;
	if ((r = (_Dcomplex*)calloc(N, sizeof(_Dcomplex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	for (unsigned i = 0; i < N; i++) {
		*(r + i) = *(z + i);
	}
	return r;
}

_Dcomplex ccos(_Dcomplex z) {
	_Dcomplex r = { cos(z.real) * cosh(z.imag), -sin(z.real) * sinh(z.imag) };
	return r;
}

_Dcomplex ccosh(_Dcomplex z) {
	_Dcomplex r = { cosh(z.real) * cos(z.imag), sinh(z.real) * sin(z.imag) };
	return r;
}

_Dcomplex cdiv(_Dcomplex x, _Dcomplex y) {
	return p2c(pdiv(c2p(x), c2p(y)));
}

_Dcomplex cexp(_Dcomplex z) {
	_Dcomplex r = { exp(z.real) * cos(z.imag), exp(z.real) * sin(z.imag) };
	return r;
}

_Dcomplex clog(_Dcomplex z) {
	_Dpolar p = c2p(z);
	_Dcomplex r = { log(p.radial), p.angular };
	return r;
}

_Dcomplex clog10(_Dcomplex z) {
	_Dcomplex t = clog(z);
	_Dcomplex r = { t.real / log(10), t.imag / log(10) };
	return r;
}

_Dcomplex clogn(_Dcomplex z, double n) {
	_Dcomplex t = clog(z);
	_Dcomplex r = { t.real / log(n), t.imag / log(n) };
	return r;
}

_Dcomplex cmul(_Dcomplex x, _Dcomplex y) {
	return p2c(pmul(c2p(x), c2p(y)));
}

_Dcomplex conj(_Dcomplex z) {
	_Dcomplex r = { z.real, -z.imag };
	return r;
}

_Dcomplex cpow(_Dcomplex x, _Dcomplex y) {
	return cexp(cmul(y, clog(x)));
}

_Dcomplex csin(_Dcomplex z) {
	_Dcomplex r = { sin(z.real) * cosh(z.imag), cos(z.real) * sinh(z.imag) };
	return r;
}

_Dcomplex csinh(_Dcomplex z) {
	_Dcomplex r = { sinh(z.real) * cos(z.imag), cosh(z.real) * sin(z.imag) };
	return r;
}

_Dcomplex csub(_Dcomplex x, _Dcomplex y) {
	_Dcomplex r = { x.real - y.real, x.imag - y.imag };
	return r;
}

_Dcomplex csqrt(_Dcomplex z) {
	if (z.imag != 0) {
		_Dcomplex r = {
			sqrt((cabs(z) + z.real) / 2),
			z.imag / fabs(z.imag) * sqrt((cabs(z) - z.real) / 2)
		};
		return r;
	}
	else {
		_Dcomplex r = { sqrt(z.real), 0. };
		return r;
	}
}

void cswap(_Dcomplex* a, _Dcomplex* b) {
	_Dcomplex tmp = { 0., 0. };
	tmp = *a;
	*a = *b;
	*b = tmp;
	return;
}

_Dcomplex ctan(_Dcomplex z) {
	double v = cos(2 * z.real) + cosh(2 * z.imag);
	_Dcomplex t = { sin(2 * z.real), sinh(2 * z.imag) };
	_Dcomplex r = { t.real / v, t.imag / v };
	return r;
}

_Dcomplex ctanh(_Dcomplex z) {
	double v = cosh(2 * z.real) + cos(2 * z.imag);
	_Dcomplex t = { sinh(2 * z.real), sin(2 * z.imag) };
	_Dcomplex r = { t.real / v, t.imag / v };
	return r;
}

double  norm(_Dcomplex z) {
	return z.real * z.real + z.imag * z.imag;
}

_Dcomplex p2c(_Dpolar p) {
	_Dcomplex r = { p.radial * cos(p.angular), p.radial * sin(p.angular) };
	return r;
}

_Dpolar	padd(_Dpolar a, _Dpolar b) {
	return c2p(cadd(p2c(a), p2c(b)));
}

_Dpolar	pdiv(_Dpolar a, _Dpolar b) {
	_Dpolar r = { a.radial / b.radial, a.angular - b.angular };
	return r;
}

_Dpolar	pmul(_Dpolar a, _Dpolar b) {
	_Dpolar r = { a.radial * b.radial, a.angular + b.angular };
	return r;
}

_Dpolar	psub(_Dpolar a, _Dpolar b) {
	return c2p(csub(p2c(a), p2c(b)));
}


_Dcomplex _Cbuild(double _Re, double _Im) {
    _Dcomplex r = { _Re, _Im };
    return r;
}

_Dcomplex _Cmulcc(_Dcomplex _X, _Dcomplex _Y) {
    return _Cbuild(_X.real * _Y.real - _X.imag * _Y.imag, _X.real * _Y.imag + _X.imag * _Y.real);
}

_Dcomplex _Cmulcr(_Dcomplex _X, double _Y) {
    return _Cbuild(_X.real * _Y, _X.imag * _Y);
}

_Fcomplex _FCbuild(float _Re, float _Im) {
    _Fcomplex r = { _Re, _Im };
    return r;
}

_Fcomplex _FCmulcc(_Fcomplex _X, _Fcomplex _Y) {
    return _FCbuild(_X.real * _Y.real - _X.imag * _Y.imag, _X.real * _Y.imag + _X.imag * _Y.real);
}

_Fcomplex _FCmulcr(_Fcomplex _X, float _Y) {
    return _FCbuild(_X.real * _Y, _X.imag * _Y);
}

_Lcomplex _LCbuild(long double _Re, long double _Im) {
    _Lcomplex r = { _Re, _Im };
    return r;
}

_Lcomplex _LCmulcc(_Lcomplex _X, _Lcomplex _Y) {
    return _LCbuild(_X.real * _Y.real - _X.imag * _Y.imag, _X.real * _Y.imag + _X.imag * _Y.real);
}

_Lcomplex _LCmulcr(_Lcomplex _X, long double _Y) {
    return _LCbuild(_X.real * _Y, _X.imag * _Y);
}
