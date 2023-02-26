// 
// complex.h
// 

#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Types
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef _C_POLAR_T
    #define _C_POLAR_T
    typedef struct _C_double_polar {
        double radial;
        double angular;
    } _C_double_polar;

    typedef struct _C_float_polar {
        float _Val[2];
    } _C_float_polar;

    typedef struct _C_ldouble_polar {
        long double _Val[2];
    } _C_ldouble_polar;
#endif

typedef _C_double_polar  _Dpolar;
typedef _C_float_polar   _Fpolar;
typedef _C_ldouble_polar _Lpolar;

#ifndef _C_complex_T
    #define _C_complex_T
    typedef struct _C_double_complex {
        double real;
        double imag;
    } _C_double_complex;

    typedef struct _C_float_complex {
        float real;
        float imag;
    } _C_float_complex;

    typedef struct _C_ldouble_complex {
        long double real;
        long double imag;
    } _C_ldouble_complex;
#endif

typedef _C_double_complex  _Dcomplex;
typedef _C_float_complex   _Fcomplex;
typedef _C_ldouble_complex _Lcomplex;

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Macros
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#define _DCOMPLEX_(re, im)  _Cbuild(re, im)
#define _FCOMPLEX_(re, im)  _FCbuild(re, im)
#define _LCOMPLEX_(re, im)  _LCbuild(re, im)

#define _Complex_J _Cbuild(0.0, 1.0)
#define J          _Complex_J

#ifndef PI
#define PI (3.1415926535897932384626433832795028841971693993751058209749445923)
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Functions
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

_Dpolar	    c2p(_Dcomplex z);
double      cabs(_Dcomplex z);
_Dcomplex   cacos(_Dcomplex z);
_Dcomplex   cacosh(_Dcomplex z);
_Dcomplex   cadd(_Dcomplex x, _Dcomplex y);
double      carg(_Dcomplex z);
_Dcomplex   casin(_Dcomplex z);
_Dcomplex   casinh(_Dcomplex z);
_Dcomplex   catan(_Dcomplex z);
_Dcomplex   catanh(_Dcomplex z);
_Dcomplex*  ccopy(_Dcomplex* z, unsigned N);
_Dcomplex   ccos(_Dcomplex z);
_Dcomplex   ccosh(_Dcomplex z);
_Dcomplex   cdiv(_Dcomplex x, _Dcomplex y);
_Dcomplex   cexp(_Dcomplex z);
_Dcomplex   clog(_Dcomplex z);
_Dcomplex   clog10(_Dcomplex z);
_Dcomplex   clogn(_Dcomplex z, double n);
_Dcomplex   cmul(_Dcomplex x, _Dcomplex y);
_Dcomplex   conj(_Dcomplex z);
_Dcomplex   cpow(_Dcomplex x, _Dcomplex y);
_Dcomplex   csin(_Dcomplex z);
_Dcomplex   csinh(_Dcomplex z);
_Dcomplex   csub(_Dcomplex x, _Dcomplex y);
_Dcomplex   csqrt(_Dcomplex z);
void	      cswap(_Dcomplex* a, _Dcomplex* b);
_Dcomplex   ctan(_Dcomplex z);
_Dcomplex   ctanh(_Dcomplex z);
double      norm(_Dcomplex z);
_Dcomplex   p2c(_Dpolar p);
_Dpolar	    padd(_Dpolar a, _Dpolar b);
_Dpolar	    pdiv(_Dpolar a, _Dpolar b);
_Dpolar	    pmul(_Dpolar a, _Dpolar b);
_Dpolar	    psub(_Dpolar a, _Dpolar b);

_Dcomplex _Cbuild(double _Re, double _Im);
_Dcomplex _Cmulcc(_Dcomplex _X, _Dcomplex _Y);
_Dcomplex _Cmulcr(_Dcomplex _X, double _Y);

_Fcomplex _FCbuild(float _Re, float _Im);
_Fcomplex _FCmulcc(_Fcomplex _X, _Fcomplex _Y);
_Fcomplex _FCmulcr(_Fcomplex _X, float _Y);

_Lcomplex _LCbuild(long double _Re, long double _Im);
_Lcomplex _LCmulcc(_Lcomplex _X, _Lcomplex _Y);
_Lcomplex _LCmulcr(_Lcomplex _X, long double _Y);
