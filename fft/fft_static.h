// 
// fft_static.h
// 
//		Static function declaration for fft.c
// 

#pragma once

#include "fft.h"

static double _get_forward_norm(unsigned N, norm_mode mode);

static double _get_backward_norm(unsigned N, norm_mode mode);

static _Dcomplex* _execute_dft(_Dcomplex* data, unsigned int N, bool is_forward);

static _Dcomplex* _raw_dft(_Dcomplex* data, unsigned int N, bool is_forward, double fct);

static _Dcomplex* _execute_fft(_Dcomplex* data, unsigned int N, bool is_forward);

static _Dcomplex* _raw_fft(_Dcomplex* data, unsigned int N, bool is_forward, double fct);

static _Dcomplex* _raw_radix2fft(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input);

static _Dcomplex* _raw_radarfft(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input);

static _Dcomplex* _raw_bluesteinfft(_Dcomplex* data, unsigned N, bool is_forward, bool keep_input);
