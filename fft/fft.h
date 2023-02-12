// 
// fft.h
// 
// pheno
// Feb 04, 2023
// Feb 11, 2023
// 

#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "complex.h"
#include "numeric.h"

typedef enum {
	backward,
	ortho,
	forward
} norm_mode;

static double _get_forward_norm(unsigned N, norm_mode mode);

static double _get_backward_norm(unsigned N, norm_mode mode);

static complex* _execute_dft(complex* data, unsigned int N, bool is_forward);

static complex* _raw_dft(complex* data, unsigned int N, bool is_forward, complex fct);

/// <summary>
/// Calculate Decrete Fourier Transform
/// <para>.</para>
/// <para> N must be same as data size </para>
/// <para> mode have { backward / ortho / forward } </para>
/// </summary>
/// <param name="data  "> complex* </param>
/// <param name="N       "> unsigned int </param>
/// <param name="mode"> norm_mode </param>
/// <returns> 
/// <para> mode = backward, 1 * Sum { x[n] * pow(e, -j2PI nk/N) } </para>
/// <para> mode = ortho, (1 / sqrt(N)) * Sum { x[n] * pow(e, -j2PI nk/N) } </para>
/// <para> mode = forward, (1 / N) * Sum { x[n] * pow(e, -j2PI nk/N) } </para>
/// </returns>
complex* dft(complex* data, unsigned int N, norm_mode mode);

/// <summary>
/// Calculate Inverse Decrete Fourier Transform
/// <para>.</para>
/// <para> N must be same as data size </para>
/// <para> mode have { backward / ortho / forward } </para>
/// </summary>
/// <param name="data  "> complex* </param>
/// <param name="N       "> unsigned int </param>
/// <param name="mode"> norm_mode </param>
/// <returns> 
/// <para> mode = backward, (1 / N) * Sum { x[n] * pow(e, j2PI nk/N) } </para>
/// <para> mode = ortho, (1 / sqrt(N)) * Sum { x[n] * pow(e, j2PI nk/N) } </para>
/// <para> mode = forward, 1 * Sum { x[n] * pow(e, j2PI nk/N) } </para>
/// </returns>
complex* idft(complex* data, unsigned int N, norm_mode mode);

static complex* _raw_fft(complex* data, unsigned int N, bool is_forward, complex fct);

/// <summary>
/// Calculate Fast Fourier Transform
/// <para>.</para>
/// <para> N must be same as data size </para>
/// <para> mode have { backward / ortho / forward } </para>
/// </summary>
/// <param name="data  "> complex* </param>
/// <param name="N       "> unsigned int </param>
/// <param name="mode"> norm_mode </param>
/// <returns> 
/// <para> mode = backward, 1 * Sum { x[n] * pow(e, -j2PI nk/N) } </para>
/// <para> mode = ortho, (1 / sqrt(N)) * Sum { x[n] * pow(e, -j2PI nk/N) } </para>
/// <para> mode = forward, (1 / N) * Sum { x[n] * pow(e, -j2PI nk/N) } </para>
/// </returns>
complex* fft(complex* data, unsigned int N, norm_mode mode);

/// <summary>
/// Calculate Inverse Fast Fourier Transform
/// <para>.</para>
/// <para> N must be same as data size </para>
/// <para> mode have { backward / ortho / forward } </para>
/// </summary>
/// <param name="data  "> complex* </param>
/// <param name="N       "> unsigned int </param>
/// <param name="mode"> norm_mode </param>
/// <returns> 
/// <para> mode = backward, (1 / N) * Sum { x[n] * pow(e, j2PI nk/N) } </para>
/// <para> mode = ortho, (1 / sqrt(N)) * Sum { x[n] * pow(e, j2PI nk/N) } </para>
/// <para> mode = forward, 1 * Sum { x[n] * pow(e, j2PI nk/N) } </para>
/// </returns>
complex* ifft(complex* data, unsigned int N, norm_mode mode);

static complex* _raw_radix8fft(complex* data, bool is_forward);

static complex* radix8fft(complex* data);

static complex* radix8ifft(complex* data);

static complex* _raw_radix2fft(complex* data, unsigned N, bool is_forward);

static complex* radix2fft(complex* data, unsigned N);

static complex* radix2ifft(complex* data, unsigned N);

static complex* _raw_radarfft(complex* data, unsigned N, bool is_forward);

static complex* radarfft(complex* data, unsigned N);

static complex* radarifft(complex* data, unsigned N);

static complex* _raw_hybrid_radixfft(complex* data, unsigned N, bool is_forward);

static complex* hybrid_radixfft(complex* data, unsigned N);

static complex* hybrid_radixifft(complex* data, unsigned N);
