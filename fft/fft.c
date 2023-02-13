// 
// fft.c
// 
// pheno
// Feb 04, 2023
// Feb 11, 2023
// 

#include "fft.h"

static double _get_forward_norm(unsigned N, norm_mode mode) {
	switch (mode) {
	case backward:
		return 1.;
	case ortho:
		return sqrt(N);
	case forward:
		return 1. * N;
	default:
		return 1.;
	}
}

static double _get_backward_norm(unsigned N, norm_mode mode) {
	switch (mode) {
	case backward:
		return 1. * N;
	case ortho:
		return sqrt(N);
	case forward:
		return 1.;
	default:
		return 1. * N;
	}
}

static complex* _execute_dft(complex* data, unsigned int N, bool is_forward) {

	complex* r = 0;
	if ((r = (complex*)calloc(N, sizeof(complex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	int b = (is_forward) ? 1 : -1;
	for (int k = 0; k < (signed)N; k++) {
		for (int n = 0; n < (signed)N; n++) {
			complex e = { 0., -2. * PI * k * n * b / N };
			*(r + k) = cadd(*(r + k),
				cmul(*(data + n), cexp(e))
			);
		}
	}

	return r;
}

static complex* _raw_dft(complex* data, unsigned int N, bool is_forward, complex fct) {
	complex* result = _execute_dft(data, N, is_forward);
	for (unsigned i = 0; i < N; i++) {
		*(result + i) = cmul(*(result + i), fct);
	}
	return result;
}

complex* dft(complex* data, unsigned int N, norm_mode mode) {
	complex fct = { 1. / _get_forward_norm(N, mode), 0. };
	return _raw_dft(data, N, true, fct);
}

complex* idft(complex* data, unsigned int N, norm_mode mode) {
	complex fct = { 1. / _get_backward_norm(N, mode), 0. };
	return _raw_dft(data, N, false, fct);
}

static complex* _execute_fft(complex* data, unsigned int N, bool is_forward) {
	if (N == 1) return data;
	if (N == 8) return _raw_radix8fft(data, is_forward, true);
	if (N == 2) return _raw_radix2fft(data, N, is_forward, true);
	if (isprime(N)) return _raw_radarfft(data, N, is_forward, true);
	if (isPowerOf2(N)) return _raw_radix2fft(data, N, is_forward, true);
	else return _raw_hybrid_radixfft(data, N, is_forward, true);
}

static complex* _raw_fft(complex* data, unsigned int N, bool is_forward, complex fct) {
	complex* result = _execute_fft(data, N, is_forward);
	for (unsigned i = 0; i < N; i++) {
		*(result + i) = cmul(*(result + i), fct);
	}
	return result;
}

complex* fft(complex* data, unsigned int N, norm_mode mode) {
	complex fct = { 1. / _get_forward_norm(N, mode), 0. };
	return _raw_fft(data, N, true, fct);
}

complex* ifft(complex* data, unsigned int N, norm_mode mode) {
	complex fct = { 1. / _get_backward_norm(N, mode), 0. };
	return _raw_fft(data, N, false, fct);
}

static complex* _raw_radix8fft(complex* data, bool is_forward, bool keep_input) {

	int k = (is_forward) ? 1 : -1;

	complex wn[4] = {
		{ 1., 0. },											// exp(-i0PI/8)
		{ cos(2. * PI / 8.), -sin(2. * PI / 8.) * k },		// exp(-i2PI/8)
		{ cos(4. * PI / 8.), -sin(4. * PI / 8.) * k },		// exp(-i4PI/8)
		{ cos(6. * PI / 8.), -sin(6. * PI / 8.) * k }		// exp(-i6PI/8)
	};

	complex tmp1[8] = { 0 };
	complex tmp2[8] = { 0 };

	// Decimation In Time - Fast Fourier Transform

	tmp1[0] = cadd(data[0], data[4]);
	tmp1[1] = csub(data[0], data[4]);
	tmp1[2] = cadd(data[2], data[6]);
	tmp1[3] = cmul(csub(data[2], data[6]), wn[2]);
	tmp1[4] = cadd(data[1], data[5]);
	tmp1[5] = csub(data[1], data[5]);
	tmp1[6] = cadd(data[3], data[7]);
	tmp1[7] = cmul(csub(data[3], data[7]), wn[2]);
	if (!keep_input) free(data);

	tmp2[0] = cadd(tmp1[0], tmp1[2]);
	tmp2[1] = cadd(tmp1[1], tmp1[3]);
	tmp2[2] = csub(tmp1[0], tmp1[2]);
	tmp2[3] = csub(tmp1[1], tmp1[3]);
	tmp2[4] = cadd(tmp1[4], tmp1[6]);
	tmp2[5] = cmul(cadd(tmp1[5], tmp1[7]), wn[1]);
	tmp2[6] = cmul(csub(tmp1[4], tmp1[6]), wn[2]);
	tmp2[7] = cmul(csub(tmp1[5], tmp1[7]), wn[3]);

	complex* result = 0;
	if ((result = (complex*)calloc(8, sizeof(complex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	*(result + 0) = cadd(tmp2[0], tmp2[4]);
	*(result + 1) = cadd(tmp2[1], tmp2[5]);
	*(result + 2) = cadd(tmp2[2], tmp2[6]);
	*(result + 3) = cadd(tmp2[3], tmp2[7]);
	*(result + 4) = csub(tmp2[0], tmp2[4]);
	*(result + 5) = csub(tmp2[1], tmp2[5]);
	*(result + 6) = csub(tmp2[2], tmp2[6]);
	*(result + 7) = csub(tmp2[3], tmp2[7]);

	return result;
}

static complex* radix8fft(complex* data) {
	return _raw_radix8fft(data, true, true);
}

static complex* radix8ifft(complex* data) {
	return _raw_radix8fft(data, false, true);
}

static complex* _raw_radix2fft(complex* data, unsigned N, bool is_forward, bool keep_input) {

	if (N == 8) return _raw_radix8fft(data, is_forward, false);
	if (N == 1) return data;

	complex* Peven = 0;
	if ((Peven = (complex*)calloc((N / 2), sizeof(complex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	complex* Podd = 0;
	if ((Podd = (complex*)calloc((N / 2), sizeof(complex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	for (int i = 0; i < (signed)N / 2; i++) {
		*(Peven + i) = *(data + i * 2);
		*(Podd + i) = *(data + i * 2 + 1);
	}
	if (!keep_input) free(data);
	
	Peven = _raw_radix2fft(Peven, N / 2, is_forward, false);
	Podd = _raw_radix2fft(Podd, N / 2, is_forward, false);

	complex* result = 0;
	if ((result = (complex*)calloc(N, sizeof(complex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	complex w = { 0 };
	int k = (is_forward) ? 1 : -1;
	for (int i = 0; i < (signed)N / 2; i++) {
		complex e = { 0., -2. * PI * i * k / N };
		w = cexp(e);
		*(result + i) = cadd(*(Peven + i), cmul(w, *(Podd + i)));
		*(result + i + N / 2) = csub(*(Peven + i), cmul(w, *(Podd + i)));
	}

	free(Peven);
	free(Podd);

	return result;
}

static complex* radix2fft(complex* data, unsigned N) {
	return _raw_radix2fft(data, N, true, true);
}

static complex* radix2ifft(complex* data, unsigned N) {
	return _raw_radix2fft(data, N, false, true);
}

static complex* _raw_radarfft(complex* data, unsigned N, bool is_forward, bool keep_input) {

	if (N == 1) return data;

	unsigned M = nextPowerOf2(2 * (N - 1) - 1);
	unsigned g = findPrimitiveRoot(N);

	// a[q] = x_g^q = x[ pow(g, q) (mod N) ]
	complex* a = 0;
	if ((a = (complex*)calloc((size_t)N - 1, sizeof(complex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	// b[q] = x_g^-q = x[ pow(g, -q) (mod N) ]
	complex* b = 0;
	if ((b = (complex*)calloc((size_t)N - 1, sizeof(complex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	unsigned* g_mod_minus_q = 0;
	if ((g_mod_minus_q = (unsigned*)calloc((size_t)N - 1, sizeof(unsigned))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}

	int k = (is_forward) ? 1 : -1;
	*g_mod_minus_q = expModuloNInverse(g, N);
	*a = *(data + expModuloN(g, 0, N));
	complex e = { 0., -2 * PI * k / N };
	*b = cexp(e);
	for (unsigned q = 1; q < N - 1; q++) {
		*(a + q) = *(data + expModuloN(g, q, N));
		complex e = { 0., -2 * PI * *(g_mod_minus_q + q - 1) * k / N };
		*(b + q) = cexp(e);
		// pow(g, -2) (mod N)
		//      = { [pow(g, -1) (mod N)] * [pow(g, -1) (mod N)] } % N
		//		= [ pow(g, -1) * pow(g, -1) ] (mod N)
		*(g_mod_minus_q + q) = (*(g_mod_minus_q + q - 1) * *g_mod_minus_q) % N;
	}

	// Padding zero between the 0th element and the 1st element 
	// until the length of a is extended to M.
	complex* new_a = 0;
	complex zero = { 0., 0. };
	if ((new_a = (complex*)realloc(a, M * sizeof(complex))) == NULL) {
		printf("Allocation failed!\n");
		free(a);
		return 0;
	}
	a = NULL;
	// Repeat the array until the length increases to M
	complex* new_b = 0;
	if ((new_b = (complex*)realloc(b, M * sizeof(complex))) == NULL) {
		printf("Allocation failed!\n");
		free(b);
		return 0;
	}
	b = NULL;

	for (unsigned i = 1; i < N - 1; i++) {
		*(new_a + i + M - N + 1) = zero;
		cswap(new_a + i, new_a + i + M - N + 1);
	}
	for (unsigned i = N - 1; i < M; i++) {
		if (i < (M - N + 1) + 1) *(new_a + i) = zero;
		*(new_b + i) = *(new_b + i - N + 1);
	}

	complex* fft_a = (M == 8) ? radix8fft(new_a) : radix2fft(new_a, M);
	free(new_a);
	complex* fft_b = (M == 8) ? radix8fft(new_b) : radix2fft(new_b, M);
	free(new_b);

	complex* fft_ab = 0;
	if ((fft_ab = (complex*)calloc(M, sizeof(complex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	for (unsigned i = 0; i < M; i++) {
		*(fft_ab + i) = cmul(*(fft_a + i), *(fft_b + i));
	}
	free(fft_a);
	free(fft_b);

	complex* conv_ab = (M == 8) ? radix8ifft(fft_ab) : radix2ifft(fft_ab, M);
	free(fft_ab);

	complex* result = 0;
	if ((result = (complex*)calloc(N, sizeof(complex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	for (unsigned n = 0; n < N; n++) {
		*result = cadd(*result, *(data + n));
	}
	complex cM = { (double)M, 0. };
	*(result + 1) = cadd(*data, cdiv(*conv_ab, cM));
	for (unsigned q = 1; q < N - 1; q++) {
		*(result + *(g_mod_minus_q + q - 1)) = cadd(*data, cdiv(*(conv_ab + q), cM));
	}
	if (!keep_input) free(data);
	free(g_mod_minus_q);
	free(conv_ab);
	return result;
}

static complex* radarfft(complex* data, unsigned N) {
	return _raw_radarfft(data, N, true, true);
}

static complex* radarifft(complex* data, unsigned N) {
	return _raw_radarfft(data, N, false, true);
}

static complex* _raw_hybrid_radixfft(complex* data, unsigned N, bool is_forward, bool keep_input) {

	if (N == 1) return data;
	if (N == 8) return _raw_radix8fft(data, is_forward, false);
	if (N == 2) return _raw_radix2fft(data, N, is_forward, false);
	if (isprime(N)) return _raw_radarfft(data, N, is_forward, false);

	unsigned* factors = factor(N);
	unsigned N1 = *(factors + 1), N2 = N / N1;
	free(factors);

	complex** X = 0;
	if ((X = (complex**)calloc(N1, sizeof(complex*))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	for (unsigned i = 0; i < N1; i++) {
		if ((*(X + i) = (complex*)calloc(N2, sizeof(complex))) == NULL) {
			printf("Allocation failed!\n");
			return 0;
		}
	}

	for (unsigned n1 = 0; n1 < N1; n1++) {
		for (unsigned n2 = 0; n2 < N2; n2++) {
			X[n1][n2] = data[N1 * n2 + n1];
		}
	}
	if (!keep_input) free(data);

	int k = (is_forward) ? 1 : -1;
	for (unsigned n1 = 0; n1 < N1; n1++) {
		complex* inner = 0;
		if ((inner = (complex*)calloc(N2, sizeof(complex))) == NULL) {
			printf("Allocation failed!\n");
			return 0;
		}
		for (unsigned i = 0; i < N2; i++) {
			*(inner + i) = X[n1][i];
		}
		complex* fft_inner = _raw_hybrid_radixfft(inner, N2, is_forward, false);
		for (unsigned n2 = 0; n2 < N2; n2++) {
			complex e = { 0., -2 * PI * n1 * n2 * k / N };
			X[n1][n2] = cmul(fft_inner[n2], cexp(e));
		}
		free(fft_inner);
	}

	for (unsigned n2 = 0; n2 < N2; n2++) {
		complex* outer = 0;
		if ((outer = (complex*)calloc(N1, sizeof(complex))) == NULL) {
			printf("Allocation failed!\n");
			return 0;
		}
		for (unsigned i = 0; i < N1; i++) {
			*(outer + i) = X[i][n2];
		}
		complex* fft_outer = _raw_hybrid_radixfft(outer, N1, is_forward, false);
		for (unsigned n1 = 0; n1 < N1; n1++) {
			X[n1][n2] = fft_outer[n1];
		}
		free(fft_outer);
	}

	complex* result = 0;
	if ((result = (complex*)calloc(N, sizeof(complex))) == NULL) {
		printf("Allocation failed!\n");
		return 0;
	}
	for (unsigned n1 = 0; n1 < N1; n1++) {
		for (unsigned n2 = 0; n2 < N2; n2++) {
			*(result + N2 * n1 + n2) = X[n1][n2];
		}
	}

	for (unsigned i = 0; i < N1; i++) {
		free(*(X + i));
	}
	free(X);

	return result;
}

static complex* hybrid_radixfft(complex* data, unsigned N) {
	return _raw_hybrid_radixfft(data, N, true, true);
}

static complex* hybrid_radixifft(complex* data, unsigned N) {
	return _raw_hybrid_radixfft(data, N, false, true);
}
