#pragma once
#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <windows.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"

extern const double PI;

namespace Reconstruction {

	__global__ void product(float* a, float* b, float* res, const int nxa, const int nya, const int nxb);
	__global__ void init(float* a, const int nx, const int ny, const float value);
	__global__ void transpose_col(float* a, float* res, const int nx, const int ny);
	__global__ void intg_source(float* source, float* res, const int ni, const int nb, const int ne);
	__global__ void lintg(float* sysmat, float* matfrac, float* matatn, float* res, int nd, int ni, int ne, int nm);
	//__global__ void hesseq(float* source, float* sysmat, float* n, float* nbar, float* matatn, float* lintg, float* suma, float* res, int nd, int ni, int nb, int ne, int nm);

	__global__ void gradq(float* sysmat, float* n, float* nbar, float* grad_n, float* res, const int nd, const int ni, const int nb, const int nm);

	__global__ void hesseq1(float* source, float* matatn, float* lintg, float* res, const int ni, const int nb, const int ne, const int nm);
	__global__ void hesseq2(float* sysmat, float* n, float* nbar, float* suma, float* tmp, float* res, const int nd, const int ni, const int nb, const int nm);

	void calc_product(float* a, float* b, float* res, int nxa, int nya, int nxb);
	float** calc_product(float** a, float** b, int nxa, int nya, int nxb);

	void init_matrix(float* a, int nx, int ny, float value);
	void init_matrix(float** a, int nx, int ny, float value);

	float** calc_T(float** a, int nx, int ny);

	float** calc_n(float*** source, int ni, int nb, int ne, int v_begin);

	void calc_lintg(float** res, float** sysmat, float** matfrac, float** matatn, int nd, int ni, int ne, int nm, float pixsize);

	void calc_gradq(float** res, float** sysmat, float* n, float** nbar, float*** grad_n, int nd, int ni, int nb, int nm, int v_begin, float pixsize);

	void calc_hesseq(float*** res, float*** source, float** sysmat, float* n, float** nbar, float** matatn, float** lintg, float* suma, int nd, int ni, int nb, int ne, int nm, int v_begin, float pixsize);

		//h[j][k][m] += smr[i][j] * (n[i][b] / nbar[i][b]) * suma[i] * source[v_begin + i][b][e] * matatn[k][e] * matatn[m][e] * lintg[i][e];

}