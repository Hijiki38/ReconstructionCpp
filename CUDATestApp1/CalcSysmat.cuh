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

	__global__ void calc_coeff_cbct(float* result, const int nd, const int center,
		const int w, const float theta, const float sdd, const int rotcount);

	__global__ void calc_coeff(float* result, const int nd, const int center,
		const int w, const float theta, const float sdd, const int rotcount);

<<<<<<< HEAD
<<<<<<< Updated upstream
	int calc_sysmat2(float* elem, int* rowptr, int* colind, const int nv,
		const int nd, const int center, const int sdd);
=======
=======
>>>>>>> e4107061f770b4d61d36e1cdab14719abc7504d0
	__global__ void calc_coeff_test(float* result, const int nd, const int center,
		const int w, const float theta, const float sdd, const int rotcount);

	void devicereset();

	int calc_sysmat(float* elem, int* rowptr, int* colind, const int vn, const int nd, 
		const int center, const float sdd, const bool write_sysmat = false);

	int calc_sysmat2(float* elem, int* rowptr, int* colind, const int v_start, const int v_size, 
		const int nv, const int nd, const int center, const float sdd, const bool init = false, const bool write_sysmat = false);
<<<<<<< HEAD
>>>>>>> Stashed changes
=======
>>>>>>> e4107061f770b4d61d36e1cdab14719abc7504d0

	double cpuSecond();
}