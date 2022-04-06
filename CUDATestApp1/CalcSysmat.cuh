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
		const int w, const int theta, const int sdd, const int rotcount);

	void calc_sysmat(float* result, const int nd, const int center,
		const int w, const int theta, const int sdd, const int rotcount);

	int calc_sysmat2(float* elem, int* rowptr, int* colind, const int nv,
		const int nd, const int center, const int sdd);

	double cpuSecond();
}