#pragma once
#include <cstdio>
//#include <cstdlib>
#include <time.h>
#include <windows.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"

namespace Reconstruction {

	__global__ void calc_coeff_cbct(float* result, const int nd, const int center,
		const int w, const int theta, const int sdd, const int rotcount);

	void calc_sysmat(float* result, const int nd, const int center,
		const int w, const int theta, const int sdd, const int rotcount);

	double cpuSecond();
}