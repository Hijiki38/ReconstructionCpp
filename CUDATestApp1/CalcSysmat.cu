#include "CalcSysmat.cuh"


#define CHECK(call)\
{\
    const cudaError_t error = call;\
    if (error != cudaSuccess)\
        {\
            printf("Error: %s:%d",__FILE__,__LINE__);\
            printf("code:%d, reason: %s\n", error, cudaGetErrorString(error));\
            exit(1);\
        }\
}

namespace Reconstruction {

	__global__ void calc_coeff_cbct(float* result, const int nd, const int center,
		const int w, const int theta, const int sdd, const int rotcount) {

		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nd + ix;

		float relx = ix;
		float rely = iy;
		float tmp_x;
		float offset_detector, phi, intercept;

		float la1, la2, lb1, lb2, sa, sb;
		float a = 0.5; // = pixsize / 2
		float tan_angle;
		//float tan_delta = a / sdd;
		float tan_delta = 0;

		if (ix < nd && iy < nd) {
			for (int i = 0; i < rotcount; i++) {
				tmp_x = relx;
				relx = 2 * center - rely - 1;
				rely = tmp_x;
			}
			//float point::get_relative(float _num) { return _num - center + 0.5; }

			relx = relx - center + 0.5;
			rely = center - rely - 0.5;

			offset_detector = (nd - w - 1 - center + 0.5) / cosf(theta); //offset of the detector
			//phi = atan2f(offset_detector, sdd);
			phi = 0;

			tan_angle = tanf(theta + phi);
			intercept = offset_detector - rely + relx * tan_angle;

			la1 = a - (-a * (tan_angle + tan_delta) + intercept + a / sqrt(1 / (1 + (tan_angle + tan_delta) * (tan_angle + tan_delta))));
			la2 = a - (a * (tan_angle + tan_delta) + intercept + a / sqrt(1 / (1 + (tan_angle + tan_delta) * (tan_angle + tan_delta))));
			lb1 = a + (-a * (tan_angle - tan_delta) + intercept - a / sqrt(1 / (1 + (tan_angle - tan_delta) * (tan_angle - tan_delta))));
			lb2 = a + (a * (tan_angle - tan_delta) + intercept - a / sqrt(1 / (1 + (tan_angle - tan_delta) * (tan_angle - tan_delta))));

			if (la1 < 0) {
				if (la2 < 0) {
					sa = 0;
				}
				else {
					sa = a * la2 * la2 / (-la1 + la2);
				}
			}
			else if (la1 < 2 * a) {
				if (la2 < 0) {
					sa = a * la1 * la1 / (la1 - la2);
				}
				else if (la2 < 2 * a) {
					sa = a * (la1 + la2);
				}
			}
			else {
				if (la2 < 2 * a) {
					sa = a * (la1 + la2) - (la1 - 2 * a) * (la1 - 2 * a) / (2 * (la1 - la2));
				}
				else {
					sa = (2 * a) * (2 * a);
				}
			}

			if (lb1 < 0) {
				if (lb2 < 0) {
					sb = 0;
				}
				else {
					sb = a * lb2 * lb2 / (-lb1 + lb2);
				}
			}
			else if (lb1 < 2 * a) {
				if (lb2 < 0) {
					sb = a * lb1 * lb1 / (lb1 -lb2);
				}
				else if (lb2 < 2 * a) {
					sb = a * (lb1 + lb2);
				}
			}
			else {
				if (lb2 < 2 * a) {
					sb = a * (lb1 + lb2) - (lb1 - 2 * a) * (lb1 - 2 * a) / (2 * (lb1 - lb2));
				}
				else {
					sb = (2 * a) * (2 * a);
				}
			}

			result[idx] = (2 * a) * (2 * a) - (sa + sb);
		}
	}

	void calc_sysmat(float* result, const int nd, const int center,
		const int w, const int theta, const int sdd, const int rotcount) {

		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));

		int nxy = nd * nd;
		int nBytes = nxy * sizeof(float);

		float* d_res;
		CHECK(cudaMalloc((void**)&d_res, nBytes));

		CHECK(cudaMemcpy(d_res, result, nBytes, cudaMemcpyHostToDevice));

		int dimx = 32;
		int dimy = 32;
		dim3 block(dimx, dimy);
		dim3 grid((nd + block.x - 1) / block.x, (nd + block.y - 1) / block.y);

		double iStart = cpuSecond();
		calc_coeff_cbct << < grid, block >> > (d_res, nd, center, w, theta, sdd, rotcount);
		CHECK(cudaDeviceSynchronize());
		double iElaps = cpuSecond() - iStart;

		CHECK(cudaGetLastError());
		CHECK(cudaMemcpy(result, d_res, nBytes, cudaMemcpyDeviceToHost));

		CHECK(cudaFree(d_res));

		//printf("Elapsed: %lf [s] \n", iElaps);

		cudaDeviceReset();

		return;
	}

	int calc_sysmat2(float* elem, int* rowptr, int* colind, const int nv, const int nd, const int center, const int sdd) {

		float area = 0;
		float theta = 0;
		int rotatecount = 0;
		bool firstelem = true;

		int nonzero = 0;
		float* tmpmat = (float*)malloc(sizeof(float) * nd * nd);

		FILE* fp;

		fp = fopen("C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sysmatgpu.csv", "w");


		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));

		int nxy = nd * nd;
		int nBytes = nxy * sizeof(float);

		float* d_res;
		CHECK(cudaMalloc((void**)&d_res, nBytes));

		//CHECK(cudaMemcpy(d_res, tmpmat, nBytes, cudaMemcpyHostToDevice));

		int dimx = 32;
		int dimy = 32;
		dim3 block(dimx, dimy);
		dim3 grid((nd + block.x - 1) / block.x, (nd + block.y - 1) / block.y);

		char str[100000];
		char buf[24];

		str[0] = '\0';

		for (int v = 0; v < nv; v++)
		{

			printf("\r%d / %d", v, nv);

			if (theta >= PI / 4)
			{ //“Š‰eŠp‚ª45“x‚ð’´‚¦‚½‚ç‰æ‘œ‚ð90“x‰E‰ñ“]‚³‚¹“Š‰eŠp‚ð - 45“x‚É
				rotatecount++;
				theta -= PI / 2;
			}

			for (int w = 0; w < nd; w++)
			{

				double iStart = cpuSecond();
				calc_coeff_cbct << < grid, block >> > (d_res, nd, center, w, theta, sdd, rotatecount);
				CHECK(cudaDeviceSynchronize());
				double iElaps = cpuSecond() - iStart;

				CHECK(cudaGetLastError());
				CHECK(cudaMemcpy(tmpmat, d_res, nBytes, cudaMemcpyDeviceToHost));

				//printf("Elapsed: %lf [s] \n", iElaps);

				firstelem = true;

				//printf("hoge");
				str[0] = '\0';

				for (int y = 0; y < nd; y++)
				{
					for (int x = 0; x < nd; x++)
					{
						area = tmpmat[y * nd + x];
						if (area != 0) {
							//if (nonzero == MAXMATERIALS - 1) {
							//	exit(1);
							//}
							elem[nonzero] = area;
							colind[nonzero] = nd * y + x;
							if (firstelem) {
								rowptr[nd * v + w] = nonzero;
								firstelem = false;
							}
							nonzero++;
						}
						snprintf(buf, 24, "%f,", area);
						strcat(str, buf);
						//fprintf(fp, "%f", area);
					}
				}

				strcat(str, "\n");
				fprintf(fp, str);

			}

			theta += 2 * PI / nv;
		}

		CHECK(cudaFree(d_res));
		cudaDeviceReset();

		fclose(fp);

		return nonzero;
	}

	double cpuSecond() {

		SYSTEMTIME st;
		GetLocalTime(&st);

		return ((double)st.wSecond + (double)st.wMilliseconds * 1.e-3);

	}
}
