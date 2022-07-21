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

	__global__ void calc_coeff_test(float* result, const int nd, const int center,
		const int w, const float theta, const float sdd, const int rotcount) {

		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nd + ix;

		if (ix < nd && iy < nd) {
			result[idx] = 0;
		}

	}


	__global__ void calc_coeff(float* result, const int nd, const int center,
		const int w, const float theta, const float sdd, const int rotcount) {

		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nd + ix;

		float relx = (float)ix;
		float rely = (float)iy;
		float tmp_x;
		float offset_detector, intercept;

		float la1, la2, lb1, lb2, sa, sb;
		float a = 0.5; // = pixsize / 2
		float tan_angle = tanf(theta);
		float cos_angle = cosf(theta);


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
			intercept = offset_detector - rely + relx * tan_angle;

			la1 = a - (-a * tan_angle + intercept + a / cos_angle);
			la2 = a - (a * tan_angle + intercept + a / cos_angle);
			lb1 = a + (-a * tan_angle + intercept - a / cos_angle);
			lb2 = a + (a * tan_angle + intercept - a / cos_angle);

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
				else {
					sa = a * (la1 + la2) - (la2 - 2 * a) * (la2 - 2 * a) / (2 * (la2 - la1));
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
				else {
					sb = a * (lb1 + lb2) - (lb2 - 2 * a) * (lb2 - 2 * a) / (2 * (lb2 - lb1));
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


	__global__ void calc_coeff_cbct(float* result, const int nd, const int center,
		const int w, const float theta, const float sdd, const int rotcount) {

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
		float tan_delta = a / sdd;

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
			phi = atan2f(offset_detector, sdd);

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
				else {
					sa = a * (la1 + la2) - (la2 - 2 * a) * (la2 - 2 * a) / (2 * (la2 - la1));
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
					sb = a * lb1 * lb1 / (lb1 - lb2);
				}
				else if (lb2 < 2 * a) {
					sb = a * (lb1 + lb2);
				}
				else {
					sb = a * (lb1 + lb2) - (lb2 - 2 * a) * (lb2 - 2 * a) / (2 * (lb2 - lb1));
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

	__global__ void calc_coeff_cbct_l(float* result, const int nd, const int center,
		const float theta, const float phi, const float sod_norm, const float sdd_norm, const int rotcount) {

		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nd + ix;

		float relx = ix;
		float rely = iy;
		float tmp_x;
		float offset_detector;
		float intercept;

		float l0, l1, k;
		float a = 0.5; //pixsize / 2
		float tan_angle;
		float tan_delta = a / sdd_norm;

		float pi = 3.141592;

		if (ix < nd && iy < nd) {
			for (int i = 0; i < rotcount; i++) {
				tmp_x = relx;
				relx = 2 * center - rely - 1;
				rely = tmp_x;
			}
			//float point::get_relative(float _num) { return _num - center + 0.5; }

			relx = relx - center + 0.5;
			rely = center - rely - 0.5;

			tan_angle = tanf(theta + phi);
			//intercept = offset_detector - rely + relx * tan_angle; //���Y�s�N�Z����ł̂��ؕ�
			intercept = sod_norm * tanf(phi) / (cosf(phi) * cosf(theta + phi)) - rely + relx * tan_angle; //���Y�s�N�Z����ł̂��ؕ�

			k = 0;

			if (theta + phi < pi / 2 && theta + phi >= 0) {  //0~45
				l0 = (-a - intercept) / tan_angle;
				l1 = (a - intercept) / tan_angle;
				if (l1 >= a && l0 < -a) {
					k = 1;
				}
				else if (l1 < a && l0 < -a) {
					k = 0.5 + l1 / (2 * a); //might be negative
				}
				else {
					k = 0.5 - l0 / (2 * a); //might be negative
				}
			}
			else if (theta + phi > -pi / 2 && theta + phi < 0) { //-45~0
				l0 = (a - intercept) / tan_angle;
				l1 = (-a - intercept) / tan_angle;
				if (l1 >= a && l0 < -a) {
					k = 1;
				}
				else if (l1 < a && l0 < -a) {
					k = 0.5 + l1 / (2 * a); //might be negative
				}
				else {
					k = 0.5 - l0 / (2 * a); //might be negative
				}
			}
			//else if (theta + phi > pi / 2) { //45~
			//	l0 = -a * tan_angle + intercept;
			//	l1 = a * tan_angle + intercept;
			//	if (l1 >= a && l0 < -a) {
			//		k = 1;
			//	}
			//	else if (l1 < a && l0 < -a) {
			//		k = 0.5 + l1 / (2 * a); //might be negative
			//	}
			//	else {
			//		k = 0.5 - l0 / (2 * a); //might be negative
			//	}
			//}
			//else { //~-45
			//	l0 = -a * tan_angle + intercept;
			//	l1 = a * tan_angle + intercept;
			//	if (l0 >= a && l1 < -a) {
			//		k = 1;
			//	}
			//	else if (l0 < a && l1 < -a) {
			//		k = 0.5 + l1 / (2 * a); //might be negative
			//	}
			//	else {
			//		k = 0.5 - l0 / (2 * a); //might be negative
			//	}
			//}

			result[idx] = k / cosf(theta + phi);
		}
	}

	void deviceinit() {

	}

	void devicereset() {
		CHECK(cudaDeviceReset());
	}

	int calc_sysmat(float* elem, int* rowptr, int* colind, const int nv, const int nd, const int center, const float sdd, const bool write_sysmat) {

		float area = 0;
		float theta = 0;
		int rotatecount = 0;
		bool firstelem = true;

		int nonzero = 0;
		float* tmpmat = (float*)malloc(sizeof(float) * nd * nd);


		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));

		int nxy = nd * nd;
		int nBytes = nxy * sizeof(float);

		float* d_res;
		printf("\ncudamalloc");
		CHECK(cudaMalloc((void**)&d_res, nBytes));
		printf("\ncompleted");
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
			{ //���e�p��45�x�𒴂�����摜��90�x�E��]�������e�p�� - 45�x��
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

					}
				}

			}

			theta += 2 * PI / nv;
		}

		CHECK(cudaFree(d_res));
		CHECK(cudaDeviceReset());

		//printf("check elem[100]: %f\n", elem[100]);

		return nonzero;
	}

	int calc_sysmat_l(float* elem, int* rowptr, int* colind, const int nv, const int nd, const int center, const float sod, const float sdd, const bool write_sysmat) {

		float area = 0;
		float theta = 0;
		float thetatmp = 0;
		float phi = 0;
		int rotatecount = 0;
		int additionalrot = 0;
		bool firstelem = true;

		float offset_detector;

		int nonzero = 0;
		float* tmpmat = (float*)malloc(sizeof(float) * nd * nd);


		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));

		int nxy = nd * nd;
		int nBytes = nxy * sizeof(float);

		float* d_res;
		printf("\ncudamalloc");
		CHECK(cudaMalloc((void**)&d_res, nBytes));
		printf("\ncompleted");
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
			{ //���e�p��45�x�𒴂�����摜��90�x�E��]�������e�p�� - 45�x��
				rotatecount++;
				theta -= PI / 2;
			}


			for (int w = 0; w < nd; w++)
			{

				offset_detector = (nd - w - 1 - center + 0.5); //offset of the detector
				phi = atan2f(offset_detector, sdd);

				additionalrot = 0;
				thetatmp = theta;
				if (theta + phi >= PI / 4) {
					additionalrot = 1;
					thetatmp -= PI / 2;
				}
				else if (theta + phi < -PI / 4) {
					additionalrot = -1;
					thetatmp += PI / 2;
				}

				double iStart = cpuSecond();
				calc_coeff_cbct_l << < grid, block >> > (d_res, nd, center, thetatmp, phi, sod, sdd, rotatecount + additionalrot);
				CHECK(cudaDeviceSynchronize());
				double iElaps = cpuSecond() - iStart;


				//	__global__ void calc_coeff_cbct_l(float* result, const int nd, const int center,
				//const float theta, const float phi, const float sod_norm, const float sdd_norm, const int rotcount);

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
						if (area > 0) {
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

					}
				}

			}

			theta += 2 * PI / nv;
		}

		CHECK(cudaFree(d_res));
		CHECK(cudaDeviceReset());

		//printf("check elem[100]: %f\n", elem[100]);

		return nonzero;
	}

	int calc_sysmat2(float* elem, int* rowptr, int* colind, const int v_start, const int v_size, const int nv, const int nd, const int center, const float sdd, const bool init, const bool write_sysmat) {

		float area = 0;
		float theta = 0;
		int rotatecount = 0;
		bool firstelem = true;

		int nonzero = 0;
		float* tmpmat = (float*)malloc(sizeof(float) * nd * nd);



		if (write_sysmat) {
			FILE* fp;
			fp = fopen("C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sysmatgpu.csv", "w");
		}

		if (init) {
			int dev = 0;
			cudaDeviceProp deviceprop;
			CHECK(cudaGetDeviceProperties(&deviceprop, dev));
			CHECK(cudaSetDevice(dev));
		}


		int nxy = nd * nd;
		int nBytes = nxy * sizeof(float);

		float* d_res;
		CHECK(cudaMalloc((void**)&d_res, nBytes));

		int dimx = 32;
		int dimy = 32;
		dim3 block(dimx, dimy);
		dim3 grid((nd + block.x - 1) / block.x, (nd + block.y - 1) / block.y);

		//char str[100000];
		//char buf[24];

		//str[0] = '\0';



		for (int i = 0; i < v_start; i++) {
			if (theta >= PI / 4)
			{ //���e�p��45�x�𒴂�����摜��90�x�E��]�������e�p�� - 45�x��
				rotatecount++;
				theta -= PI / 2;
			}
			theta += 2 * PI / nv;
		}

		for (int v = v_start; v < v_start + v_size; v++)
		{

			printf("\r%d / %d to %d  (theta: %f)", v, v_start, v_start + v_size, theta);

			if (theta >= PI / 4)
			{ //���e�p��45�x�𒴂�����摜��90�x�E��]�������e�p�� - 45�x��
				rotatecount++;
				theta -= PI / 2;
			}

			for (int w = 0; w < nd; w++)
			{

				//double iStart = cpuSecond();
				//calc_coeff << < grid, block >> > (d_res, nd, center, w, theta, sdd, rotatecount);
				//printf("\ncalling kernel");
				calc_coeff << < grid, block >> > (d_res, nd, center, w, theta, sdd, rotatecount);
				//printf("... finished\n");
				CHECK(cudaDeviceSynchronize());
				//printf("... finished\n");
				//double iElaps = cpuSecond() - iStart;

				CHECK(cudaGetLastError());

				//printf("\n cudaMemcpy: %d", w);
				CHECK(cudaMemcpy(tmpmat, d_res, nBytes, cudaMemcpyDeviceToHost));
				//printf("... finished\n");
				
				//printf("Elapsed: %lf [s] \n", iElaps);

				firstelem = true;

				//printf("hoge");
				//if (write_sysmat) {
				//	str[0] = '\0';
				//}

				//printf("\n convert to spmat");

				for (int y = 0; y < nd; y++)
				{
					for (int x = 0; x < nd; x++)
					{
						area = tmpmat[y * nd + x];
						//if (area != 0 && (x - center) * (x - center) + (y - center) * (y - center) < center * center) {
						if (area != 0) {
							elem[nonzero] = area;
							colind[nonzero] = nd * y + x;
							if (firstelem) {
								rowptr[nd * (v - v_start) + w] = nonzero;
								firstelem = false;
							}
							nonzero++;
						}
						//if (write_sysmat) {
						//	snprintf(buf, 24, "%f,", area);
						//	strcat(str, buf);
						//}
					}
				}

				//printf("... finished\n");

				//if (write_sysmat) {
				//	strcat(str, "\n");
				//	fprintf(fp, str);
				//}


			}

			theta += 2 * PI / nv;
		}

		printf("check elem[100]: %f\n", elem[100]);

		//printf("\ncudafree\n");
		CHECK(cudaFree(d_res));
		//printf("... finished\n");
		//printf("\ncudadevicereset\n");
		CHECK(cudaDeviceReset());
		//printf("... finished\n");

		//if (write_sysmat) {
		//	fclose(fp);
		//}

		// for debug
		//exit(0);



		//printf("\nfree tmpmat\n");
		free(tmpmat);
		//printf("... finished\n");

		////for debug
		//elem[0] = 0;
		//rowptr[0] = 0;
		//colind[0] = 0;
		//return 0;
		////


		//printf("\n %d !!!!!!!!!!!!!!!!!!!!!!!! \n", nonzero);

		return nonzero;
	}

	double cpuSecond() {

		SYSTEMTIME st;
		GetLocalTime(&st);

		return ((double)st.wSecond + (double)st.wMilliseconds * 1.e-3);

	}
}
