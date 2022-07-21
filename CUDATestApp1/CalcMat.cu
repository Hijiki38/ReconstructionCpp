#include "CalcMat.cuh"


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

	__global__ void product(float* a, float* b, float* res, const int nxa, const int nya, const int nxb) {

		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nxb + ix;

		if (ix < nxb && iy < nya) {
			res[idx] = 0;
			for (int i = 0; i < nxa; i++) {
				res[idx] += a[iy * nxa + i] * b[i * nxb + ix];
			}
		}
	}

	__global__ void init(float* a, const int nx, const int ny, const float value) {

		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nx + ix;

		if (ix < nx && iy < ny) {
			a[idx] = value;
		}
	}

	__global__ void transpose_col(float* a, float* res, const int nx, const int ny) {

		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;

		if (ix < nx && iy < ny) {
			res[iy * nx + ix] = a[ix * ny + iy];
		}
	}

	__global__ void intg_source(float* source, float* res, const int ni, const int nb, const int ne){ //n[i][b] += source[v_begin + i][b][e];

		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nb + ix;

		if (ix < nb && iy < ni) {
			res[idx] = 0;
			for (int i = 0; i < ne; i++) {
				res[idx] += source[(iy * nb * ne) + (ix * ne) + i];
			}
		}
	}

	__global__ void gradq(float* sysmat, float* n, float* nbar, float* grad_n, float* res, const int nd, const int ni, const int nb, const int nm) {

		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nm + ix; // nd*nd X nm

		if (ix < nm && iy < nd * nd) {
			res[idx] = 0;
			for (int i = 0; i < ni; i++) {
				for (int b = 0; b < nb; b++) {
					//res[idx] += -sysmat[i * nd * nd + iy] * (1 - (n[i * nb + b] / nbar[i * nb + b])) * grad_n[b * nm + ix];
					// 
					res[idx] += -sysmat[i * nd * nd + iy] * (1 - (n[i * nb + b] / nbar[i * nb + b])) * grad_n[(i * nb * nm) + (b * nm) + ix];
					//res[idx] += (1 - (n[i * nb + b] / nbar[i * nb + b])) * grad_n[(i * nb * nm) + (b * nm) + ix];
				}
			}
		}

	}

	__global__ void gradq2_1(float* source, float* matatn, float* lintg, float* res, const int ni, const int nb, const int ne, const int nm) {
		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nm + ix; // ni*nb X nm

		if (ix < nm && iy < ni * nb) {
			res[idx] = 0;
			for (int e = 0; e < ne; e++) {
				res[idx] += source[(iy / nb) * nb * ne + (iy % nb) * ne + e] * lintg[(iy / nb) * ne + e] * matatn[ix * ne + e];
			}
			//res[idx] /= 10;

		}
	}

	__global__ void gradq2_2(float* sysmat, float* n, float* nbar, float* suma, float* tmp, float* res, const int nd, const int ni, const int nb, const int nm) {
		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nm + ix; // nd*nd X nm

		if (ix < nm && iy < nd * nd) {
			res[idx] = 0;

			for (int i = 0; i < ni; i++) {
				for (int b = 0; b < nb; b++) {
					res[idx] += sysmat[i * nd * nd + iy] * ((n[i * nb + b] / nbar[i * nb + b]) - 1) * tmp[i * nb * nm + b * nm + ix];  // i b m
					//res[idx] += sysmat[i * nd * nd + iy] * ((n[i * nb + b] * tmp[i * nb * nm + b * nm + ix] / nbar[i * nb + b]) - tmp[i * nb * nm + b * nm + ix]);  // i b m
					//res[idx] += suma[i] * ((n[i * nb + b] * tmp[i * nb * nm + b * nm + ix] / nbar[i * nb + b]) - tmp[i * nb * nm + b * nm + ix]);  // i b m
					//res[idx] += (1 - (n[i * nb + b] / nbar[i * nb + b])) * grad_n[(i * nb * nm) + (b * nm) + ix];
				}
			}
		}
	}


	__global__ void nbard_k(float* source, float* lintg, float* res, float* bnmean, const int ni, const int nb, const int ne, const int nm) {
		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nb + ix; // ni X nb
		float tmp;

		if (ix < nb && iy < ni) {
			res[idx] = 0;
			for (int e = 0; e < ne; e++) {
				tmp = ((e * 0.1 + 15) - bnmean[iy % ni]);
				res[idx] += -source[iy * nb * ne + ix * ne + e] * tmp * tmp * lintg[(iy / nb) * ne + e] * 0.5;
			}

		}
	}

	__global__ void gradp(float* n, float* nbar, float* nbard, float* res, const int nd, const int ni, const int nb) {
		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nd + ix; // nd X nd

		if (ix < nd && iy < nd) {
			res[idx] = 0;

			for (int i = 0; i < ni; i++) {
				for (int b = 0; b < nb; b++) {
					res[idx] += (1 - (n[i * nb + b] / nbar[i * nb + b])) * -nbard[i * nb + b];  // i b m
				}
			}
		}
	}


	//for (int j = 0; j < nd * nd; j++) {
	//	//std::cout << "\r" << j << "/" << nd * nd;
	//	for (int k = 0; k < nm; k++) {
	//		for (int i = 0; i < block_size; i++) {
	//			for (int b = 0; b < nb; b++) {
	//				grad_q[j][k] += -smr[i][j] * (1 - (n[i][b] / nbar[i][b])) * grad_n[b][k];
	//			}
	//		}
	//	}
	//}

	__global__ void lintg(float* sysmat, float* matfrac, float* matatn, float* res, int nd, int ni, int ne, int nm) {

		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * ne + ix;

		if (ix < ne && iy < ni) {
			res[idx] = 0;

			for (int j = 0; j < nd * nd; j++) {
				for (int k = 0; k < nm; k++) {
					res[idx] += sysmat[iy * nd * nd + j] * matfrac[j * nm + k] * matatn[k * ne + ix];
				}
			}
		}
		//lintg[i][e] += smr[i][j] * matfrac[j][k] * matatn[k][e];
	}

	//h[j][k][m] += smr[i][j] * (n[i][b] / nbar[i][b]) * suma[i] * source[v_begin + i][b][e] * matatn[k][e] * matatn[m][e] * lintg[i][e];

//res1[i][b][k][m]

	__global__ void hesseq1(float* source, float* matatn, float* lintg, float* res, const int ni, const int nb, const int ne, const int nm) {
		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nm * nm + ix; // ni*nb X nm*nm

		if (ix < nm * nm && iy < ni * nb) {
			res[idx] = 0;
			for (int e = 0; e < ne; e++) {
				//res[idx] += source[(iy / ni) * nb * ne + (iy % nb) * ne + e] * matatn[(ix % nm) * ne + e] * matatn[(ix / nm) * ne + e] * lintg[(iy / ni) * ne + e];
				res[idx] += source[(iy / nb) * nb * ne + (iy % nb) * ne + e] * matatn[(ix / nm) * ne + e] * matatn[(ix % nm) * ne + e] * lintg[(iy / nb) * ne + e];
			}
			//res[idx] /= 10;
		}
	}
	__global__ void hesseq2(float* sysmat, float* n, float* nbar, float* suma, float* tmp, float* res, const int nd, const int ni, const int nb, const int nm) {
		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nm * nm + ix; // nd*nd X nm*nm

		if (ix < nm * nm && iy < nd * nd) {
			res[idx] = 0;
			for (int i = 0; i < ni; i++) {
				for (int b = 0; b < nb; b++) {
					//res[idx] += sysmat[i * nd * nd + iy] * (n[i * nb + b] / nbar[i * nb + b]) * suma[i] * tmp[i * nb * nm * nm + b * nm * nm + (ix / nm) * nm + (ix % nm)];
					
					//res[idx] += sysmat[i * nd * nd + iy] * (n[i * nb + b] / nbar[i * nb + b]) * suma[i] * tmp[i * nb * nm * nm + b * nm * nm + ix];
					res[idx] += sysmat[i * nd * nd + iy] * suma[i] * tmp[i * nb * nm * nm + b * nm * nm + ix] * 0.5;
				}
			}
		}
	}

	__global__ void hesseqp(float* source, float* matatn, float* lintg, float* res, const int ni, const int nb, const int ne, const int nm) {
		unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int idx = iy * nm + ix; // ni*nb X nm

		if (ix < nm && iy < ni * nb) {
			res[idx] = 0;
			for (int e = 0; e < ne; e++) {
				//res[idx] += source[(iy / ni) * nb * ne + (iy % nb) * ne + e] * matatn[(ix % nm) * ne + e] * matatn[(ix / nm) * ne + e] * lintg[(iy / ni) * ne + e];
				res[idx] += source[(iy / nb) * nb * ne + (iy % nb) * ne + e] * matatn[ix * ne + e] * lintg[(iy / nb) * ne + e];
			}
		}
	}
	//__global__ void hesseqp2(float* nbar, float* nbard, float* suma, float* tmp, float* res, const int nd, const int ni, const int nb, const int nm) {
	//	unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
	//	unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
	//	unsigned int idx = iy * nm + ix; // nd*nd X nm

	//	if (ix < nm && iy < nd * nd) {
	//		res[idx] = 0;
	//		for (int i = 0; i < ni; i++) {
	//			for (int b = 0; b < nb; b++) {
	//				res[idx] += suma[i] * (nbard[i * nb + b] / nbar[i * nb + b]) * tmp[i * nb * nm + b * nm + ix] * 0.5;
	//			}
	//		}
	//	}
	//}

	//__global__ void hessep(float* n, float* nbar, float* nbard, float res, const int nd, const int ni, const int nb, const int nm) {
	//	unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
	//	unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
	//	unsigned int idx = iy * nm + ix; // nd*nd X nm

	//	if (ix < nm && iy < nd * nd) {
	//		res[idx] = 0;
	//		float tmp;
	//		for (int i = 0; i < ni; i++) {
	//			for (int b = 0; b < nb; b++) {
	//				tmp = nbard[i * nb + b] / nbar[i * nb + b];
	//				res[idx] += n[i * nb + b] * tmp * tmp + (1 - n[i * nb + b] / nbar[i * nb + b]) * nbard[i * nb + b] * 0.5;
	//			}
	//		}
	//	}
	//}

	//hesseqp1 << <grid1, block >> > (d_source, d_matatn, d_lintg, d_res1, ni, nb, ne, nm); // res1: ni*nb X nm
	//CHECK(cudaDeviceSynchronize());
	//CHECK(cudaGetLastError());

	//int count = 0;

	//hesseqp2 << <grid2, block >> > (d_nbar, d_suma, d_res1, d_res2, nd, ni, nb, nm); // res2: nd*nd X nm
	//CHECK(cudaDeviceSynchronize());
	//CHECK(cudaGetLastError());
	//CHECK(cudaMemcpy(res2, d_res2, nBytes_res2, cudaMemcpyDeviceToHost));


	//	hesseq2 << <grid2, block >> > (d_sysmat, d_n, d_nbar, d_suma, d_res1, d_res2, nd, ni, nb, nm);

	//for (int j = 0; j < nd * nd; j++) {
	//	std::cout << "\r" << j << "/" << nd * nd;
	//	for (int k = 0; k < nm; k++) {
	//		for (int m = 0; m < nm; m++) {
	//			for (int i = 0; i < block_size; i++) {
	//				for (int b = 0; b < nb; b++) {
	//					for (int e = 0; e < ne; e++) {
	//						h[j][k][m] += smr[i][j] * (n[i][b] / nbar[i][b]) * suma[i] * source[v_begin + i][b][e] * matatn[k][e] * matatn[m][e] * lintg[i][e];
	//					}
	//				}
	//			}
	//		}
	//	}
	//}



	void calc_product(float* a, float* b, float* res, int nxa, int nya, int nxb) {

		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));

		int nBytes_a = nya * nxa * sizeof(float);
		int nBytes_b = nxa * nxb * sizeof(float);
		int nBytes_r = nya * nxb * sizeof(float);

		float *d_a, *d_b, *d_res;
		CHECK(cudaMalloc((void**)&d_a, nBytes_a));
		CHECK(cudaMalloc((void**)&d_b, nBytes_b));
		CHECK(cudaMalloc((void**)&d_res, nBytes_r));
		CHECK(cudaMemcpy(d_a, a, nBytes_a, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_b, b, nBytes_b, cudaMemcpyHostToDevice));
		//CHECK(cudaMemcpy(d_res, res, nBytes_r, cudaMemcpyHostToDevice));

		int dimx = 32;
		int dimy = 32;
		dim3 block(dimx, dimy);
		dim3 grid((nxb + block.x - 1) / block.x, (nya + block.y - 1) / block.y);


		//double iStart = cpuSecond();
		product <<< grid, block >>> (d_a, d_b, d_res, nxa, nya, nxb);
		CHECK(cudaDeviceSynchronize());
		//double iElaps = cpuSecond() - iStart;

		CHECK(cudaGetLastError());
		CHECK(cudaMemcpy(res, d_res, nBytes_r, cudaMemcpyDeviceToHost));

		CHECK(cudaFree(d_a));
		CHECK(cudaFree(d_b));
		CHECK(cudaFree(d_res));
		CHECK(cudaDeviceReset());
	}

	float** calc_product(float** a, float** b, int nxa, int nya, int nxb) {
		
		float* a1 = (float*)malloc(sizeof(float) * nya * nxa);
		float* b1 = (float*)malloc(sizeof(float) * nxa * nxb);
		float* r1 = (float*)malloc(sizeof(float) * nya * nxb);

		float** res = (float**)malloc(sizeof(float*) * nya);
		for (int i = 0; i < nya; i++) {
			res[i] = (float*)malloc(sizeof(float) * nxb);
		}

		for (int i = 0; i < nya; i++) {
			for (int j = 0; j < nxa; j++) {
				a1[i * nxa + j] = a[i][j];
			}
		}

		for (int i = 0; i < nxa; i++) {
			for (int j = 0; j < nxb; j++) {
				b1[i * nxb + j] = b[i][j];
			}
		}

		for (int i = 0; i < nya; i++) {
			for (int j = 0; j < nxb; j++) {
				r1[i * nxb + j] = res[i][j];
			}
		}

		calc_product(a1, b1, r1, nxa, nya, nxb);

		for (int i = 0; i < nya; i++) {
			for (int j = 0; j < nxb; j++) {
				res[i][j] = r1[i * nxb + j];
			}
		}

		free(a1);
		free(b1);
		free(r1);

		return res;
	}

	void init_matrix(float* a, int nx, int ny, float value) {

		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));

		int nBytes = nx * ny * sizeof(float);

		float* d_a;
		CHECK(cudaMalloc((void**)&d_a, nBytes));
		CHECK(cudaMemcpy(d_a, a, nBytes, cudaMemcpyHostToDevice));

		int dimx = 32;
		int dimy = 32;
		dim3 block(dimx, dimy);
		dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);


		//double iStart = cpuSecond();
		init << < grid, block >> > (d_a, nx, ny, value);
		CHECK(cudaDeviceSynchronize());
		//double iElaps = cpuSecond() - iStart;

		CHECK(cudaGetLastError());
		CHECK(cudaMemcpy(a, d_a, nBytes, cudaMemcpyDeviceToHost));

		CHECK(cudaFree(d_a));
		CHECK(cudaDeviceReset());

	}

	void init_matrix(float** a, int nx, int ny, float value) {

		float* a1 = (float*)malloc(sizeof(float) * nx * ny);

		for (int i = 0; i < ny; i++) {
			for (int j = 0; j < nx; j++) {
				a1[i * nx + j] = a[i][j];
			}
		}

		init_matrix(a1, nx, ny, value);

		for (int i = 0; i < ny; i++) {
			for (int j = 0; j < nx; j++) {
				a[i][j] = a1[i * nx + j];
			}
		}

		free(a1);
	}

	float** calc_T(float** a, int nx, int ny) {

		float* a1 = (float*)malloc(sizeof(float) * nx * ny);

		float** res = (float**)malloc(sizeof(float*) * nx);
		for (int i = 0; i < nx; i++) {
			res[i] = (float*)malloc(sizeof(float) * ny);
		}

		for (int i = 0; i < ny; i++) {
			for (int j = 0; j < nx; j++) {
				a1[i * nx + j] = a[i][j];
			}
		}

		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));

		int nBytes = nx * ny * sizeof(float);

		float* d_a;
		float* d_res;
		CHECK(cudaMalloc((void**)&d_a, nBytes));
		CHECK(cudaMalloc((void**)&d_res, nBytes));
		CHECK(cudaMemcpy(d_a, a, nBytes, cudaMemcpyHostToDevice));

		int dimx = 32;
		int dimy = 32;
		dim3 block(dimx, dimy);
		dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);


		//double iStart = cpuSecond();
		transpose_col << < grid, block >> > (d_a, d_res, nx, ny);
		CHECK(cudaDeviceSynchronize());
		//double iElaps = cpuSecond() - iStart;

		CHECK(cudaGetLastError());
		CHECK(cudaMemcpy(res, d_res, nBytes, cudaMemcpyDeviceToHost));

		CHECK(cudaFree(d_a));
		CHECK(cudaFree(d_res));
		CHECK(cudaDeviceReset());

		for (int i = 0; i < ny; i++) {
			for (int j = 0; j < nx; j++) {
				a[i][j] = a1[i * nx + j];
			}
		}

		free(a1);

		return a;

	}

	float** calc_n(float*** source, int ni, int nb, int ne, int v_begin) {  //n[i][b] += source[v_begin + i][b][e];

		float* source_block = (float*)malloc(sizeof(float) * ni * nb * ne);
		float* res1 = (float*)malloc(sizeof(float) * ni * nb);
		float** res = (float**)malloc(sizeof(float*) * ni);
		for (int i = 0; i < ni; i++) {
			res[i] = (float*)malloc(sizeof(float) * nb);
		}


		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nb; j++) {
				for (int k = 0; k < ne; k++) {
					source_block[(i * nb * ne) + (j * ne) + k] = source[i + v_begin][j][k];
				}
			}
		}

		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));

		int nBytes_s = ni * nb * ne * sizeof(float);
		int nBytes_r = ni * nb * sizeof(float);

		float* d_source;
		float* d_res;
		CHECK(cudaMalloc((void**)&d_source, nBytes_s));
		CHECK(cudaMalloc((void**)&d_res, nBytes_r));
		CHECK(cudaMemcpy(d_source, source_block, nBytes_s, cudaMemcpyHostToDevice));

		int dimx = 32;
		int dimy = 32;
		dim3 block(dimx, dimy);
		dim3 grid((nb + block.x - 1) / block.x, (ni + block.y - 1) / block.y);


		//double iStart = cpuSecond();
		intg_source << < grid, block >> > (d_source, d_res, ni, nb, ne);  //(float* source, float* res, int ni, int nb, int ne){
		CHECK(cudaDeviceSynchronize());
		//double iElaps = cpuSecond() - iStart;

		CHECK(cudaGetLastError());
		CHECK(cudaMemcpy(res1, d_res, nBytes_r, cudaMemcpyDeviceToHost));

		CHECK(cudaFree(d_source));
		CHECK(cudaFree(d_res));
		CHECK(cudaDeviceReset());

		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nb; j++) {
				res[i][j] = res1[i * nb + j];
			}
		}

		free(source_block);
		free(res1);

		return res;
	
	}

	void calc_nbard(float** res, float*** source, float** lintg, float* bnmean, int nd, int ni, int nb, int ne, int nm, int v_begin) {

		int nBytes_source = sizeof(float) * ni * nb * ne;
		int nBytes_lintg = sizeof(float) * ni * ne;
		int nBytes_bnmean = sizeof(float) * nb;
		int nBytes_res1 = sizeof(float) * ni * nb; // ni*nb

		float* _source = (float*)malloc(nBytes_source);
		float* _lintg = (float*)malloc(nBytes_lintg);
		float* _bnmean = (float*)malloc(nBytes_bnmean);
		float* res1 = (float*)malloc(nBytes_res1);

		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nb; j++) {
				for (int k = 0; k < ne; k++) {
					_source[(i * nb * ne) + (j * ne) + k] = source[i + v_begin][j][k];
				}
			}
			for (int j = 0; j < ne; j++) {
				_lintg[(i * ne) + j] = lintg[i][j];
			}
		}

		for (int i = 0; i < nb; i++) {
			_bnmean[i] = bnmean[i];
		}

		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));

		float* d_source;
		float* d_lintg;
		float* d_bnmean;
		float* d_res1;

		CHECK(cudaMalloc((void**)&d_source, nBytes_source));
		CHECK(cudaMalloc((void**)&d_lintg, nBytes_lintg));
		CHECK(cudaMalloc((void**)&d_bnmean, nBytes_bnmean));
		CHECK(cudaMalloc((void**)&d_res1, nBytes_res1));

		CHECK(cudaMemcpy(d_source, _source, nBytes_source, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_lintg, _lintg, nBytes_lintg, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_bnmean, bnmean, nBytes_bnmean, cudaMemcpyHostToDevice));

		int dimx = 32;
		int dimy = 32;
		dim3 block(dimx, dimy);
		dim3 grid((nb + block.x - 1) / block.x, (ni + block.y - 1) / block.y);

		nbard_k << <grid, block >> > (d_source, d_lintg, d_res1, d_bnmean, ni, nb, ne, nm);
		CHECK(cudaDeviceSynchronize());
		CHECK(cudaGetLastError());
		CHECK(cudaMemcpy(res1, d_res1, nBytes_res1, cudaMemcpyDeviceToHost));


		CHECK(cudaFree(d_source));
		CHECK(cudaFree(d_lintg));
		CHECK(cudaFree(d_bnmean));
		CHECK(cudaFree(d_res1));
		CHECK(cudaDeviceReset());

		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nb; j++) {
				res[i][j] += res1[i * nb + j];
			}
		}

		free(_source);
		free(_lintg);
		free(_bnmean);
		free(res1);

	}

	void calc_lintg(float** res, float** sysmat, float** matfrac, float** matatn, int nd, int ni, int ne, int nm, float pixsize) {

		int nBytes_sysmat = sizeof(float) * ni * nd * nd;
		int nBytes_matfrac = sizeof(float) * nd * nd * nm;
		int nBytes_matatn = sizeof(float) * nm * ne;
		int nBytes_res1 = sizeof(float) * nd * nd * ne;
		int nBytes_res2 = sizeof(float) * ni * ne;
		int nBytes_res = sizeof(float) * ni * ne;

		float* _sysmat = (float*)malloc(nBytes_sysmat);
		float* _matfrac = (float*)malloc(nBytes_matfrac);
		float* _matatn = (float*)malloc(nBytes_matatn);
		float* res2 = (float*)malloc(nBytes_res2);

		float* res_tmp = (float*)malloc(nBytes_res1);

		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nd * nd; j++) {
				_sysmat[(i * nd * nd) + j] = sysmat[i][j] * pixsize;
			}
		}
		for (int i = 0; i < nm; i++) {
			for (int j = 0; j < ne; j++) {
				_matatn[(i * ne) + j] = matatn[i][j];
			}
		}
		for (int i = 0; i < nd * nd; i++) {
			for (int j = 0; j < nm; j++) {
				_matfrac[(i * nm) + j] = matfrac[i][j];
			}
		}

		//float** res = (float**)malloc(sizeof(float*) * ni); //ni * ne
		//for (int i = 0; i < ni; i++) {
		//	res[i] = (float*)malloc(sizeof(float) * ne);
		//}

		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));

		float* d_sysmat;
		float* d_matfrac;
		float* d_matatn;
		////float* d_tmp;
		//float* d_res1;
		//float* d_res2;
		float* d_res;


		CHECK(cudaMalloc((void**)&d_sysmat, nBytes_sysmat));
		CHECK(cudaMalloc((void**)&d_matfrac, nBytes_matfrac));
		CHECK(cudaMalloc((void**)&d_matatn, nBytes_matatn));
		CHECK(cudaMalloc((void**)&d_res, nBytes_res));
		//CHECK(cudaMalloc((void**)&d_res1, nBytes_res1));
		//CHECK(cudaMalloc((void**)&d_res2, nBytes_res2));

		CHECK(cudaMemcpy(d_sysmat, _sysmat, nBytes_sysmat, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_matfrac, _matfrac, nBytes_matfrac, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_matatn, _matatn, nBytes_matatn, cudaMemcpyHostToDevice));

		int dimx = 32;
		int dimy = 32;
		dim3 block(dimx, dimy);
		dim3 grid((ne + block.x - 1) / block.x, (ni + block.y - 1) / block.y);
		//im3 grid1((ne + block.x - 1) / block.x, ((nd * nd) + block.y - 1) / block.y);
		//dim3 grid2((ne + block.x - 1) / block.x, (ni + block.y - 1) / block.y);

		//lintg[i][e] += smr[i][j] * matfrac[j][k] * matatn[k][e];

		lintg <<<grid, block>>> (d_sysmat, d_matfrac, d_matatn, d_res, nd, ni, ne, nm);


		//product << < grid1, block >> > (d_matfrac, d_matatn, d_res1, nm, nd*nd, ne);
		//

		////for (int i = 0; i < nd * nd * ne; i+=1) {
		////	printf("%f,", res_test[i]);
		////}

		////lintg << <grid, block >> > (d_sysmat, d_matfrac, d_matatn, d_res, nd, ni, ne, nm);
		//CHECK(cudaDeviceSynchronize());
		//CHECK(cudaGetLastError());
		////CHECK(cudaMemcpy(res_tmp, d_res1, nBytes_res1, cudaMemcpyDeviceToHost));

		//product << < grid2, block >> > (d_sysmat, d_res1, d_res2, nd * nd, ni, ne);
		CHECK(cudaDeviceSynchronize());
		CHECK(cudaGetLastError());
		CHECK(cudaMemcpy(res2, d_res, nBytes_res, cudaMemcpyDeviceToHost));

		CHECK(cudaFree(d_sysmat));
		CHECK(cudaFree(d_matfrac));
		CHECK(cudaFree(d_matatn));
		CHECK(cudaFree(d_res));
		//CHECK(cudaFree(d_res1));
		//CHECK(cudaFree(d_res2));
		CHECK(cudaDeviceReset());

		printf("\n 2nd \n");

		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < ne; j++) {
				res[i][j] = res2[i * ne + j];
				//if (res2[i] != 0) printf("%f,", res2[i]);
			}
		}

		free(_sysmat);
		free(_matfrac);
		free(_matatn);
		free(res2);
		free(res_tmp);

		//return res;
	}

	void calc_gradq(float** res, float** sysmat, float* n, float** nbar, float*** grad_n, int nd, int ni, int nb, int nm, int v_begin, float pixsize) {

		int nBytes_sysmat = sizeof(float) * ni * nd * nd;
		int nBytes_n = sizeof(float) * ni * nb;
		int nBytes_gradn = sizeof(float) * ni * nb * nm;
		int nBytes_res1 = sizeof(float) * nd * nd * nm; // nd*nd X nm

		//for (int j = 0; j < nd * nd; j++) {
		//	//std::cout << "\r" << j << "/" << nd * nd;
		//	for (int k = 0; k < nm; k++) {
		//		for (int i = 0; i < block_size; i++) {
		//			for (int b = 0; b < nb; b++) {
		//				grad_q[j][k] += -smr[i][j] * (1 - (n[i][b] / nbar[i][b])) * grad_n[b][k];
		//			}
		//		}
		//	}
		//}

		float* _sysmat = (float*)malloc(nBytes_sysmat);
		float* _n = (float*)malloc(nBytes_n);
		float* _nbar = (float*)malloc(nBytes_n);
		float* _gradn = (float*)malloc(nBytes_gradn);
		float* res1 = (float*)malloc(nBytes_res1);

		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nb; j++) {
				_n[(i * nb) + j] = n[(i + v_begin) + j * ni];
				_nbar[(i * nb) + j] = nbar[i][j];
				for (int k = 0; k < nm; k++) {
					_gradn[(i * nb * nm) + (j * nm) + k] = grad_n[i][j][k];
				}
			}
			for (int j = 0; j < nd * nd; j++) {
				_sysmat[(i * nd * nd) + j] = sysmat[i][j] * pixsize;
				/*if (sysmat[i][j] != 0) printf("\nhoge");*/
			}

		}
		//for (int i = 0; i < nb; i++) {
		//	for (int j = 0; j < nm; j++) {

		//		_gradn[(i * nm) + j] = grad_n[i][j];
		//	}
		//}


		//float** res = (float**)malloc(sizeof(float*) * nd * nd); //nd*nd X nm 
		//for (int i = 0; i < nd * nd; i++) {
		//	res[i] = (float*)malloc(sizeof(float) * nm);
		//}

		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));

		float* d_sysmat;
		float* d_n;
		float* d_nbar;
		float* d_gradn;
		float* d_res1;


		CHECK(cudaMalloc((void**)&d_sysmat, nBytes_sysmat));
		CHECK(cudaMalloc((void**)&d_n, nBytes_n));
		CHECK(cudaMalloc((void**)&d_nbar, nBytes_n));
		CHECK(cudaMalloc((void**)&d_gradn, nBytes_gradn));
		CHECK(cudaMalloc((void**)&d_res1, nBytes_res1));

		CHECK(cudaMemcpy(d_sysmat, _sysmat, nBytes_sysmat, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_n, _n, nBytes_n, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_nbar, _nbar, nBytes_n, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_gradn, _gradn, nBytes_gradn, cudaMemcpyHostToDevice));

		int dimx = 32;
		int dimy = 32;
		dim3 block(dimx, dimy);
		dim3 grid((nm + block.x - 1) / block.x, ((nd * nd) + block.y - 1) / block.y);

		//h[j][k][m] += smr[i][j] * (n[i][b] / nbar[i][b]) * suma[i] * source[v_begin + i][b][e] * matatn[k][e] * matatn[m][e] * lintg[i][e];

		//res1[i][b][k][m]


		gradq << <grid, block >> > (d_sysmat, d_n, d_nbar, d_gradn, d_res1, nd, ni, nb, nm);
		CHECK(cudaDeviceSynchronize());
		CHECK(cudaGetLastError());
		CHECK(cudaMemcpy(res1, d_res1, nBytes_res1, cudaMemcpyDeviceToHost));

		CHECK(cudaFree(d_sysmat));
		CHECK(cudaFree(d_n));
		CHECK(cudaFree(d_nbar));
		CHECK(cudaFree(d_gradn));
		CHECK(cudaFree(d_res1));
		CHECK(cudaDeviceReset());

		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nb; j++) {
				if(i % 100 == 0) printf("\nn,nbar, (n/nbar) = %f, %f, %f", _n[i * nb + j], nbar[i][j], _n[i * nb + j] / nbar[i][j]);
			}
		}

		//for (int i = 0; i < ni; i++) {
		//	for (int j = 0; j < nd * nd; j++) {
		//		if(sysmat[i][j] != 0) printf("\nsmr(%d, %d) = %f", i, j, sysmat[i][j]);
		//	}
		//}

		for (int i = 0; i < nd * nd; i++) {
			for (int j = 0; j < nm; j++) {
				res[i][j] = res1[i * nm + j];
				//if (res[i][j] != 0 && i % 1000 == 0) {
				//	printf("at(j,k) = %d, %d, value = %f \n", i, j, res[i][j]);
				//}
				//printf("%f,", res[i][j]);
			}
		}

		free(_sysmat);
		free(_n);
		free(_nbar);
		free(_gradn);
		free(res1);

		//return res;

	}

	void calc_gradq2(float** res, float** sysmat, float*** source, float* n, float** nbar, float** matatn, float** lintg, float* suma, int nd, int ni, int nb, int ne, int nm, int v_begin, float pixsize) {

		int nBytes_source = sizeof(float) * ni * nb * ne;
		int nBytes_sysmat = sizeof(float) * ni * nd * nd;
		int nBytes_n = sizeof(float) * ni * nb;
		int nBytes_matatn = sizeof(float) * nm * ne;
		int nBytes_lintg = sizeof(float) * ni * ne;
		int nBytes_suma = sizeof(float) * ni;
		int nBytes_res1 = sizeof(float) * ni * nb * nm; // ni*nb X nm
		int nBytes_res2 = sizeof(float) * nd * nd * nm; // nd*nd X nm

		float* _source = (float*)malloc(nBytes_source);
		float* _sysmat = (float*)malloc(nBytes_sysmat);
		float* _n = (float*)malloc(nBytes_n);
		float* _nbar = (float*)malloc(nBytes_n);
		float* _matatn = (float*)malloc(nBytes_matatn);
		float* _lintg = (float*)malloc(nBytes_lintg);
		float* res1 = (float*)malloc(nBytes_res1);
		float* res2 = (float*)malloc(nBytes_res2);

		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nb; j++) {
				for (int k = 0; k < ne; k++) {
					_source[(i * nb * ne) + (j * ne) + k] = source[i + v_begin][j][k];
				}
				_n[(i * nb) + j] = n[(i + v_begin) + j * ni];
				_nbar[(i * nb) + j] = nbar[i][j];
			}
			for (int j = 0; j < ne; j++) {
				_lintg[(i * ne) + j] = lintg[i][j];
			}
			for (int j = 0; j < nd * nd; j++) {
				_sysmat[(i * nd * nd) + j] = sysmat[i][j] * pixsize;
			}
		}
		for (int i = 0; i < nm; i++) {
			for (int j = 0; j < ne; j++) {
				_matatn[(i * ne) + j] = matatn[i][j];
			}
		}

		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));

		float* d_source;
		float* d_sysmat;
		float* d_n;
		float* d_nbar;
		float* d_matatn;
		float* d_lintg;
		float* d_suma;
		float* d_res1;
		float* d_res2;


		CHECK(cudaMalloc((void**)&d_source, nBytes_source));
		CHECK(cudaMalloc((void**)&d_sysmat, nBytes_sysmat));
		CHECK(cudaMalloc((void**)&d_n, nBytes_n));
		CHECK(cudaMalloc((void**)&d_nbar, nBytes_n));
		CHECK(cudaMalloc((void**)&d_matatn, nBytes_matatn));
		CHECK(cudaMalloc((void**)&d_lintg, nBytes_lintg));
		CHECK(cudaMalloc((void**)&d_suma, nBytes_suma));
		CHECK(cudaMalloc((void**)&d_res1, nBytes_res1));
		CHECK(cudaMalloc((void**)&d_res2, nBytes_res2));

		CHECK(cudaMemcpy(d_source, _source, nBytes_source, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_sysmat, _sysmat, nBytes_sysmat, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_n, _n, nBytes_n, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_nbar, _nbar, nBytes_n, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_matatn, _matatn, nBytes_matatn, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_lintg, _lintg, nBytes_lintg, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_suma, suma, nBytes_suma, cudaMemcpyHostToDevice));

		int dimx = 32;
		int dimy = 32;
		dim3 block(dimx, dimy);
		dim3 grid1((nm + block.x - 1) / block.x, ((ni * nb) + block.y - 1) / block.y);
		dim3 grid2((nm + block.x - 1) / block.x, ((nd * nd) + block.y - 1) / block.y);

		gradq2_1 << <grid1, block >> > (d_source, d_matatn, d_lintg, d_res1, ni, nb, ne, nm);
		CHECK(cudaDeviceSynchronize());
		CHECK(cudaGetLastError());

		printf("hogehoge");

		gradq2_2 << <grid2, block >> > (d_sysmat, d_n, d_nbar, d_suma, d_res1, d_res2, nd, ni, nb, nm);
		CHECK(cudaDeviceSynchronize());
		CHECK(cudaGetLastError());
		CHECK(cudaMemcpy(res2, d_res2, nBytes_res2, cudaMemcpyDeviceToHost));

		printf("piyo");

		CHECK(cudaFree(d_source));
		CHECK(cudaFree(d_sysmat));
		CHECK(cudaFree(d_n));
		CHECK(cudaFree(d_nbar));
		CHECK(cudaFree(d_matatn));
		CHECK(cudaFree(d_lintg));
		CHECK(cudaFree(d_suma));
		CHECK(cudaFree(d_res1));
		CHECK(cudaFree(d_res2));
		CHECK(cudaDeviceReset());

		printf("hagehage");

		for (int i = 0; i < nd * nd; i++) {
			for (int j = 0; j < nm; j++) {
				res[i][j] = res2[i * nm + j];
			}
			//if (res[i][0] < -1000) {
			//	printf("\ngradq, n/nbar  %f  %f", res[i][0], n[]);
			//}
		}

		free(_source);
		free(_sysmat);
		free(_n);
		free(_nbar);
		free(_matatn);
		free(_lintg);
		free(res1);
		free(res2);

	}

	float calc_gradp(float* n, float** nbar, float** nbard, int nd, int ni, int nb, int ne, int v_begin) {

		int nBytes_n = sizeof(float) * ni * nb;
		int nBytes_res = sizeof(float) * nd * nd; // nd X nd

		float* _n = (float*)malloc(nBytes_n);
		float* _nbar = (float*)malloc(nBytes_n);
		float* _nbard = (float*)malloc(nBytes_n);
		float* res1 = (float*)malloc(nBytes_res);

		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nb; j++) {
				_n[(i * nb) + j] = n[(i + v_begin) + j * ni];
				_nbar[(i * nb) + j] = nbar[i][j];
				_nbard[(i * nb) + j] = nbard[i][j];
			}
		}

		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));

		float* d_n;
		float* d_nbar;
		float* d_nbard;
		float* d_res1;

		CHECK(cudaMalloc((void**)&d_n, nBytes_n));
		CHECK(cudaMalloc((void**)&d_nbar, nBytes_n));
		CHECK(cudaMalloc((void**)&d_nbard, nBytes_n));

		CHECK(cudaMalloc((void**)&d_res1, nBytes_res));

		CHECK(cudaMemcpy(d_n, _n, nBytes_n, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_nbar, _nbar, nBytes_n, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_nbard, _nbard, nBytes_n, cudaMemcpyHostToDevice));

		int dimx = 32;
		int dimy = 32;
		dim3 block(dimx, dimy);
		dim3 grid((nd + block.x - 1) / block.x, (nd + block.y - 1) / block.y);

		gradp << <grid, block >> > (d_n, d_nbar, d_nbard, d_res1, nd, ni, nb);
		CHECK(cudaDeviceSynchronize());
		CHECK(cudaGetLastError());
		CHECK(cudaMemcpy(res1, d_res1, nBytes_res, cudaMemcpyDeviceToHost));

		printf("piyo");


		CHECK(cudaFree(d_n));
		CHECK(cudaFree(d_nbar));		
		CHECK(cudaFree(d_nbard));
		CHECK(cudaFree(d_res1));
		CHECK(cudaDeviceReset());

		printf("hagehage");


		float res = 0;// res2[0];

		for (int j = 0; j < nd * nd; j++) {
			res += res1[j];
		}

		free(_n);
		free(_nbar);
		free(_nbard);
		free(res1);

		return res;

	}

	void calc_hesseq(float*** res, float*** source, float** sysmat, float* n, float** nbar, float** matatn, float** lintg, float* suma, int nd, int ni, int nb, int ne, int nm, int v_begin, float pixsize) {

		int nBytes_source = sizeof(float) * ni * nb * ne;
		int nBytes_sysmat = sizeof(float) * ni * nd * nd;
		int nBytes_n = sizeof(float) * ni * nb;
		int nBytes_matatn = sizeof(float) * nm * ne;
		int nBytes_lintg = sizeof(float) * ni * ne;
		int nBytes_suma = sizeof(float) * ni;
		int nBytes_res1 = sizeof(float) * ni * nb * nm * nm; // ni*nb X nm*nm
		int nBytes_res2 = sizeof(float) * nd * nd * nm * nm;

		float* _source = (float*)malloc(nBytes_source);
		float* _sysmat = (float*)malloc(nBytes_sysmat);
		float* _n      = (float*)malloc(nBytes_n);
		float* _nbar   = (float*)malloc(nBytes_n);
		float* _matatn = (float*)malloc(nBytes_matatn);
		float* _lintg  = (float*)malloc(nBytes_lintg);
		//float* res1    = (float*)malloc(nBytes_res1);
		float* res2    = (float*)malloc(nBytes_res2);

		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nb; j++) {
				for (int k = 0; k < ne; k++) {
					_source[(i * nb * ne) + (j * ne) + k] = source[i + v_begin][j][k];
				}
				_n[(i * nb) + j] = n[(i + v_begin) + j * ni];
				_nbar[(i * nb) + j] = nbar[i][j];
			}
			for (int j = 0; j < nd * nd; j++) {
				_sysmat[(i * nd * nd) + j] = sysmat[i][j] * pixsize;
			}
			for (int j = 0; j < ne; j++) {
				_lintg[(i * ne) + j] = lintg[i][j];
			}
		}
		for (int i = 0; i < nm; i++) {
			for (int j = 0; j < ne; j++) {
				_matatn[(i * ne) + j] = matatn[i][j];
			}
		}


		//float*** res_tmp = (float***)malloc(sizeof(float**) * ni * nb); //ni*nb X nm X nm
		//for (int i = 0; i < ni * nb; i++) {
		//	res_tmp[i] = (float**)malloc(sizeof(float*) * nm);
		//	for (int j = 0; j < nm; j++) {
		//		res_tmp[i][j] = (float*)malloc(sizeof(float) * nm);
		//	}
		//}


		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));


		float* d_source;
		float* d_sysmat;
		float* d_n;
		float* d_nbar;
		float* d_matatn;
		float* d_lintg;
		float* d_suma;
		float* d_res1;
		float* d_res2;


		CHECK(cudaMalloc((void**)&d_source, nBytes_source));
		CHECK(cudaMalloc((void**)&d_sysmat, nBytes_sysmat));
		CHECK(cudaMalloc((void**)&d_n, nBytes_n));
		CHECK(cudaMalloc((void**)&d_nbar, nBytes_n));
		CHECK(cudaMalloc((void**)&d_matatn, nBytes_matatn));
		CHECK(cudaMalloc((void**)&d_lintg, nBytes_lintg));
		CHECK(cudaMalloc((void**)&d_suma, nBytes_suma));
		CHECK(cudaMalloc((void**)&d_res1, nBytes_res1));
		CHECK(cudaMalloc((void**)&d_res2, nBytes_res2));

		CHECK(cudaMemcpy(d_source, _source, nBytes_source, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_sysmat, _sysmat, nBytes_sysmat, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_n, _n, nBytes_n, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_nbar, _nbar, nBytes_n, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_matatn, _matatn, nBytes_matatn, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_lintg, _lintg, nBytes_lintg, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_suma, suma, nBytes_suma, cudaMemcpyHostToDevice));

		int dimx = 32;
		int dimy = 32;
		dim3 block(dimx, dimy);
		dim3 grid1(((nm * nm) + block.x - 1) / block.x, ((ni * nb) + block.y - 1) / block.y);
		dim3 grid2(((nm * nm) + block.x - 1) / block.x, ((nd * nd) + block.y - 1) / block.y);

		//h[j][k][m] += smr[i][j] * (n[i][b] / nbar[i][b]) * suma[i] * source[v_begin + i][b][e] * matatn[k][e] * matatn[m][e] * lintg[i][e];

		//res1[i][b][k][m]

		hesseq1 << <grid1, block >> > (d_source, d_matatn, d_lintg, d_res1, ni, nb, ne, nm); // res1: ni*nb X nm*nm
		CHECK(cudaDeviceSynchronize());
		CHECK(cudaGetLastError());

		//res2 = h[j][k][m]
		//CHECK(cudaMemcpy(res1, d_res1, nBytes_res1, cudaMemcpyDeviceToHost));

		int count = 0;
		//for (int i = 0; i < ni * nb; i++) {
		//	for (int j = 0; j < nm; j++) {
		//		for (int k = 0; k < nm; k++) {
		//			res_tmp[i][j][k] = res1[i * nm * nm + j * nm + k];
		//			if (j != k && res_tmp[i][j][k] != 0) {
		//				count++;
		//			}
		//		}
		//	}
		//}
		//printf("not 0 count: %d / %d", count, ni * nb * nm * nm);


		hesseq2 << <grid2, block >> > (d_sysmat, d_n, d_nbar, d_suma, d_res1, d_res2, nd, ni, nb, nm);
		CHECK(cudaDeviceSynchronize());
		CHECK(cudaGetLastError());
		CHECK(cudaMemcpy(res2, d_res2, nBytes_res2, cudaMemcpyDeviceToHost));

		//res[idx] += sysmat[i * nd * nd + iy] * (n[i * nb + b] / nbar[i * nb + b]) * suma[i] * tmp[i * nb * nm * nm + b * nm * nm + ix];

		/*hesseq << <grid, block >> > (d_source, d_sysmat, d_n, d_nbar, d_matatn, d_lintg, d_suma, d_res, nd, ni, nb, ne, nm); 
		CHECK(cudaDeviceSynchronize());
		CHECK(cudaGetLastError());
		CHECK(cudaMemcpy(res1, d_res, nBytes_res, cudaMemcpyDeviceToHost));*/

		CHECK(cudaFree(d_source));
		CHECK(cudaFree(d_sysmat));
		CHECK(cudaFree(d_n));
		CHECK(cudaFree(d_nbar));
		CHECK(cudaFree(d_matatn));
		CHECK(cudaFree(d_lintg));
		CHECK(cudaFree(d_suma));
		CHECK(cudaFree(d_res1));
		CHECK(cudaFree(d_res2));
		CHECK(cudaDeviceReset());

		count = 0;

		for (int i = 0; i < nd * nd; i++) {
			for (int j = 0; j < nm; j++) {
				for (int k = 0; k < nm; k++) {
					res[i][j][k] = res2[i * nm * nm + j * nm + k];
					if (j != k && res[i][j][k] != 0) {
						count++;
					}
				}
			}
		}

		printf("(2)not 0 count: %d / %d", count, nd * nd * nm * nm);

		free(_source);
		free(_sysmat);
		free(_n);
		free(_nbar);
		free(_matatn);
		free(_lintg);
		//free(res1);
		free(res2);

		//return res;

	}

	void calc_hesseqp(float* res, float*** source, float** nbar, float** nbard, float** matatn, float** lintg, float* suma, int nd, int ni, int nb, int ne, int nm, int v_begin) {

		int nBytes_source = sizeof(float) * ni * nb * ne;
		int nBytes_n = sizeof(float) * ni * nb;
		int nBytes_matatn = sizeof(float) * nm * ne;
		int nBytes_lintg = sizeof(float) * ni * ne;
		int nBytes_suma = sizeof(float) * ni;
		int nBytes_res1 = sizeof(float) * ni * nb * nm; // ni*nb X nm*nm
		//int nBytes_res2 = sizeof(float) * nd * nd * nm;

		float* _source = (float*)malloc(nBytes_source);
		float* _nbar = (float*)malloc(nBytes_n);
		float* _nbard = (float*)malloc(nBytes_n);
		float* _matatn = (float*)malloc(nBytes_matatn);
		float* _lintg = (float*)malloc(nBytes_lintg);
		float* res1 = (float*)malloc(nBytes_res1);
		//float* res2 = (float*)malloc(nBytes_res2);

		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nb; j++) {
				for (int k = 0; k < ne; k++) {
					_source[(i * nb * ne) + (j * ne) + k] = source[i + v_begin][j][k];
				}
				_nbar[(i * nb) + j] = nbar[i][j];
				_nbard[(i * nb) + j] = nbard[i][j];
			}
			for (int j = 0; j < ne; j++) {
				_lintg[(i * ne) + j] = lintg[i][j];
			}
		}
		for (int i = 0; i < nm; i++) {
			for (int j = 0; j < ne; j++) {
				_matatn[(i * ne) + j] = matatn[i][j];
			}
		}


		int dev = 0;
		cudaDeviceProp deviceprop;
		CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		CHECK(cudaSetDevice(dev));


		float* d_source;
		float* d_nbar;
		float* d_nbard;
		float* d_matatn;
		float* d_lintg;
		float* d_suma;
		float* d_res1;
		//float* d_res2;


		CHECK(cudaMalloc((void**)&d_source, nBytes_source));
		CHECK(cudaMalloc((void**)&d_nbar, nBytes_n));
		CHECK(cudaMalloc((void**)&d_nbard, nBytes_n));
		CHECK(cudaMalloc((void**)&d_matatn, nBytes_matatn));
		CHECK(cudaMalloc((void**)&d_lintg, nBytes_lintg));
		CHECK(cudaMalloc((void**)&d_suma, nBytes_suma));
		CHECK(cudaMalloc((void**)&d_res1, nBytes_res1));
		//CHECK(cudaMalloc((void**)&d_res2, nBytes_res2));

		CHECK(cudaMemcpy(d_source, _source, nBytes_source, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_nbar, _nbar, nBytes_n, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_nbard, _nbar, nBytes_n, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_matatn, _matatn, nBytes_matatn, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_lintg, _lintg, nBytes_lintg, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_suma, suma, nBytes_suma, cudaMemcpyHostToDevice));

		int dimx = 32;
		int dimy = 32;
		dim3 block(dimx, dimy);
		dim3 grid((nm + block.x - 1) / block.x, ((ni * nb) + block.y - 1) / block.y);
		//dim3 grid2((nm + block.x - 1) / block.x, ((nd * nd) + block.y - 1) / block.y);


		hesseqp << <grid, block >> > (d_source, d_matatn, d_lintg, d_res1, ni, nb, ne, nm); // res1: ni*nb X nm
		CHECK(cudaDeviceSynchronize());
		CHECK(cudaGetLastError());
		CHECK(cudaMemcpy(res1, d_res1, nBytes_res1, cudaMemcpyDeviceToHost));

		for (int m = 0; m < nm; m++) {
			for (int i = 0; i < ni; i++) {
				for (int b = 0; b < nb; b++) {
					res[m] += suma[i] * (nbard[i][b] / nbar[i][b]) * res1[i * nb * nm + b * nm + m] * 0.5;
					//res[m] = 0;
				}
			}
		}



		//int count = 0;

		//hesseqp2 << <grid2, block >> > (d_nbar, d_nbard, d_suma, d_res1, d_res2, nd, ni, nb, nm); // res2: nd*nd X nm
		//CHECK(cudaDeviceSynchronize());
		//CHECK(cudaGetLastError());
		//CHECK(cudaMemcpy(res2, d_res2, nBytes_res2, cudaMemcpyDeviceToHost));

		//__global__ void hesseqp2(float* nbar, float* nbard, float* suma, float* tmp, float* res, const int nd, const int ni, const int nb, const int nm) {
		//	unsigned int ix = blockDim.x * blockIdx.x + threadIdx.x;
		//	unsigned int iy = blockDim.y * blockIdx.y + threadIdx.y;
		//	unsigned int idx = iy * nm + ix; // nd*nd X nm

		//	if (ix < nm && iy < nd * nd) {
		//		res[idx] = 0;
		//		for (int i = 0; i < ni; i++) {
		//			for (int b = 0; b < nb; b++) {
		//				res[idx] += suma[i] * (nbard[i * nb + b] / nbar[i * nb + b]) * tmp[i * nb * nm + b * nm + ix] * 0.5;
		//			}
		//		}
		//	}
		//}


		CHECK(cudaFree(d_source));
		CHECK(cudaFree(d_nbar));
		CHECK(cudaFree(d_nbard));
		CHECK(cudaFree(d_matatn));
		CHECK(cudaFree(d_lintg));
		CHECK(cudaFree(d_suma));
		CHECK(cudaFree(d_res1));
		//CHECK(cudaFree(d_res2));
		CHECK(cudaDeviceReset());

		//count = 0;

		//for (int i = 0; i < nd * nd; i++) {
		//	for (int j = 0; j < nm; j++) {
		//		res[i][j] = res2[i * nm + j];
		//	}
		//}

		//for (int i = 0; i < nd * nd; i++) {
		//	for (int j = 0; j < nm; j++) {
		//		res[i] += res2[i * nm + j];
		//	}
		//}


		free(_source);
		free(_nbar);
		free(_nbard);
		free(_matatn);
		free(_lintg);
		free(res1);
		//free(res2);

	}

	float calc_hessep(float* n, float** nbar, float** nbard, int ni, int nb) {

		float res = 0;
		float tmp;
		for (int i = 0; i < ni; i++) {
			for (int b = 0; b < nb; b++) {
				tmp = nbard[i][b] / nbar[i][b];
				res += n[i * nb + b] * tmp * tmp + (1 - n[i * nb + b] / nbar[i][b]) * nbard[i][b] * 0.5;
			}
		}

		return res;

		//int nBytes_source = sizeof(float) * ni * nb * ne;
		//int nBytes_n = sizeof(float) * ni * nb;
		//int nBytes_matatn = sizeof(float) * nm * ne;
		//int nBytes_lintg = sizeof(float) * ni * ne;
		//int nBytes_suma = sizeof(float) * ni;
		//int nBytes_res1 = sizeof(float) * ni * nb * nm; // ni*nb X nm*nm
		//int nBytes_res2 = sizeof(float) * nd * nd * nm;

		//float* _source = (float*)malloc(nBytes_source);
		//float* _nbar = (float*)malloc(nBytes_n);
		//float* _nbard = (float*)malloc(nBytes_n);
		//float* _matatn = (float*)malloc(nBytes_matatn);
		//float* _lintg = (float*)malloc(nBytes_lintg);
		//float* res2 = (float*)malloc(nBytes_res2);

		//for (int i = 0; i < ni; i++) {
		//	for (int j = 0; j < nb; j++) {
		//		for (int k = 0; k < ne; k++) {
		//			_source[(i * nb * ne) + (j * ne) + k] = source[i + v_begin][j][k];
		//		}
		//		_nbar[(i * nb) + j] = nbar[i][j];
		//		_nbard[(i * nb) + j] = nbard[i][j];
		//	}
		//	for (int j = 0; j < ne; j++) {
		//		_lintg[(i * ne) + j] = lintg[i][j];
		//	}
		//}
		//for (int i = 0; i < nm; i++) {
		//	for (int j = 0; j < ne; j++) {
		//		_matatn[(i * ne) + j] = matatn[i][j];
		//	}
		//}


		//int dev = 0;
		//cudaDeviceProp deviceprop;
		//CHECK(cudaGetDeviceProperties(&deviceprop, dev));
		//CHECK(cudaSetDevice(dev));


		//float* d_source;
		//float* d_nbar;
		//float* d_nbard;
		//float* d_matatn;
		//float* d_lintg;
		//float* d_suma;
		//float* d_res1;
		//float* d_res2;


		//CHECK(cudaMalloc((void**)&d_source, nBytes_source));
		//CHECK(cudaMalloc((void**)&d_nbar, nBytes_n));
		//CHECK(cudaMalloc((void**)&d_nbard, nBytes_n));
		//CHECK(cudaMalloc((void**)&d_matatn, nBytes_matatn));
		//CHECK(cudaMalloc((void**)&d_lintg, nBytes_lintg));
		//CHECK(cudaMalloc((void**)&d_suma, nBytes_suma));
		//CHECK(cudaMalloc((void**)&d_res1, nBytes_res1));
		//CHECK(cudaMalloc((void**)&d_res2, nBytes_res2));

		//CHECK(cudaMemcpy(d_source, _source, nBytes_source, cudaMemcpyHostToDevice));
		//CHECK(cudaMemcpy(d_nbar, _nbar, nBytes_n, cudaMemcpyHostToDevice));
		//CHECK(cudaMemcpy(d_nbard, _nbar, nBytes_n, cudaMemcpyHostToDevice));
		//CHECK(cudaMemcpy(d_matatn, _matatn, nBytes_matatn, cudaMemcpyHostToDevice));
		//CHECK(cudaMemcpy(d_lintg, _lintg, nBytes_lintg, cudaMemcpyHostToDevice));
		//CHECK(cudaMemcpy(d_suma, suma, nBytes_suma, cudaMemcpyHostToDevice));

		//int dimx = 32;
		//int dimy = 32;
		//dim3 block(dimx, dimy);
		//dim3 grid1((nm + block.x - 1) / block.x, ((ni * nb) + block.y - 1) / block.y);
		//dim3 grid2((nm + block.x - 1) / block.x, ((nd * nd) + block.y - 1) / block.y);


		//hessep << <grid, block >> > (d_nbar, d_nbard, d_suma, d_res1, d_res2, nd, ni, nb, nm); // res2: nd*nd X nm
		//CHECK(cudaDeviceSynchronize());
		//CHECK(cudaGetLastError());
		//CHECK(cudaMemcpy(res2, d_res2, nBytes_res2, cudaMemcpyDeviceToHost));


		//CHECK(cudaFree(d_source));
		//CHECK(cudaFree(d_nbar));
		//CHECK(cudaFree(d_nbard));
		//CHECK(cudaFree(d_matatn));
		//CHECK(cudaFree(d_lintg));
		//CHECK(cudaFree(d_suma));
		//CHECK(cudaFree(d_res1));
		//CHECK(cudaFree(d_res2));
		//CHECK(cudaDeviceReset());

		//count = 0;

		//for (int i = 0; i < nd * nd; i++) {
		//	for (int j = 0; j < nm; j++) {
		//		res[i][j] = res2[i * nm + j];
		//	}
		//}


		//free(_source);
		//free(_nbar);
		//free(_nbard);
		//free(_matatn);
		//free(_lintg);
		////free(res1);
		//free(res2);
	}

}
