#include "CalcImgdiff.cuh"


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

/*

namespace Reconstruction {

	__global__ void calc_imgdiff(float int row, int num_col,) {

	}

	void CalcImgdiff(float* attenu, float* imgdiff, const float* sino, const float* elem, const float* rowptr, const float* colind, int blocksize, int nd, int nv, int ne)
	{


		float* _sysmatrow = (float*)malloc(sizeof(float) * nd * nd);

		bool block = true;
		bool init = true;

		bool eachblockmode = false;

		for (int j = 0; j < nd * nd; j++) {
			imgdiff[j] = 0;
		}

		calc_imgdiff();

		for (int j = 0; j < blocksize; j++) {
			calc_imgdiff();
			//(*_sysmatblock).Extract_row_dense(j, nd * nd, _sysmatrow);
			calc_imgdiff(imgdiff, _sysmatrow, attenu, sino[i * blocksize + j], nd * nd);
			Reconstruction::mul_array1(imgdiff, nd, nd * nd);
		}
	}

	void Extract_row_dense(int row, int num_col, float* vec) {

		for (int i = 0; i < num_col; i++) {
			vec[i] = 0;
		}

		for (int i = rowptr.get()[row]; i < rowptr.get()[row + 1]; i++) {
			int tmp = colind.get()[i];
			vec[tmp] = elements.get()[i];
		}

	}

	void IterationRec::calc_imgdiff(float* idiff, float* smr, float* atn, float sn, int size) const {
		float sys_atn = Reconstruction::dot_array(smr, atn, size);
		float sys_sys = Reconstruction::dot_array(smr, smr, size);

		Reconstruction::mul_array1(smr, ((sys_atn - sn) / sys_sys), size);
		Reconstruction::add_array(idiff, smr, size);
	}
}

*/