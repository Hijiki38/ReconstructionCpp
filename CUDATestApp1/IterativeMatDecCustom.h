#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <numeric>
#include <math.h>
#include "PCCTsinogram.h"
#include "TV.h"
#include "point.h"
#include "geometry.h"
#include "SparseMatrix.h"
#include "methods.h"
#include "materials.h"
#include "CalcSysmat.cuh"
#include "CalcMat.cuh"
#include "DenseMatrix.h"
#include "IterativeMatDec.h"




extern const int MAXMATERIALS_MD;// 1800000000;
extern const float MAX_LINTG;
extern const float MIN_LINTG;
extern const float MAX_TANH;
extern const float MIN_TANH;
extern const double PI;

namespace Reconstruction {

	//class square_matrix_error : public std::exception {
	//	virtual const char* what() const noexcept { return "Error: The matrix is not square."; }
	//};

	//class energy_resolution_error : public std::exception {
	//	virtual const char* what() const noexcept { return "Error: The energy resolutions of materials and spectrum are different."; }
	//};

	class IterativeMatDecCustom : public IterativeMatDec {
	protected:
		float prec;
		float* bnmean;
		float* xsource;
		float* precdiff, * precdiffsum;

	public:
		IterativeMatDecCustom(Reconstruction::PCCTsinogram* s, Reconstruction::geometry* geometry, Reconstruction::materials* m, std::vector<Reconstruction::spectrum>& ss, int size) : IterativeMatDec(s, geometry, m, ss, size) {

			bnmean = (float*)malloc(sizeof(float) * s->get_nb());
			precdiff = (float*)malloc(sizeof(float) * nd * nd);
			precdiffsum = (float*)malloc(sizeof(float) * nd * nd);

			for (int i = 0; i < s->get_nb(); i++) {
				bnmean[i] = 0;
			}

			for (int i = 0; i < nd * nd; i++) {
				precdiff[i] = 0;
				precdiffsum[i] = 0;
			}

			xsource = source_spectrum[0].get_data(); //sourceフォルダにはBINで分けないスペクトルを１ファイルだけ用意


			//for debug
			bnmean[0] = 0;
			bnmean[1] = 0;
			bnmean[2] = 0;
			bnmean[3] = 0;
			bnmean[4] = 0;
			bnmean[5] = 0;

		}

		void calc_imgdiff_gpu(float** smr, int v_begin) override;

	};
}
