#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "PCCTsinogram.h"
#include "TV.h"
#include "point.h"
#include "geometry.h"
#include "SparseMatrix.h"
#include "methods.h"


static const int MAXMATERIALS = 500000;
extern const double PI;

namespace Reconstruction {



	class MLEM {
	private:
		//int n_detector, n_view, n_bin;
		Reconstruction::PCCTsinogram* sino;

		SparseMatrix sysmat;
		float *attenu, *imgdiff, *imgdiff_tv;
		float* img_result;
		float* proj_original;

		float relpar = 1.05;
		float relpard = 0.01;

		int block_num;

		Reconstruction::geometry* geometry_normalized = new geometry();


	public:
		MLEM(Reconstruction::PCCTsinogram* s, Reconstruction::geometry* geometry) {

			sino = s;

			int n_detector = (*s).get_nd();
			int n_view = (*s).get_nv();
			int n_bin = 1;// (*s).get_ne();  //‚Æ‚è‚ ‚¦‚¸single bin

			std::cout << "header, nd:" << n_detector << " nv:" << n_view << "\n";

			geometry_normalized->is_conebeam = geometry->is_conebeam;

			geometry_normalized->pixelsize = geometry->pixelsize;
			geometry_normalized->sod = geometry->sod / geometry->pixelsize;
			geometry_normalized->sdd = geometry->sdd / geometry->pixelsize;
			geometry_normalized->axis_correction = geometry->axis_correction;

			block_num = n_detector;
			//block_num = 1;

			sysmat = *(new SparseMatrix());
			attenu = (float*)malloc(n_detector * n_detector * n_bin * sizeof(float));
			imgdiff_tv = (float*)malloc(n_detector * n_detector * n_bin * sizeof(float));
			imgdiff = (float*)malloc(n_detector * n_detector * n_bin * sizeof(float));

			img_result = (float*)malloc((unsigned long)n_detector * n_detector * sizeof(float));
			proj_original = (float*)malloc((unsigned long)n_view * n_detector * sizeof(float));
		}

		float* reconstruction(int itr, int tvitr = -1);

		void generate_sysmat(bool is_conebeam);
		//void forwardprojection(float* proj, float* img, bool is_conebeam);
		//void backprojection(float* proj, float* img, bool is_conebeam);

		float calc_area(float intercept, float detectoroffset, float angle);
		float calc_area_cbct(float intercept, float detectoroffset, float angle);
		float _calc_subarea(float l1, float l2, float a);

		Reconstruction::geometry* get_geometry();

	};
}

