#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <math.h>
#include "PCCTsinogram.h"
#include "TV.h"
#include "point.h"
#include "geometry.h"
#include "SparseMatrix.h"
#include "methods.h"
#include "CalcSysmat.cuh"


static const int MAXMATERIALS = 450000000;
extern const double PI;

namespace Reconstruction {



	class IterationRec {
	protected:
		Reconstruction::PCCTsinogram* sino;

		std::unique_ptr<SparseMatrix> sysmat;
		//std::vector<float*> attenu, imgdiff;   //material images
		float* attenu, * imgdiff;

		float relpar;// = 1.05;
		//float relpard = 0.01;

		int block_size;

		

		Reconstruction::geometry* geometry_normalized = new geometry();


	public:
		IterationRec(Reconstruction::PCCTsinogram* s, Reconstruction::geometry* geometry, float par, int size) {

			relpar = par;

			sino = s;

			int n_detector = (*s).get_nd();
			int n_view = (*s).get_nv();
			int n_bin = (*s).get_nb();  

			std::cout << "header, nd:" << n_detector << " nv:" << n_view << "\n";

			geometry_normalized->is_conebeam = geometry->is_conebeam;

			geometry_normalized->pixelsize = geometry->pixelsize;
			geometry_normalized->sod = geometry->sod / geometry->pixelsize;
			geometry_normalized->sdd = geometry->sdd / geometry->pixelsize;
			geometry_normalized->axiscor_pixels = geometry->axiscor_pixels;

			block_size = n_detector * size;
			//block_num = 1;

			//sysmat = *(new SparseMatrix());
			attenu = (float*)malloc(n_detector * n_detector * n_bin * sizeof(float));
			imgdiff = (float*)malloc(n_detector * n_detector * n_bin * sizeof(float));

			for (int i = 0; i < n_detector * n_detector * n_bin; i++) {
				attenu[i] = 1;
				imgdiff[i] = 0;
			}
		}

		float* reconstruction(int itr, bool usegpu);

		std::unique_ptr<SparseMatrix> generate_sysmat_gpu(int begin, int size, bool init);
		std::unique_ptr<SparseMatrix> generate_sysmat(bool use_gpu = false, bool write_sysmat = false);

		virtual void calc_imgdiff(float* idiff, float* smr, float* atn, float sn, int size) const;
		virtual void calc_attenu(float* atn, float* idiff, int size) const = 0;

		float calc_area(float intercept, float detectoroffset, float angle);
		float calc_area_cbct(float intercept, float detectoroffset, float angle);
		float _calc_subarea(float l1, float l2, float a);

		Reconstruction::geometry* get_geometry();

	};
}

