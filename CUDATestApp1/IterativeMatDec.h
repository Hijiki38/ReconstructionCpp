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




static const int MAXMATERIALS_MD = 45000000;
static const float MAX_LINTG = 1e+15;
static const float MIN_LINTG = 1e-10;
static const float MAX_TANH = 0.9999999;
static const float MIN_TANH = 0.0000001;
extern const double PI;

namespace Reconstruction {

	class square_matrix_error : public std::exception {
		virtual const char* what() const noexcept { return "Error: The matrix is not square."; }
	};

	class energy_resolution_error : public std::exception {
		virtual const char* what() const noexcept { return "Error: The energy resolutions of materials and spectrum are different."; }
	};

	class IterativeMatDec {
	protected:
		Reconstruction::PCCTsinogram* sino;

		std::unique_ptr<SparseMatrix> sysmat;
		//DenseMatrix matfrac, matfracprev, imgdiff, matatn, * source;
		float **matfrac, **matfracprev, **imgdiff, **matatn, ***source;


		//float relpar;// = 1.05;
		//float relpard = 0.01;

		int block_size;

		float pixsize;

		Reconstruction::geometry* geometry_normalized = new geometry();
		Reconstruction::materials* mat;
		std::vector<Reconstruction::spectrum> source_spectrum;

		int nd, nv, nb, ne, nm;


	public:
		IterativeMatDec(Reconstruction::PCCTsinogram* s, Reconstruction::geometry* geometry, Reconstruction::materials* m, std::vector<Reconstruction::spectrum>& ss, int size) {

			//relpar = par;
			sino = s;
			mat = m;
			source_spectrum = ss;

			//std::cout << "spectrum energy size: " << ss[0].get_size() << ", matlist energy size:" << m->get_matlist()[0].attenu.get_size();
			if (ss[0].get_size() != m->get_matlist()[0].attenu.get_size()) throw energy_resolution_error();

			nd = (*s).get_nd();
			nv = (*s).get_nv();
			nb = (*s).get_nb();
			nm = (*m).get_matlist().size();
			ne = ss[0].get_size(); //‚Æ‚è‚ ‚¦‚¸0.1keV(1500)


			//float** _matfrac, ** _matfracprev, ** _imgdiff, ** _matatn, *** _source; 

			std::cout << "header, nd:" << nd << " nv:" << nv << "\n";

			geometry_normalized->is_conebeam = geometry->is_conebeam;

			geometry_normalized->pixelsize = geometry->pixelsize;
			geometry_normalized->sod = geometry->sod / geometry->pixelsize;
			geometry_normalized->sdd = geometry->sdd / geometry->pixelsize;
			geometry_normalized->axiscor_pixels = geometry->axiscor_pixels;

			pixsize = geometry->pixelsize;

			block_size = nd * size;


			matfrac = (float**)malloc(nd * nd * sizeof(float*));
			matfracprev = (float**)malloc(nd * nd * sizeof(float*));
			imgdiff = (float**)malloc(nd * nd * sizeof(float*));
			matatn = (float**)malloc(nm * sizeof(float*));
			source = (float***)malloc(nd * nv * sizeof(float**)); // ndxnv X nb X ne matrix

			//_matfrac = (float**)malloc(nd*nd * sizeof(float*));
			//_matfracprev = (float**)malloc(nd * nd * sizeof(float*));
			//_imgdiff = (float**)malloc(nd*nd * sizeof(float*));
			//_matatn = (float**)malloc(nm * sizeof(float*));
			//_source = (float***)malloc(nd * nv * sizeof(float**)); // ndxnv X nb X ne matrix


			//std::cout << "nb:" << nb << std::endl;

			for (int i = 0; i < nd * nv; i++) {
				source[i] = (float**)malloc(sizeof(float*) * nb);
				for (int b = 0; b < nb; b++) {
					source[i][b] = (float*)malloc(sizeof(float) * ne);
					float* _tmp = source_spectrum[b].get_data();
					for (int e = 0; e < ne; e++) {
						//std::cout << "i, b, e:" << i << ", " << b << ", " << e << std::endl;
						source[i][b][e] = _tmp[e] * sino->get_sinobgvec()[b * nd * nv + i];
					}
				}
			}

			for (int i = 0; i < nd * nd; i++) {
				matfrac[i] = (float*)malloc(nm * sizeof(float));
				matfracprev[i] = (float*)malloc(nm * sizeof(float));
				imgdiff[i] = (float*)malloc(nm * sizeof(float));
				for (int j = 0; j < nm; j++) {
					matfrac[i][j] = 1;// 1 / nm;
					matfracprev[i][j] = 1;// 1 / nm;
					imgdiff[i][j] = 0;
				}
			}

			for (int i = 0; i < nm; i++) {
				matatn[i] = (float*)malloc(ne * sizeof(float));
				float* _tmp = mat->get_matlist()[i].attenu.get_data();
				for (int j = 0; j < ne; j++) {
					matatn[i][j] = _tmp[j];
				}
			}

			//matfrac = DenseMatrix(_matfrac, nd * nd, nm);
			//matfracprev = DenseMatrix(_matfracprev, nd * nd, nm);
			//imgdiff = DenseMatrix(_imgdiff, nd * nd, nm);
			//matatn = DenseMatrix(_matatn, nm, ne);
			//source = (DenseMatrix*)malloc(sizeof(DenseMatrix) * nd * nv);
			//for (int i = 0; i < nd * nv; i++) {
			//	//std::cout << i << std::endl;
			//	source[i] = DenseMatrix(_source[i], nb, ne);
			//}

		}

		float** reconstruction(int itr);

		std::unique_ptr<SparseMatrix> generate_sysmat_gpu(int begin, int size, bool init);
		std::unique_ptr<SparseMatrix> generate_sysmat(bool use_gpu = false, bool write_sysmat = false);

		std::vector<int>& calc_neighbor(int j, int size);
		virtual void calc_imgdiff_gpu(float** smr, int v_begin);
		virtual void calc_imgdiff(float** smr, int v_begin);
		//virtual void calc_attenu(float** atn, float** idiff, int size) const = 0;

		//float** gaussian_elim(float** input, int size, bool log = false);
		int gaussian_elim(float** input, float** output, int size, bool log = false);
		std::vector<std::vector<float>>& gaussian_elim(std::vector<std::vector<float>>& input);

		float calc_area(float intercept, float detectoroffset, float angle);
		float calc_area_cbct(float intercept, float detectoroffset, float angle);
		float _calc_subarea(float l1, float l2, float a);

		Reconstruction::geometry* get_geometry();

	};
}

