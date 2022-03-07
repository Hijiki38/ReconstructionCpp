#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "Eigen/Sparse"
#include "sinogram.h"
#include "TV.h"
#include "point.h"
#include "geometry.h"

//using namespace std;
//using namespace Eigen;

using Trip = Eigen::Triplet<float>;
using SpMat = Eigen::SparseMatrix<float>;

static const int MAXMATERIALS = 500000;
extern const double PI;

namespace Reconstruction {



	class ART {
	private:
		int n_detector, n_view;
		Reconstruction::sinogram* sino;

		SpMat sysmat;
		vector<Trip> materials;
		Eigen::VectorXf attenu, imgdiff_art, imgdiff_tv;

		//bool is_conebeam = false;

		float relpar = 1.05;
		float relpard = 0.01;

		int block_num;
		
		Reconstruction::geometry* geometry_normalized = new geometry();


	public:
		ART(Reconstruction::sinogram* s, Reconstruction::geometry* geometry) {

			sino = s;
			n_detector = (*s).get_nd();
			n_view = (*s).get_nv();

			std::cout << "header, nd:" << n_detector << " nv:" << n_view << "\n";

			geometry_normalized->is_conebeam = geometry->is_conebeam;

			geometry_normalized->pixelsize = geometry->pixelsize;
			geometry_normalized->sod = geometry->sod / geometry->pixelsize;
			geometry_normalized->sdd = geometry->sdd / geometry->pixelsize;
			geometry_normalized->axis_correction = geometry->axis_correction;

			block_num = n_detector;
			//block_num = 1;

			sysmat = *(new SpMat(n_view * n_detector, n_detector * n_detector));
			materials = *(new vector<Trip>(MAXMATERIALS));
			attenu = Eigen::VectorXf::Zero(n_detector * n_detector);
			imgdiff_tv = Eigen::VectorXf::Zero(n_detector * n_detector);
			imgdiff_art = Eigen::VectorXf::Zero(n_detector * n_detector);
		}

		Eigen::VectorXf* reconstruction(int itr, int tvitr);

		void generate_sysmat(bool is_conebeam);

		float calc_area(float intercept, float detectoroffset, float angle);
		float calc_area_cbct(float intercept, float detectoroffset, float angle);
		float _calc_subarea(float l1, float l2, float a);

		Reconstruction::geometry* get_geometry();

	};
}

