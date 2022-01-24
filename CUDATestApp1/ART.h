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
		Eigen::VectorXf attenu, imgdiff;
		float relpar = 1.05;
		float relpard = 0.01;
		float sod, sdd;
		float axis_correlation;
		float pixelsize;


	public:
		ART(Reconstruction::sinogram* s, float sod, float sdd, float axis_correlation, float pixelsize) {

			sino = s;
			n_detector = (*s).get_nd();
			n_view = (*s).get_nv();

			std::cout << "header, nd:" << n_detector << " nv:" << n_view << "\n";

			this->pixelsize = pixelsize;
			this->sod = sod / pixelsize;
			this->sdd = sdd / pixelsize;
			this->axis_correlation = axis_correlation;


			sysmat = *(new SpMat(n_view * n_detector, n_detector * n_detector));
			materials = *(new vector<Trip>(MAXMATERIALS));
			attenu = Eigen::VectorXf::Zero(n_detector * n_detector);
			imgdiff = Eigen::VectorXf::Zero(n_detector * n_detector);
		}

		Eigen::VectorXf* reconstruction(int itr, int tvitr);

	};
}

