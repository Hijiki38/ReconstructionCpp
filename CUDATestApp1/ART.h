#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "Eigen/Sparse"
//#include "sinogram.h"

using namespace std;
using namespace Eigen;


namespace Reconstruction {


typedef Eigen::SparseMatrix<float> SpMat;
typedef Eigen::Triplet<float> Trip;

static const int MAXMATERIALS = 500000;
extern const double PI;

	class ART {
	private:
		int n_detector, n_view;
		Reconstruction::sinogram* sino;

		SpMat sysmat;
		vector<Trip> materials;
		VectorXf attenu, imgdiff;
		float relpar = 1.05;


	public:
		ART(sinogram* s) {

			sino = s;
			n_detector = (*s).get_nd();
			n_view = (*s).get_nv();

			cout << "header, nd:" << n_detector << " nv:" << n_view << "\n";

			sysmat = *(new SpMat(n_view * n_detector, n_detector * n_detector));
			materials = *(new vector<Trip>(MAXMATERIALS));
			attenu = VectorXf::Zero(n_detector * n_detector);
			imgdiff = VectorXf::Zero(n_detector * n_detector);
		}

		VectorXf* reconstruction();

	};
}

