#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "Eigen/Sparse"
#include "sinogram.h"

using namespace std;
using namespace Eigen;

typedef Eigen::SparseMatrix<float> SpMat;
typedef Eigen::Triplet<float> Trip;

static const int MAXMATERIALS = 500000;
extern const double PI;

class ART {
private:
	int n_detector, n_view;
	sinogram* sino;

	SpMat sysmat;
	vector<Trip> materials;
	VectorXf attenu, attenutmp; // , proj;
	float relpar = 1.05;


public:
	/*ART(vector<float> sino, float nd, float nv) {
		sinogram = Map<VectorXf, Unaligned>(sino.data(), sino.size());
		n_detector = nd;
		n_view = nv;

		sysmat = *(new SpMat(n_view * n_detector, n_detector * n_detector));
		materials = *(new vector<Trip>(MAXMATERIALS));
		attenu = *(new VectorXf(n_detector * n_detector));
	}*/


	ART(sinogram* s) {

		sino = s;
		//n_detector = (*sino).get_nd();
		//n_view = (*sino).get_nv();
		n_detector = (*s).get_nd();
		n_view = (*s).get_nv();

		cout << "header, nd:" << n_detector << " nv:" << n_view << "\n";

		sysmat = *(new SpMat(n_view * n_detector, n_detector * n_detector));
		materials = *(new vector<Trip>(MAXMATERIALS));
		//attenu = *(new VectorXf(n_detector * n_detector));
		attenu = VectorXf::Zero(n_detector * n_detector);
	}

	VectorXf* reconstruction();

};