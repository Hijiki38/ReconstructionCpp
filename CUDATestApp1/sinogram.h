#pragma once

#include <stdio.h>
#include <iostream>
#include <string>
#include "Eigen/Sparse"

using namespace std;
using namespace Eigen;

class sinogram {
private:
	VectorXf eigen_sinovec;
	int n_d;
	int n_v;

public:

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	sinogram(VectorXf svec, int d, int v) {
		eigen_sinovec = svec;
		n_d = d;
		n_v = v;
	}

	sinogram(vector<float> svec, int d, int v) {
		VectorXf tmp = Map<VectorXf, Unaligned>(svec.data(), svec.size());
		eigen_sinovec = tmp; // Map<VectorXf, Unaligned>(svec.data(), svec.size());
		n_d = d;
		n_v = v;
	}

	sinogram() {
		//eigen_sinovec = *(new VectorXf());
		n_d = 1;
		n_v = 1;
	}

	VectorXf get_sinovec() { return eigen_sinovec; }
	int get_nd() { return n_d; }
	int get_nv() { return n_v; }

};