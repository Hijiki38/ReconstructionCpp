#pragma once

#include <stdio.h>
#include <iostream>
#include <string>
#include "Eigen/Sparse"

//using namespace std;
using namespace Eigen;

namespace Reconstruction {
	class sinogram {
	private:
		VectorXf eigen_sinovec;
		int n_d;
		int n_v;

	public:

		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		sinogram(VectorXf* svec, int d, int v) {
			eigen_sinovec = *svec;
			n_d = d;
			n_v = v;
		}

		sinogram(std::vector<float>* svec, int d, int v) {
			VectorXf tmp = Map<VectorXf, Unaligned>(svec->data(), svec->size());
			eigen_sinovec = tmp; // Map<VectorXf, Unaligned>(svec.data(), svec.size());
			n_d = d;
			n_v = v;
		}

		sinogram() {
			//eigen_sinovec = *(new VectorXf());
			n_d = 1;
			n_v = 1;
		}

		VectorXf get_sinovec();
		int get_nd();
		int get_nv();


		static sinogram* read_sinogram(std::string inpath);

	};
}
