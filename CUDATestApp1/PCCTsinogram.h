#pragma once

#include <stdio.h>
#include <iostream>
#include "Eigen/Sparse"

using namespace std;
using namespace Eigen;

namespace Reconstruction {
	class PCCTsinogram {
	private:
		float* sinovec;
		int n_d; //number of the detector
		int n_v; //number of view
		int n_e; //number of energy bin

	public:

		PCCTsinogram(float* s, int d, int v, int e) {
			sinovec = s;
			n_d  = d;
			n_v  = v;
			n_e  = e;
		}

		float* get_sinovec() { return sinovec; }
		int get_nd() { return n_d; }
		int get_nv() { return n_v; }
		int get_ne() { return n_e; }

		static PCCTsinogram* read_PCCTsinogram(std::vector<std::string> inpaths);

	};
}
