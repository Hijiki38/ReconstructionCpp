#pragma once

#include <stdio.h>
#include <iostream>
#include "Eigen/Sparse"

//using namespace std;
using namespace Eigen;

namespace Reconstruction {
	class PCCTsinogram {
	private:
		float* sinovec;
		float* sinobgvec;
		int n_d; //number of the detector
		int n_v; //number of view
		int n_b; //number of energy bin

	public:

		//PCCTsinogram(float* s, int d, int v, int b) {
		//	sinovec = s;
		//	n_d  = d;
		//	n_v  = v;
		//	n_b  = b;
		//}

		PCCTsinogram(float* obj, float* bg, int d, int v, int b) {
			sinovec = obj;
			sinobgvec = bg;
			n_d = d;
			n_v = v;
			n_b = b;
		}

		float* get_sinovec() { return sinovec; }
		float* get_sinobgvec() { return sinobgvec; }
		int get_nd() { return n_d; }
		int get_nv() { return n_v; }
		int get_nb() { return n_b; }

		static PCCTsinogram* read_PCCTsinogram(std::vector<std::string> inpaths_obj, std::vector<std::string> inpaths_bg);

		void convert_negativelog();

	};
}
