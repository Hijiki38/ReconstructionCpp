#pragma once
#include "Eigen/Sparse"
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>

namespace Reconstruction {
	class TV {
	private:
		float totalvariation;
		float xyn, xpyn, xny, xy, xpy, xnyp, xyp;
		float v, v1, v2, v3;

	public:
		TV(float tv) {
			totalvariation = tv;
		}
		TV() {
			totalvariation = 0;
		}
		void calc_tv(Eigen::VectorXf* img, Eigen::VectorXf* imgdiff, float param, bool updatetv);
		float get_tv();
		void set_tv(float tv);
	};
}
