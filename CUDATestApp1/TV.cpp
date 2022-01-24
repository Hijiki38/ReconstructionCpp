#include "TV.h"

namespace Reconstruction{

	void TV::calc_tv(Eigen::VectorXf* img, Eigen::VectorXf* imgdiff, float param, bool updatetv) {

		int rowsize = static_cast<int>(sqrt((*img).rows()));
		float tv = 0.0;

		for (int i = 1; i < rowsize - 1; i++) {
			for (int j = 1; j < rowsize - 1; j++) {
				xyn = (*img)((i - 1) * rowsize + j);
				xpyn = (*img)((i - 1) * rowsize + j + 1);
				xny = (*img)(i * rowsize + j - 1);
				xy = (*img)(i * rowsize + j);
				xpy = (*img)(i * rowsize + j + 1);
				xnyp = (*img)((i + 1) * rowsize + j - 1);
				xyp = (*img)((i + 1) * rowsize + j);
				v = sqrt(param + (xy - xny) * (xy - xny) + (xy - xyn) * (xy - xyn));
				v1 = (2 * xy - xny - xyn) / v;
				v2 = (xpy - xy) / sqrt(param + (xpy - xy) * (xpy - xy) + (xpy - xpyn) * (xpy - xpyn));
				v3 = (xyp - xy) / sqrt(param + (xyp - xy) * (xyp - xy) + (xyp - xnyp) * (xyp - xnyp));

				//std::cout << "\n" << sqrt(param + (xpy - xy) * (xpy - xy) + (xpy - xpyn) * (xpy - xpyn));
				//std::cout << "\n" << sqrt(param + (xyp - xy) * (xyp - xy) + (xyp - xnyp) * (xyp - xnyp));

				//std::cout << "v1, v2, v3: " << v1 << "," << v2 << "," << v3 << "\n";

				//std::cout << i * rowsize + j << " / " << imgdiff->size() << "\n";
				(*imgdiff)(i * rowsize + j) = v1 - v2 - v3;

				tv += v;
			}
		}

		if (updatetv) {
			totalvariation = tv;
		}
	}

	float TV::get_tv() { return totalvariation; }

	void TV::set_tv(float tv) { totalvariation = tv; }

}
