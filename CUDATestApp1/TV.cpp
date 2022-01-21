#include <math.h>
#include "Eigen/Sparse"

using namespace Eigen;

namespace Reconstruction{
	class TV {
	private:


	public:
		static float calc_tv(VectorXf* img, VectorXf* imgdiff, float param) {

			int imgsize = (*img).rows();
			float tv = 0.0;

			for (int i = 1; i < imgsize - 1; i++) {
				for (int j = 1; j < imgsize - 1; j++) {
					float xyn = (*img)(i - 1, j);
					float xpyn = (*img)(i - 1, j + 1);
					float xny = (*img)(i, j - 1);
					float xy = (*img)(i, j);
					float xpy = (*img)(i, j + 1);
					float xnyp = (*img)(i + 1, j - 1);
					float xyp = (*img)(i + 1, j);

					double v = sqrt(param + (xy - xny) * (xy - xny) + (xy - xyn) * (xy - xyn));
					double v1 = (2 * xy - xny - xyn) / v;
					double v2 = (xpy - xy) / sqrt(param + (xpy - xy) * (xpy - xy) + (xpy - xpyn) * (xpy - xpyn));
					double v3 = (xyp - xy) / sqrt(param + (xyp - xy) * (xyp - xy) + (xyp - xnyp) * (xyp - xnyp));
					(*imgdiff)(i * imgsize + j) = v1 - v2 - v3;
					tv += v;
				}
			}

			return tv;

		}

	};
}
