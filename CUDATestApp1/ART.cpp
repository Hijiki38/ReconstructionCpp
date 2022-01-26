#include "ART.h"

using std::cout;
using namespace Eigen;

namespace Reconstruction {

	VectorXf* ART::reconstruction(int itr, int tvitr) {

		int x, y, v, w;
		//int origx, origy, tmpx; // relx, rely, tmpx;
		int center;
		float relx, rely;

		point point_abs(0, 0, 0);

		float theta, phi;
		float offset_detector, offset_detector_relative, intercept_Y;
		float center_relative_x, center_relative_y;
		Trip* material = NULL;

		center = n_detector / 2;
		center_relative_x = 0;
		center_relative_y = axis_correlation;

		cout << "\nStart Generationg System Matrix " << n_detector << " " << n_view << "\n";
		point_abs.set_center(center);
		for (y = 0; y < n_detector; y++) {
			cout << "\rGenerating System Matrix:" << y << " / " << n_detector << "\n";
			for (x = 0; x < n_detector; x++) {
				theta = 0;
				point_abs.set_xy(x, y);
				for (v = 0; v < n_view; v++) {
					phi = 0;
					if (theta >= PI / 4) { //“Š‰eŠp‚ª45“x‚ð’´‚¦‚½‚ç‰æ‘œ‚ð90“x‰E‰ñ“]‚³‚¹“Š‰eŠp‚ð - 45“x‚É
						theta -= PI / 2;
						point_abs.rotate90(true);
						relx = point_abs.get_relative(point_abs.get_x());
						rely = point_abs.get_relative(point_abs.get_inverted(point_abs.get_y()));
					}
					for (w = 0; w < n_detector; w++) {
						offset_detector = (n_detector - w - 1 - center + 0.5) / cos(theta); //offset of the detector
						//offset_detector_relative = offset_detector * (sod + relx) / sdd;
						//phi = atan2f(offset_detector, sdd);
						intercept_Y = offset_detector - rely + relx * tan(theta + phi);

						if (intercept_Y <= 0.5 * (1 - tan(theta)) and intercept_Y >= -0.5 * (1 - tan(theta))) {
							material = new Trip(static_cast<float>(n_detector * v + w), static_cast<float>(n_detector * y + x), 1);
							materials.push_back(*material);
						}
						else if (intercept_Y <= 0.5 * (1 + tan(theta)) and intercept_Y >= 0.5 * (1 - tan(theta))) {
							material = new Trip(static_cast<float>(n_detector * v + w), static_cast<float>(n_detector * y + x), 1);
							materials.push_back(*material);
						}

					}
					theta += 2 * PI / n_view;
				}
			}
		}
		cout << "End Generating System Matrix";

		sysmat.setFromTriplets(materials.begin(), materials.end());


		float prevsum, diff;
		int itrcount;
		float sys_atn;// , sys_sys;
		float* sys_sys = (float*)malloc(sizeof(float) * n_detector * n_view);
		VectorXf _sino = (*sino).get_sinovec();
		//VectorXf* _tmp = static_cast<VectorXf*>(malloc(sizeof(VectorXf) * n_detector * n_detector));
		VectorXf* _tmp = new VectorXf(n_detector * n_detector);

		TV* tv = new TV();
		TV* tvd = new TV();

		for (int i = 0; i < n_detector * n_view; i++) {
			*(_tmp) = sysmat.row(i);
			sys_sys[i] = VectorXf(*(_tmp)).dot(VectorXf(*(_tmp)));
			cout << "\rGenerating (sys * sys) :" << i << " / " << n_detector * n_view << "\n";
		}

		itrcount = 0;
		while (itrcount < itr) { //getchar() != '\n') {
			//Iteration(ART)
			//attenutmp = attenu;
			prevsum = attenu.sum();
			for (int i = 0; i < n_detector * n_view; i++) {
				*(_tmp) = sysmat.row(i);
				sys_atn = VectorXf(*(_tmp)).dot(attenu);
				attenu = attenu - relpar * ((sys_atn - _sino[i]) / sys_sys[i]) * *(_tmp);
			}
			diff = attenu.sum() - prevsum;
			cout << "\rIteration:" << itrcount << ", diff:" << diff << string(10, ' ');
			//cout << "\rIteration:" << attenu[0] << "," << attenu[1] << "," << attenu[2];



			for (int i = 0; i < tvitr; i++) {

				cout << "\rIteration:" << itrcount << ", TVitration:" << i;

				*(_tmp) = attenu;
				tv->calc_tv(&attenu, &imgdiff, 0.0001, true);

				attenu = attenu - imgdiff;
				tvd->calc_tv(&attenu, &imgdiff, 0.0001, true);

				if (tv < tvd) {
					attenu = *(_tmp);
				}

				tv = tvd;
			}

			itrcount++;
		}

		cout << "\n end iteration!";

		return &attenu;
	}

}

