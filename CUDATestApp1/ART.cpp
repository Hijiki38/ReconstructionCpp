#include "ART.h"

using std::cout;
using namespace Eigen;

namespace Reconstruction {

	VectorXf* ART::reconstruction() {

		int x, y, v, w;
		int origx, origy, tmpx; // relx, rely, tmpx;
		int center;
		float relx, rely;
		float theta;
		float offsetD, offsetY;
		Trip* material = NULL;

		center = n_detector / 2;


		cout << "\nStart Generationg System Matrix " << n_detector << " " << n_view << "\n";
		for (y = 0; y < n_detector; y++) {
			cout << "\rGenerating System Matrix:" << y << " / " << n_detector << "\n";
			for (x = 0; x < n_detector; x++) {
				theta = 0;
				origx = x;
				origy = y;
				relx = origx - center + 0.5;
				rely = n_detector - origy - 1 - center + 0.5;
				for (v = 0; v < n_view; v++) {
					if (theta >= PI / 4) { //“Š‰eŠp‚ª45“x‚ð’´‚¦‚½‚ç‰æ‘œ‚ð90“x‰E‰ñ“]‚³‚¹“Š‰eŠp‚ð - 45“x‚É
						theta -= PI / 2;
						tmpx = origx;
						origx = n_detector - origy - 1;
						origy = tmpx;
						relx = origx - center + 0.5;
						rely = n_detector - origy - 1 - center + 0.5;
					}
					for (w = 0; w < n_detector; w++) {
						offsetD = (n_detector - w - 1 - center + 0.5) / cos(theta); //offset of the detector
						offsetY = offsetD - rely + relx * tan(theta);


						if (offsetY <= 0.5 * (1 - tan(theta)) and offsetY >= -0.5 * (1 - tan(theta))) {
							material = new Trip(static_cast<float>(n_detector * v + w), static_cast<float>(n_detector * y + x), 1);
							materials.push_back(*material);
						}
						else if (offsetY <= 0.5 * (1 + tan(theta)) and offsetY >= 0.5 * (1 - tan(theta))) {
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

		float tv = 0;
		float tvd = 0;

		for (int i = 0; i < n_detector * n_view; i++) {
			*(_tmp) = sysmat.row(i);
			sys_sys[i] = VectorXf(*(_tmp)).dot(VectorXf(*(_tmp)));
			cout << "\nsyssys" << sys_sys[i];
		}

		itrcount = 0;
		while (itrcount < 5) { //getchar() != '\n') {
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

			//Iteration(TV)

			for (int i = 0; i < 3; i++) {

				cout << "\rIteration:" << itrcount << ", TVitration:" << i;

				*(_tmp) = attenu;
				tv = Reconstruction::TV.calc_tv(attenu, imgdiff, 0.0001); //VectorXf* img, VectorXf* imgdiff, float param

				attenu = attenu - imgdiff;
				tvd = Reconstruction::TV.calc_tv(attenu, imgdiff, 0.0001);

				if (tv < tvd) {
					attenu = *(_tmp);
				}

				tv = tvd;
			}



			//for (int i = 0; 1 < 5; 1++) {
			//	Buffer.BlockCopy(img, 0, temp, 0, width * height);
			//	for (int y = 0; y < height; ++y) {
			//		for (int x = 0; x < width; ++x) img[x, y] -= D * diff[x, y];
			//	}
			//	double cvd = CalcTv(img, diff, E);
			//	if (vd < cvd) {
			//		Buffer.BlockCopy(temp, 0, img, 0, width * height);
			//		break;
			//	}
			//	vd = cvd;
			//}
			/*return vd;*/

			itrcount++;
		}

		cout << "\n end iteration!";

		return &attenu;
	}


	/*
	class ART {
	private:
		int n_detector, n_view, center;
		float sinogram;

		SpMat sysmat;
		vector<Trip> materials;
		VectorXf attenu, attenutmp, proj;
		float relpar;


	public:
		ART(float* sino, float nd, float nv) {
			sinogram = *sino;
			n_detector = nd;
			n_view = nv;

			relpar = 1.05;

			if (n_detector % 2 == 0) {
				center = n_detector / 2;
			}
			else {
				cerr << "WARNING: The number of detectors is uneven. This might leads unexpected results.";
				center = (n_detector - 1) / 2;
			}

			sysmat = *(new SpMat(n_view * n_detector, n_detector * n_detector));
			materials = *(new vector<Trip>(MAXMATERIALS));
			attenu = *(new VectorXf(n_detector * n_detector));
			proj = *(new VectorXf(n_view * n_detector));

		}

		VectorXf *reconstruction() {

			float diff;
			int itrcount;

			int x, y, v, w, i;
			int origx, origy, relx, rely, tmpx;
			int center;
			float theta;
			float offsetD, offsetY;

			Trip *material = NULL;
			float sys_atn, sys_sys;
			MatrixXf output;



			for (y = 0; y < n_detector; y++) {
				for (x = 0; x < n_detector; x++) {
					theta = 0;
					origx = x;
					origy = y;
					relx = origx - center + 0.5;
					rely = n_detector - origy - 1 - center + 0.5;
					for (v = 0; v < n_view; v++) {
						if (theta >= PI / 4) { //“Š‰eŠp‚ª45“x‚ð’´‚¦‚½‚ç‰æ‘œ‚ð90“x‰E‰ñ“]‚³‚¹“Š‰eŠp‚ð - 45“x‚É
							theta -= PI / 2;
							tmpx = origx;
							origx = n_detector - origy - 1;
							origy = tmpx;
							relx = origx - center + 0.5;
							rely = n_detector - origy - 1 - center + 0.5;
						}
						for (w = 0; w < n_detector; w++) {
							offsetD = (n_detector - w - 1 - center + 0.5) / cos(theta); //offset of the detector
								offsetY = offsetD - rely + relx * tan(theta);
								if (offsetY <= 0.5 * (1 - tan(theta)) and offsetY >= -0.5 * (1 - tan(theta))) {
									material = new Trip(static_cast<float>(n_detector * v + w), static_cast<float>(n_detector * y + x), 1);
									materials.push_back(*material);
								}
								else if (offsetY <= 0.5 * (1 + tan(theta)) and offsetY >= 0.5 * (1 - tan(theta))) {
									material = new Trip(static_cast<float>(n_detector * v + w), static_cast<float>(n_detector * y + x), 1);
									materials.push_back(*material);
								}

						}
						theta += 2 * PI / n_view;
					}
				}
			}

			sysmat.setFromTriplets(materials.begin(), materials.end());

			while (getchar() != '\n') {
				attenutmp = attenu;
				for (i = 0; i < n_detector * n_view; i++) {
					sys_atn = VectorXf(sysmat.row(i)).dot(attenu);
					sys_sys = VectorXf(sysmat.row(i)).dot(VectorXf(sysmat.row(i)));
					attenu = attenu - relpar * ((sys_atn - proj[i]) / sys_sys) * sysmat.row(i);
				}
				diff = attenu.sum() - attenutmp.sum();
				cout << "\rIteration:" << itrcount << ", diff:" << diff << string(10, ' ');
				itrcount++;
			}

			return &attenu;
			}

		}

		//float* genSysmat() {

		//}

		//float* iteration() {

		//}
	};*/

}

