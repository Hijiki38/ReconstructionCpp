#include "MLEM.h"

using std::cout;
using namespace Eigen;

namespace Reconstruction {


	float* MLEM::reconstruction(int itr, int tvitr)
	{
		int nd = (*sino).get_nd();
		int nv = (*sino).get_nv();
		int ne = (*sino).get_ne();
		float prevsum, diff;
		int itrcount;
		float sys_atn, sys_sys;
		float* syssys = (float*)malloc(sizeof(float) * nd * nv);
		float* sysrow_dev;
		float* syssys_dev;
		float* _sino = (*sino).get_sinovec();
		float* _sysmatblock;
		VectorXf* _tmp = new VectorXf(nd * nd);

		TV* tv = new TV();
		TV* tvd = new TV();

		bool block = true;

		generate_sysmat(geometry_normalized->is_conebeam);

		//cudaMalloc(&sysrow_dev, sizeof(float) * nd * nv);
		//cudaMalloc(&syssys_dev, sizeof(float) * nd * nv);

		//cout << "\n";
		//for (int i = 0; i < nd * nv; i++)
		//{
		//	*(_tmp) = sysmat.row(i);
		//	cout << "\rGenerating (sys * sys) :" << i << " / " << nd * nv;
		//	//cudaMemcpy(sysrow_dev, _tmp, sizeof(float) * nd * nv, cudaMemcpyHostToDevice);
		//	//cuda

		//	//syssys[i] = VectorXf(*(_tmp)).dot(VectorXf(*(_tmp)));
		//	//cout << "\rGenerating (sys * sys) :" << i << " / " << nd * nv;
		//}



		//// GPUで計算
		//gpu_function << <(N + 255) / 256, 256 >> > (dev_x, dev_y);

		//// GPU⇒CPUのデータコピー
		//cudaMemcpy(host_y, dev_y, N * sizeof(float), cudaMemcpyDeviceToHost);


		cout << "\n";
		//for (int i = 0; i < nd * nv; i++) 
		//{
		//	*(_tmp) = sysmat.row(i);
		//	syssys[i] = VectorXf(*(_tmp)).dot(VectorXf(*(_tmp)));
		//	//cout << "\rGenerating (sys * sys) :" << i << " / " << nd * nv;
		//}

		itrcount = 0;
		while (itrcount < itr)
		{
			//prevsum = attenu.sum();

			if (block == true) { //block-ART
				for (int i = 0; i < nv * nd / block_num; i++) {

					_sysmatblock = sysmat.Extract_row(sysmat, i * block_num, block_num);

					//_sysmatblock = new MatrixXf(sysmat.middleRows(i * block_num, block_num).eval());
					//imgdiff_art = std::vector<float>(nd * nd, 0);
					//imgdiff = (float*)malloc(nd * nd * sizeof(float));
					for (int j = 0; j < nd * nd * ne; j++) {
						imgdiff[j] = 0;
					}

					//imgdiff_art = Eigen::VectorXf::Zero(nd * nd);
					for (int j = 0; j < block_num; j++) {

						sys_atn = _sysmatblock.row(i * block_num + j).dot(attenu);
						sys_sys = _sysmatblock.row(i * block_num + j).dot(_sysmatblock->row(i * block_num + j));
						imgdiff = imgdiff + ((sys_atn - _sino[i * block_num + j]) / sys_sys) * _sysmatblock->row(j);

						//sys_atn = VectorXf(_sysmatblock->row(j)).dot(attenu);
						//sys_sys = VectorXf(_sysmatblock->row(j)).dot(VectorXf(_sysmatblock->row(j)));
						//imgdiff = imgdiff + ((sys_atn - _sino[i * block_num + j]) / sys_sys) * VectorXf(_sysmatblock->row(j));
					}
					attenu = attenu - (relpar / block_num) * imgdiff;
					delete _sysmatblock;
				}
				diff = attenu.sum() - prevsum;
				cout << "\rIteration:" << itrcount << ", diff:" << diff << string(10, ' ');
			}
			else {
				for (int i = 0; i < nd * nv; i++) // ray-by-ray ART
				{
					*(_tmp) = sysmat.row(i);
					sys_atn = VectorXf(*(_tmp)).dot(attenu);
					syssys[i] = VectorXf(*(_tmp)).dot(VectorXf(*(_tmp)));
					attenu = attenu - relpar * ((sys_atn - _sino[i]) / syssys[i]) * *(_tmp);
				}
				diff = attenu.sum() - prevsum;
				cout << "\rIteration:" << itrcount << ", diff:" << diff << string(10, ' ');
			}


			if (tvitr != -1) {
				for (int i = 0; i < tvitr; i++)
				{
					cout << "\nIteration:" << itrcount << ", TVitration:" << i;

					*(_tmp) = attenu;
					tv->calc_tv(&attenu, &imgdiff_tv, 0.0001, true);

					attenu = attenu - imgdiff_tv;
					tvd->calc_tv(&attenu, &imgdiff_tv, 0.0001, true);

					cout << "tv:" << tv << ", tvd:" << tvd << "\n";
					if (tv < tvd)
					{
						attenu = *(_tmp);
					}

					tv = tvd;
				}
			}


			itrcount++;
		}

		cout << "\n end iteration!";

		return &attenu;
	}


	void MLEM::generate_sysmat(bool is_conebeam) {
		int center;
		float relx = 0;
		float rely = 0;

		point point_abs(0, 0, 0);

		float theta, phi;
		float offset_detector, offset_detector_relative, intercept_Y;
		float center_relative_x, center_relative_y;
		float area = 0;
		Trip* material = NULL;

		center = nd / 2;
		center_relative_x = 0;
		center_relative_y = geometry_normalized->axis_correction;

		cout << "\nStart Generationg System Matrix " << nd << " " << nv << "\n";
		point_abs.set_center(center);
		for (int y = 0; y < nd; y++)
		{
			cout << "\rGenerating System Matrix:" << y << " / " << nd << "\n";
			for (int x = 0; x < nd; x++)
			{
				theta = 0;
				point_abs.set_xy(x, y);
				for (int v = 0; v < nv; v++)
				{
					phi = 0;
					if (theta >= PI / 4)
					{ //投影角が45度を超えたら画像を90度右回転させ投影角を - 45度に
						theta -= PI / 2;
						point_abs.rotate90();
						relx = point_abs.get_relative(point_abs.get_x());
						rely = point_abs.get_relative(point_abs.get_inverted(point_abs.get_y()));
					}
					for (int w = 0; w < nd; w++)
					{
						offset_detector = (nd - w - 1 - center + 0.5) / cos(theta); //offset of the detector
						if (is_conebeam)
						{
							//offset_detector_relative = offset_detector * (geometry_normalized->sod + relx) / geometry_normalized->sdd;
							phi = atan2f(offset_detector, geometry_normalized->sdd);
						}
						intercept_Y = offset_detector - rely + relx * tan(theta + phi);

						area = calc_area(intercept_Y, offset_detector, theta);
						if (area != 0) {
							material = new Trip(static_cast<float>(nd * v + w), static_cast<float>(nd * y + x), area);
							materials.push_back(*material);
						}



						//if (intercept_Y <= 0.5 * (1 - tan(theta)) and intercept_Y >= -0.5 * (1 - tan(theta))) 
						//{
						//	material = new Trip(static_cast<float>(nd * v + w), static_cast<float>(nd * y + x), 1);
						//	materials.push_back(*material);
						//}
						//else if (intercept_Y <= 0.5 * (1 + tan(theta)) and intercept_Y >= 0.5 * (1 - tan(theta))) 
						//{
						//	material = new Trip(static_cast<float>(nd * v + w), static_cast<float>(nd * y + x), 1);
						//	materials.push_back(*material);
						//}

					}
					theta += 2 * PI / nv;
				}
			}
		}
		cout << "End Generating System Matrix";

		sysmat.setFromTriplets(materials.begin(), materials.end());
	}

	float ART::calc_area(float intercept, float detectoroffset, float angle) {
		float la1, la2, lb1, lb2, sa, sb;
		float a = 0.5;

		la1 = a - (-a * tan(angle) + intercept + a / cos(angle));
		la2 = a - (a * tan(angle) + intercept + a / cos(angle));
		lb1 = a + (-a * tan(angle) + intercept - a / cos(angle));
		lb2 = a + (a * tan(angle) + intercept - a / cos(angle));



		sa = _calc_subarea(la1, la2, a);
		sb = _calc_subarea(lb1, lb2, a);

		return (2 * a) * (2 * a) - (sa + sb);
	}

	float ART::calc_area_cbct(float intercept, float detectoroffset, float angle) {
		float la1, la2, lb1, lb2, sa, sb;
		float a = 0.5;
		float tan_angle = tan(angle);
		float tan_delta = a / geometry_normalized->sdd;

		la1 = a - (-a * (tan_angle + tan_delta) + intercept + a / sqrt(1 / (1 + (tan_angle + tan_delta) * (tan_angle + tan_delta))));
		la2 = a - (a * (tan_angle + tan_delta) + intercept + a / sqrt(1 / (1 + (tan_angle + tan_delta) * (tan_angle + tan_delta))));
		lb1 = a + (-a * (tan_angle - tan_delta) + intercept - a / sqrt(1 / (1 + (tan_angle - tan_delta) * (tan_angle - tan_delta))));
		lb2 = a + (a * (tan_angle - tan_delta) + intercept - a / sqrt(1 / (1 + (tan_angle - tan_delta) * (tan_angle - tan_delta))));

		sa = _calc_subarea(la1, la2, a);
		sb = _calc_subarea(lb1, lb2, a);

		return (2 * a) * (2 * a) - (sa + sb);
	}

	float ART::_calc_subarea(float l1, float l2, float a) {

		if (l1 >= 2 * a) {
			if (l2 >= 2 * a) {
				return (2 * a) * (2 * a);
			}
			else {
				return a * (l1 + l2) - (l1 - 2 * a) * (l1 - 2 * a) / (2 * (l1 - l2));
			}
		}
		else if (l1 >= 0) {
			if (l2 >= 2 * a) {
				return a * (l1 + l2) - (l2 - 2 * a) * (l2 - 2 * a) / (2 * (l2 - l1));
			}
			else if (l2 >= 0) {
				return a * (l1 + l2);
			}
			else {
				return a * l1 * l1 / (l1 + abs(l2));
			}
		}
		else {
			if (l2 >= 0) {
				return a * l2 * l2 / (abs(l1) + l2);
			}
			else {
				return 0;
			}
		}
	}

	geometry* ART::get_geometry() {
		return geometry_normalized;
	}
}

