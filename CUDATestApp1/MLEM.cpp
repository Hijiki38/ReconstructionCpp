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
		//float* syssys = (float*)malloc(sizeof(float) * nd * nv);
		//float* sysrow_dev;
		//float* syssys_dev;
		float* _sino = (*sino).get_sinovec();

		SparseMatrix _sysmatblock;
		float* _sysmatrow = (float*)malloc(sizeof(float) * nd * nd);

		TV* tv = new TV();
		TV* tvd = new TV();

		bool block = true;

		//for (int i = 0; i < nd * nd; i++) {
		//	_sysmatrow[i] = 0;
		//}

		sysmat = *(generate_sysmat());

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



		//// GPU‚ÅŒvZ
		//gpu_function << <(N + 255) / 256, 256 >> > (dev_x, dev_y);

		//// GPUËCPU‚Ìƒf[ƒ^ƒRƒs[
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
			prevsum = Reconstruction::sum_array(attenu);

			if (block == true) { //block-ART
				for (int i = 0; i < nv * nd / block_num; i++) {

					_sysmatblock = *(sysmat.Create_blockmat(i * block_num, block_num));

					//_sysmatblock = new MatrixXf(sysmat.middleRows(i * block_num, block_num).eval());
					//imgdiff_art = std::vector<float>(nd * nd, 0);
					//imgdiff = (float*)malloc(nd * nd * sizeof(float));
					for (int j = 0; j < nd * nd; j++) {
						imgdiff[j] = 0;
					}

					//imgdiff_art = Eigen::VectorXf::Zero(nd * nd);
					for (int j = 0; j < block_num; j++) {

						std::cout << "size:" << sizeof(_sysmatrow) << "nd: " << nd << std::endl;
						_sysmatblock.Extract_row_dense(j, nd * nd, _sysmatrow);

						sys_atn = Reconstruction::dot_array(_sysmatrow, attenu);
						sys_sys = Reconstruction::dot_array(_sysmatrow, _sysmatrow);
						imgdiff = Reconstruction::add_array(imgdiff, Reconstruction::mul_array(((sys_atn - _sino[i * block_num + j]) / sys_sys), _sysmatrow));

						//sys_atn = VectorXf(_sysmatblock->row(j)).dot(attenu);
						//sys_sys = VectorXf(_sysmatblock->row(j)).dot(VectorXf(_sysmatblock->row(j)));
						//imgdiff = imgdiff + ((sys_atn - _sino[i * block_num + j]) / sys_sys) * VectorXf(_sysmatblock->row(j));
					}
					//attenu = attenu - (relpar / block_num) * imgdiff;
					attenu = Reconstruction::sub_array(attenu, Reconstruction::mul_array((relpar / block_num), imgdiff));
					//delete _sysmatblock;
				}
				diff = Reconstruction::sum_array(attenu) - prevsum;
				cout << "\rIteration:" << itrcount << ", diff:" << diff << string(10, ' ');
			}

			itrcount++;
		}

		cout << "\n end iteration!";

		//free(_sysmatrow);

		return attenu;
	}


	std::unique_ptr<SparseMatrix> MLEM::generate_sysmat() {
		int center;
		float relx = 0;
		float rely = 0;

		int nd = sino->get_nd();
		int nv = sino->get_nv();

		int nonzerocount = 0;
		bool firstelem = true;

		point point_abs(0, 0, 0);

		int rotatecount = 0;

		float theta, phi;
		float offset_detector, offset_detector_relative, intercept_Y;
		float center_relative_x, center_relative_y;
		float area = 0;
	/*	Trip* material = NULL;*/

		float* elements = (float*)malloc(sizeof(float) * 5000000);	//all nonzero values
		int* rowptr = (int*)malloc(sizeof(int) * (nd * nv + 1));		//indices of the first nonzero element in each row
		int* colind = (int*)malloc(sizeof(int) * 5000000);		//the column indices of the corresponding elements

		int nonzero = 0;		//the number of nonzero elements

		center = nd / 2;
		center_relative_x = 0;
		center_relative_y = geometry_normalized->axis_correction;

		cout << "\nStart Generationg System Matrix " << nd << " " << nv << "\n";
		point_abs.set_center(center);

		theta = 0;
		for (int v = 0; v < nv; v++) 
		{
			cout << "\rGenerating System Matrix:" << y << " / " << nd << "\n";
			std::cout << "nonzero: " << nonzero << std::endl;
			std::cout << "colind: " << colind[nonzero - 1] << std::endl;
			std::cout << "elem: " << elements[nonzero - 1] << std::endl;
			std::cout << "rptr: " << rowptr[nd * y - 1] << std::endl;

			if (theta >= PI / 4)
			{ //“Š‰eŠp‚ª45“x‚ğ’´‚¦‚½‚ç‰æ‘œ‚ğ90“x‰E‰ñ“]‚³‚¹“Š‰eŠp‚ğ - 45“x‚É
				rotatecount++;
				theta -= PI / 2;
				point_abs.rotate90();
				relx = point_abs.get_relative(point_abs.get_x());
				rely = point_abs.get_relative(point_abs.get_inverted(point_abs.get_y()));
			}

			for (int w = 0; w < nd; w++) 
			{
				for (int y = 0; y < nd; y++)
				{

					for (int x = 0; x < nd; x++)
					{
						//theta = 0;
						point_abs.set_xy(x, y);
						for (int v = 0; v < nv; v++)
						{
							firstelem = true;
							phi = 0;
							if (theta >= PI / 4)
							{ //“Š‰eŠp‚ª45“x‚ğ’´‚¦‚½‚ç‰æ‘œ‚ğ90“x‰E‰ñ“]‚³‚¹“Š‰eŠp‚ğ - 45“x‚É
								theta -= PI / 2;
								point_abs.rotate90();
								relx = point_abs.get_relative(point_abs.get_x());
								rely = point_abs.get_relative(point_abs.get_inverted(point_abs.get_y()));
							}
							for (int w = 0; w < nd; w++)
							{
								offset_detector = (nd - w - 1 - center + 0.5) / cos(theta); //offset of the detector
								if (geometry_normalized->is_conebeam)
								{
									//offset_detector_relative = offset_detector * (geometry_normalized->sod + relx) / geometry_normalized->sdd;
									phi = atan2f(offset_detector, geometry_normalized->sdd);
								}
								intercept_Y = offset_detector - rely + relx * tan(theta + phi);

								area = calc_area(intercept_Y, offset_detector, theta);
								if (area != 0) {
									elements[nonzero] = area;
									colind[nonzero] = nd * y + x;
									if (firstelem) {
										rowptr[nd * v + w] = nonzero;
										firstelem = false;
									}

									//material = new Trip(static_cast<float>(nd * v + w), static_cast<float>(nd * y + x), area);
									//materials.push_back(*material);
									nonzero++;
								}
							}
							theta += 2 * PI / nv;
						}
					}
				}
			}
		}

		


		for (int y = 0; y < nd; y++)
		{
			cout << "\rGenerating System Matrix:" << y << " / " << nd << "\n";
			std::cout << "nonzero: " << nonzero << std::endl;
			std::cout << "colind: " << colind[nonzero - 1] << std::endl;
			std::cout << "elem: " << elements[nonzero - 1] << std::endl;
			std::cout << "rptr: " << rowptr[nd * y - 1] << std::endl;

			for (int x = 0; x < nd; x++)
			{
				theta = 0;
				point_abs.set_xy(x, y);
				for (int v = 0; v < nv; v++)
				{
					firstelem = true;
					phi = 0;
					if (theta >= PI / 4)
					{ //“Š‰eŠp‚ª45“x‚ğ’´‚¦‚½‚ç‰æ‘œ‚ğ90“x‰E‰ñ“]‚³‚¹“Š‰eŠp‚ğ - 45“x‚É
						theta -= PI / 2;
						point_abs.rotate90();
						relx = point_abs.get_relative(point_abs.get_x());
						rely = point_abs.get_relative(point_abs.get_inverted(point_abs.get_y())); //y²‚Í‹t‚È‚Ì‚Å”½“]
					}
					for (int w = 0; w < nd; w++)
					{
						offset_detector = (nd - w - 1 - center + 0.5) / cos(theta); //offset of the detector
						if (geometry_normalized->is_conebeam)
						{
							//offset_detector_relative = offset_detector * (geometry_normalized->sod + relx) / geometry_normalized->sdd;
							phi = atan2f(offset_detector, geometry_normalized->sdd);
						}
						intercept_Y = offset_detector - rely + relx * tan(theta + phi);

						area = calc_area(intercept_Y, offset_detector, theta);
						if (area != 0) {
							elements[nonzero] = area;
							colind[nonzero] = nd * y + x;
							if (firstelem) {
								rowptr[nd * v + w] = nonzero;
								firstelem = false;
							}

							//material = new Trip(static_cast<float>(nd * v + w), static_cast<float>(nd * y + x), area);
							//materials.push_back(*material);
							nonzero++;
						}
					}
					theta += 2 * PI / nv;
				}
			}
		}

		cout << "End Generating System Matrix";

		//sysmat.setFromTriplets(materials.begin(), materials.end());
		//sysmat = SparseMatrix(elements, rowptr, colind, nonzero);

		std::unique_ptr<SparseMatrix> sysmatptr(new SparseMatrix(elements, rowptr, colind, nonzero));
		
		return sysmatptr;
		//sysmat = *(sysmatptr);

	}

	float MLEM::calc_area(float intercept, float detectoroffset, float angle) {
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

	float MLEM::calc_area_cbct(float intercept, float detectoroffset, float angle) {
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

	float MLEM::_calc_subarea(float l1, float l2, float a) {

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

	geometry* MLEM::get_geometry() {
		return geometry_normalized;
	}
}

