#include "IterationRec.h"

//using std::cout;
//using namespace Eigen;

namespace Reconstruction {

	float* IterationRec::reconstruction(int itr, int tvitr)
	{
		int nd = (*sino).get_nd();
		int nv = (*sino).get_nv();
		int ne = (*sino).get_ne();
		float prevsum, diff;
		int itrcount;
		float sys_atn, sys_sys;
		float* _sino = (*sino).get_sinovec();

		std::unique_ptr<SparseMatrix> _sysmatblock;
		//std::shared_ptr<SparseMatrix> sysmathoge;
		float* _sysmatrow = (float*)malloc(sizeof(float) * nd * nd);

		bool block = true;
		bool init = true;

		//itrcount = 0;  //blockÇ≤Ç∆Ç…sysmatÇ¬Ç≠ÇÈÉîÉ@Å[ÉWÉáÉì
		//while (itrcount < itr)
		//{
		//	prevsum = Reconstruction::sum_array(attenu, nd * nd);

		//	if (block == true) { //block-ART
		//		for (int i = 0; i < nv * nd / block_size; i++) {

		//			std::cout << "gen sysmat:" << i << "/" << nv * nd / block_size << std::endl;
		//			sysmat = move(generate_sysmat_gpu(i,1,init));
		//			//std::unique_ptr<SparseMatrix> sysmathoge = generate_sysmat_gpu(i, 1, init);


		//			if (init) { init = false; }

		//			for (int j = 0; j < nd * nd; j++) {
		//				imgdiff[j] = 0;
		//			}

		//			for (int j = 0; j < block_size; j++) {
		//				sysmat->Extract_row_dense(j, nd * nd, _sysmatrow);
		//				calc_imgdiff(imgdiff, _sysmatrow, attenu, _sino[i * block_size + j], nd * nd);
		//			}

		//			calc_attenu(attenu, imgdiff, nd * nd);
		//		}

		//		diff = Reconstruction::sum_array(attenu, nd * nd) - prevsum;
		//		std::cout << "\rIteration:" << itrcount << ", diff:" << diff << string(10, ' ');

		//	}

		//	itrcount++;
		//}

		//Reconstruction::devicereset();

		//return attenu;

		sysmat = generate_sysmat(false, false);
		//sysmat = *hoge;

		std::cout << "\nStart iteration";

		itrcount = 0;
		while (itrcount < itr)
		{
			prevsum = Reconstruction::sum_array(attenu, nd * nd);

			if (block == true) { //block-ART
				for (int i = 0; i < nv * nd / block_size; i++) {

					_sysmatblock = (*sysmat).Create_blockmat(i * block_size, block_size);
					//_sysmatblock = *(sysmat.Create_blockmat(i * block_size, block_size));

					for (int j = 0; j < nd * nd; j++) {
						imgdiff[j] = 0;
					}

					//imgdiff_art = Eigen::VectorXf::Zero(nd * nd);
					for (int j = 0; j < block_size; j++) {

						(*_sysmatblock).Extract_row_dense(j, nd * nd, _sysmatrow);
						calc_imgdiff(imgdiff, _sysmatrow, attenu, _sino[i * block_size + j], nd * nd);
						//imgdiff = imgdiff + ((sys_atn - _sino[i * block_num + j]) / sys_sys) * _sysmatrow;
					}

					calc_attenu(attenu, imgdiff, nd * nd);
					//attenu = attenu - (relpar / block_num) * imgdiff;



				}
				diff = Reconstruction::sum_array(attenu, nd * nd) - prevsum;
				std::cout << "\rIteration:" << itrcount << ", diff:" << diff << string(10, ' ');
			}

			itrcount++;


		}

		std::cout << "\n end iteration!";

		free(_sysmatrow);

		Reconstruction::devicereset();

		return attenu;
	}

	void IterationRec::calc_imgdiff(float* idiff, float* smr, float* atn, float sn, int size) const {
		float sys_atn = Reconstruction::dot_array(smr, atn, size);
		float sys_sys = Reconstruction::dot_array(smr, smr, size);

		Reconstruction::mul_array1(smr, ((sys_atn - sn) / sys_sys), size);
		Reconstruction::add_array(idiff, smr, size);
	}

	std::unique_ptr<SparseMatrix> IterationRec::generate_sysmat_gpu(int begin, int size, bool init) {

		int nd = sino->get_nd();
		int nv = sino->get_nv();
		int center = nd / 2;

		//float* elements = (float*)malloc(sizeof(float) * MAXMATERIALS);	//all nonzero values
		//int* rowptr = (int*)malloc(sizeof(int) * (nd * size + 1));		//indices of the first nonzero element in each row
		//int* colind = (int*)malloc(sizeof(int) * MAXMATERIALS);		//the column indices of the corresponding elements

		std::unique_ptr<float[]> elements = std::make_unique<float[]>(MAXMATERIALS);
		std::unique_ptr<int[]> rowptr = std::make_unique<int[]>((nd * size + 1));	
		std::unique_ptr<int[]> colind = std::make_unique<int[]>(MAXMATERIALS);	

		int nonzero = 0;		//the number of nonzero elements

		std::cout << "\nStart Generating System Matrix(GPU) " << nd << " " << nv << "\n";
		nonzero = Reconstruction::calc_sysmat2(elements.get(), rowptr.get(), colind.get(), begin, size, nv, nd, center, geometry_normalized->sdd);
		rowptr[nd * size] = nonzero;
		std::cout << "\nSystem matrix generated!(GPU), nonzero = " << nonzero << "\n";

		std::unique_ptr<SparseMatrix> sysmatptr(new SparseMatrix(elements, rowptr, colind, nonzero));
		return sysmatptr;

	}

	std::unique_ptr<SparseMatrix> IterationRec::generate_sysmat(bool use_gpu, bool write_sysmat) {
		int center;
		float relx = 0;
		float rely = 0;

		int nd = sino->get_nd();
		int nv = sino->get_nv();

		bool firstelem = true;

		point point_abs(0, 0, 0);

		int rotatecount = 0;

		float theta, phi;
		float offset_detector, offset_detector_relative, intercept_Y;
		float center_relative_x, center_relative_y;
		float area = 0;

		//float* elements = (float*)malloc(sizeof(float) * MAXMATERIALS);	//all nonzero values
		//int* rowptr = (int*)malloc(sizeof(int) * (nd * nv + 1));		//indices of the first nonzero element in each row
		//int* colind = (int*)malloc(sizeof(int) * MAXMATERIALS);		//the column indices of the corresponding elements

		std::unique_ptr<float[]> elements = std::make_unique<float[]>(MAXMATERIALS);
		std::unique_ptr<int[]> rowptr = std::make_unique<int[]>((nd * nv + 1));
		std::unique_ptr<int[]> colind = std::make_unique<int[]>(MAXMATERIALS);


		int nonzero = 0;		//the number of nonzero elements


		ofstream ofs("C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sysmat.csv"); //for debug


		center = nd / 2;
		center_relative_x = 0;
		center_relative_y = geometry_normalized->axiscor_pixels;


		if (use_gpu) {
			std::cout << "\nStart Generating System Matrix(GPU) " << nd << " " << nv << "\n";
			//nonzero = Reconstruction::calc_sysmat2(elements, rowptr, colind, 0, nv, nd, center, geometry_normalized->sdd);
			nonzero = Reconstruction::calc_sysmat(elements.get(), rowptr.get(), colind.get(), nv, nd, center, geometry_normalized->sdd);
			rowptr[nd * nv] = nonzero;
			std::cout << "\nSystem matrix generated!(GPU), nonzero = " << nonzero << "\n";

			std::unique_ptr<SparseMatrix> sysmatptr(new SparseMatrix(elements, rowptr, colind, nonzero));
			return sysmatptr;
		}
		else {
			point_abs.set_center(center);
			theta = 0;

			for (int v = 0; v < nv; v++)
			{
				std::cout << "\rGenerating System Matrix(CPU):" << v << " / " << nv;

				if (theta >= PI / 4)
				{ //ìäâeäpÇ™45ìxÇí¥Ç¶ÇΩÇÁâÊëúÇ90ìxâEâÒì]Ç≥ÇπìäâeäpÇ - 45ìxÇ…
					rotatecount++;
					theta -= PI / 2;
				}

				for (int w = 0; w < nd; w++)
				{
					firstelem = true;
					for (int y = 0; y < nd; y++)
					{
						for (int x = 0; x < nd; x++)
						{
							//theta = 0;
							point_abs.set_xy(x, y);
							for (int i = 0; i < rotatecount; i++) {
								point_abs.rotate90();
							}
							//point_abs.rotate90(rotatecount);
							relx = point_abs.get_relative(point_abs.get_x());
							rely = point_abs.get_relative(point_abs.get_inverted(point_abs.get_y()));

							offset_detector = (nd - w - 1 - center + 0.5) / cos(theta); //offset of the detector
							if (geometry_normalized->is_conebeam)
							{
								//offset_detector_relative = offset_detector * (geometry_normalized->sod + relx) / geometry_normalized->sdd;
								phi = atan2f(offset_detector, geometry_normalized->sdd);
								intercept_Y = offset_detector - rely + relx * tan(theta + phi);
								area = calc_area_cbct(intercept_Y, offset_detector, theta);
							}
							else {
								phi = 0;
								intercept_Y = offset_detector - rely + relx * tan(theta + phi);
								area = calc_area(intercept_Y, offset_detector, theta);
							}

							if (area != 0) {
								if (nonzero == MAXMATERIALS - 1) {
									std::cout << "System matrix is too big!";
									exit(1);
								}
								elements[nonzero] = area;
								colind[nonzero] = nd * y + x;
								if (firstelem) {
									rowptr[nd * v + w] = nonzero;
									firstelem = false;
								}
								nonzero++;
							}

							if (write_sysmat) {
								ofs << area << ", ";
							}
						}
					}

					if (write_sysmat) {
						ofs << "\n";
					}
				}

				theta += 2 * PI / nv;
			}

			rowptr[nd * nv] = nonzero;

			std::cout << "End Generating System Matrix";

			std::unique_ptr<SparseMatrix> sysmatptr(new SparseMatrix(elements, rowptr, colind, nonzero));
			return sysmatptr;
		}
	}

	float IterationRec::calc_area(float intercept, float detectoroffset, float angle) {
		float la1, la2, lb1, lb2, sa, sb;
		float a = 0.5; //pixelsize / 2

		la1 = a - (-a * tan(angle) + intercept + a / cos(angle));
		la2 = a - (a * tan(angle) + intercept + a / cos(angle));
		lb1 = a + (-a * tan(angle) + intercept - a / cos(angle));
		lb2 = a + (a * tan(angle) + intercept - a / cos(angle));

		sa = _calc_subarea(la1, la2, a);
		sb = _calc_subarea(lb1, lb2, a);

		return (2 * a) * (2 * a) - (sa + sb);
	}

	float IterationRec::calc_area_cbct(float intercept, float detectoroffset, float angle) {
		float la1, la2, lb1, lb2, sa, sb;
		float a = 0.5; //pixelsize / 2
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

	float IterationRec::_calc_subarea(float l1, float l2, float a) {

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

	geometry* IterationRec::get_geometry() {
		return geometry_normalized;
	}
}

