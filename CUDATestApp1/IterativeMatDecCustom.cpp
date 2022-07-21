#include "IterativeMatDecCustom.h"

//using std::cout;
//using namespace Eigen;

namespace Reconstruction {

	//float** IterativeMatDecCustom::reconstruction(int itr)
	//{
	//	float prevsum, diff;
	//	int itrcount;
	//	float sys_atn, sys_sys;
	//	float* _sino = sino->get_sinovec();
	//	float* _sinobg = sino->get_sinobgvec();

	//	bool block = true;
	//	bool init = true;

	//	std::unique_ptr<SparseMatrix> _sysmatblock;

	//	float** sysmat1proj = (float**)malloc(sizeof(float*) * block_size);
	//	for (int i = 0; i < block_size; i++) {
	//		sysmat1proj[i] = (float*)malloc(sizeof(float) * nd * nd);
	//	}

	//	std::cout << "start gen sysmat" << std::endl;

	//	//sysmat = generate_sysmat(true, true);
	//	sysmat = generate_sysmat(true, false);


	//	std::cout << "preparing source spectrum..." << std::endl;
	//	for (int i = 0; i < nd * nv; i++) { //正規化された線源スペクトルにbgの光子数をかける
	//		for (int j = 0; j < nb; j++) {
	//			for (int k = 0; k < ne; k++) {
	//				source[i][j][k] *= _sinobg[i];
	//			}
	//		}
	//	}

	//	std::cout << "\nStart iteration";

	//	itrcount = 0;
	//	while (itrcount < itr)
	//	{
	//		//prevsum = std::accumulate(matfrac[0].begin(), matfrac[0].end(), 0);
	//		prevsum = Reconstruction::sum_array(matfrac[0], nd * nd);



	//		if (block == true) { //block-ART

	//			diff = 0;

	//			for (int j = 0; j < nd * nd; j++) {
	//				for (int k = 0; k < nm; k++) {
	//					matfractmp[j][k] = matfrac[j][k];
	//				}
	//			}

	//			for (int i = 0; i < nv * nd / block_size; i += 10) {
	//				//for (int i = 0; i < nv * nd / block_size; i++) {

	//					//for (int j = 0; j < nd * nd; j++) {
	//					//	for (int k = 0; k < nm; k++) {
	//					//		matfractmp[j][k] = matfrac[j][k];
	//					//	}
	//					//}




	//				_sysmatblock = (*sysmat).Create_blockmat(i * block_size, block_size);

	//				std::cout << "\ninitializing imgdiff...";
	//				init_matrix(imgdiff, nm, nd * nd, 0);
	//				std::cout << "finished" << std::endl;


	//				for (int j = 0; j < block_size; j++) {
	//					(*_sysmatblock).Extract_row_dense(j, nd * nd, sysmat1proj[j]);
	//				}

	//				std::cout << "### update image (" << i << "/" << nv * nd / block_size << ") ###" << std::endl;
	//				calc_imgdiff_gpu(sysmat1proj, i * block_size);

	//				for (int j = 0; j < nd * nd; j++) {
	//					for (int k = 0; k < nm; k++) {
	//						diff += imgdiff[j][k];

	//						imgdiffsum[j][k] += imgdiff[j][k] / nv;

	//						//matfrac[j][k] -= imgdiff[j][k];
	//						////if (matfrac[j][k] < 0 || (j / nd - nd / 2) * (j / nd - nd / 2) + (j % nd - nd / 2) * (j % nd - nd / 2) > (nd / 2) * (nd / 2)) matfrac[j][k] = 0;
	//						//if (matfrac[j][k] < 0) matfrac[j][k] = 0;
	//						////if (matfrac[j][k] < -10000) matfrac[j][k] = -10000;
	//						//if (matfrac[j][k] > 1000) matfrac[j][k] = 1000;
	//						//matfracprev[j][k] = matfractmp[j][k];

	//					}
	//				}

	//				int count;
	//				//for debug
	//				for (int m = 0; m < nm; m++) {

	//					std::string outfilename = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\material_decomposition\\output1";
	//					outfilename += mat->get_matlist()[m].name;
	//					outfilename += std::to_string(i);
	//					outfilename += ".csv";
	//					std::ofstream ofs(outfilename);

	//					count = 0;
	//					for (int k = 0; k < nd; k++)
	//					{
	//						for (int j = 0; j < nd; j++)
	//						{
	//							ofs << matfrac[count][m] - imgdiffsum[count][m] << ',';
	//							count++;
	//						}
	//						ofs << '\n';
	//					}
	//				}

	//				//if (i > 19) break;
	//				//break;
	//			}

	//			for (int j = 0; j < nd * nd; j++) {
	//				for (int k = 0; k < nm; k++) {
	//					matfrac[j][k] -= imgdiffsum[j][k];
	//					if (matfrac[j][k] < 0) matfrac[j][k] = 0;
	//					//if (matfrac[j][k] < -10000) matfrac[j][k] = -10000;
	//					if (matfrac[j][k] > 1000) matfrac[j][k] = 1000;
	//					matfracprev[j][k] = matfractmp[j][k];
	//				}
	//			}


	//			std::cout << "\nIteration:" << itrcount << ", diff:" << diff << std::string(10, ' ');
	//		}

	//		itrcount++;


	//	}

	//	std::cout << "\n end iteration!";

	//	for (int i = 0; i < block_size; i++) {
	//		free(sysmat1proj[i]);
	//	}
	//	free(sysmat1proj);

	//	//Reconstruction::devicereset();

	//	return matfrac;

	//}

	void IterativeMatDecCustom::calc_imgdiff_gpu(float** smr, int v_begin) {

		////frac: nd*nd X nm matrix,  mat_atn : nm X ne matrix,  source_spectrum : nb X ne (source[j])

		int nm = mat->get_matlist().size();
		int ne = source_spectrum[0].get_size();

		float w = 1;

		//

		//float** testmat = (float**)malloc(sizeof(float*) * 3);
		//for (int i = 0; i < 3; i++) {
		//	testmat[i] = (float*)malloc(sizeof(float) * 3);
		//}

		//testmat[0][0] = 0;
		//testmat[0][1] = 1;
		//testmat[0][2] = 1;
		//testmat[1][0] = 5;
		//testmat[1][1] = 1;
		//testmat[1][2] = 3;
		//testmat[2][0] = 2;
		//testmat[2][1] = 0;
		//testmat[2][2] = 1;

		//float** res = gaussian_elim(testmat, 3);

		//for (int i = 0; i < 3; i++) {
		//	for (int j = 0; j < 3; j++) {
		//		std::cout << res[i][j] << ",";
		//	}
		//	std::cout << std::endl;
		//}

		//return;

		float** bdist = (float**)malloc(sizeof(float*) * nb);
		for (int i = 0; i < nb; i++) {
			bdist[i] = (float*)malloc(sizeof(float) * ne);
		}


		float** nbar = (float**)malloc(sizeof(float*) * block_size); //nd*nv X nb
		float** nbard = (float**)malloc(sizeof(float*) * block_size); //nd*nv X nb
		float** lintg = (float**)malloc(sizeof(float*) * block_size); //nd*nv X ne
		for (int i = 0; i < block_size; i++) {
			nbar[i] = (float*)malloc(sizeof(float) * nb);
			nbard[i] = (float*)malloc(sizeof(float) * nb);
			lintg[i] = (float*)malloc(sizeof(float) * ne);
		}

		float* suma = (float*)malloc(sizeof(float) * block_size); //nd*nv
		float* eta = (float*)malloc(sizeof(float) * nm); //nm (0.00001)
		float* gamma = (float*)malloc(sizeof(float) * nm); //nm (0.00001)


		float*** grad_n = (float***)malloc(sizeof(float**) * block_size); //ni X nb X nm

		float grad_p;

		float** grad_q = (float**)malloc(sizeof(float*) * nd * nd); //nd*nd X nm
		float** grad_s = (float**)malloc(sizeof(float*) * nd * nd); //nd*nd X nm

		for (int i = 0; i < block_size; i++) {
			grad_n[i] = (float**)malloc(sizeof(float*) * nb);
			for (int j = 0; j < nb; j++) {
				grad_n[i][j] = (float*)malloc(sizeof(float) * nm);
			}
		}

		for (int i = 0; i < nd * nd; i++) {
			grad_q[i] = (float*)malloc(sizeof(float) * nm);
			grad_s[i] = (float*)malloc(sizeof(float) * nm);
		}

		float*** h = (float***)malloc(sizeof(float**) * nd * nd); //nd*nd X nm X nm
		float* h_qp = (float*)malloc(sizeof(float) * nm);
		float h_p;
		float*** hs = (float***)malloc(sizeof(float**) * nd * nd); //nd*nd X nm X nm
		float*** h_inv = (float***)malloc(sizeof(float**) * nd * nd); //nd*nd X nm X nm
		float*** hs_inv = (float***)malloc(sizeof(float**) * nd * nd); //nd*nd X nm X nm
		for (int i = 0; i < nd * nd; i++) {
			h[i] = (float**)malloc(sizeof(float*) * (nm + 1));
			h_inv[i] = (float**)malloc(sizeof(float*) * (nm + 1));
			hs[i] = (float**)malloc(sizeof(float*) * (nm + 1));
			hs_inv[i] = (float**)malloc(sizeof(float*) * (nm + 1));
			for (int j = 0; j < nm + 1; j++) {
				h[i][j] = (float*)malloc(sizeof(float) * (nm + 1));
				h_inv[i][j] = (float*)malloc(sizeof(float) * (nm + 1));
				hs[i][j] = (float*)malloc(sizeof(float) * (nm + 1));
				hs_inv[i][j] = (float*)malloc(sizeof(float) * (nm + 1));
				for (int k = 0; k < nm + 1; k++) {
					h[i][j][k] = 0;
					h_inv[i][j][k] = 0;
					hs[i][j][k] = 0;
					hs_inv[i][j][k] = 0;
				}
			}
		}

		for (int i = 0; i < nm; i++) {
			h_qp[i] = 0;

			eta[i] = 0.001;// 0.00001;
			gamma[i] = 0.01;
		}

		eta[0] = 0.001; //air
		eta[1] = 0.001; //al
		eta[2] = 0.001; //c
		eta[3] = 0.001; //ti
		//eta[4] = 0.001; //cu

		gamma[0] = 0.01;
		gamma[1] = 0.01;
		gamma[2] = 0.01;
		gamma[3] = 0.01;
		//gamma[4] = 0.01;

		float epsilon = 1;

		//eta[0] = 0.0001; //air
		//eta[1] = 0.0001; //al
		//eta[2] = 0.0001; //c
		//eta[3] = 0.0001; //ti
		//eta[4] = 0.0001; //cu

		//gamma[0] = 0.1;
		//gamma[1] = 0.1;
		//gamma[2] = 0.1;
		//gamma[3] = 0.1;
		//gamma[4] = 0.1;

		//eta[0] = 0.001; //air
		//eta[1] = 0.001; //al
		//eta[2] = 0.001; //c
		//eta[3] = 0.01; //ti

		//gamma[0] = 0.01;
		//gamma[1] = 0.01;
		//gamma[2] = 0.01;
		//gamma[3] = 0.001;



		//std::cout << "\npreparation... " << std::endl;
		//n = calc_n(source, block_size, nb, ne, v_begin);
		//for (int i = 0; i < block_size; i++) {
		//	for (int b = 0; b < nb; b++) {
		//		//n0[i][b] = sino->get_sinobgvec()[i * nb + nb];
		//		//n1[i][b] = sino->get_sinovec()[i * nb + nb];
		//		//for (int e = 0; e < ne; e++) {
		//		//	n[i][b] += source[v_begin + i][b][e];
		//		//}
		//	}
		//}

		//std::cout << "n[100][2] = " << n[100][2] << std::endl;

		float tmpf = 0;


		std::cout << "preparation... " << std::endl;
		calc_lintg(lintg, smr, matfrac, matatn, nd, block_size, ne, nm, pixsize);

		//for (int i = 0; i < block_size; i++) {
		//	std::cout << "\r" << i << "/" << block_size;
		//	for (int e = 0; e < ne; e++) {
		//		lintg[i][e] = 0;
		//		for (int j = 0; j < nd * nd; j++) {
		//			if (smr[i][j] != 0) {
		//				for (int k = 0; k < nm; k++) {
		//					lintg[i][e] += smr[i][j] * matfrac[j][k] * matatn[k][e];
		//				}
		//			}
		//		}
		//		//std::cout << lintg[i][e] << " at " << i << "," << e << std::endl;
		//	}
		//}

		//take negative and exp of lintg
		std::cout << "lintg[1][8] = " << lintg[1][8] << std::endl;


		for (int i = 0; i < block_size; i++) {
			for (int e = 0; e < ne; e++) {
				if (i % 50 == 0 && e % 100 == 0) { //lintg[i][e] != 0 &&
					std::cout << "lintg: " << lintg[i][e] << " at " << i << "," << e << std::endl;
				}
				lintg[i][e] = std::exp(-lintg[i][e]);
				if (lintg[i][e] > 0) {
					if (lintg[i][e] < MIN_LINTG) lintg[i][e] = MIN_LINTG;
					if (lintg[i][e] > MAX_LINTG) lintg[i][e] = MAX_LINTG;
				}
				else {
					if (lintg[i][e] < -MAX_LINTG) lintg[i][e] = -MAX_LINTG;
					if (lintg[i][e] > -MIN_LINTG) lintg[i][e] = -MIN_LINTG;
				}

			}
		}
		std::cout << "lintg[100][800] = " << lintg[100][800] << std::endl;

		std::cout << "preparation(2)... ";
		float evalue = 0;
		for (int b = 0; b < nb; b++) {
			for (int e = 0; e < ne; e++) {
				evalue = (e + 150) / 10;
				if (b == 0 && e < 15) {
					std::cout << "e, evalue: " << e << " " << evalue;
				}
				bdist[b][e] = std::exp(-0.5 * prec * (evalue - bnmean[b]));

				for (int i = 0; i < block_size; i++) {
					source[v_begin + i][b][e] = xsource[e] * bdist[b][e];
				}
			}
		}


		for (int i = 0; i < block_size; i++) {
			for (int b = 0; b < nb; b++) {
				nbar[i][b] = 0;
				for (int e = 0; e < ne; e++) {
					nbar[i][b] += source[v_begin + i][b][e] * lintg[i][e];
					
					//nbar[i][b] += xsource[e] * bdist[b][e] * lintg[i][e];

					//if (i % 10 == 0 && e % 100 == 0) {
					//	std::cout << "i,b,e, source, lintg = " << i << " " << b << " " << e << ", " << source[v_begin + i][b][e] << ", " << lintg[i][e] << std::endl;
					//}
				}
				//nbar[i][b] /= ne;
				//nbar[i][b] /= 10;

				//if (nbar[i][b] == 0) {
				//	std::cout << "zero at (i,b)= " << i << ", " << b << std::endl;
				//	for (int e = 0; e < ne; e++) {
				//		std::cout << "e, source, lintg = " << e << ", " << source[v_begin + i][b][e] << ", " << lintg[i][e] << std::endl;
				//	}
				//}
			}
		}
		//std::cout << "nbar[100][2] = " << nbar[100][2] << std::endl;

		std::cout << "preparation(3)... ";
		calc_nbard(nbard, source, lintg, bnmean, nd, block_size, nb, ne, nm, v_begin);

		std::cout << " completed." << std::endl;



		std::cout << "calc suma: ";
		for (int i = 0; i < block_size; i++) {
			suma[i] = 0;
			for (int j = 0; j < nd * nd; j++) {
				suma[i] += smr[i][j];
			}
			suma[i] *= pixsize;
		}

		std::cout << "calc grad n: ";
		// grad n
		for (int b = 0; b < nb; b++) {
			for (int k = 0; k < nm; k++) {
				for (int i = 0; i < block_size; i++) {
					grad_n[i][b][k] = 0;
					for (int e = 0; e < ne; e++) {
						//if(source[v_begin + i][b][e] * matatn[k][e] * lintg[i][e] != 0 && i == 10 && e % 100 == 0)	std::cout << "source, matatn, lintg, all" << source[v_begin + i][b][e] << " " << matatn[k][e] << " " << lintg[i][e] << " " << source[v_begin + i][b][e] * matatn[k][e] * lintg[i][e] << std::endl;

						grad_n[i][b][k] += source[v_begin + i][b][e] * pixsize * matatn[k][e] * lintg[i][e];

						//std::cout << grad_n[b][k] << std::endl;
					}
					//grad_n[i][b][k] /= ne;
					//grad_n[i][b][k] /= 10;
				}

			}
		}
		std::cout << "grad_n:" << std::endl;
		//for (int b = 0; b < nb; b++) {
		//	for (int k = 0; k < nm; k++) {
		//		std::cout << grad_n[5][b][k] << " ";
		//	}
		//	std::cout << std::endl;
		//}
		std::cout << "completed." << std::endl;

		//grad q
		std::cout << "calc grad q: ";
		//calc_gradq(grad_q, smr, sino->get_sinovec(), nbar, grad_n, nd, block_size, nb, nm, v_begin, pixsize);
		calc_gradq2(grad_q, smr, source, sino->get_sinovec(), nbar, matatn, lintg, suma, nd, block_size, nb, ne, nm, v_begin, pixsize);

		//void calc_gradq2(float** res, float*** source, float* n, float** nbar, float** matatn, float** lintg, float* suma, int nd, int ni, int nb, int ne, int nm, int v_begin, float pixsize);

		for (int j = 0; j < nd * nd; j++) {
			for (int k = 0; k < nm; k++) {
				if (grad_q[j][k] != 0 && j % 50000 == 0) {
					std::cout << "grad_q " << j << " " << k << ":" << grad_q[j][k] << " ";
					std::cout << std::endl;
				}

			}
		}
		std::cout << "completed." << std::endl;

		std::cout << "calc grad q: ";
		grad_p = calc_gradp(sino->get_sinovec(), nbar, nbard, nd, block_size, nb, ne, v_begin);

		// grad s
		//bool sfrag = false;

		std::cout << "calc grad s: ";
		for (int j = 0; j < nd * nd; j++) {
			//std::cout << "\r" << j << "/" << nd * nd;
			std::vector<int> neighbors = calc_neighbor(j, nd);
			for (int k = 0; k < nm; k++) {
				grad_s[j][k] = 0;
				for (int l : neighbors) {
					float tmp = std::tanh((2 * matfrac[j][k] - matfracprev[j][k] - matfracprev[l][k]) / gamma[k]);
					if (tmp > MAX_TANH) tmp = MAX_TANH;
					if (tmp < MIN_TANH) tmp = MIN_TANH;
					//grad_s[j][k] += eta[k] * (w / gamma[k]) * std::tanh((2 * matfrac[j][k] - matfracprev[j][k] - matfracprev[l][k]) / gamma[k]);
					grad_s[j][k] += eta[k] * (w / gamma[k]) * tmp;
					//if (!sfrag && std::abs((2 * matfrac[j][k] - matfracprev[j][k] - matfracprev[l][k]) / gamma[k]) < 10) {
					//	sfrag = true;
					//	//std::cout << "s OK. first:" << j << " " << k << std::endl;
					//}
				}
			}
		}
		for (int j = 0; j < nd * nd; j += 1) {
			for (int k = 0; k < nm; k++) {
				if (grad_s[j][k] != 0 && j % 50000 == 0) {
					std::cout << "grad_s " << j << " " << k << ":" << grad_s[j][k] << " ";
					std::cout << std::endl;
				}
			}
		}
		std::cout << "completed." << std::endl;


		//derivative of regularization term   eps 2 r
		for (int j = 0; j < nd * nd; j++) {
			grad_s[j][nm] += epsilon * 2 * prec;
		}





		//hessian q
		std::cout << "calc hessian q: ";
		calc_hesseq(h, source, smr, sino->get_sinovec(), nbar, matatn, lintg, suma, nd, block_size, nb, ne, nm, v_begin, pixsize);

		for (int j = 0; j < block_size; j++) {
			if (smr[j][0] != 0) 		std::cout << "nonzero smr[0]:" << smr[j][0] << " at " << j << std::endl;
		}


		std::cout << "h[0]" << std::endl;
		for (int j = 0; j < nm; j++) {
			for (int k = 0; k < nm; k++) {
				std::cout << h[0][j][k] << ", ";
			}
			std::cout << std::endl;
		}
		std::cout << "completed." << std::endl;

		//hessian qp
		calc_hesseqp(h_qp, source, nbar, nbard, matatn, lintg, suma, nd, block_size, nb, ne, nm, v_begin);


		//hessian p
		//h_p = calc_hessep(sino->get_sinovec(), nbar, nbard, block_size, nb);
		h_p = 0;
		float tmp;
		for (int i = 0; i < block_size; i++) {
			for (int b = 0; b < nb; b++) {
				tmp = nbard[i][b] / nbar[i][b];
				h_p += sino->get_sinovec()[i * nb + b] * tmp * tmp + (1 - sino->get_sinovec()[i * nb + b] / nbar[i][b]) * nbard[i][b] * 0.5;
			}
		}


		//merge hessian
		for (int j = 0; j < nd * nd; j++) {
			for (int k = 0; k < nm; k++) {
				h[j][nm][k] = h_qp[k];
				h[j][k][nm] = h_qp[k];
			}
			h[j][nm][nm] = h_p;
		}

		std::cout << "completed." << std::endl;


		//hessian s
		//sfrag = false;
		std::cout << "calc hessian s: ";
		for (int j = 0; j < nd * nd; j++) {
			std::vector<int> neighbors = calc_neighbor(j, nd);
			for (int k = 0; k < nm; k++) {
				//for (int m = 0; m < nm; m++) {
				//	h[j][k][m] = 0;
				//}
				for (int l : neighbors) {
					float tmp = std::tanh((2 * matfrac[j][k] - matfracprev[j][k] - matfracprev[l][k]) / gamma[k]);
					if (tmp > MAX_TANH) tmp = MAX_TANH;
					if (tmp < MIN_TANH) tmp = MIN_TANH;
					//h[j][k][k] += eta[k] * (2 * w / std::sqrt(gamma[k])) * (1 - std::pow(std::tanh((2 * matfrac[j][k] - matfracprev[j][k] - matfracprev[l][k]) / gamma[k]), 2));
					h[j][k][k] += eta[k] * (2 * w / std::sqrt(gamma[k])) * (1 - std::pow(tmp, 2));
					//hs[j][k][k] += eta[k] * (2 * w / std::sqrt(gamma[k])) * (1 - std::pow(tmp, 2));
					//if (!sfrag && std::abs((2 * matfrac[j][k] - matfracprev[j][k] - matfracprev[l][k]) / gamma[k]) < 10) {
					//	sfrag = true;
					//	std::cout << "s OK. first:" << j << " " << k << std::endl;
					//}
				}
			}
		}

		//hessian(2nd derivative) regularization p
		for (int j = 0; j < nd * nd; j++) {
			h[j][nm][nm] = epsilon * 0.5;
		}


		//for (int j = 0; j < nd * nd; j++) {
		//	for (int k = 0; k < nm; k++) {
		//		if (h[j][k][k] != 0) {
		//			std::cout << h[j][k][2] << " ";
		//			std::cout << std::endl;
		//		}

		//	}

		//}
		std::cout << "h[0](q & s)" << std::endl;
		for (int j = 0; j < nm; j++) {
			for (int k = 0; k < nm; k++) {
				std::cout << h[0][j][k] << ", ";
			}
			std::cout << std::endl;
		}
		std::cout << "completed." << std::endl;

		//calc inv of hessian
		std::cout << "invert hessian: ";
		std::vector<int> hesselist;
		for (int j = 0; j < nd * nd; j++) {
			//if (gaussian_elim(h[j], h_inv[j], nm) == 0) {
			//	hesselist.push_back(j);
			//}
			//gaussian_elim(hs[j], hs_inv[j], nm);
			gaussian_elim(h[j], h_inv[j], nm);

			//if(j < 100) std::cout << "inv:" << j << std::endl;
			//if (j > 261118) {
			//	std::cout << "gaussian elimination: " << j << std::endl;
			//	h_inv[j] = gaussian_elim(h[j], nm, true);
			//}
			//else {
			//	h_inv[j] = gaussian_elim(h[j], nm);
			//}

		}
		//for (int j = 0; j < nd * nd; j++) {
		//	for (int k = 0; k < nm; k++) {
		//		if (h_inv[j][k][k] != 0) {
		//			std::cout << j << ", " << k << ": " << h_inv[j][k][nm-1-k] << " ";
		//			std::cout << std::endl;
		//		}

		//	}
		//}
		//std::cout << "h[51456][2][2]: " << h[51456][2][2] << std::endl;



		//update material images(by Newton method)
		std::cout << "update images: ";

		for (int j = 0; j < nd * nd; j++) {
			for (int k = 0; k < nm; k++) {
				imgdiff[j][k] = 0;
			}
		}

		for (const auto& j : hesselist) {
			for (int k = 0; k < nm; k++) {
				for (int m = 0; m < nm; m++) {
					imgdiff[j][k] += grad_q[j][k] * h_inv[j][k][m];
				}
			}
		}

		for (int j = 0; j < nd * nd; j++) {
			for (int k = 0; k < nm; k++) {
				//imgdiff[j][k] = 0;
				for (int m = 0; m < nm+1; m++) {
					//imgdiff[j][k] += (grad_q[j][k] + grad_s[j][k]) * h_inv[j][k][m];
					imgdiff[j][k] += (grad_q[j][m] + grad_s[j][m]) * h_inv[j][k][m];
					////imgdiff[j][k] += grad_q[j][k] * h_inv[j][k][m] + grad_s[j][k] * hs_inv[j][k][m];
					//imgdiff[j][k] += grad_s[j][k] * hs_inv[j][k][m];
					////if((grad_q[j][k] + grad_s[j][k]) != 0) 
					if (j < 10) std::cout << "j,k,m, imgdiff, grad, hinv:" << j << ", " << k << ", " << m << ", " << imgdiff[j][k] << ", " << grad_q[j][k] + grad_s[j][k] << ", " << h_inv[j][k][m] << std::endl;
				}
			}
			for (int m = 0; m < nm + 1; m++) {
				precdiff[j] += (grad_q[j][m] + grad_s[j][m]) * h_inv[j][nm][m];
			}
		}

		//for (int j = 0; j < nd * nd; j += 1) {
		//	for (int k = 0; k < nm; k++) {
		//		if (imgdiff[j][k] != 0) {
		//			std::cout << j << ", " << k << ": " << imgdiff[j][k] << " ";
		//			std::cout << std::endl;
		//		}
		//	}
		//}
		std::cout << "completed." << std::endl;



		free(suma);
		free(eta);
		free(gamma);

		for (int i = 0; i < block_size; i++) {
			free(nbar[i]);
			free(lintg[i]);
		}

		free(nbar);
		free(lintg);



		for (int i = 0; i < block_size; i++) {
			for (int j = 0; j < nb; j++) {
				free(grad_n[i][j]);
			}
			free(grad_n[i]);
		}
		free(grad_n);
		for (int i = 0; i < nd * nd; i++) {
			free(grad_q[i]);
			free(grad_s[i]);
		}
		free(grad_q);
		free(grad_s);

		for (int i = 0; i < nd * nd; i++) {
			for (int j = 0; j < nm; j++) {
				free(h[i][j]);
				free(h_inv[i][j]);
				free(hs[i][j]);
				free(hs_inv[i][j]);
			}
			free(h[i]);
			free(h_inv[i]);
			free(hs[i]);
			free(hs_inv[i]);
		}
		free(h);
		free(h_inv);
		free(hs);
		free(hs_inv);

	}
}

