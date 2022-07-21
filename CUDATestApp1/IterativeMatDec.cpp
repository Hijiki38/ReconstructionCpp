#include "IterativeMatDec.h"

//using std::cout;
//using namespace Eigen;

namespace Reconstruction {

	float** IterativeMatDec::reconstruction(int itr)
	{
		float prevsum, diff;
		int itrcount;
		float sys_atn, sys_sys;
		float* _sino = sino->get_sinovec();
		float* _sinobg = sino->get_sinobgvec();

		bool block = true;
		bool init = true;

		std::unique_ptr<SparseMatrix> _sysmatblock;

		float** sysmat1proj = (float**)malloc(sizeof(float*) * block_size);
		for (int i = 0; i < block_size; i++) {
			sysmat1proj[i] = (float*)malloc(sizeof(float) * nd * nd);
		}

		//test
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

		//std::cout << std::endl;
		//for (int i = 0; i < 3; i++) {
		//	for (int j = 0; j < 3; j++) {
		//		std::cout << res[i][j] << ",";
		//	}
		//	std::cout << std::endl;
		//}

		//exit(0);
		//test

		//for initializing image
		//std::string initimgname = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\input\\fourmetals_simulation\\proj\\initialimg2.csv";
		//std::ofstream iifs(initimgname);

		//std::ifstream stream(initimgname);
		//std::string line;

		//int counttmp = 0;

		//while (getline(stream, line)) {
		//	std::vector<std::string> strs = Reconstruction::splitstring(line, ',');
		//	for (int i = 0; i < strs.size(); i++) {
		//		for (int j = 0; j < nm; j++) {
		//			matfrac[counttmp % (nd * nd)][j] = stof(strs[i]);
		//		}	
		//		//std::cout << "j k; " << counttmp % (nd * nd) << "  " << counttmp / (nd * nd) << std::endl;
		//		counttmp++;
		//	}
		//}



		std::cout << "start gen sysmat" << std::endl;

		//sysmat = generate_sysmat(true, true);
		sysmat = generate_sysmat(true, false);


		std::cout << "preparing source spectrum..." << std::endl;
		for (int i = 0; i < nd * nv; i++) { //³‹K‰»‚³‚ê‚½üŒ¹ƒXƒyƒNƒgƒ‹‚Ébg‚ÌŒõŽq”‚ð‚©‚¯‚é
			for (int j = 0; j < nb; j++) {
				for (int k = 0; k < ne; k++) {
					source[i][j][k] *= _sinobg[i];
				}
			}
		}

		std::cout << "\nStart iteration";

		bool eachproj = false; //update image for each projection

		itrcount = 0;
		while (itrcount < itr)
		{
			//prevsum = std::accumulate(matfrac[0].begin(), matfrac[0].end(), 0);
			prevsum = Reconstruction::sum_array(matfrac[0], nd * nd);

			

			if (block == true) { //block-ART

				diff = 0;

				if (!eachproj) {
					for (int j = 0; j < nd * nd; j++) {
						for (int k = 0; k < nm; k++) {
							matfractmp[j][k] = matfrac[j][k];
						}
					}
				}


				for (int i = 0; i < nv * nd / block_size; i+=10) {
				//for (int i = 0; i < nv * nd / block_size; i++) {

					//each proj ver
					if (eachproj) {
						for (int j = 0; j < nd * nd; j++) {
							for (int k = 0; k < nm; k++) {
								matfractmp[j][k] = matfrac[j][k];
							}
						}
					}


					_sysmatblock = (*sysmat).Create_blockmat(i * block_size, block_size);
						
					std::cout << "\ninitializing imgdiff...";
					init_matrix(imgdiff, nm, nd*nd, 0);
					std::cout << "finished" << std::endl;


					for (int j = 0; j < block_size; j++) {
						(*_sysmatblock).Extract_row_dense(j, nd * nd, sysmat1proj[j]);
					}

					std::cout << "### update image (" << i << "/" << nv * nd / block_size << ") ###" << std::endl;
					calc_imgdiff_gpu(sysmat1proj, i * block_size);

					for (int j = 0; j < nd*nd; j++) {
						for (int k = 0; k < nm; k++) {
							diff += imgdiff[j][k];

							//all proj ver
							if (eachproj) {
								matfrac[j][k] -= imgdiff[j][k] / nv;
								//if (matfrac[j][k] < 0 || (j / nd - nd / 2) * (j / nd - nd / 2) + (j % nd - nd / 2) * (j % nd - nd / 2) > (nd / 2) * (nd / 2)) matfrac[j][k] = 0;
								if (matfrac[j][k] < 0) matfrac[j][k] = 0;
								//if (matfrac[j][k] < -10000) matfrac[j][k] = -10000;
								if (matfrac[j][k] > 1000) matfrac[j][k] = 1000;
								matfracprev[j][k] = matfractmp[j][k];
							}
							else {
								imgdiffsum[j][k] += imgdiff[j][k] / nv;
							}




						}	
					}

					int count;
					//for debug
					for (int m = 0; m < nm; m++) {

						std::string outfilename = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\material_decomposition\\output1";
						outfilename += mat->get_matlist()[m].name;
						outfilename += std::to_string(i);
						outfilename += ".csv";
						std::ofstream ofs(outfilename);

						count = 0;
						for (int k = 0; k < nd; k++)
						{
							for (int j = 0; j < nd; j++)
							{
								if (eachproj) {
									ofs << matfrac[count][m] << ',';
								}
								else {
									ofs << matfrac[count][m] - imgdiffsum[count][m] << ',';
								}
								
								count++;
							}
							ofs << '\n';
						}
					}

					//if (i > 19) break;
					//break;
				}

				float disx, disy;

				if (eachproj) {

				}
				else {
					for (int j = 0; j < nd * nd; j++) {
						for (int k = 0; k < nm; k++) {
							disx = (j / nd) - (nd / 2);
							disy = (j % nd) - (nd / 2);
							if (disx * disx + disy * disy <= (nd / 2) * (nd / 2)) {
								matfrac[j][k] -= imgdiffsum[j][k];
							}
							//matfrac[j][k] -= imgdiffsum[j][k];
							if (matfrac[j][k] < 0) matfrac[j][k] = 0;
							//if (matfrac[j][k] < -10000) matfrac[j][k] = -10000;
							if (matfrac[j][k] > 10000) matfrac[j][k] = 10000;
							matfracprev[j][k] = matfractmp[j][k];
						}
					}
				}



				std::cout << "\nIteration:" << itrcount << ", diff:" << diff << std::string(10, ' ');
			}

			itrcount++;


		}

		std::cout << "\n end iteration!";

		for (int i = 0; i < block_size; i++) {
			free(sysmat1proj[i]);
		}
		free(sysmat1proj);

		//Reconstruction::devicereset();

		return matfrac;

	}

	std::vector<int>& IterativeMatDec::calc_neighbor(int j, int size) {
		std::vector<int> result;

		if (j >= size) {
			result.push_back(j - size);
		}
		if (j + size < size * size) {
			result.push_back(j + size);
		}
		if (j % size != 0) {
			result.push_back(j - 1);
		}
		if (j % size != size - 1) {
			result.push_back(j + 1);
		}

		return result;
	}


	void IterativeMatDec::calc_imgdiff_gpu(float** smr, int v_begin) {

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


		float** nbar	= (float**)malloc(sizeof(float*) * block_size); //nd*nv X nb
		float** lintg   = (float**)malloc(sizeof(float*) * block_size); //nd*nv X ne
		for (int i = 0; i < block_size; i++) {
			nbar[i] = (float*)malloc(sizeof(float) * nb);
			lintg[i] = (float*)malloc(sizeof(float) * ne);
		}

		float* suma		= (float*)malloc(sizeof(float) * block_size); //nd*nv
		float* eta		= (float*)malloc(sizeof(float) * nm); //nm (0.00001)
		float* gamma	= (float*)malloc(sizeof(float) * nm); //nm (0.00001)


		float*** grad_n = (float***)malloc(sizeof(float**) * block_size); //ni X nb X nm

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
		float*** hs = (float***)malloc(sizeof(float**) * nd * nd); //nd*nd X nm X nm
		float*** h_inv = (float***)malloc(sizeof(float**) * nd * nd); //nd*nd X nm X nm
		float*** hs_inv = (float***)malloc(sizeof(float**) * nd * nd); //nd*nd X nm X nm
		for (int i = 0; i < nd * nd; i++) {
			h[i] = (float**)malloc(sizeof(float*) * nm);
			h_inv[i] = (float**)malloc(sizeof(float*) * nm);
			hs[i] = (float**)malloc(sizeof(float*) * nm);
			hs_inv[i] = (float**)malloc(sizeof(float*) * nm);
			for (int j = 0; j < nm; j++) {
				h[i][j] = (float*)malloc(sizeof(float) * nm);
				h_inv[i][j] = (float*)malloc(sizeof(float) * nm);
				hs[i][j] = (float*)malloc(sizeof(float) * nm);
				hs_inv[i][j] = (float*)malloc(sizeof(float) * nm);
				for (int k = 0; k < nm; k++) {
					h[i][j][k] = 0;
					h_inv[i][j][k] = 0;
					hs[i][j][k] = 0;
					hs_inv[i][j][k] = 0;
				}
			}
		}

		for (int i = 0; i < nm; i++) {
			eta[i] = 0.001;// 0.00001;
			gamma[i] = 0.01;
		}

		//eta[0] = 0.01; //air
		//eta[1] = 0.001; //al
		//eta[2] = 1; //c
		//eta[3] = 1; //ti
		////eta[4] = 0.001; //cu

		//gamma[0] = 0.001;
		//gamma[1] = 0.01;
		//gamma[2] = 0.00001;
		//gamma[3] = 0.00001;
		////gamma[4] = 0.01;

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

		eta[0] = 0.001; //air
		eta[1] = 0.001; //al
		eta[2] = 0.001; //c
		eta[3] = 0.01; //ti

		gamma[0] = 0.01;
		gamma[1] = 0.01;
		gamma[2] = 0.01;
		gamma[3] = 0.001;



		std::cout << "\npreparation... " << std::endl;
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


		std::cout << "preparation(2)... " << std::endl;
		calc_lintg(lintg, smr, matfrac, matatn, nd, block_size, ne, nm, pixsize);

		//for (int i = 0; i < block_size; i++) {
		//	std::cout << "\r" << i << "/" << block_size;
		//	for (int e = 0; e < ne; e++) {
		//		lintg[i][e] = 0;
		//		for (int j = 0; j < nd * nd; j++) {
		//			if (smr[i][j] != 0) {
		//				for (int k = 0; k < nm; k++) {
		//					lintg[i][e] += smr[i][j] * matfrac[j][k] * matatn[k][e] * pixsize;
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

		std::cout << "preparation(3)... ";
		for (int i = 0; i < block_size; i++) {
			for (int b = 0; b < nb; b++) {
				nbar[i][b] = 0;
				for (int e = 0; e < ne; e++) {
					nbar[i][b] += source[v_begin + i][b][e] * lintg[i][e];
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

		//h[j][k][m] += smr[i][j] * (n[i][b] / nbar[i][b]) * suma[i] * source[v_begin + i][b][e] * matatn[k][e] * matatn[m][e] * lintg[i][e];

		//for (int j = 0; j < 2; j++) {
		//	for (int k = 0; k < nm; k++) {
		//		for (int m = 0; m < nm; m++) {
		//			h[j][k][m] = 0;
		//			for (int i = 0; i < block_size; i++) {
		//				for (int b = 0; b < nb; b++) {
		//					for (int e = 0; e < ne; e++) {
		//						h[j][k][m] += pixsize * smr[i][j] * (sino->get_sinovec()[i + b * block_size] / nbar[i][b]) * suma[i] * source[v_begin + i][b][e] * matatn[k][e] * matatn[m][e] * lintg[i][e];
		//					}
		//				}
		//			}
		//		}
		//	}
		//}

		//std::cout << "h[0](cpu)" << std::endl;
		//for (int j = 0; j < nm; j++) {
		//	for (int k = 0; k < nm; k++) {
		//		std::cout << h[0][j][k] << ", ";
		//	}
		//	std::cout << std::endl;
		//}
		std::cout << "completed." << std::endl;

		//for (int j = 0; j < nd * nd; j++) {
		//	for (int k = 0; k < nm; k++) {
		//		for (int m = 0; m < nm; m++) {
		//			h[j][k][m] /= ne;
		//		}
		//	}
		//}
		/*for (int j = 0; j < nd * nd; j++) {
			std::cout << "\r" << j << "/" << nd * nd;
			for (int k = 0; k < nm; k++) {
				for (int m = 0; m < nm; m++) {
					for (int i = 0; i < block_size; i++) {
						for (int b = 0; b < nb; b++) {
							for (int e = 0; e < ne; e++) {
								h[j][k][m] += smr[i][j] * (n[i][b] / nbar[i][b]) * suma[i] * source[v_begin + i][b][e] * matatn[k][e] * matatn[m][e] * lintg[i][e];
							}
						}
					}
				}
			}
		}*/
		//std::cout << "h[51456][2][2]: " << h[51456][2][2] << std::endl;
		//for (int j = 0; j < nd * nd; j++) {
		//	for (int k = 0; k < nm; k++) {
		//		if (h[j][k][k] != 0) {
		//			std::cout << h[j][k][2] << " ";
		//			std::cout << std::endl;
		//		}

		//	}

		//}


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
		std::cout << "update images: " ;

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
				for (int m = 0; m < nm; m++) {
					//imgdiff[j][k] += (grad_q[j][k] + grad_s[j][k]) * h_inv[j][k][m];
					imgdiff[j][k] += (grad_q[j][m] + grad_s[j][m]) * h_inv[j][k][m];
					////imgdiff[j][k] += grad_q[j][k] * h_inv[j][k][m] + grad_s[j][k] * hs_inv[j][k][m];
					//imgdiff[j][k] += grad_s[j][k] * hs_inv[j][k][m];
					////if((grad_q[j][k] + grad_s[j][k]) != 0) 
					//if(imgdiff[j][k] < -100000000) 
					if(j < 10) std::cout << "j,k,m, imgdiff, grad, hinv:" << j << ", " << k << ", " << m << ", " << imgdiff[j][k] << ", " << grad_q[j][m] + grad_s[j][m] << ", " << h_inv[j][k][m] << std::endl;
				}
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

		for (int i = 0; i < nd*nd; i++) {
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



	int IterativeMatDec::gaussian_elim(float** input, float** result, int size, bool log) {

		float tmp = 0;
		float* tmpv;

		//float** result = (float**)malloc(sizeof(float*) * size);
		//for (int i = 0; i < size; i++) {
		//	result[i] = (float*)malloc(sizeof(float) * size);
		//}

		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				result[i][j] = 0;
			}
			result[i][i] = 1;
		}

		for (int i = 0; i < size; i++) {

			if (log) {
				std::cout << "step:" << i << std::endl;
				for (int k = 0; k < size; k++) {
					for (int l = 0; l < size; l++) {
						std::cout << input[k][l] << ", ";
					}
					for (int l = 0; l < size; l++) {
						std::cout << result[k][l] << ", ";
					}
					std::cout << std::endl;
				}
			}

			if (input[i][i] != 0) {
				tmp = 1 / input[i][i];
			}
			else {
				for (int j = i + 1; j < size; j++) {
					if (input[j][i] != 0) {
						tmpv = input[j];
						input[j] = input[i];
						input[i] = tmpv;

						tmpv = result[j];
						result[j] = result[i];
						result[i] = tmpv;

						tmp = 1 / input[i][i];

						if (log) {
							std::cout << "swap:" << i << "," << j << std::endl;
							for (int k = 0; k < size; k++) {
								for (int l = 0; l < size; l++) {
									std::cout << input[k][l] << ", ";
								}
								for (int l = 0; l < size; l++) {
									std::cout << result[k][l] << ", ";
								}
								std::cout << std::endl;
							}
						}

						break;
					}
					if (j == size - 1) {
						//std::cout << "invert error!!" << std::endl;
						//for (int k = 0; k < size; k++) {
						//	for (int l = 0; l < size; l++) {
						//		std::cout << input[k][l] << ", ";
						//	}
						//	for (int l = 0; l < size; l++) {
						//		std::cout << result[k][l] << ", ";
						//	}
						//	std::cout << std::endl;
						//}
						return -1;
						//exit(1);
					}
				}
			}
			for (int j = 0; j < size; j++) {
				input[i][j] *= tmp;
				result[i][j] *= tmp;
			}
			for (int j = 0; j < size; j++) {
				if (i != j) {
					tmp = input[j][i];
					for (int k = 0; k < size; k++) {
						input[j][k] -= input[i][k] * tmp;
						result[j][k] -= result[i][k] * tmp;
					}
				}
			}
		}

		return 0;// result;
	}

	std::vector<std::vector<float>>& IterativeMatDec::gaussian_elim(std::vector<std::vector<float>>& input) {

		if (input[0].size() != input.size()) throw square_matrix_error();
		int size = input.size();
		int tmp = 0;
		std::vector<float> tmpv;

		std::vector<std::vector<float>> result(std::vector<std::vector<float>>(size, std::vector<float>(size, 0)));

		for (int i = 0; i < size; i++) {
			result[i][i] = 1;
		}

		for (int i = 0; i < size; i++) {
			if (input[i][i] != 0) {
				tmp = 1 / input[i][i];
			}
			else {
				for (int j = i+1; j < size; j++) {
					if (input[j][i] != 0) {
						tmpv = input[j];
						input[j] = input[i];
						input[i] = tmpv;

						tmpv = result[j];
						result[j] = result[i];
						result[i] = tmpv;

						tmp = 1 / input[i][i];
						break;
					}
					if (j == size - 1) {
						std::cout << "invert error!!" << std::endl;
						exit(1);
					}
				}
			}

			for (int j = 0; j < size; j++) {
				input[i][j] *= tmp;
				result[i][j] *= tmp;
			}
			for (int j = 0; j < size; j++) {
				if (i != j) {
					tmp = input[j][i];
					for (int k = 0; k < size; k++) {
						input[j][k] -= input[i][k] * tmp;
						result[j][k] -= result[i][k] * tmp;
					}
				}
			}
		}

		return result;
	}

	std::unique_ptr<SparseMatrix> IterativeMatDec::generate_sysmat_gpu(int begin, int size, bool init) {

		int nd = sino->get_nd();
		int nv = sino->get_nv();
		int center = nd / 2;

		//float* elements = (float*)malloc(sizeof(float) * MAXMATERIALS);	//all nonzero values
		//int* rowptr = (int*)malloc(sizeof(int) * (nd * size + 1));		//indices of the first nonzero element in each row
		//int* colind = (int*)malloc(sizeof(int) * MAXMATERIALS);		//the column indices of the corresponding elements

		std::unique_ptr<float[]> elements = std::make_unique<float[]>(MAXMATERIALS_MD);
		std::unique_ptr<int[]> rowptr = std::make_unique<int[]>((nd * size + 1));
		std::unique_ptr<int[]> colind = std::make_unique<int[]>(MAXMATERIALS_MD);

		int nonzero = 0;		//the number of nonzero elements

		//std::cout << "\nStart Generating System Matrix(GPU) " << begin << " / " << nv << "\n";
		nonzero = Reconstruction::calc_sysmat2(elements.get(), rowptr.get(), colind.get(), begin, size, nv, nd, center, geometry_normalized->sdd);
		rowptr[nd * size] = nonzero;
		//std::cout << "\nSystem matrix generated!(GPU), nonzero = " << nonzero << "\n";

		std::unique_ptr<SparseMatrix> sysmatptr(new SparseMatrix(elements, rowptr, colind, nonzero));
		return sysmatptr;

	}

	std::unique_ptr<SparseMatrix> IterativeMatDec::generate_sysmat(bool use_gpu, bool write_sysmat) {
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

		std::unique_ptr<float[]> elements = std::make_unique<float[]>(MAXMATERIALS_MD);
		std::unique_ptr<int[]> rowptr = std::make_unique<int[]>((nd * nv + 1));
		std::unique_ptr<int[]> colind = std::make_unique<int[]>(MAXMATERIALS_MD);


		int nonzero = 0;		//the number of nonzero elements


		std::ofstream ofse("C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sysmat\\elem.csv"); //for debug
		std::ofstream ofsr("C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sysmat\\rowptr.csv"); //for debug
		std::ofstream ofsc("C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sysmat\\colind.csv"); //for debug


		center = nd / 2;
		center_relative_x = 0;
		center_relative_y = geometry_normalized->axiscor_pixels;


		if (use_gpu) {
			std::cout << "\nStart Generating System Matrix(GPU) " << nd << " " << nv << "\n";
			//nonzero = Reconstruction::calc_sysmat2(elements.get(), rowptr.get(), colind.get(), 0, nv, nd, center, geometry_normalized->sdd);
			nonzero = Reconstruction::calc_sysmat(elements.get(), rowptr.get(), colind.get(), nv, nd, center, geometry_normalized->sdd);
			rowptr[nd * nv] = nonzero;
			std::cout << "\nSystem matrix generated!(GPU), nonzero = " << nonzero << "\n";

			
			if (write_sysmat) {
				int nlcount = 0;
				for (int i = 0; i < nonzero; i++, nlcount++) {
					ofse << elements[i];
					ofsc << colind[i];
					if (nlcount >= 10000) {
						ofse << ",";
						ofsc << ",";
					}
					else {
						ofse << std::endl;
						ofsc << std::endl;
					}
				}
				nlcount = 0;
				for (int i = 0; i < nd * nv + 1; i++, nlcount++) {
					ofsr << rowptr[i];
					if (nlcount >= 10000) {
						ofsr << ",";
					}
					else {
						ofsr << std::endl;
					}
				}
			}

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
				{ //“Š‰eŠp‚ª45“x‚ð’´‚¦‚½‚ç‰æ‘œ‚ð90“x‰E‰ñ“]‚³‚¹“Š‰eŠp‚ð - 45“x‚É
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
								if (nonzero == MAXMATERIALS_MD - 1) {
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

							//if (write_sysmat) {
							//	ofs << area << ", ";
							//}
						}
					}

					//if (write_sysmat) {
					//	ofs << "\n";
					//}
				}

				theta += 2 * PI / nv;
			}

			rowptr[nd * nv] = nonzero;

			std::cout << "End Generating System Matrix";

			std::unique_ptr<SparseMatrix> sysmatptr(new SparseMatrix(elements, rowptr, colind, nonzero));
			return sysmatptr;
		}
	}

	float IterativeMatDec::calc_area(float intercept, float detectoroffset, float angle) {
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

	float IterativeMatDec::calc_area_cbct(float intercept, float detectoroffset, float angle) {
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

	float IterativeMatDec::_calc_subarea(float l1, float l2, float a) {

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

	geometry* IterativeMatDec::get_geometry() {
		return geometry_normalized;
	}
}

