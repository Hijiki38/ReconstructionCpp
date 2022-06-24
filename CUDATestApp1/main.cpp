#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "ART.h"
#include "MLEM.h"
#include "IterationRec.h"
#include "sinogram.h"
#include "PCCTsinogram.h"
#include "methods.h"
#include "geometry.h"
#include "crtdbg.h"
#include "IterativeMatDec.h"

#define malloc(X) _malloc_dbg(X,_NORMAL_BLOCK,__FILE__,__LINE__)
#define new ::new(_NORMAL_BLOCK, __FILE__, __LINE__)

const double PI = 3.14159265358979;

enum class rec_name {
	FBP,
	ART,
	MLEM,
	MATDEC,
};

int main(int argc, char *argv[]) {
	
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_DELAY_FREE_MEM_DF | _CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF); //detects heap colluption

	int mode;
	int count=0;
	int blockproj = 1; //projection per single block

	std::string inpath, projdir, bgdir, matdir, sourcedir;
	//Reconstruction::geometry geo = {true, 205.7, 1100.0, -0.9345, 0.01869}; //sod, sdd, cor(pixels), pixsize
	//Reconstruction::geometry geo = { true, 100.0, 610.0, 63, 0.01639 };
	Reconstruction::geometry geo = { true, 400.0, 800.0, 0, 0.005 };
	Reconstruction::sinogram* sg;
	Reconstruction::PCCTsinogram* pcsg;

	Reconstruction::ART* art;
	Reconstruction::MLEM* mlem;
	Reconstruction::IterativeMatDec* matdec;
	VectorXf* result = nullptr;
	float* resultvec;
	float** resultfracvec;

	std::vector<float> outvec;
	int ressize;
	int nd;
	int nv;
	
	int debugmode = 3; //-1:disable, 0,1,2...:reconstruction mode

	std::cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+\n|R|E|C|O|N|S|T|R|U|C|T|O|R|\n+-+-+-+-+-+-+-+-+-+-+-+-+-+";
	//if (argc == 1) 
	//{
	//	//inpath = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\margedsinogram10.csv"a;
	//	inpath = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\fourmetal10_180f32x1000.csv";
	//	//inpath = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sinogram_total_f32.csv";
	//	//inpath = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sino_ball2.csv"; 
	//	//inpath = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sino_star360.csv";
	//	
	//	//inpath = "D:\\tmp\\10-30_f32_2.csv";
	//	//inpath = "E:\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sino_star.csv";
	//}
	//else 
	//{
	//	inpath = argv[1];
	//}

	//sg = Reconstruction::sinogram::read_sinogram(inpath);
	//if (sg == nullptr) 
	//{
	//	"Can not read sinogram file!";
	//	return 1;
	//}

	if (argc == 1) {
		//projdir   = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\input\\fourmetals_simulation\\proj\\ideal6bin\\obj";
		//bgdir     = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\input\\fourmetals_simulation\\proj\\ideal6bin\\bg";
		//matdir    = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\input\\fourmetals_simulation\\material\\1keV";
		//sourcedir = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\input\\fourmetals_simulation\\source\\ideal6bin\\1keV";
		projdir = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\input\\fourmetals_simulation\\proj\\eq5bin\\obj";
		bgdir = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\input\\fourmetals_simulation\\proj\\eq5bin\\bg";
		matdir = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\input\\fourmetals_simulation\\material\\0.1keV";
		sourcedir = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\input\\fourmetals_simulation\\source\\s5bin\\0.1keV";
	}

	std::vector<std::string> objpaths, bgpaths;
	for (const auto& file : std::filesystem::directory_iterator(projdir)) {
		std::cout << "Loading: " << file.path().string() << std::endl;
		objpaths.push_back(file.path().string());
	}
	for (const auto& file : std::filesystem::directory_iterator(bgdir)) {
		std::cout << "Loading: " << file.path().string() << std::endl;
		bgpaths.push_back(file.path().string());
	}



	pcsg = Reconstruction::PCCTsinogram::read_PCCTsinogram(objpaths, bgpaths);
	if (pcsg == nullptr)
	{
		"Can not read sinogram file!";
		return 1;
	}

	Reconstruction::materials mat = Reconstruction::materials::read_materials(matdir);
	std::cout << "loaded material files: " << mat.get_matlist().size();
	std::vector<Reconstruction::spectrum> source;

	std::cout << "Load source files" << std::endl;
	for (const auto& file : std::filesystem::directory_iterator(sourcedir)) {
		std::cout << "Loading: " << file.path().string() << std::endl;
		source.push_back(Reconstruction::spectrum::read_data(file.path().string()));
	}
	std::cout << "finished" << std::endl;

	//pcsg->convert_negativelog();  //convert sinogram to negative log

	std::cout << "\nread sinogram completed. d, v =" << (*pcsg).get_nd() << ", " << (*pcsg).get_nv();
	std::cout << "\nInput reconstruction mode: (0:FBP, 1:ART, 2;MLEM)";
	if (debugmode == -1) {
		while (1)
		{
			std::cin >> mode;
			if (mode == static_cast<int>(rec_name::FBP)
				|| mode == static_cast<int>(rec_name::ART)
				|| mode == static_cast<int>(rec_name::MLEM))
			{
				break;
			}
			else
			{
				std::cout << "Undefined value.";
			}
		}
	}
	else {
		mode = debugmode;
	}


	if (mode == static_cast<int>(rec_name::FBP)) 
	{
		//art = new Reconstruction::ART(sg);
		//cout << "activate FBP";
		//result = (*art).reconstruction();
	}
	else if (mode == static_cast<int>(rec_name::ART)) 
	{
		std::ofstream ofs("C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\vsARToutput.csv");

		art = new Reconstruction::ART(pcsg, &geo, blockproj);
		std::cout << "\nactivate ART";
		resultvec = (*art).reconstruction(5, false);

		nd = (*pcsg).get_nd();
		nv = (*pcsg).get_nv();

		for (int i = 0; i < nd; i++)
		{
			std::cout << "\r writing..  " << i << " / " << nd - 1;
			for (int j = 0; j < nd; j++)
			{
				ofs << resultvec[count] << ',';
				count++;
			}
			ofs << '\n';
		}
	}
	else if (mode == static_cast<int>(rec_name::MLEM))
	{
		std::ofstream ofs("C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\vsMLEMoutput.csv");

		mlem = new Reconstruction::MLEM(pcsg, &geo, blockproj);
		std::cout << "\nactivate MLEM";
		resultvec = (*mlem).reconstruction(5, false);

		nd = (*pcsg).get_nd();
		nv = (*pcsg).get_nv();

		for (int i = 0; i < nd; i++)
		{
			std::cout << "\r writing..  " << i << " / " << nd - 1;
			for (int j = 0; j < nd; j++)
			{
				ofs << resultvec[count] << ',';
				count++;
			}
			ofs << '\n';
		}
	}
	else if (mode == static_cast<int>(rec_name::MATDEC))
	{
		matdec = new Reconstruction::IterativeMatDec(pcsg, &geo, &mat, source, blockproj);
		std::cout << "\nactivate Material decomposition";
		resultfracvec = (*matdec).reconstruction(1);

		nd = (*pcsg).get_nd();
		nv = (*pcsg).get_nv();

		for (int m = 0; m < mat.get_matlist().size(); m++) {

			std::string outfilename = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\material_decomposition\\output1";
			outfilename += mat.get_matlist()[m].name;
			outfilename += ".csv";
			std::ofstream ofs(outfilename);

			count = 0;
			for (int i = 0; i < nd; i++)
			{
				std::cout << "\r writing..  " << i << " / " << nd - 1;
				for (int j = 0; j < nd; j++)
				{
					ofs << resultfracvec[count][m] << ',';
					count++;
				}
				ofs << '\n';
			}
		}
	}





}





