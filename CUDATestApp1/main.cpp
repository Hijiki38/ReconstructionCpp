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

#define malloc(X) _malloc_dbg(X,_NORMAL_BLOCK,__FILE__,__LINE__)
#define new ::new(_NORMAL_BLOCK, __FILE__, __LINE__)

const double PI = 3.14159265358979;

enum class rec_name {
	FBP,
	ART,
	MLEM,
};

int main(int argc, char *argv[]) {
	
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_DELAY_FREE_MEM_DF | _CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF);

	int mode;
	int count=0;
	string inpath;
	//Reconstruction::geometry geo = {true, 205.7, 1100.0, -0.9345, 0.01869}; //sod, sdd, cor, pixsize
	Reconstruction::geometry geo = { false, 400.0, 800.0, 0, 0.05 };
	Reconstruction::sinogram* sg;
	Reconstruction::PCCTsinogram* pcsg;
	Reconstruction::ART* art;
	Reconstruction::MLEM* mlem;
	VectorXf* result = nullptr;
	float* resultvec;

	vector<float> outvec;
	int ressize;
	int nd;
	int nv;
	
	int debugmode = 1; //-1:disable, 0,1,2...:reconstruction mode

	std::cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+\n|R|E|C|O|N|S|T|R|U|C|T|O|R|\n+-+-+-+-+-+-+-+-+-+-+-+-+-+";
	if (argc == 1) 
	{
		//inpath = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sinogram_total_f32.csv";
		// 
		//inpath = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sino_ball2.csv"; 
		inpath = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sino_star360.csv";
		
		//inpath = "D:\\tmp\\10-30_f32_2.csv";
		//inpath = "E:\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sino_star.csv";
	}
	else 
	{
		inpath = argv[1];
	}

	sg = Reconstruction::sinogram::read_sinogram(inpath);
	if (sg == nullptr) 
	{
		"Can not read sinogram file!";
		return 1;
	}

	std::vector<std::string> inpaths;
	inpaths.push_back(inpath);

	pcsg = Reconstruction::PCCTsinogram::read_PCCTsinogram(inpaths);
	if (pcsg == nullptr)
	{
		"Can not read sinogram file!";
		return 1;
	}

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
		ofstream ofs("C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\vsARToutput.csv");

		art = new Reconstruction::ART(pcsg, &geo);
		std::cout << "\nactivate ART";
		resultvec = (*art).reconstruction(10, 5);

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
		ofstream ofs("C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\vsMLEMoutput.csv");

		mlem = new Reconstruction::MLEM(pcsg, &geo);
		std::cout << "\nactivate MLEM";
		resultvec = (*mlem).reconstruction(10, 5);

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





}





