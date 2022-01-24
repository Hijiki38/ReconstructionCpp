#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "ART.h"
#include "sinogram.h"
#include "methods.h"

const double PI = 3.14159265358979;

enum class rec_name {
	FBP,
	ART,
};

int main(int argc, char *argv[]) {

	int mode;
	int count=0;
	string inpath;
	Reconstruction::sinogram* sg;
	Reconstruction::ART* art;
	VectorXf* result = nullptr;

	vector<float> outvec;
	int ressize;
	int nd;
	int nv;
	ofstream ofs("C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\vsARToutput.csv");


	std::cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+\n|R|E|C|O|N|S|T|R|U|C|T|O|R|\n+-+-+-+-+-+-+-+-+-+-+-+-+-+";
	if (argc == 1) {
		//inpath = ".";
		inpath = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sino_star.csv";
		//inpath = "E:\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sino_star.csv";
	}
	else {
		inpath = argv[1];
	}

	std::cout << "\nInput reconstruction mode: (0:FBP, 1:ART)";
	while (1) {
		std::cin >> mode;
		sg = (*sg).read_sinogram(inpath);
		std::cout << "\nread sinogram completed. d, v =" << (*sg).get_nd() << ", " << (*sg).get_nv();
		
		if (mode == static_cast<int>(rec_name::FBP)) {
			break;
		}
		else if (mode == static_cast<int>(rec_name::ART)) {
			break;
		}
		else {
			std::cout << "Undefined value.";
		}
	}

	if (mode == static_cast<int>(rec_name::FBP)) {
		//art = new Reconstruction::ART(sg);
		//cout << "activate FBP";
		//result = (*art).reconstruction();
	}
	else if (mode == static_cast<int>(rec_name::ART)) {
		art = new Reconstruction::ART(sg, 205.7, 1100.0, -0.9345, 0.01869);
		std::cout << "\nactivate ART";
		result = (*art).reconstruction(3,5);
	}


	outvec.resize((*result).size());
	VectorXf::Map(&outvec[0], (*result).size()) = (*result);

	std::cout << "received vector:" << (*result)[0] << "," << (*result)[1] << "," << (*result)[2];

	
	//ofstream ofs("E:\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\vsARToutput.csv");
	ressize = (*result).size();
	nd = (*sg).get_nd();
	nv = (*sg).get_nv();

	//std::cout << "outvec_len:" << outvec.size();

	for (int i = 0; i < ressize / nd; i++) {
		std::cout << "\r writing..  " << i << " / " << (ressize / nd) - 1;
		for (int j = 0; j < nd; j++) {
			//std::cout << count << "\n";
			ofs << outvec[count] << ',';
			count++;
		}
		ofs << '\n';
	}
}





