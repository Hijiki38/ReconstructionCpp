#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "ART.h"
#include "sinogram.h"
#include "util.h"
using namespace std;

const double PI = 3.14159265358979;

enum class rec_name {
	FBP,
	ART,
};

int main(int argc, char *argv[]) {

	int mode;
	int count=0;
	string inpath;
	sinogram* sg;
	ART* art;
	VectorXf* result;

	cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+\n|R|E|C|O|N|S|T|R|U|C|T|O|R|\n+-+-+-+-+-+-+-+-+-+-+-+-+-+";
	if (argc == 1) {
		//inpath = ".";
		inpath = "C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\sino_star360.csv";
	}
	else {
		inpath = argv[1];
	}

	cout << "\nInput reconstruction mode: (0:FBP, 1:ART)";
	while (1) {
		cin >> mode;
		if (mode == static_cast<int>(rec_name::FBP)) {
			sg = read_sinogram(inpath);
			art = new ART(sg);
			cout << "activate FBP";
			result = (*art).reconstruction();
			break;
		}
		else if (mode == static_cast<int>(rec_name::ART)) {
			sg = read_sinogram(inpath);
			cout << "\nread sinogram completed. d, v =" << (*sg).get_nd() << ", " << (*sg).get_nv();
			art = new ART(sg);
			cout << "\nactivate ART";
			result = (*art).reconstruction();
			break;
		}
		else {
			cout << "Undefined value.";
		}
	}

	vector<float> outvec;
	outvec.resize((*result).size());
	VectorXf::Map(&outvec[0], (*result).size()) = (*result);

	cout << "received vector:" << (*result)[0] << "," << (*result)[1] << "," << (*result)[2];

	ofstream ofs("C:\\Users\\takum\\Dropbox\\Aoki_Lab\\util\\Reconstructor\\output\\vsARToutput.csv");
	int ressize = (*result).size();
	int nd = (*sg).get_nd();
	cout << "outvec_len:" << outvec.size();
	for (int i = 0; i < ressize / nd; i++) {
		cout << "\r writing..  " << i << " / " << ressize;
		for (int j = 0; j < nd; j++) {
			cout << count << "\n";
			ofs << outvec[count] << ',';
			count++;
		}
		ofs << '\n';
	}
}

/*sinogram* read_sinogram(string inpath) {
	ifstream stream(inpath);
	string line;
	vector<float> sinovec;
	int count_row = 0;
	int count_all = 0;

	while (getline(stream, line)) {
		vector<string> strs = split(line, ',');
		for (int i = 0; i < strs.size(); i++) {
			sinovec.push_back(stof(strs.at(i)));
			count_all++;
		}
		count_row++;
	}

	sinogram sino(sinovec, count_all / count_row, count_row);
	return &sino;
}

vector<string> split(string& input, char delimiter) {
	istringstream stream(input);
	string field;
	vector<string> output;

	while (getline(stream, field, delimiter)) {
		output.push_back(field);
	}

	return output;
}*/





