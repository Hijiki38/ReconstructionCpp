#pragma once

#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "sinogram.h"

using std::cout;

sinogram* read_sinogram(string inpath);
vector<string> split(string& input, char delimiter);

sinogram* read_sinogram(string inpath) {
	ifstream stream(inpath);
	string line;
	vector<float> sinovec;
	int count_row = 0;
	int count_all = 0;

	cout << "read sinogram";

	while (getline(stream, line)) {
		cout << "readrow:" << count_row << "\n";
		vector<string> strs = split(line, ',');
		for (int i = 0; i < strs.size(); i++) {
			sinovec.push_back(stof(strs.at(i)));
			count_all++;
			cout << "readelem:" << count_all << "\n";
		}
		count_row++;
	}

	cout << "d: " << (count_all / count_row) << ", v: " << count_row << "\n";

	sinogram *sino = new sinogram(sinovec, count_all / count_row, count_row);

	cout << "d: " << (*sino).get_nd() << ", v: " << (*sino).get_nv() << "\n";

	return sino;
}

vector<string> split(string& input, char delimiter) {
	istringstream stream(input);
	string field;
	vector<string> output;

	while (getline(stream, field, delimiter)) {
		output.push_back(field);
	}

	return output;
}

