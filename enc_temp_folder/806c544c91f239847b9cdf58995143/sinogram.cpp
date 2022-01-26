#include"sinogram.h"
#include"methods.h"


namespace Reconstruction {

	VectorXf sinogram::get_sinovec() { return eigen_sinovec; }
	
	int sinogram::get_nd() { return n_d; }
	
	int sinogram::get_nv() { return n_v; }

	sinogram* sinogram::read_sinogram(string inpath) {
		ifstream stream(inpath);
		string line;
		vector<float> sinovec;
		int count_row = 0;
		int count_all = 0;

		cout << "read sinogram";

		while (getline(stream, line)) {
			cout << "\rreadrow:" << count_row;
			vector<string> strs = Reconstruction::splitstring(line, ',');
			for (int i = 0; i < strs.size(); i++) {
				sinovec.push_back(stof(strs.at(i)));
				count_all++;
				cout << "\rreadelem:" << count_all;
			}
			count_row++;
		}

		cout << "\nd: " << (count_all / count_row) << ", v: " << count_row;

		sinogram* sino = new sinogram(sinovec, count_all / count_row, count_row);

		cout << "\nd: " << (*sino).get_nd() << ", v: " << (*sino).get_nv();

		return sino;
	}
}