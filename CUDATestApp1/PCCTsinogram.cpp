#include"PCCTsinogram.h"
#include"methods.h"
#include <thread>
#include <time.h>


namespace Reconstruction {

	PCCTsinogram* PCCTsinogram::read_PCCTsinogram(std::vector<std::string> inpaths) {

		float* tmp_sino;
		vector<float> tmp_sinovec;
		bool is_first_file = true;
		int count_row = 0;
		int count_all = 0;
		int count_bins = 0;

		for (const auto& inpath : inpaths) {
			ifstream stream(inpath);
			string line;

			cout << "readbin:" << count_bins << "\n";

			while (getline(stream, line)) {
				vector<string> strs = Reconstruction::splitstring(line, ',');
				for (int i = 0; i < strs.size(); i++) {
					tmp_sinovec.push_back(stof(strs.at(i)));
					if(is_first_file){ count_all++; }
				}

				if(is_first_file){ count_row++; }
			}

			count_bins++;
			is_first_file = false;
		}

		tmp_sino = (float*)malloc(sizeof(float) * tmp_sinovec.size());
		std::copy(tmp_sinovec.begin(), tmp_sinovec.end(), tmp_sino);

		cout << "d: " << (count_all / count_row) << ", v: " << count_row << ", e: " << count_bins << "\n";

		PCCTsinogram* sino = new PCCTsinogram(tmp_sino, count_all / count_row, count_row, count_bins);

		cout << "d: " << (*sino).get_nd() << ", v: " << (*sino).get_nv() << ", e: " << (*sino).get_ne() << "\n";

		return sino;
	}
}