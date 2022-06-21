#include"PCCTsinogram.h"
#include"methods.h"
#include <thread>
#include <time.h>


namespace Reconstruction {

	PCCTsinogram* PCCTsinogram::read_PCCTsinogram(std::vector<std::string> inpaths_obj, std::vector<std::string> inpaths_bg) {

		float* tmp_sino, * tmp_bgsino;
		std::vector<float> tmp_sinovec, tmp_bgsinovec;
		bool is_first_file = true;
		int count_row = 0;
		int count_all = 0;
		int count_bins = 0;

		//OBJ
		for (const auto& inpath : inpaths_obj) {
			std::ifstream stream(inpath);
			std::string line;

			while (getline(stream, line)) {
				std::vector<std::string> strs = Reconstruction::splitstring(line, ',');
				for (int i = 0; i < strs.size(); i++) {
					tmp_sinovec.push_back(stof(strs.at(i)));
					if(is_first_file){ count_all++; }
				}

				if(is_first_file){ count_row++; }
			}

			count_bins++;
			is_first_file = false;

			std::cout << "\nreadbin:" << count_bins << "\n";
		}

		tmp_sino = (float*)malloc(sizeof(float) * tmp_sinovec.size());
		std::copy(tmp_sinovec.begin(), tmp_sinovec.end(), tmp_sino);

		//BG
		std::cout << "read bg" << std::endl;
		for (const auto& inpath : inpaths_bg) {
			std::ifstream stream(inpath);
			std::string line;

			while (getline(stream, line)) {
				std::vector<std::string> strs = Reconstruction::splitstring(line, ',');
				for (int i = 0; i < strs.size(); i++) {
					tmp_bgsinovec.push_back(stof(strs.at(i)));
				}
			}
		}

		tmp_bgsino = (float*)malloc(sizeof(float) * tmp_bgsinovec.size());
		std::copy(tmp_bgsinovec.begin(), tmp_bgsinovec.end(), tmp_bgsino);


		//generate instance
		std::cout << "d: " << (count_all / count_row) << ", v: " << count_row << ", e: " << count_bins << "\n";

		PCCTsinogram* sino = new PCCTsinogram(tmp_sino, tmp_bgsino, count_all / count_row, count_row, count_bins);

		std::cout << "d: " << (*sino).get_nd() << ", v: " << (*sino).get_nv() << ", e: " << (*sino).get_nb() << "\n";

		return sino;
	}

	void PCCTsinogram::convert_negativelog() {
		int size = get_nd() * get_nv() * get_nb();
		float* a = get_sinovec();
		for (int i = 0; i < size; i++) {
			a[i] = -1 * std::log(a[i]);
		}
	}
}