#pragma once
#include<map>
#include<algorithm>
#include<functional>
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <assert.h>
#include "methods.h"
#include "spectrum.h"

namespace Reconstruction {

	struct single_material {
		std::string name;
		float density;
		Reconstruction::spectrum attenu;
	};

	class materials {
	private:
		std::vector<single_material> mat_list;

	public:
		materials(){}

		materials(std::string inpath) {
			read_materials(inpath);
		}

		static materials read_materials(std::string inpath) {
			materials result;

			std::cout << "\nread materials: " << inpath << std::endl;

			for (const auto& file : std::filesystem::directory_iterator(inpath)) {

				std::cout << "\nread " << file.path().string() << std::endl;

				std::string materialname, filename;
				std::vector<std::string> splited_filename;
				float dens;

				Reconstruction::spectrum coeff = spectrum::read_data(file.path().string());

				filename = file.path().stem().string();
				splited_filename = Reconstruction::splitstring(filename, '_');  //(Material name)_(density).csv
				assert(splited_filename.size() == 2);
				materialname = splited_filename[0];
				dens = stof(splited_filename[1]);

				std::cout << "registering " << materialname << std::endl;
				//coeff.read_data(file.path().string());

				single_material mat = { materialname, dens, coeff };
				result.mat_list.push_back(mat);

			}

			std::cout << "materials are successfully read." << std::endl;
			return result;
		}

		std::vector<single_material> get_matlist() { return mat_list; }
	};
}