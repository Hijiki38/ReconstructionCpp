#pragma once
#include<map>
#include<algorithm>
#include<functional>
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include "methods.h"


namespace Reconstruction {
	class spectrum {
	private:
		float* data;
		int size;

	public:
		//spectrum(std::string p) {
		//	data = (float*)malloc(s * sizeof(float));
		//	size = s;
		//	read_data(p, data, size);
		//}

		spectrum(float* d, int s) {
			data = d;
			size = s;
		}

		float* get_data();

		int get_size();

		static spectrum read_data(std::string path);

	};
}