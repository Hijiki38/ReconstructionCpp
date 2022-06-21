#include "spectrum.h"

namespace Reconstruction {

	float* spectrum::get_data() {
		return data;
	}

	int spectrum::get_size() { return size; }

	spectrum spectrum::read_data(std::string path) {

		std::cout << "start reading spectrum data\n";
		float* tmp;

		std::ifstream stream = std::ifstream(path);
		std::string line;

		int count = 0;
		while (getline(stream, line)) {
			count++;
		}

		stream.clear();
		stream.seekg(0, std::ios_base::beg);

		std::cout << "count:" << count << std::endl;

		tmp = (float*)malloc(sizeof(float) * count);
		count = 0;

		while (getline(stream, line)) {
			std::vector<std::string> strs = Reconstruction::splitstring(line, ' ');
			tmp[count] = stof(strs[1]);

			//std::cout << "added: " << stof(strs[0]) << ", " << tmp[count];

			count++;
		}

		std::cout << "count:" << count << std::endl;

		spectrum result(tmp, count);

		return result;
	}
}