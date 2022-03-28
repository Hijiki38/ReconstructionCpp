#include "methods.h"

namespace Reconstruction {
	vector<string> splitstring(string& input, char delimiter) {
		istringstream stream(input);
		string field;
		vector<string> output;

		while (getline(stream, field, delimiter)) {
			output.push_back(field);
		}

		return output;
	}

	float* add_array(float* a, float* b) {
		if (sizeof(a) != sizeof(b)) {
			std::cerr << "size of arrays must be same." << std::endl;
			return NULL;
		}
		float* result = (float*)malloc(sizeof(a));
		for (int i = 0; i < (sizeof(result) / sizeof(float)); i++) {
			result[i] = a[i] + b[i];
		}
		return result;
	}

	float* sub_array(float* a, float* b) {
		if (sizeof(a) != sizeof(b)) {
			std::cerr << "size of arrays must be same." << std::endl;
			return NULL;
		}
		float* result = (float*)malloc(sizeof(a));
		for (int i = 0; i < (sizeof(result) / sizeof(float)); i++) {
			result[i] = a[i] - b[i];
		}
		return result;
	}

	float* mul_array(float* a, float* b) {
		if (sizeof(a) != sizeof(b)) {
			std::cerr << "size of arrays must be same." << std::endl;
				return NULL;
		}
		float* result = (float*)malloc(sizeof(a));
		for (int i = 0; i < (sizeof(result) / sizeof(float)); i++) {
			result[i] = a[i] * b[i];
		}
		return result;
	}

	float* mul_array(float attenu, float* a) {
		float* result = (float*)malloc(sizeof(a));
		for (int i = 0; i < (sizeof(result) / sizeof(float)); i++) {
			result[i] = a[i] * attenu;
		}
		return result;
	}

	float dot_array(float* a, float* b) {
		if (sizeof(a) != sizeof(b)) {
			std::cerr << "size of arrays must be same." << std::endl;
				return NULL;
		}
		float result = 0;
		for (int i = 0; i < (sizeof(a) / sizeof(float)); i++) {
			result += a[i] * b[i];
		}
		return result;
	}

	float sum_array(float *a) {
		float result = 0;
		for (int i = 0; i < (sizeof(a) / sizeof(float)); i++) {
			result += a[i];
		}
		return result;
	}
}