#include "methods.h"
#include <stdexcept>

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

	void add_array(float* a, float* b, std::size_t length) {

		//float* result = (float*)malloc(sizeof(float) * length);

		for (int i = 0; i < length; i++) {
			a[i] = a[i] + b[i];
		}
	}

	void sub_array(float* a, float* b, std::size_t length) {
		//float* result = (float*)malloc(sizeof(float) * length);
		for (int i = 0; i < length; i++) {
			a[i] = a[i] - b[i];
		}
	}

	void mul_array1(float*a, float attenu, std::size_t length) {
		//float* result = (float*)malloc(sizeof(float) * length);
		for (int i = 0; i < length; i++) {
			a[i] = a[i] * attenu;
		}
	}

	void mul_array2(float* a, float* b, std::size_t length) {
		//float* result = (float*)malloc(sizeof(float) * length);
		for (int i = 0; i < length; i++) {
			a[i] = a[i] * b[i];
		}
	}

	float dot_array(float* a, float* b, std::size_t length) {
		float result = 0;
		for (int i = 0; i < length; i++) {
			result += a[i] * b[i];
		}
		return result;
	}

	float sum_array(float *a, std::size_t length) {
		float result = 0;
		for (int i = 0; i < length; i++) {
			result += a[i];
		}
		return result;
	}
}