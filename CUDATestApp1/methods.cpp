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
}