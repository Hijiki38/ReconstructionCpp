#pragma once
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "sinogram.h"


namespace Reconstruction {
	vector<string> splitstring(string& input, char delimiter);
	float* add_array(float* a, float* b);
	float* sub_array(float* a, float* b);
	float* mul_array(float* a, float* b);
	float* mul_array(float attenu, float* a);
	float dot_array(float* a, float* b);
	float sum_array(float* a);
}


