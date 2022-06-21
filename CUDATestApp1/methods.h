#pragma once
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "sinogram.h"


namespace Reconstruction {
	std::vector<std::string> splitstring(std::string& input, char delimiter);
	void add_array(float* a, float* b, std::size_t length);
	void sub_array(float* a, float* b, std::size_t length);
	void mul_array1(float* a, float attenu, std::size_t length);
	void mul_array2(float* a, float* b, std::size_t length);
	float dot_array(float* a, float* b, std::size_t length);
	float sum_array(float* a, std::size_t length);
}


