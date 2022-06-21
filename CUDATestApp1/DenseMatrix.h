#pragma once
#pragma once

#include <stdio.h>
#include <iostream>
#include <string>
//#include "Eigen/Sparse"


namespace Reconstruction {
	class matrix_product_error : public std::exception {
		virtual const char* what() const noexcept { return "Error: Invalid number of columns and rows."; }
	};

	//template <typename T>
	class DenseMatrix {
	private:
		std::vector<std::vector<float>> elements;

	public:

		int row;
		int col;

		DenseMatrix() {

		}

		DenseMatrix(int y, int x) {
			elements = *(new std::vector<std::vector<float>>(y, std::vector<float>(x, 0)));
			row = y;
			col = x;
		}

		DenseMatrix(float**const elem, int y, int x) {
			elements = *(new std::vector<std::vector<float>>(y, std::vector<float>(x, 0)));
			row = y;
			col = x;
			for (int i = 0; i < y; i++) {
				for (int j = 0; j < x; j++) {
					//std::cout << "(y,x) = " << i << ", " << j << std::endl;
					elements[i][j] = elem[i][j];
				}
			}
		}

		DenseMatrix(float* const elem, int size, int axis) { //0:vertical  1:horizontal
			if (axis == 0) {
				elements = *(new std::vector<std::vector<float>>(size, std::vector<float>(1, 0)));
				row = size;
				col = 1;
			}
			else {
				elements = *(new std::vector<std::vector<float>>(1, std::vector<float>(size, 0)));
				row = 1;
				col = size;
			}
		}

		DenseMatrix(std::vector<float> elem) {
			elements.push_back(elem);
			row = 1;
			col = elem.size();
		}

		~DenseMatrix() {
			
		}

		std::vector<float>& operator[](const int i) {
			return elements[i];
		}

		DenseMatrix operator-(DenseMatrix other) const {
			if (row != other.row || col != other.col) throw matrix_product_error();

			DenseMatrix result(row, col);
			for (int i = 0; i < row; i++) {
				for (int j = 0; j < col; j++) {
					result[i][j] -= other[i][j];
				}
			}

			return result;
		}

		DenseMatrix operator*(DenseMatrix other) const {
			if (col != other.row) throw matrix_product_error();

			DenseMatrix result(row, other.col);
			for (int i = 0; i < row; i++) {
				for (int j = 0; j < other.col; j++) {
					for (int k = 0; k < col; k++) {
						result[i][j] += elements[i][k] * other[k][j];
					}
				}
			}

			return result;
		}

		DenseMatrix& operator*=(const float x) {
			for (int i = 0; i < row; i++) {
				for (int j = 0; j < col; j++) {
					elements[i][j] *= x;
				}
			}

			return *this;
		}

		DenseMatrix operator*(const float x) {
			return DenseMatrix(*this) *= x;
		}

		static DenseMatrix& T(DenseMatrix mat) {
			DenseMatrix result(mat.col, mat.row);
			for (int i = 0; i < mat.row; i++) {
				for (int j = 0; j < mat.col; j++) {
					result[j][i] = mat[i][j];
				}
			}

			return result;
		}

		static DenseMatrix hadamard(DenseMatrix mat1, DenseMatrix mat2) {
			if (mat1.row != mat2.row || mat1.col != mat2.col) throw matrix_product_error();

			DenseMatrix result(mat1.row, mat1.col);
			for (int i = 0; i < mat1.row; i++) {
				for (int j = 0; j < mat1.col; j++) {
					result[i][j] = mat1[i][j] * mat2[i][j];
				}
			}

			return result;
		}

		DenseMatrix& exp() {
			for (int i = 0; i < row; i++) {
				for (int j = 0; j < col; j++) {
					elements[i][j] = std::exp(elements[i][j]);
				}
			}

			return *this;
		}

		std::vector<std::vector<float>>& get_elem() {
			return elements;
		}


		

	};
}