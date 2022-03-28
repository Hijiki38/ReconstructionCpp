#pragma once

#include <stdio.h>
#include <iostream>
#include <string>
#include "Eigen/Sparse"

using namespace std;
using namespace Eigen;

namespace Reconstruction {
	class SparseMatrix {
	private:
		//CSR format
		float* elements;	//all nonzero values
		int* rowptr;		//indices of the first nonzero element in each row
		int* colind;		//the column indices of the corresponding elements
		int nonzero;		//the number of nonzero elements

	public:
		SparseMatrix() {

		}

		SparseMatrix(float* elem, int* rptr, int* cind, int nz) {
			elements = elem;
			rowptr = rptr;
			colind = cind;
			nonzero = nz;
		}

		float* get_elements() { return elements; }
		int* get_rowptr() { return rowptr; }
		int* get_colind() { return colind; }
		int get_nonzero() { return nonzero; }

		SparseMatrix Extract_blockmat(int begin, int rows) {

			float* elements_new;	
			int* rowptr_new;	
			int* colind_new;	

			float* elements_orig = elements;	
			int* rowptr_orig = rowptr;		
			int* colind_orig = colind;

			int elem_begin, elem_num;


			elem_begin = colind_orig[begin];
			elem_num = colind_orig[begin + rows] - elem_begin;

			elements_new = (float*)malloc(elem_num * sizeof(float));
			rowptr_new = (int*)malloc((rows + 1) * sizeof(int));
			colind_new = (int*)malloc(elem_num * sizeof(float));



			if (colind_new != NULL && elements_new != NULL && rowptr_new != NULL) {
				for (int i = 0; i < rows; i++) {
					rowptr_new[i] = rowptr_orig[begin + i] - elem_begin;
				}
				rowptr_new[rows] = elem_num;

				for (int i = 0; i < elem_num; i++) {
					elements_new[i] = elements_orig[i + elem_begin];
					colind_new[i] = colind_orig[i + elem_begin];
				}
			}
			else {
				std::cerr << "Refering Null Array!";
			}

			return SparseMatrix(elements_new, rowptr_new, colind_new, elem_num);

		}

		float* Extract_row_dense(int row, int num_col) {

		}

	};
}