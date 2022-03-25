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

		SparseMatrix Extract_blockmat(SparseMatrix mat, int begin, int rows) {

		}

		float* Extract_row(SparseMatrix mat, int row) {

		}

	};
}