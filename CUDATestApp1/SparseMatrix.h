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
<<<<<<< Updated upstream
		float* elements;	//all nonzero values
		int* rowptr;		//indices of the first nonzero element in each row
		int* colind;		//the column indices of the corresponding elements
=======
		std::unique_ptr<float[]> elements;	//all nonzero values
		std::unique_ptr<int[]> rowptr;		//indices of the first nonzero element in each row
		std::unique_ptr<int[]> colind;		//the column indices of the corresponding elements
>>>>>>> Stashed changes
		int nonzero;		//the number of nonzero elements

	public:
<<<<<<< Updated upstream
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
=======
		SparseMatrix(std::unique_ptr<float[]>& elem, std::unique_ptr<int[]>& rptr, std::unique_ptr<int[]>& cind, int nz)
			: elements(move(elem))
			, rowptr(move(rptr))
			, colind(move(cind))
			, nonzero(nz){}


		//SparseMatrix(float* elem, int* rptr, int* cind, int nz) {
		//	elements = elem;
		//	rowptr = rptr;
		//	colind = cind;
		//	nonzero = nz;
		//}

		//~SparseMatrix() {
		//	std::cout << "destruction!";
		//	if (elements) 
		//		elements.reset();
		//	if (rowptr) 
		//		rowptr.reset();
		//	if (colind) 
		//		colind.reset();
		//	std::cout << " ...finished!\n";
		//}

		~SparseMatrix(){}

		//std::unique_ptr<float[]> get_elements() { return elements; }
		//std::unique_ptr<int[]> get_rowptr() { return rowptr; }
		//std::unique_ptr<int[]> get_colind() { return colind; }
		//int get_nonzero() { return nonzero; }
>>>>>>> Stashed changes

		std::unique_ptr<SparseMatrix> Create_blockmat(int begin, int rows) {

			int elem_num = rowptr[begin + rows] - rowptr[begin];

			float* elements_new = (float*)malloc(elem_num * sizeof(float));
			int* rowptr_new = (int*)malloc((rows + 1) * sizeof(int));
			int* colind_new = (int*)malloc(elem_num * sizeof(float));

			//std::cout << "\nstart generating blockmat\n";

			for (int i = 0; i < rows; i++) {
				rowptr_new[i] = rowptr[begin + i] - rowptr[begin];
				//std::cout << rowptr_new[i] << std::endl;
			}
			rowptr_new[rows] = elem_num;

			for (int i = 0; i < elem_num; i++) {
				elements_new[i] = elements[i + rowptr[begin]];
				colind_new[i] = colind[i + rowptr[begin]];
				//std:cout << "new colind: " << colind_new[i] << std::endl;
			}

<<<<<<< Updated upstream


=======
>>>>>>> Stashed changes
			std::unique_ptr<SparseMatrix> spmat(new SparseMatrix(elements_new, rowptr_new, colind_new, elem_num));

			return spmat;
		}

		void Extract_row_dense(int row, int num_col, float* vec) {

			for (int i = 0; i < num_col; i++) {
				vec[i] = 0;
			}

			for (int i = rowptr[row]; i < rowptr[row + 1]; i++) {
				int tmp = colind[i];
				vec[colind[i]] = elements[i];
			}

		}

	};
}