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
		std::shared_ptr<float[]> elements;	//all nonzero values
		std::shared_ptr<int[]> rowptr;		//indices of the first nonzero element in each row
		std::shared_ptr<int[]> colind;		//the column indices of the corresponding elements
		int nonzero;		//the number of nonzero elements

		//float* elements;	//all nonzero values
		//int* rowptr;		//indices of the first nonzero element in each row
		//int* colind;		//the column indices of the corresponding elements
		//int nonzero;		//the number of nonzero elements

	public:
		SparseMatrix(const std::shared_ptr<float[]>& elem, const std::shared_ptr<int[]>& rptr, const std::shared_ptr<int[]>& cind, int nz) {
			elements = elem;
			rowptr = rptr;
			colind = cind;
			nonzero = nz;
		}

		//SparseMatrix(float* elem, int* rptr, int* cind, int nz) {
		//	elements = elem;
		//	rowptr = rptr;
		//	colind = cind;
		//	nonzero = nz;
		//}

		//~SparseMatrix() {

		//	std::cout << "uoo\n";

		//	free(elements);
		//	free(rowptr);
		//	free(colind);
		//}

		std::shared_ptr<float[]> get_elements() { return elements; }
		std::shared_ptr<int[]> get_rowptr() { return rowptr; }
		std::shared_ptr<int[]> get_colind() { return colind; }
		int get_nonzero() { return nonzero; }

		std::unique_ptr<SparseMatrix> Create_blockmat(int begin, int rows) {

			int elem_num = rowptr.get()[begin + rows] - rowptr.get()[begin];

			//float* elements_new = (float*)malloc(elem_num * sizeof(float));
			//int* rowptr_new = (int*)malloc((rows + 1) * sizeof(int));
			//int* colind_new = (int*)malloc(elem_num * sizeof(float));

			std::unique_ptr<float[]> elements_new = std::make_unique<float[]>(elem_num);	//all nonzero values
			std::unique_ptr<int[]> rowptr_new = std::make_unique<int[]>((rows + 1));		//indices of the first nonzero element in each row
			std::unique_ptr<int[]> colind_new = std::make_unique<int[]>(elem_num);		//the column indices of the corresponding elements

			//std::cout << "\nstart generating blockmat\n";

			for (int i = 0; i < rows; i++) {
				rowptr_new.get()[i] = rowptr.get()[begin + i] - rowptr.get()[begin];
				//std::cout << rowptr_new[i] << std::endl;
			}
			rowptr_new.get()[rows] = elem_num;

			for (int i = 0; i < elem_num; i++) {
				elements_new.get()[i] = elements.get()[i + rowptr.get()[begin]];
				colind_new.get()[i] = colind.get()[i + rowptr.get()[begin]];
				//std:cout << "new colind: " << colind_new[i] << std::endl;
			}

			std::unique_ptr<SparseMatrix> spmat(new SparseMatrix(std::move(elements_new), std::move(rowptr_new), std::move(colind_new), std::move(elem_num)));

			return spmat;
		}

		void Extract_row_dense(int row, int num_col, float* vec) {

			for (int i = 0; i < num_col; i++) {
				vec[i] = 0;
			}

			for (int i = rowptr.get()[row]; i < rowptr.get()[row + 1]; i++) {
				int tmp = colind.get()[i];
				vec[tmp] = elements.get()[i];
			}

		}

	};
}