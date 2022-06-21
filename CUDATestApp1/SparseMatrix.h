#pragma once

#include <stdio.h>
#include <iostream>
#include <string>
#include "Eigen/Sparse"

//using namespace std;
using namespace Eigen;

namespace Reconstruction {
	class SparseMatrix {
	private:
		//CSR format
		std::unique_ptr<float[]> elements;	//all nonzero values
		std::unique_ptr<int[]> rowptr;		//indices of the first nonzero element in each row
		std::unique_ptr<int[]> colind;		//the column indices of the corresponding elements
		int nonzero;		//the number of nonzero elements

	public:
		SparseMatrix(std::unique_ptr<float[]>& elem, std::unique_ptr<int[]>& rptr, std::unique_ptr<int[]>& cind, int nz)
			: elements(move(elem))
			, rowptr(move(rptr))
			, colind(move(cind))
			, nonzero(nz){}


		~SparseMatrix(){}

		std::unique_ptr<SparseMatrix> Create_blockmat(int begin, int rows) {

			int elem_num = rowptr.get()[begin + rows] - rowptr.get()[begin];

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

			std::unique_ptr<SparseMatrix> spmat(new SparseMatrix(elements_new, rowptr_new, colind_new, elem_num));

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