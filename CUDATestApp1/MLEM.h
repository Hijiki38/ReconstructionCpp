#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "IterationRec.h"
#include "PCCTsinogram.h"
#include "point.h"
#include "geometry.h"
#include "SparseMatrix.h"
#include "methods.h"

extern const double PI;

namespace Reconstruction {

	class MLEM : public Reconstruction::IterationRec {

	public:

		MLEM(Reconstruction::PCCTsinogram* s, Reconstruction::geometry* geometry)
			: IterationRec(s, geometry) {};

		void calc_imgdiff(float* idiff, float* smr, float* atn, float sn, int size) const override;
		void calc_attenu(float* atn, float* idiff, float par, int size) const override;

	};
}

