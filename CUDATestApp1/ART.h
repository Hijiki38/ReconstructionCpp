#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "IterationRec.h"
#include "PCCTsinogram.h"
#include "TV.h"
#include "point.h"
#include "geometry.h"
#include "SparseMatrix.h"
#include "methods.h"

extern const double PI;

namespace Reconstruction {

	class ART : public Reconstruction::IterationRec {

	public:

		ART(Reconstruction::PCCTsinogram* s, Reconstruction::geometry* geometry, int blocksize)
			: IterationRec(s, geometry, 1.05, blocksize) {};

		void calc_imgdiff(float* idiff, float* smr, float* atn, float sn, int size) const override;
		void calc_attenu(float* atn, float* idiff, int size) const override;

	};
}

