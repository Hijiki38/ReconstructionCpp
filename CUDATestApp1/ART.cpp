#include "ART.h"

using std::cout;
using namespace Eigen;

namespace Reconstruction {

	void ART::calc_imgdiff(float* idiff, float* smr, float* atn, float sn, int size) const {
		float sys_atn = Reconstruction::dot_array(smr, atn, size);
		float sys_sys = Reconstruction::dot_array(smr, smr, size);

		Reconstruction::mul_array1(smr, ((sys_atn - sn) / sys_sys), size);
		Reconstruction::add_array(idiff, smr, size);
	}

	void ART::calc_attenu(float* atn, float* idiff, int size) const {
		Reconstruction::mul_array1(idiff, (relpar / block_size), size);
		Reconstruction::sub_array(atn, idiff, size); //ART:  attenu = attenu - (relpar / block_num) * imgdiff;
	}
}

