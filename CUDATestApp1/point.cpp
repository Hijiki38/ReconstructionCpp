#include "point.h"

namespace Reconstruction {

	float point::get_x() { return x; }
	float point::get_y() { return y; }
	float point::get_relative(float _num) { return _num - center + 0.5; }
	float point::get_inverted(float _num) { return 2 * center - _num - 1; }

	void point::set_x(float _x) { x = _x; }
	void point::set_y(float _y) { y = _y; }
	void point::set_xy(float _x, float _y) { x = _x; y = _y; }
	void point::set_center(float _center) { center = _center; }

	void point::rotate90(bool is_clockwise) {
		float tmp_x = x;
		
		x = get_inverted(y);
		y = tmp_x;


		//tmpx = origx;
		//origx = n_detector - origy - 1;
		//origy = tmpx;
		//relx = origx - center + 0.5;
		//rely = n_detector - origy - 1 - center + 0.5;
	}

}