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

	void point::rotate90() {

		float tmp_x = x;
		x = get_inverted(y);
		y = tmp_x;
	}

	void point::rotate90(int rotatecount) {

		float tmp_x = x;

		switch (rotatecount % 4) {
		case 0:
			break;
		case 1:
			x = get_inverted(y);
			y = tmp_x;
			break;
		case 2:
			x = get_inverted(x);
			y = get_inverted(y);
			break;
		case 3:
			x = get_inverted(y);
			y = get_inverted(tmp_x);
			break;
		}
	}

}