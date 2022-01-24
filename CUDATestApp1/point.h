#pragma once


namespace Reconstruction {
	//template <class float> class point {
	class point{
	private:
		float x;
		float y;
		float center;
	public:
		point(float x, float y, float center) {
			this->x = x;
			this->y = y;
			this->center = center;
		}
		float get_x();
		float get_y();
		float get_relative(float _num);
		float get_inverted(float _num);
		void set_x(float _x);
		void set_y(float _y);
		void set_xy(float _x, float _y);
		void set_center(float _center);
		void rotate90(bool is_clockwise);
	};
}