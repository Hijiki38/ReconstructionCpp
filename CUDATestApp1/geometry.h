#pragma once

namespace Reconstruction {
	struct geometry {
		bool is_conebeam = false;
		float sod = 0;
		float sdd = 0;
		float axis_correction = 0;
		float pixelsize = 0;
	};
}
