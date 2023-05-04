/* MPL License: https://github.com/nickcastel50/conway-geom/blob/typescript_api/LICENSE.md */
//Curve Implementation of a Curve

#pragma once
#include <vector>
#include <glm/glm.hpp>

namespace webifc::geometry {

	struct IfcCurve
	{
		std::vector<glm::dvec3> points;
		std::vector<uint16_t> indices;
		void Add(glm::dvec3 pt);
		void Add(glm::dvec2 pt);
		glm::dvec2 Get2d(size_t i) const;
		glm::dvec3 Get3d(size_t i) const;
		void Invert();
		bool IsCCW() const;

		private:
			static constexpr double EPS_TINY = 1e-9;
	};

}
