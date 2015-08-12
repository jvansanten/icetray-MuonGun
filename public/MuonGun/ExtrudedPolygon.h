/** $Id$
 * @file
 * @author Jakob van Santen <jakob.van.santen@desy.de>
 *
 * $Revision$
 * $Date$
 */

#ifndef I3MUONGUN_EXTRUDEDPOLYGON_H_INCLUDED
#define I3MUONGUN_EXTRUDEDPOLYGON_H_INCLUDED

#include <MuonGun/Surface.h>
#include <dataclasses/I3Direction.h>

namespace I3MuonGun {

	class ExtrudedPolygon : public simclasses::Surface {
	public:
		ExtrudedPolygon(const std::vector<I3Position> &points, double padding=0.);

		virtual std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const;
		virtual bool operator==(const Surface&) const;
		
		std::vector<double> GetX() const;
		std::vector<double> GetY() const;
		std::vector<double> GetZ() const;
		
		struct vec2 {
			vec2(double xi, double yi);
			static vec2 from_I3Position(const I3Position &pos);
			static vec2 normalized(double xi, double yi);
			
			double x, y;
		};
	
	private:

		struct side {
			side(const vec2 &point, const vec2 &next_point);
			vec2 origin;
			vec2 vector;
			double length;
			I3Direction normal;
		};
		std::vector<side> sides_;
		std::pair<double, double> z_range_;
		double cap_area_;
		
		std::pair<double, double> GetDistanceToCaps(const I3Position &p, const I3Direction &dir) const;
		std::pair<double, double> GetDistanceToHull(const I3Position &p, const I3Direction &dir) const;

		bool PointInHull(const I3Position &p) const;
		double GetDistanceToCap(const I3Position &p, const I3Direction &dir, double cap_z) const;
	};
	
	I3_POINTER_TYPEDEFS(ExtrudedPolygon);

}

#endif // I3MUONGUN_EXTRUDEDPOLYGON_H_INCLUDED
