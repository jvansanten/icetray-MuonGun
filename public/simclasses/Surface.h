/** $Id: Surface.h 135924 2015-08-06 14:17:50Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 135924 $
 * $Date: 2015-08-06 16:17:50 +0200 (Thu, 06 Aug 2015) $
 */

#ifndef SIMCLASSES_SURFACE_H_INCLUDED
#define SIMCLASSES_SURFACE_H_INCLUDED

#include <icetray/I3PointerTypedefs.h>
#include <icetray/serialization.h>

class I3Direction;
class I3Position;

namespace simclasses {

/**
 * @brief A closed surface
 *
 * Surface knows how to find the intersections of a ray with itself.
 */
class Surface {
public:
	virtual ~Surface();
	/**
	 * Find the points where a ray intersects the surface
	 *
	 * @param[in] p   The origin of the ray
	 * @param[in] dir The direction of the ray
	 * @returns a pair of distances from the origin of the ray to the
	 *          intersections with the surface. A distance of NAN means
	 *          that the ray never intersects the surface, and negative
	 *          distances mean that the intersection in question is
	 *          "behind" the origin.
	 */
	virtual std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const = 0;
	std::pair<double, double> no_intersection() const
	{
		return std::make_pair(NAN, NAN);
	}

private:
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(Surface);

}

BOOST_CLASS_VERSION(simclasses::Surface, 0);

#endif
