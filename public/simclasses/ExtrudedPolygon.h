/** $Id$
 * @file
 * @author Jakob van Santen <jakob.van.santen@desy.de>
 *
 * $Revision$
 * $Date$
 */

#ifndef SIMCLASSES_EXTRUDEDPOLYGON_H_INCLUDED
#define SIMCLASSES_EXTRUDEDPOLYGON_H_INCLUDED

#include <simclasses/SamplingSurface.h>
#include <simclasses/detail/ExtrudedPolygonBase.h>

namespace simclasses {

class ExtrudedPolygon : public ExtrudedPolygonBase<SamplingSurface> {
private:
	typedef ExtrudedPolygonBase<SamplingSurface> Base;
public:
	ExtrudedPolygon(const std::vector<I3Position> &points, double padding=0.) : Base(points, padding) {};
	~ExtrudedPolygon();
private:
	ExtrudedPolygon() {}
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(ExtrudedPolygon);

}

BOOST_CLASS_VERSION(simclasses::ExtrudedPolygon, 0);

#endif // SIMCLASSES_EXTRUDEDPOLYGON_H_INCLUDED
