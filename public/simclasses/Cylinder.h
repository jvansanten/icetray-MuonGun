/** $Id$
 * @file
 * @author Jakob van Santen <jakob.van.santen@desy.de>
 *
 * $Revision$
 * $Date$
 */

#ifndef SIMCLASSES_CYLINDER_H_INCLUDED
#define SIMCLASSES_CYLINDER_H_INCLUDED

#include <simclasses/SamplingSurface.h>
#include <simclasses/detail/CylinderBase.h>

namespace simclasses {

/**
 * @brief A cylinder aligned with the z axis
 */
class Cylinder : public detail::CylinderBase<SamplingSurface> {
private:
	typedef detail::CylinderBase<SamplingSurface> Base;
public:
	virtual ~Cylinder();
	Cylinder(double length, double radius, I3Position center=I3Position(0,0,0)) : CylinderBase(length, radius, center) {};

private:
	Cylinder() {}
	
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(Cylinder);

}

BOOST_CLASS_VERSION(simclasses::Cylinder, 0);

#endif // SIMCLASSES_CYLINDER_H_INCLUDED