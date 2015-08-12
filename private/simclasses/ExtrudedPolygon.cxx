/** $Id$
 * @file
 * @author Jakob van Santen <jakob.van.santen@desy.de>
 *
 * $Revision$
 * $Date$
 */

#include <simclasses/ExtrudedPolygon.h>

namespace simclasses {

ExtrudedPolygon::~ExtrudedPolygon() {}

template <typename Archive>
void
ExtrudedPolygon::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("Base", base_object<Base>(*this));
}

}

I3_SERIALIZABLE(simclasses::ExtrudedPolygon);
