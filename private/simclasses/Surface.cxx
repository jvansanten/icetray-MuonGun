/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <simclasses/Surface.h>

namespace simclasses {

Surface::~Surface() {}

template <typename Archive>
void
Surface::serialize(Archive &ar __attribute__ ((unused)), unsigned version __attribute__ ((unused)))
{}

}

I3_SERIALIZABLE(simclasses::Surface);
