
#include <MuonGun/Distribution.h>

namespace I3MuonGun {

template <typename Archive>
void
Distribution::serialize(Archive &ar, unsigned version)
{
	ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
}

Distribution::~Distribution() {};

}

I3_SERIALIZABLE(I3MuonGun::Distribution);

