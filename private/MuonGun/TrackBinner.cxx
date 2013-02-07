/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/TrackBinner.h>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <icetray/I3Units.h>

namespace I3MuonGun {

TrackBinner::TrackBinner(double mindepth, double maxdepth, unsigned steps)
{	
	std::vector<double> multbins, rbins;
	{
		using namespace boost::assign;
		multbins += 0.5, 1.5, 2.5, 3.5, 9.5, 19.5, 39.5, 99.5;
		rbins += 0, 5, 10, 15, 25, 45;
	}
	
	using namespace histogram::binning;
	histogram::histogram<5>::bin_specification mspecs;
	mspecs[0] = uniform<cosine>::create(0, M_PI/2, 11);
	double dh = (maxdepth-mindepth)/(2*(steps-1));
	mspecs[1] = uniform<>::create(mindepth-dh, maxdepth-dh, steps);
	mspecs[2] = multbins;
	mspecs[3] = rbins;
	mspecs[4] = uniform<histogram::binning::log10>::create(1, 1e7, 101);
	
	{
		histogram::histogram<2>::bin_specification specs = {{mspecs[0], uniform<histogram::binning::log10>::create(1e2, 1e11, 101)}};
		primary_ = boost::make_shared<histogram::histogram<2> >(specs);
	}
	
	{
		histogram::histogram<3>::bin_specification specs = {{
		    uniform<cosine>::create(0, M_PI/2, 101), mspecs[1], uniform<>::create(0.5, 99.5, 100)}};
		multiplicity_ = boost::make_shared<histogram::histogram<3> >(specs);
	}
	
	{
		histogram::histogram<4>::bin_specification specs = {{mspecs[0], mspecs[1], mspecs[2],
			uniform<power<2> >::create(0, 250, 101)}};
		radius_ = boost::make_shared<histogram::histogram<4> >(specs);
	}
	
	energy_ = boost::make_shared<histogram::histogram<5> >(mspecs);
	
}

void TrackBinner::Consume(boost::shared_ptr<const TrackBundle> tracks,
    double energy, double zenith, double weight)
{
	{
		boost::array<double, 2> values = {{zenith, energy}};
		primary_->fill(values, weight);

	}
	for (std::map<double, std::vector<CompactTrack> >::const_iterator
	    i = tracks->begin(); i != tracks->end(); i++) {

		double mult = i->second.size();
		boost::array<double, 3> mult_values = {{zenith, i->first/I3Units::km, mult}};
		multiplicity_->fill(mult_values, weight);
		boost::array<double, 4> radius_values = {{zenith, i->first/I3Units::km, mult, 0}};
		boost::array<double, 5> energy_values = {{zenith, i->first/I3Units::km, mult, 0, 0}};
		BOOST_FOREACH(const CompactTrack &track, i->second) {
			radius_values[3] = energy_values[3] = track.GetRadius();
			energy_values[4] = track.GetEnergy();
			radius_->fill(radius_values, weight/mult);
			energy_->fill(energy_values, weight/mult);
		}
	}
}



}