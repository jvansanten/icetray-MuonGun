/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/histogram.h>
#include <boost/foreach.hpp>
#include <boost/python/slice.hpp>

namespace bp = boost::python;

static bp::object
to_dashi(boost::shared_ptr<I3MuonGun::histogram::histogram_base> h)
{
	bp::object histmodule = bp::import("dashi.histogram");
	bp::object numpy = bp::import("numpy");
	bp::list edges;
	BOOST_FOREACH(const std::vector<double> &dim, h->binedges())
		edges.append(bp::list(dim));
	
	bp::object dhist = histmodule.attr("histogram")(h->ndim(), edges);
	
	bp::list shape;
	size_t ndim = h->ndim();
	size_t tsize = 1;
	BOOST_FOREACH(size_t i, h->shape()) {
		tsize *= i;
		shape.append(i);
	}
	bp::tuple shapet(shape);
	bp::object dtype = numpy.attr("float64");
	
	bp::object buffer = bp::object(bp::handle<>(PyBuffer_FromMemory(
	    (void*)h->raw_bincontent(), sizeof(double)*tsize)));
	dhist.attr("_h_bincontent").attr("__setitem__")(bp::slice(),
	    numpy.attr("ndarray")(shapet, dtype, buffer));
	
	buffer =  bp::object(bp::handle<>(PyBuffer_FromMemory(
	    (void*)h->raw_squaredweights(), sizeof(double)*tsize)));
	dhist.attr("_h_squaredweights").attr("__setitem__")(bp::slice(),
	    numpy.attr("ndarray")(shapet, dtype, buffer));
	
	return dhist;
}

static boost::shared_ptr<I3MuonGun::histogram::histogram_base> test_histogram()
{
	using namespace I3MuonGun::histogram;
	
	histogram<3>::bin_specification edges;
	{
		using namespace binning;
		// edges[0] = uniform<>::create(0, 1, 11);
		edges[0] = uniform<cosine>::create(0, M_PI/2., 11);
		// edges[0] = boost::get<boost::shared_ptr<binning::scheme> >(edges[0])->edges();
		edges[1] = uniform<binning::log10>::create(1e3, 1e5, 11);
		edges[2] = uniform<power<2> >::create(1, 10, 11);
	}
	
	boost::shared_ptr<histogram<3> > h(new histogram<3>(edges));
	
	boost::array<double, 3> values = {{0.5, 1e4, 8}};
	
	for (int i=0; i < 25; i++) {
		values[0] = 1e-1*i;
		h->fill(values);
	}
	
	return h;
}


void register_histogram()
{
	using namespace I3MuonGun::histogram;
	namespace bp = boost::python;
	
	bp::def("test_histogram", &test_histogram);
	
	bp::class_<histogram_base, shared_ptr<histogram_base>, boost::noncopyable>("histogram", bp::no_init)
		.def("to_dashi", &to_dashi)
			;
	try {
		bp::import("dashi.histogram");
		bp::import("numpy");
	} catch (bp::error_already_set) {
		
		PyErr_Clear();
	}
	
}