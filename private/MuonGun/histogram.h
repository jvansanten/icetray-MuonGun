
#ifndef I3MUONGUN_HISTOGRAM_H_INCLUDED
#define I3MUONGUN_HISTOGRAM_H_INCLUDED

#include <boost/function.hpp>
#include <boost/make_shared.hpp>
#include <boost/variant.hpp>

#ifdef NDEBUG
#define BOOST_DISABLE_ASSERTS
#endif
#include <boost/multi_array.hpp>
#include <boost/array.hpp>

// An n-dimensional histogram class, loosely inspired by dashi.histogram

namespace I3MuonGun {

namespace binning {

struct scheme {
	virtual const std::vector<double>& edges() const = 0;
	virtual const size_t index(double value) const = 0;
	
	virtual ~scheme() {};
};

// general case: arbitrarily-spaced edges, findable by binary search
class general : public scheme {
public:
	general(const std::vector<double> &edges)
	{
		if (edges.front() > -std::numeric_limits<double>::infinity())
			edges_.push_back(-std::numeric_limits<double>::infinity());
		std::copy(edges.begin(), edges.end(), std::back_inserter(edges_));
		if (edges.back() < std::numeric_limits<double>::infinity())
			edges_.push_back(std::numeric_limits<double>::infinity());
	}
	
	const std::vector<double>& edges() const
	{ return edges_; }
	
	const size_t index(double value) const
	{
		size_t j = std::distance(edges_.begin(),
		    std::upper_bound(edges_.begin(),
		    edges_.end(), value));
		assert(j > 0);
		return j-1;
	}
private:
	std::vector<double> edges_;
};

struct identity {
	static inline double map(double v) { return v; }
	static inline double imap(double v) { return v; }
};

struct log10 {
	static inline double map(double v) { return std::pow(10, v); }
	static inline double imap(double v) { return std::log10(v); }
};

template <int N>
struct power {
	static inline double map(double v) { return std::pow(v, N); }
	static inline double imap(double v) { return std::pow(v, 1./N); }
};

template <>
struct power<2> {
	static inline double map(double v) { return v*v; }
	static inline double imap(double v) { return std::sqrt(v); }
};


// Optimal case: bin edges uniform under some transformation
// between set limits
template <typename Transformation = identity >
class uniform : public scheme {
public:
	uniform(double low, double high, size_t nsteps)
	    : offset_(Transformation::imap(low)),
	    range_(Transformation::imap(high)-Transformation::imap(low)),
	    min_(map(0)), max_(map(1)), nsteps_(nsteps)
	{
		edges_.reserve(nsteps+2);
		edges_.push_back(-std::numeric_limits<double>::infinity());
		for (size_t i = 0; i < nsteps; i++)
			edges_.push_back(map(i/double(nsteps_-1)));
		edges_.push_back(std::numeric_limits<double>::infinity());
	}
	
	const std::vector<double>& edges() const
	{ return edges_; }
	
	const size_t index(double value) const
	{
		if (value < min_)
			return 0;
		else if (value >= max_)
			return edges_.size()-2;
		else {
			return size_t(floor((nsteps_-1)*imap(value)))+1;
		}
	}
private:
	inline double map(double value) const
	{
		return Transformation::map(range_*value + offset_);
	}
	inline double imap(double value) const
	{
		return (Transformation::imap(value)-offset_)/range_;
	}
	
	double offset_, range_, min_, max_;
	size_t nsteps_;
	std::vector<double> edges_;
};

};

template <size_t N, typename T = double>
class histogram {
public:
	typedef boost::array<boost::variant< std::vector<double>, boost::shared_ptr<binning::scheme> >, N> bin_specification;
public:
	// Construct with non-uniform bins in all dimensions
	histogram(const boost::array<std::vector<double>, N> &edges)
	{
		for (size_t i=0; i < N; i++)
			binners_[i] = boost::make_shared<binning::general>(edges[i]);
	
		make_datacube();
	}
	
	// Construct with uniform bins in all dimensions
	histogram(const boost::array<boost::shared_ptr<binning::scheme>, N> &schemes)
	    : binners_(schemes)
	{
		make_datacube();
	}
	
	// Construct with a mix of uniform and non-uniform bins in different dimensions
	histogram(const boost::array<boost::variant< std::vector<double>, boost::shared_ptr<binning::scheme> >, N> &schemes)
	{
		for (size_t i=0; i < N; i++)
			binners_[i] = boost::apply_visitor(bin_visitor(), schemes[i]);
		make_datacube();
	}
	
	void fill(const boost::array<double, N> &values, T weight=1)
	{
		boost::array<typename datacube_type::index, N> idx;
		for (size_t i=0; i < N; i++) {
			if (std::isnan(values[i]))
				return;
			idx[i] = binners_[i]->index(values[i]);
		}
		
		bincontent_(idx) += weight;
		squaredweights_(idx) += weight*weight;
	}
	
	const boost::array<std::vector<double>, N> & edges() const { return edges_; }
	
private:
	void make_datacube()
	{
		boost::array<size_t, N> dims;
		for (size_t i=0; i < N; i++) {
			edges_[i] = binners_[i]->edges();
			dims[i] = edges_[i].size()-1u;			
		}
		
		bincontent_.resize(dims);
		squaredweights_.resize(dims);
		
		std::fill(bincontent_.data(), bincontent_.data()+bincontent_.size(), 0);
		std::fill(squaredweights_.data(), squaredweights_.data()+squaredweights_.size(), 0);
	}
	
	typedef boost::shared_ptr<binning::scheme > scheme_ptr;
	struct bin_visitor : public boost::static_visitor<scheme_ptr> {
		scheme_ptr operator()(const scheme_ptr & v) const
		{ return v; }
		
		scheme_ptr operator()(const std::vector<double> & v) const
		{ return boost::make_shared<binning::general>(v); }
	};

private:
	boost::array<scheme_ptr, N> binners_;
	boost::array<std::vector<double>, N> edges_;
	
	typedef boost::multi_array<T, N> datacube_type;
	datacube_type bincontent_;
	datacube_type squaredweights_;

};

}

#endif
