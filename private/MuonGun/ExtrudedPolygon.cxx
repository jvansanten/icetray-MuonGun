/** $Id$
 * @file
 * @author Jakob van Santen <jakob.van.santen@desy.de>
 *
 * $Revision$
 * $Date$
 */

#include <MuonGun/ExtrudedPolygon.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>
#include <phys-services/I3RandomService.h>

#include <boost/foreach.hpp>
#include <cmath>

namespace {

typedef simclasses::polygon::vec2 vec2;

/// A counterclockwise curve is the basic building block of a convex hull
class ccw_curve : public std::vector<vec2> {
public:
	// Add a point to the curve
	void operator()(const vec2 &p)
	{
		// Remove points until the curve will be counterclockwise
		while (size() >= 2 && !ccw((*this)[size()-2], (*this)[size()-1], p))
			pop_back();
		push_back(p);
	}
private:
	static bool
	ccw(const vec2 &o, const vec2 &a, const vec2 &b)
	{
		// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
		// positive, if OAB makes a counter-clockwise turn,
		// negative for clockwise turn, and zero if the points are collinear.
		return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x) > 0;
	}
};

/// @brief Computes the convex hull of a set of 2D points.
/// Input: an iterable sequence of (x, y) pairs representing the points.
/// Output: a list of vertices of the convex hull in counter-clockwise order,
///   starting from the vertex with the lexicographically smallest coordinates.
/// Implements Andrew's monotone chain algorithm. O(n log n) complexity.
/// 
/// Lifted from http://code.icecube.wisc.edu/svn/sandbox/ckopper/eventinjector/python/util/__init__.py
std::vector<vec2>
convex_hull(const std::vector<I3Position> &positions)
{
	std::vector<vec2> hull;
	
	// Build a set of unique points, sorted lexicographically
	std::set<vec2> points;
	std::transform(positions.begin(), positions.end(),
	    std::inserter(points, points.end()), vec2::from_I3Position);
	
	// Boring case: 1 point (perhaps repeated)
	if (points.size() <= 1) {
		std::copy(points.begin(), points.end(), std::back_inserter(hull));
		return hull;
	}
	
	// Build lower and upper hulls
	std::vector<vec2> lower = std::for_each(points.begin(), points.end(), ccw_curve());
	std::vector<vec2> upper = std::for_each(points.rbegin(), points.rend(), ccw_curve());
	
	// Concatenation of the lower and upper hulls gives the convex hull.
	// Last point of each list is omitted because it is repeated at the
	// beginning of the other list.
	std::copy(lower.begin(), lower.end()-1, std::back_inserter(hull));
	std::copy(upper.begin(), upper.end()-1, std::back_inserter(hull));
	
	return hull;
}


/// @brief Expand the x-y footprint by moving each edge out by a distance
/// *padding*.
/// 
/// A convex polygon can be offset by moving each vertex parallel to the
/// edges by a distance that is inversely proportional to the sine of the
/// counterclockwise angle between the edges that meet at each vertex.
/// This breaks down for edges that are [anti]parallel or, but neither
/// case should occur for maximally simplified polygons.
std::vector<vec2>
expand_polygon(const std::vector<vec2> &hull, double padding)
{
	std::vector<vec2> points;
	for (std::vector<vec2>::const_iterator p = hull.begin(); p != hull.end(); p++) {
		std::vector<vec2>::const_iterator next = boost::next(p);
		if (next == hull.end())
			next = hull.begin();
		std::vector<vec2>::const_iterator prev = boost::prior(
		    p == hull.begin() ? hull.end() : p);
		// normalized vector connecting this vertex to the next one
		vec2 d = vec2::normalized(next->x-p->x, next->y-p->y);
		// and the previous vertex to this one
		vec2 prev_d = vec2::normalized(p->x-prev->x, p->y-prev->y);
		// sine of the inner angle between the segments that meet here
		double det = prev_d.x*d.y - prev_d.y*d.x;
		if (det == 0.)
			log_fatal("Edges can't be [anti]parallel");
		vec2 outwards(prev_d.x-d.x, prev_d.y-d.y);
		points.push_back(vec2(p->x + outwards.x*padding/det, p->y + outwards.y*padding/det));
	}
	
	return points;
}

std::pair<double, double>
z_range(const std::vector<I3Position> &positions)
{
	std::pair<double, double> range(std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
	BOOST_FOREACH(const I3Position &pos, positions) {
		if (pos.GetZ() < range.first)
			range.first = pos.GetZ();
		if (pos.GetZ() > range.second)
			range.second = pos.GetZ();
	}
	
	return range;
}

std::pair<double, double>
make_ordered_pair(double a, double b)
{
	if (std::isnan(a) || (a > b))
		return std::make_pair(b, a);
	else
		return std::make_pair(a, b);
}

}

namespace simclasses {

namespace polygon {

vec2::vec2(double xi, double yi) : x(xi), y(yi)
{}

template <typename Archive>
void vec2::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("X", x);
	ar & make_nvp("Y", y);
}

vec2
vec2::from_I3Position(const I3Position &p)
{
	return vec2(p.GetX(), p.GetY());
}

vec2
vec2::normalized(double xi, double yi)
{
	double l = hypot(xi, yi);
	return vec2(xi/l, yi/l);
}

bool
operator<(const vec2 &a, const vec2 &b)
{
	if (a.x < b.x)
		return true;
	else if (a.x > b.x)
		return false;
	else if (a.y < b.y)
		return true;
	else
		return false;
}

side::side(const vec2 &p, const vec2 &np) : origin(p),
    vector(np.x-p.x, np.y-p.y), length(hypot(vector.x, vector.y)),
	normal(vector.y/length, -vector.x/length, 0.)
{}

}

template <typename Base>
ExtrudedPolygonBase<Base>::ExtrudedPolygonBase(const std::vector<I3Position> &points, double padding)
{
	using simclasses::polygon::vec2;
	
	std::pair<double, double> zrange = z_range(points);
	std::vector<vec2> hull = convex_hull(points);
	if (padding != 0) {
		hull = expand_polygon(hull, padding);
		zrange.first -= padding;
		zrange.second += padding;
	}
	initWithHull(hull, zrange);
}

template <typename Base>
void
ExtrudedPolygonBase<Base>::initWithHull(const std::vector<polygon::vec2> &hull,
    const std::pair<double, double> &z_range)
{
	using simclasses::polygon::vec2;
	using simclasses::polygon::side;
	
	z_range_ = z_range;
	cap_area_ = 0;
	sides_.clear();
	for (std::vector<vec2>::const_iterator p = hull.begin(); p != hull.end(); p++) {
		std::vector<vec2>::const_iterator np = boost::next(p);
		if (np == hull.end())
			np = hull.begin();
		
		sides_.push_back(side(*p, *np));
		
		// area of a simple polygon in the x-y plane
		cap_area_ += (p->x*np->y - np->x*p->y);
	}
	cap_area_ /= 2.;
}

template <typename Base>
std::vector<double>
ExtrudedPolygonBase<Base>::GetX() const
{
	std::vector<double> x;
	BOOST_FOREACH(const polygon::side &sidey, sides_)
		x.push_back(sidey.origin.x);
	return x;
}

template <typename Base>
std::vector<double>
ExtrudedPolygonBase<Base>::GetY() const
{
	std::vector<double> y;
	BOOST_FOREACH(const polygon::side &sidey, sides_)
		y.push_back(sidey.origin.y);
	return y;
}

template <typename Base>
std::vector<double>
ExtrudedPolygonBase<Base>::GetZ() const
{
	std::vector<double> z;
	z.push_back(z_range_.first);
	z.push_back(z_range_.second);
	
	return z;
}

/// Calculate the most extreme displacements to points on the 2D hull
template <typename Base>
std::pair<double, double>
ExtrudedPolygonBase<Base>::GetDistanceToHull(const I3Position &pos, const I3Direction &dir) const
{
	std::pair<double, double> offsets = Surface::no_intersection();
	
	if (dir.GetX() + dir.GetY() == 0.)
		log_fatal("Direction must have a horizontal component");
	
	BOOST_FOREACH(const polygon::side &sidey, sides_) {
		// Components of the vector connecting the test point to the
		// origin of the line segment;
		double x = sidey.origin.x - pos.GetX();
		double y = sidey.origin.y - pos.GetY();
		
		// Proportional distance along the line segment to the
		// intersection point
		double alpha = (dir.GetX()*y - dir.GetY()*x)
		    / (dir.GetY()*sidey.vector.x - dir.GetX()*sidey.vector.y);
		
		// Is there an intersection?
		if ((alpha >= 0.) && (alpha < 1.)) {
			// Distance along ray to intersection point
			double beta = dir.GetX() != 0. ?
			    (x + alpha*sidey.vector.x)/dir.GetX() :
			    (y + alpha*sidey.vector.y)/dir.GetY();
			// NB: reversed comparison is equivalent to
			// (isnan(offsets.first) || beta < offsets.first)
			if (!(beta >= offsets.first))
				offsets.first = beta;
			if (!(beta <= offsets.second))
				offsets.second = beta;
		}
	}
	double sin_theta = std::sqrt(1.-dir.GetZ()*dir.GetZ());
	offsets.first /= sin_theta;
	offsets.second /= sin_theta;
	
	return offsets;
}

/// Test whether point is inside the 2D hull by ray casting
template <typename Base>
bool
ExtrudedPolygonBase<Base>::PointInHull(const I3Position &pos) const
{
	int crossings = 0;
	for (std::vector<polygon::side>::const_iterator p = sides_.begin(); p != sides_.end(); p++) {
		std::vector<polygon::side>::const_iterator np = boost::next(p);
		if (np == sides_.end())
			np = sides_.begin();
		
		// only consider segments whose y range spans the current point
		if (!(((p->origin.y > pos.GetY()) && (np->origin.y <= pos.GetY())) ||
		      ((p->origin.y <= pos.GetY()) && (np->origin.y > pos.GetY()))))
			continue;
		// find the intersection of the segment with a horizontal line
		double xc = p->origin.x + (pos.GetY() - p->origin.y)*p->vector.x/p->vector.y;
		// is the crossing to the right of the test point?
		if (pos.GetX() < xc)
			crossings++;
	}
	
	// the point is inside iff the ray crosses an odd number of times
	return (crossings % 2) == 1;
}

template <typename Base>
double
ExtrudedPolygonBase<Base>::GetDistanceToCap(const I3Position &p, const I3Direction &dir, double cap_z) const
{
	return (cap_z-p.GetZ())/dir.GetZ();
}

template <typename Base>
std::pair<double, double>
ExtrudedPolygonBase<Base>::GetDistanceToCaps(const I3Position &p, const I3Direction &dir) const
{
	return make_ordered_pair(GetDistanceToCap(p, dir, z_range_.first), GetDistanceToCap(p, dir, z_range_.second));
}

template <typename Base>
std::pair<double, double>
ExtrudedPolygonBase<Base>::GetIntersection(const I3Position &p, const I3Direction &dir) const
{
	if (std::abs(dir.GetZ()) == 1.) {
		// perfectly vertical track: only check intersections with caps
		if (!PointInHull(p))
			return Surface::no_intersection();
		else
			return GetDistanceToCaps(p, dir);
	} else if (dir.GetZ() == 0.) {
		// perfectly horizontal track: only check intersections with sides
		if (p.GetZ() < z_range_.first || p.GetZ() > z_range_.second)
			return Surface::no_intersection();
		else
			return GetDistanceToHull(p, dir);
	} else {
		// general case: both rho and z components nonzero
		std::pair<double, double> sides = GetDistanceToHull(p, dir);
		std::pair<double, double> caps = GetDistanceToCaps(p, dir);
		
		if (caps.first >= sides.second || caps.second <= sides.first)
			return Surface::no_intersection();
		else {
			return std::make_pair(std::max(sides.first, caps.first), std::min(sides.second, caps.second));
		}
	}
}

template <typename Base>
double
ExtrudedPolygonBase<Base>::GetArea(const I3Direction &dir) const
{
	double area = 0;
	BOOST_FOREACH(const polygon::side &sidey, sides_) {
		double inner = sidey.normal*dir;
		if (inner < 0)
			area += -inner*sidey.length;
	}
	area *= (z_range_.second-z_range_.first);
	area += std::abs(dir.GetZ())*cap_area_;
	
	return area;
}

template <typename Base>
double
ExtrudedPolygonBase<Base>::GetMaximumArea() const
{
	double side_area = 0;
	BOOST_FOREACH(const polygon::side &sidey, sides_) {
		side_area += sidey.length;
	}
	// the largest possible projected area occurs for a flat square
	side_area *= (z_range_.second-z_range_.first)/2.;
	double ct_max = cos(atan(side_area/cap_area_));
	
	return cap_area_*fabs(ct_max) + side_area*sqrt(1.-ct_max*ct_max);
}

template <typename Base>
double
ExtrudedPolygonBase<Base>::GetAverageSideArea() const
{
	// the projected area of a plane, averaged over a 2\pi rotation that
	// passes through the normal, is
	// A*\int_0^\pi \Theta(\sin\alpha)\sin\alpha d\alpha / 2\pi = A/\pi
	double area = 0;
	BOOST_FOREACH(const polygon::side &sidey, sides_)
		area += sidey.length;
	area *= GetLength()/M_PI;
	
	return area;
}

template <typename Base>
I3Position
ExtrudedPolygonBase<Base>::SampleImpactPosition(const I3Direction &dir, I3RandomService &rng) const
{
	// first, pick which face it's going to hit
	double area = 0;
	double height = GetLength();
	std::vector<double> prob;
	std::pair<double, double> x_range(std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
	std::pair<double, double> y_range(std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
	
	BOOST_FOREACH(const polygon::side &sidey, sides_) {
		double inner = sidey.normal*dir;
		if (inner < 0)
			area += -inner*sidey.length*height;
		prob.push_back(area);
		// Find the bounding box of the cap
		if (sidey.origin.x < x_range.first)
			x_range.first = sidey.origin.x;
		if (sidey.origin.x > x_range.second)
			x_range.second = sidey.origin.x;
		if (sidey.origin.y < y_range.first)
			y_range.first = sidey.origin.y;
		if (sidey.origin.y > y_range.second)
			y_range.second = sidey.origin.y;
	}
	area += std::abs(dir.GetZ())*cap_area_;
	prob.push_back(area);
	std::vector<double>::iterator target =
	    std::lower_bound(prob.begin(), prob.end(), rng.Uniform(0, area));
	if (target == boost::prior(prob.end())) {
		// top or bottom face
		// triangulation would be more efficient here, but also more complicated
		I3Position pos(NAN,NAN,dir.GetZ() > 0 ? z_range_.first : z_range_.second);
		do {
			pos.SetX(rng.Uniform(x_range.first, x_range.second));
			pos.SetY(rng.Uniform(y_range.first, y_range.second));
		} while (!PointInHull(pos));
		return pos;
	} else {
		// side face
		std::vector<polygon::side>::const_iterator sidey =
		    sides_.begin() + std::distance(prob.begin(), target);
		double horizontal = rng.Uniform();
		double vertical = rng.Uniform();
		return I3Position(
		    sidey->origin.x + horizontal*sidey->vector.x,
		    sidey->origin.y + horizontal*sidey->vector.y,
		    z_range_.first + vertical*height
		);
	}
}

template <typename Base>
template <typename Archive>
void
ExtrudedPolygonBase<Base>::save(Archive &ar, unsigned version) const
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("Base", base_object<Base>(*this));
	std::vector<polygon::vec2> hull;
	BOOST_FOREACH(const polygon::side &sidey, sides_) {
		hull.push_back(sidey.origin);
	}
	ar & make_nvp("HullXY", hull);
	ar & make_nvp("HullZ", z_range_);
}

template <typename Base>
template <typename Archive>
void
ExtrudedPolygonBase<Base>::load(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("Base", base_object<Base>(*this));
	std::vector<polygon::vec2> hull;
	std::pair<double, double> z_range;
	ar & make_nvp("HullXY", hull);
	ar & make_nvp("HullZ", z_range);
	initWithHull(hull, z_range);
}

template <typename Archive>
void
ExtrudedPolygon::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("Base", base_object<Base>(*this));
}

}

namespace I3MuonGun {

ExtrudedPolygon::~ExtrudedPolygon() {}

template <typename Archive>
void
ExtrudedPolygon::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("Base", base_object<Base>(*this));
}

#if 0
bool
ExtrudedPolygon::operator==(const Surface &s) const
{
	const ExtrudedPolygon *other = dynamic_cast<const ExtrudedPolygon*>(&s);
	if (!other)
		return false;
	else {
		if (z_range_ != other->z_range_ || sides_.size() != other->sides_.size())
			return false;
		typedef std::vector<side>::const_iterator side_iter;
		for (side_iter a=sides_.begin(), b=other->sides_.begin(); a != sides_.end(); a++, b++) {
			if (a->origin.x != b->origin.x || a->origin.y != b->origin.y)
				return false;
		}
	}

	return true;
}


#endif

}

template class simclasses::ExtrudedPolygonBase<simclasses::SamplingSurface>;
template class simclasses::ExtrudedPolygonBase<I3MuonGun::SamplingSurface>;
template class I3MuonGun::detail::UprightSurface<simclasses::ExtrudedPolygonBase<I3MuonGun::SamplingSurface> >;

I3_SERIALIZABLE(simclasses::ExtrudedPolygon);
I3_SERIALIZABLE(I3MuonGun::ExtrudedPolygon);

