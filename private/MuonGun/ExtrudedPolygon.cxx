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

#include <boost/foreach.hpp>

namespace {

typedef I3MuonGun::ExtrudedPolygon::vec2 vec2;

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

// define an absolute ordering for NaN
bool
nan_less(double a, double b)
{
	if (std::isnan(b))
		return true;
	else if (std::isnan(a))
		return false;
	else
		return a < b;
}

bool
nan_greater(double a, double b)
{
	if (std::isnan(b))
		return true;
	else if (std::isnan(a))
		return false;
	else
		return a > b;
}

std::pair<double, double>
make_ordered_pair(double a, double b)
{
	if (!(a > b))
		return std::make_pair(a, b);
	else
		return std::make_pair(b, a);
}

}

namespace I3MuonGun {

ExtrudedPolygon::vec2::vec2(double xi, double yi) : x(xi), y(yi)
{}

ExtrudedPolygon::vec2
ExtrudedPolygon::vec2::from_I3Position(const I3Position &p)
{
	return vec2(p.GetX(), p.GetY());
}

ExtrudedPolygon::vec2
ExtrudedPolygon::vec2::normalized(double xi, double yi)
{
	double l = std::hypot(xi, yi);
	return vec2(xi/l, yi/l);
}

bool
operator<(const ExtrudedPolygon::vec2 &a, const ExtrudedPolygon::vec2 &b)
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

ExtrudedPolygon::side::side(const vec2 &p, const vec2 &np) : origin(p),
    vector(np.x-p.x, np.y-p.y), length(std::hypot(vector.x, vector.y)),
	normal(vector.y/length, -vector.x/length, 0.)
{}

ExtrudedPolygon::ExtrudedPolygon(const std::vector<I3Position> &points, double padding)
{
	z_range_ = z_range(points);
	cap_area_ = 0;
	
	std::vector<vec2> hull = convex_hull(points);
	if (padding != 0) {
		hull = expand_polygon(hull, padding);
		z_range_.first -= padding;
		z_range_.second += padding;
	}
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

std::vector<double>
ExtrudedPolygon::GetX() const
{
	std::vector<double> x;
	BOOST_FOREACH(const side &sidey, sides_)
		x.push_back(sidey.origin.x);
	return x;
}

std::vector<double>
ExtrudedPolygon::GetY() const
{
	std::vector<double> y;
	BOOST_FOREACH(const side &sidey, sides_)
		y.push_back(sidey.origin.y);
	return y;
}

std::vector<double>
ExtrudedPolygon::GetZ() const
{
	std::vector<double> z;
	z.push_back(z_range_.first);
	z.push_back(z_range_.second);
	
	return z;
}

/// Calculate the most extreme displacements to points on the 2D hull
std::pair<double, double>
ExtrudedPolygon::GetDistanceToHull(const I3Position &pos, const I3Direction &dir) const
{
	std::pair<double, double> offsets = no_intersection();
	
	if (dir.GetX() + dir.GetY() == 0.)
		log_fatal("Direction must have a horizontal component");
	
	BOOST_FOREACH(const side &sidey, sides_) {
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
	
	return offsets;
}

/// Test whether point is inside the 2D hull by ray casting
bool
ExtrudedPolygon::PointInHull(const I3Position &pos) const
{
	int crossings = 0;
	for (std::vector<side>::const_iterator p = sides_.begin(); p != sides_.end(); p++) {
		std::vector<side>::const_iterator np = boost::next(p);
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

double
ExtrudedPolygon::GetDistanceToCap(const I3Position &p, const I3Direction &dir, double cap_z) const
{
	double d = (p.GetZ()-cap_z)/dir.GetZ();
	if (PointInHull(p + d*dir)) {
		return d;
	} else {
		return NAN;
	}
}

std::pair<double, double>
ExtrudedPolygon::GetDistanceToCaps(const I3Position &p, const I3Direction &dir) const
{
	return make_ordered_pair(GetDistanceToCap(p, dir, z_range_.first), GetDistanceToCap(p, dir, z_range_.second));
}

std::pair<double, double>
ExtrudedPolygon::GetIntersection(const I3Position &p, const I3Direction &dir) const
{
	if (std::abs(dir.GetZ()) == 1.) {
		// perfectly vertical track: only check intersections with caps
		return GetDistanceToCaps(p, dir);
	} else if (dir.GetZ() == 0.) {
		// perfectly horizontal track: only check intersections with sides
		return GetDistanceToHull(p, dir);
	} else {
		// general case: both rho and z components nonzero
		std::pair<double, double> sides = GetDistanceToHull(p, dir);
		std::pair<double, double> caps = GetDistanceToCaps(p, dir);
		double sin_zenith = std::sqrt(1.-dir.GetZ()*dir.GetZ());
		
		boost::array<double,4> combo =
		    {{sides.first/sin_zenith, sides.second/sin_zenith,
		      caps.first, caps.second}};
		std::pair<double, double> ordered(
		    *std::min_element(combo.begin(), combo.end(), nan_less),
		    *std::min_element(combo.begin(), combo.end(), nan_greater));
		
		return ordered;
	}
}

}


