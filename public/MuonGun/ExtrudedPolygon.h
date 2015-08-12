/** $Id$
 * @file
 * @author Jakob van Santen <jakob.van.santen@desy.de>
 *
 * $Revision$
 * $Date$
 */

#ifndef I3MUONGUN_EXTRUDEDPOLYGON_H_INCLUDED
#define I3MUONGUN_EXTRUDEDPOLYGON_H_INCLUDED

#include <MuonGun/Surface.h>
#include <MuonGun/UprightSurface.h>
#include <dataclasses/I3Direction.h>

namespace simclasses {

namespace polygon {

struct vec2 {
	vec2(double xi, double yi);
	static vec2 from_I3Position(const I3Position &pos);
	static vec2 normalized(double xi, double yi);

	double x, y;

private:
	vec2() {};
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

struct side {
	side(const vec2 &point, const vec2 &next_point);
	vec2 origin;
	vec2 vector;
	double length;
	I3Direction normal;
};

}

template <typename Base>
class ExtrudedPolygonBase : public Base {
public:
	ExtrudedPolygonBase(const std::vector<I3Position> &points, double padding=0.);

	virtual std::pair<double, double> GetIntersection(const I3Position &p, const I3Direction &dir) const;
	
	double GetArea(const I3Direction &dir) const;
	double GetMaximumArea() const;
	
	I3Position SampleImpactPosition(const I3Direction &dir, I3RandomService &rng) const;
	
	std::vector<double> GetX() const;
	std::vector<double> GetY() const;
	std::vector<double> GetZ() const;

protected:
	ExtrudedPolygonBase() {};
	double GetCapArea() const { return cap_area_; }
	double GetAverageSideArea() const;
	double GetLength() const { return z_range_.second-z_range_.first; }
	std::pair<double, double> GetZRange() const { return z_range_; }

private:
	std::vector<polygon::side> sides_;
	std::pair<double, double> z_range_;
	double cap_area_;
	
	void initWithHull(const std::vector<polygon::vec2> &hull, const std::pair<double,double> &zrange);
	
	std::pair<double, double> GetDistanceToCaps(const I3Position &p, const I3Direction &dir) const;
	std::pair<double, double> GetDistanceToHull(const I3Position &p, const I3Direction &dir) const;

	bool PointInHull(const I3Position &p) const;
	double GetDistanceToCap(const I3Position &p, const I3Direction &dir, double cap_z) const;
	
	friend class boost::serialization::access;
	template <typename Archive>
	void save(Archive &, unsigned) const;
	template <typename Archive>
	void load(Archive &, unsigned);
	BOOST_SERIALIZATION_SPLIT_MEMBER();
};

class ExtrudedPolygon : public ExtrudedPolygonBase<SamplingSurface> {
private:
	typedef ExtrudedPolygonBase<SamplingSurface> Base;
public:
	ExtrudedPolygon(const std::vector<I3Position> &points, double padding=0.) : Base(points, padding) {};
private:
	ExtrudedPolygon() {}
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(ExtrudedPolygon);

}

namespace I3MuonGun {

class ExtrudedPolygon : public detail::UprightSurface<simclasses::ExtrudedPolygonBase<SamplingSurface > > {
private:
	typedef simclasses::ExtrudedPolygonBase<SamplingSurface > ExtrudedPolygonBase;
	typedef detail::UprightSurface<ExtrudedPolygonBase > Base;
public:
	virtual ~ExtrudedPolygon();
	ExtrudedPolygon(const std::vector<I3Position> &points, double padding=0.) : Base(points, padding) {};
	
	virtual bool operator==(const SamplingSurface&) const
	{
		return false;
	}

protected:
	// UprightSurface interface
	double GetTopArea() const { return ExtrudedPolygonBase::GetCapArea(); };
	double GetSideArea() const { return ExtrudedPolygonBase::GetAverageSideArea(); };
	double GetLength() const { return ExtrudedPolygonBase::GetLength(); };
	std::pair<double, double> GetZRange() const { return ExtrudedPolygonBase::GetZRange(); };

private:
	ExtrudedPolygon() {}
	friend class boost::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(ExtrudedPolygon);

}

BOOST_CLASS_VERSION(simclasses::ExtrudedPolygon, 0);
BOOST_CLASS_VERSION(I3MuonGun::ExtrudedPolygon, 0);

#endif // I3MUONGUN_EXTRUDEDPOLYGON_H_INCLUDED
