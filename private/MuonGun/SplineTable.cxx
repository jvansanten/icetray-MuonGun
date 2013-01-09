#include <stdexcept>
#include <MuonGun/SplineTable.h>
#include <boost/serialization/binary_object.hpp>

extern "C" {
	#include <photospline/bspline.h>
}

namespace I3MuonGun {

SplineTable::SplineTable() : bias_(0)
{
	memset(&table_, sizeof(struct splinetable), 0);
}

SplineTable::SplineTable(const std::string &path)
{
	if (readsplinefitstable(path.c_str(), &table_) != 0)
		throw std::runtime_error("Couldn't read spline table " + path);
	if (splinetable_read_key(&table_, SPLINETABLE_DOUBLE, "BIAS", &bias_))
		bias_ = 0;
}

SplineTable::~SplineTable()
{
	splinetable_free(&table_);
}

int
SplineTable::Eval(double *coordinates, double *result) const
{
	int centers[table_.ndim];
	
	if (tablesearchcenters(&table_, coordinates, centers) == 0)
		*result = ndsplineeval(&table_, coordinates, centers, 0);
	else
		return EINVAL;
	
	*result -= bias_;
	
	return 0;
}

std::pair<double, double>
SplineTable::GetExtents(int dim) const
{
	if (dim < 0 || dim >= table_.ndim)
		throw std::out_of_range("Dimension index out of range");
	return std::make_pair(table_.extents[dim][0], table_.extents[dim][1]);
}

template <typename Archive>
void
SplineTable::save(Archive &ar, unsigned version) const
{
	splinetable_buffer buf;
	buf.mem_alloc = &malloc;
	buf.mem_realloc = &realloc;
	writesplinefitstable_mem(&buf, &table_);
	ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
	ar & make_nvp("NBytes", buf.size);
	ar & make_nvp("FITSFile", boost::serialization::make_binary_object(buf.data, buf.size));
	free(buf.data);
}

template <typename Archive>
void
SplineTable::load(Archive &ar, unsigned version)
{
	splinetable_buffer buf;
	buf.mem_alloc = &malloc;
	buf.mem_realloc = &realloc;
	ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
	ar & make_nvp("NBytes", buf.size);
	buf.data = buf.mem_alloc(buf.size);
	ar & make_nvp("FITSFile", boost::serialization::make_binary_object(buf.data, buf.size));
	readsplinefitstable_mem(&buf, &table_);
	free(buf.data);
}

}

I3_SERIALIZABLE(I3MuonGun::SplineTable);

