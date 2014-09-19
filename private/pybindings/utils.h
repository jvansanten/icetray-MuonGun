
#ifndef MUONGUN_PYBINDINGS_UTILS_H_INCLUDED
#define MUONGUN_PYBINDINGS_UTILS_H_INCLUDED

// Bind members as universal functions if boost::numpy is available,
// otherwise fall back to a scalar function
#ifdef USE_NUMPY
#include <boost/numpy/dstream.hpp>
#define DEF(...) .def(boost::numpy::dstream::method(__VA_ARGS__))
#else
#define DEF(...) .def(__VA_ARGS__)
#endif

#endif // MUONGUN_PYBINDINGS_UTILS_H_INCLUDED
