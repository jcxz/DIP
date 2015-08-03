/**
 * A header file with an assortment of various math utilities
 */

#ifndef MATHS_H
#define MATHS_H

#include <cmath>


namespace utils {

namespace maths {

inline float cubeDiagonalLength(float a) { return std::sqrt(2.0f * a * a); }

} // End of namespace mathss

} // End of namespace utils

#endif // MATHS_H
