/**
 * Debugging functionality
 */

#ifndef DEBUG_H
#define DEBUG_H

#include <QDebug>
#include <iostream>

/** A macro for assert statements */
#include <cassert>
#define FLUIDSIM_ASSERT(x) assert(x)

/** A debugging macro */
#ifdef FLUIDSIM_DEBUG
#  define DBGM(x) std::cerr << x << std::endl
#else
#  define DBGM(x)
#endif

/** Info messages */
#ifndef FLUIDSIM_NO_INFOS
#  define INFOM(x) std::cerr << x << std::endl
#else
#  define INFOM(x)
#endif

/** A macro to print diagnostic messages */
#ifndef FLUIDSIM_NO_WARNINGS
#  define WARNM(x) std::cerr << x << std::endl
#else
#  define WARNM(x)
#endif

/** Error messages */
#define ERRORM(x) std::cerr << x << std::endl

namespace utils {

namespace debug {

template <typename T>
inline void printArray1D(const T *p_data, int n, std::ostream & os)
{
  if (n > 0) os << "[" << p_data[0];
  for (int i = 1; i < n; ++i)
  {
    os << ", " << p_data[i];
  }
  os << "]";
}

} // End of namespace debug

} // End of namespace utils

#endif // DEBUG_H
