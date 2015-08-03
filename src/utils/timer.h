/**
 * Header file defining timing utilities
 */
#ifndef TIMER_H
#define TIMER_H

#include "utils/macros.h"
#include "utils/debug.h"

#if defined(FLUIDSIM_OS_WIN)
# ifndef WIN32_LEAN_AND_MEAN
#  define WIN32_LEAN_AND_MEAN  // do not compile unnecessary bloat
# endif
# ifndef NOMINMAX
#  define NOMINMAX             // do not define the min and max macros
# endif
# include <Windows.h>
#elif defined(FLUIDSIM_OS_LINUX)
# include <ctime>
#else
# ifdef FLUIDSIM_CC_GCC
#  warning No supported timing API on this platform
# endif
#endif


namespace utils {

namespace time {

namespace detail {

/**
 * A Windows performance counter utilizing Windows's QueryPerformanceCounter() function
 */
#ifdef FLUIDSIM_OS_WIN
class TimerWindows
{
  public:
    /**
     * Default constructor
     */
    TimerWindows(void)
      : m_start(),
        m_end(),
        m_frekv(1.0f)  // 1.0f to avoid division by zero on QueryPerformanceFrequency failure
    {
      // According to MSDN (https://msdn.microsoft.com/en-us/library/windows/desktop/ms644905%28v=vs.85%29.aspx):
      // The frequency of the performance counter is fixed at system boot and
      // is consistent across all processors. Therefore, the frequency need only
      // be queried upon application initialization, and the result can be cached.
      LARGE_INTEGER f;
      if (::QueryPerformanceFrequency(&f) == FALSE)  // get the frequency in counts per second
      {
        WARNM("QueryPerformanceFrequency failed:"
              " high resolution timer not available on this platform"
              ", you might get incorrect results");
        return;
      }

      m_frekv = double(f.QuadPart) / 1000000.0f;
    }

    /**
     * A method to fire performance counter
     */
    inline bool start(void)
    {
      return (::QueryPerformanceCounter(&m_start) != FALSE);
    }

    /**
     * A method to stop performance counter
     */
    inline bool stop(void)
    {
      return (::QueryPerformanceCounter(&m_end) != FALSE);
    }

    /**
     * A function to return the time in microseconds elapsed between
     * the calls to start and stop
     */
    inline double elapsedTime(void) const
    {
      return (m_end.QuadPart - m_start.QuadPart) / m_frekv;
    }

  private:
    LARGE_INTEGER m_start;   /// the start clock cycle counter
    LARGE_INTEGER m_end;     /// end clock cycle counter
    double m_frekv;          /// performance counter frequency
};
#endif

/**
 * Linux (Unix) performance counter
 */
#ifdef FLUIDSIM_OS_LINUX
class TimerLinux
{
  public:
    /**
     * A method to fire performance counter
     */
    inline bool start(void)
    {
      return (::clock_gettime(CLOCK_REALTIME, &m_start) == 0);
    }

    /**
     * A method to stop performance counter
     */
    inline bool stop(void)
    {
      return (::clock_gettime(CLOCK_REALTIME, &m_end) == 0);
    }

    /**
     * A function to return the time in miliseconds elapsed between
     * the calls to start and stop
     */
    inline double elapsedTime(void) const
    {
      return (m_end.tv_sec - m_start.tv_sec) * 1e6 +
             (m_end.tv_nsec - m_start.tv_nsec) * 1e-3;
    }

  private:
    struct timespec m_start;   /// start time
    struct timespec m_end;     /// end time
};
#endif

/**
 * A dummy timer used on platforms that lack supported API-s
 */
class TimerDummy
{
  public:
    /**
     * Default constructor
     */
    TimerDummy(void)
    {
      WARNM("Using dummy timer: No supported timing API on this platform");
    }

    /**
     * A method to fire performance counter
     */
    inline bool start(void)
    {
      return false;
    }

    /**
     * A method to stop performance counter
     */
    inline bool stop(void)
    {
      return false;
    }

    /**
     * A function to return the time in miliseconds elapsed between
     * the calls to start and stop
     */
    inline double elapsedTime(void) const
    {
      return 0.0f;
    }
};

} // End of namespace detail

/** Defines the timer class appropriate for the currect platform */
#if defined(FLUIDSIM_OS_WIN)
typedef detail::TimerWindows Timer;
#elif defined(FLUIDSIM_OS_LINUX)
typedef detail::TimerLinux Timer;
#else
typedef detail::TimerDummy Timer;
#endif

/**
 * A function to output elapsed time nicely formated into hours, minutes, seconds
 * and miliseconds
 *
 * @param os an output stream
 * @param elapsed_time the elapsed time in milliseconds that shall be formated
 *
 * @return the output stream originally passed in
 */
std::ostream & formatTime(double elapsed_time, std::ostream & os);

/**
 * A function to output elapsed time nicely formated into hours, minutes, seconds
 * and miliseconds
 *
 * @param os an output stream
 * @param timer the class capable of meassuring time
 *
 * @return the output stream originally passed in
 */
inline std::ostream & operator<<(std::ostream & os, const Timer & timer)
{
  return formatTime(timer.elapsedTime(), os);
}

/**
 * @brief Convenience function to convert elapsed time to nanoseconds
 * @param t the timer class
 * @return the time in nanoseconds
 */
inline double elapsedNanoSeconds(const Timer & t) { return t.elapsedTime() * 1000.0f; }

/**
 * @brief Convenience function to convert elapsed time to microseconds
 * @param t the timer class
 * @return the time in microseconds
 */
inline double elapsedMicroSeconds(const Timer & t) { return t.elapsedTime(); }

/**
 * @brief Convenience function to convert elapsed time to milliseconds
 * @param t the timer class
 * @return the time in milliseconds
 */
inline double elapsedMilliSeconds(const Timer & t) { return t.elapsedTime() / 1000.0f; }

/**
 * @brief Convenience function to convert elapsed time to seconds
 * @param t the timer class
 * @return the time in seconds
 */
inline double elapsedSeconds(const Timer & t) { return t.elapsedTime() / 1000000.0f; }

} // End of time namespace

} // End of utils namespace

#endif
