/**
 * A header files with miscellaneous goodies
 */

#ifndef MISC_H
#define MISC_H

#include "utils/debug.h"

#include <functional>
#include <QDir>
#include <QDateTime>



namespace utils {

namespace misc {

class ScopedGuard
{
  public:
    typedef std::function<void()> Function;

  public:
    ScopedGuard(Function f) : m_f(f) { }
    ~ScopedGuard(void) { m_f(); }

  private:
    Function m_f;
};

template <typename T>
class Reference
{
  public:
    Reference(T & ref) : m_ref(ref) { }
    T & get(void) { return m_ref; }
    operator T & (void) { return m_ref; }

  private:
    T & m_ref;
};

template <typename T>
class ConstReference
{
  public:
    ConstReference(const T & ref) : m_ref(ref) { }
    const T & get(void) const { return m_ref; }
    operator const T & (void) const { return m_ref; }

  private:
    const T & m_ref;
};

/**
 * @brief currentDateToDirectoryName function will generate
 * a string that can be used as a folder name and that is
 * derived from date and time current at the time of calling
 * the function.
 * @return the newly generated string
 */
inline QString currentDateToDirectoryName(void)
{
  QDateTime dt(QDateTime::currentDateTime());
  QDate d(dt.date());
  QTime t(dt.time());
  return QString("%2-%3-%4_%5-%6")
      .arg(d.year(),   4, 10, QChar('0'))
      .arg(d.month(),  2, 10, QChar('0'))
      .arg(d.day(),    2, 10, QChar('0'))
      .arg(t.hour(),   2, 10, QChar('0'))
      .arg(t.minute(), 2, 10, QChar('0'));
}

/**
 * @brief currentDateToDirectoryPath function appends to base_dir a string
 * derived from actual date and time at the time of calling the function.
 * @param base_dir the base path to which a new folder name will be appended.
 * @return base_dir with the new folder name appended to it
 */
inline QString currentDateToDirectoryPath(const QString & base_dir)
{
  return QString("%1/%2").arg(base_dir).arg(currentDateToDirectoryName());
}

/**
 * @brief makePathFromCurrentDate this function will create a folder
 * under the path given by base_dir whose name will be derived from
 * current dat and time. The new path will subsequently be provided
 * as a return value of the function.
 * @param base_dir the path where to create the new folder
 * @return the input path with the new directory name appended to it or
 * a null string in case the creation failed.
 */
inline QString makePathFromCurrentDate(const QString & base_dir = QString("."))
{
  QString path(currentDateToDirectoryPath(base_dir));
  if (!QDir().mkpath(path))
  {
    WARNM("Failed to create output folder " << path.toStdString());
    return QString();
  }
  return path;
}

} // End of namespace misc

} // End of namespace utils

#endif // MISC_H
