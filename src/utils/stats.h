/** Utility functions concerning statistics */

#ifndef STATS_H
#define STATS_H

#include <limits>


namespace utils {

namespace stats {

template <class T>
class StatisticsBase
{
  public:
    StatisticsBase(void)
      : m_min(std::numeric_limits<T>::max())
      , m_max(std::numeric_limits<T>::min())
      , m_sum()
      , m_avg()
    { }

    T min(void) const { return m_min; }
    T max(void) const { return m_max; }
    T sum(void) const { return m_sum; }
    double avg(void) const { return m_avg; }

    void add(T val, unsigned int count)
    {
      if (val > m_max) m_max = val;
      if (val < m_min) m_min = val;
      m_sum += val;
      m_avg = (m_avg * count + val) / (count + 1);
    }

    void add(const StatisticsBase<T> & other, unsigned int count, unsigned int other_count)
    {
      if (m_min > other.min()) m_min = other.min();
      if (m_max < other.max()) m_max = other.max();
      m_sum += other.sum();
      m_avg = ((m_avg * count) + (other.avg() * other_count)) / (count + other_count);
    }

    // This function only sums up the individual statistic counters,
    // which can come in handy for e.g. generating summaries
    void sumCounters(const StatisticsBase<T> & other)
    {
      m_min += other.min();
      m_max += other.max();
      m_sum += other.sum();
      m_avg += other.avg();
    }

    void clear(void)
    {
      m_min = std::numeric_limits<T>::max();
      m_max = std::numeric_limits<T>::min();
      m_sum = 0.0f;
      m_avg = 0.0f;
    }

  private:
    T m_min;
    T m_max;
    T m_sum;
    double m_avg;
};


template <typename T>
class Statistics : public StatisticsBase<T>
{
  public:
    Statistics(void)
      : StatisticsBase<T>()
      , m_count(0)
    { }

    int count(void) const { return m_count; }

    void add(T val)
    {
      StatisticsBase<T>::add(val, m_count);
      m_count++;
    }

    void add(const Statistics<T> & s)
    {
      StatisticsBase<T>::add(s, m_count, s.count());
      m_count += s.count();
    }

    // This function only sums up the individual statistic counters,
    // which can come in handy for e.g. generating summaries
    void sumCounters(const Statistics<T> & s)
    {
      StatisticsBase<T>::sumCounters(s);
      m_count += s.count();
    }

    void clear(void)
    {
      StatisticsBase<T>::clear();
      m_count = 0;
    }

  protected:
    int m_count;
};

} // End of stats namespace

} // End of utils namespace

#endif // STATS_H
