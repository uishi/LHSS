#pragma once

#include <chrono> 

struct TIMER
{
  std::chrono::high_resolution_clock::time_point start, end;

  inline void StartTimer()
  {
    start = std::chrono::high_resolution_clock::now();
  }

  inline void StopTimer()
  {
    end = std::chrono::high_resolution_clock::now();
  }

  inline double GetUs()
  {
    return (long double)(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count())/1.0;
  }

  inline double GetMs()
  {
    return (long double)(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count())/1.0;
  }

  template <class T>
  double GetTimeUs(T const &s, T const &e)
  {
    auto diff = e - s;
    return (long double)(std::chrono::duration_cast<std::chrono::microseconds>(diff).count())/1.0;
  }

  template <class T>
  double GetTimeMs(T const &s, T const &e)
  {
    auto diff = e - s;
    return (long double)(std::chrono::duration_cast<std::chrono::milliseconds>(diff).count())/1.0;
  }

  inline void WriteElapsedSec(std::ostream& os)
  {
    os << static_cast<long double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()/1.0)*1000000;
  }
};
