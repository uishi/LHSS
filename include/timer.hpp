#pragma once

#include <chrono> 

struct TIMER
{
  std::chrono::high_resolution_clock::time_point start, end;

  void StartTimer()
  {
    start = std::chrono::high_resolution_clock::now();
  }

  void StopTimer()
  {
    end = std::chrono::high_resolution_clock::now();
  }

  double GetUs()
  {
    return (long double)(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count())/1.0;
  }

  double GetMs()
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
};
