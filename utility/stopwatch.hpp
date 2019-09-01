#ifndef STOPWATCH_HPP
#define STOPWATCH_HPP

#include <chrono>
#include <vector>
#include <string>
#include <iostream>

template <class precision_type>
class Stopwatch
{

typedef std::chrono::time_point<std::chrono::system_clock> StandardTime;

private:
  int n_times;
  StandardTime gs, ge;
  std::vector<StandardTime> ls, le;
  int global_elapsed;
  std::vector<int> local_elapsed;
  std::string time_unit;

public:

  Stopwatch(int n, const std::string& ps):
    n_times(n),
    ls(n),
    le(n),
    global_elapsed(0),
    local_elapsed(n,0),
    time_unit(ps)
    { }

  void global_start (void)
  {
    gs = std::chrono::system_clock::now();
  }

  void local_start (int k)
  {
    ls[k] = std::chrono::system_clock::now();
  }

  void global_stop (void)
  {
    ge = std::chrono::system_clock::now();
    global_elapsed = std::chrono::duration_cast<precision_type>(ge-gs).count();
  }

  void local_stop (int k)
  {
    le[k] = std::chrono::system_clock::now();
    local_elapsed[k] = std::chrono::duration_cast<precision_type>(le[k]-ls[k]).count();
  }

  void show_global_elapsed(void) const
  {
    std::cout << "Elapsed time: " << global_elapsed << " " << time_unit << std::endl;
  }

  void show_local_elapsed(int k) const
  {
    std::cout << "Elapsed time: " << local_elapsed[k] << " " << time_unit << std::endl;
  }

  int get_global_elapsed(void) const { return global_elapsed; }
  int get_local_elapsed(int k) const { return local_elapsed[k]; }

};

#endif /* STOPWATCH_HPP */
