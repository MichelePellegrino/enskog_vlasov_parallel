/*! \file stopwatch.hpp
 *  \brief Header containing the class implementing the stopwatch
 *
 *  Stopwatch class evaluates and stores (partial) CPU times; it can be used to
 *  asses the global CPU time of a procedure or local sub-times for each sub-routine
 */

#ifndef STOPWATCH_HPP
#define STOPWATCH_HPP

// STL chrono has all that is needed for evaluating (CPU) times
#include <chrono>
#include <vector>
#include <string>
#include <iostream>

/*
  Template argument 'precision_type' stands for the selected unit of measurement
  for time stamps
*/
template <class precision_type>
class Stopwatch
{

typedef std::chrono::time_point<std::chrono::system_clock> StandardTime;

private:

  int n_times;                      /*!< Number of partial times */

  // Vectors for global and local times (s = start, e = end)
  StandardTime gs, ge;
  std::vector<StandardTime> ls, le;

  int global_elapsed;               /*!< Global time (in [precision_type])            */
  std::vector<int> local_elapsed;   /*!< Vector of local times (in [precision_type])  */

  std::string time_unit;

public:

  /*
    NB: 'time_unit' has to be initialized coherently with 'precision_type'!
    (unfortunately it has to be initialized manually, it's not automatic)
  */
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
