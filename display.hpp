#ifndef EV_IO_HPP
#define EV_IO_HPP

#include "motherbase.hpp"
#include "parallel_environment.hpp"

#include <iostream>

class IOHandler : protected Motherbase
{

private:
  std::istream& input;
  std::ostream& output;

public:
  IOHandler(DSMC*, std::istream&, std::ostream&);
  IOHandler(DSMC*);
  ~IOHandler() = default;

  inline std::istream& get_input(void) { return input; }
  inline std::ostream& get_output(void) { return output; }

  template <class data_type>
  friend IOHandler& operator << (IOHandler&, data_type&);

  template <class data_type>
  friend DefaultPointer<IOHandler>& operator << (DefaultPointer<IOHandler>&, data_type&);

};

template<class data_type>
IOHandler& operator << (IOHandler& lhs, data_type& rhs)
{
  if (lhs.par_env->get_rank() == 0)
    lhs.output << rhs;
  return lhs;
}

template<class data_type>
DefaultPointer<IOHandler>& operator << (DefaultPointer<IOHandler>& lhs, data_type& rhs)
{
  if (lhs->par_env->get_rank() == 0)
    lhs->output << rhs;
  return lhs;
}

#endif /* EV_IO_HPP */
