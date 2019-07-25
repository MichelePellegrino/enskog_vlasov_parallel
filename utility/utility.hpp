/*! \file ev_utilities.hpp
 *  \brief Module containing some utility template functions and classes
 */

#ifndef EV_UTILITIES
#define EV_UTILITIES

/*! \namespace ev_utility
 *  \brief A namespace containing utility functions and custom all-purpose datatypes
 */
namespace ev_utility
{

/*! \fn power_n
 *  \brief Pow function (when n, the exponent, is not known at compile time)
 *
 *  This is surely not the most efficient implementation, but it is ok for our purposes
 */
inline int power_n(int x, unsigned int n)
{
  if (n == 0)
    return 1;
  else if (n%2 == 0)
    return power_n(x, n/2)*power_n(x, n/2);
  else
    return x*power_n(x, n/2)*power_n(x, n/2);
}

// Recursive template metaprogramming
template<unsigned N> struct power_impl {
  template<typename T>
  static T calc(const T &x) {
    if (N%2 == 0)
      return power_impl<N/2>::calc(x*x);
    else if (N%3 == 0)
      return power_impl<N/3>::calc(x*x*x);
    return power_impl<N-1>::calc(x)*x;
  }
};
template<> struct power_impl<0> {
  template<typename T>
  static T calc(const T &) {
    return 1;
  }
};

/*! \fn power
 *  \brief Pow function (when N, the exponent, is known at compile time)
 *
 *  This is surely not the most efficient implementation, but it is ok for our purposes
 */
template<unsigned N, typename T>
inline T power(const T &x) {
  return power_impl<N>::calc(x);
}

} /* namespace ev_utility */

#endif /* EV_UTILITIES */
