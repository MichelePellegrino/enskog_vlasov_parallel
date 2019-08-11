/*! \file integration.hpp
 *  \brief Header containing classes for numerical integration
 */

#ifndef EV_INTEGRATION_HPP
#define EV_INTEGRATION_HPP

#define DEFAULT_EPS 1e-4

#include <functional>

#include "types.hpp"
#include "romberg.hpp"

/*! \namespace ev_numeric
 *  \brief Namespace for numerical procedures
 *
 *  In this namespace are defined classes for numerical integration.
 */
namespace ev_numeric
{


enum IntegralType {Finite, Infinite};

/*! \class NumericalIntegrator
 *  \brief A class to perform numerical integration in finite and infinite 1D domains
 */
template<IntegralType dummy_integral_type>
class NumericalIntegrator
{ };

/* Template specialization for finite-domain integrals */
template<>
class NumericalIntegrator<Finite>
{
private:
  std::function<real_number(real_number)> func;
  real_number a, b;
  real_number eps;
public:
  NumericalIntegrator() = default;
  NumericalIntegrator
  (std::function<real_number(real_number)> FUN, real_number A, real_number B, real_number EPS = DEFAULT_EPS):
    func(FUN), a(A), b(B), eps(EPS) { }
  real_number integrate ( void )
  {
    return qromb(func, a, b, eps);
  }
  real_number integrate
  (std::function<real_number(real_number)> FUN, real_number A, real_number B, real_number EPS = DEFAULT_EPS)
  {
    a = A; b = B;
    return qromb(func, a, b, eps);
  }
  ~NumericalIntegrator() = default;
};

/* Template specialization for infinite-domain integrals */
template<>
class NumericalIntegrator<Infinite>
{
private:
  std::function<real_number(real_number)> func;                                 // integrand for finite integral
  std::function<real_number(real_number)> funk                                  // integrand for infinite integral
    = [this](real_number x) { return func(1.0/x)/(x*x); };
  real_number a, b;
  real_number eps;
public:
  NumericalIntegrator() = default;
  NumericalIntegrator
  (std::function<real_number(real_number)> FUN, real_number A, real_number B, real_number EPS = DEFAULT_EPS):
    func(FUN), a(1.0/B), b(1.0/A), eps(EPS) { }
  real_number integrate ( void )
  {
    return qromb(funk, a, b, eps);
  }
  real_number integrate
  (std::function<real_number(real_number)> FUN, real_number A, real_number B, real_number EPS = DEFAULT_EPS)
  {
    a = 1.0/B;
    b = 1.0/A;
    return qromb(funk, a, b, eps);
  }
  ~NumericalIntegrator() = default;
};

} /* namespace ev_numeric */

#endif /* EV_INTEGRATION_HPP */
