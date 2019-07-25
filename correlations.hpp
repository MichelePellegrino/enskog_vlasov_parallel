#ifndef EV_CORRELATIONS_HPP
#define EV_CORRELATIONS_HPP

#include "types.hpp"
#include "utility.hpp"
#include <cmath>

namespace ev_correlations {

enum CorrelationType { CarnahanStarling, Vera };

template<CorrelationType dummy_type>
class CorrelationFunction
{

};


// Specifications:

// Carnahan-Starling, 1969
template<>
class CorrelationFunction<CarnahanStarling>
{
public:
  real_number operator () (const real_number& eta)
  {
    return 0.5 * (2.0-eta) / ( (1.0-eta)*(1.0-eta)*(1.0-eta) );
  }
};

// Vera, 1997
template<>
class CorrelationFunction<Vera>
{
private:
  real_number csi = 0;
public:
  real_number operator () (const real_number& eta)
  {
    csi = 6.0 * eta / ( ev_const::pi * ev_const::sqrt2 );
    return 3.0 * (296.0 + csi * ( -340.0 + csi * ( -25.0+(csi*csi) *
      (18.0+142.0*(ev_utility::power<7>(csi))) ) ) ) / ( 200.0 * ev_const::pi *
      ev_const::sqrt2 * ((1.0-csi)*(1.0-csi)*(1.0-csi)) );
  }
};

}

#endif /* EV_CORRELATIONS_HPP */
