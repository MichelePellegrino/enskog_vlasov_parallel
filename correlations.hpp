/*! \file correlations.hpp
 *  \brief Header containing definitions for possible correlation functions
 *
 *  Correlation functions are implemented as functors; no parameter is needed for
 *  initialization
 */

#ifndef EV_CORRELATIONS_HPP
#define EV_CORRELATIONS_HPP

#include "types.hpp"
#include "utility.hpp"

// cmath is actually redundant
#include <cmath>

namespace ev_correlations {

/*! \enum CorrelationType
    \brief Enumeration of possible expressions for the correlation function

    Values within this enumeraion are used as template argument for template classes
*/
enum CorrelationType { CarnahanStarling, Vera };

/*! \class CorrelationFunction
 *  \brief Prototype template function for custom correlation functions
 */
template<CorrelationType dummy_type>
class CorrelationFunction
{

};


/* SPECIFICATIONS */

// CARNAHAN-STARLING, 1969
/*! \class CorrelationFunction<CarnahanStarling>
 *  \brief Carnahan-Starling expression
 */
template<>
class CorrelationFunction<CarnahanStarling>
{
public:
  real_number operator () (const real_number& eta)
  {
    return 0.5 * (2.0-eta) / ( (1.0-eta)*(1.0-eta)*(1.0-eta) );
  }
};


// VERA, 1997
/*! \class CorrelationFunction<Vera>
 *  \brief Vera expression
 */
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
