#ifndef EV_POTENTIAL_HPP
#define EV_POTENTIAL_HPP

#include "types.hpp"
#include <functional>
#include <cmath>

class NondirectionalPairPotential
{
  protected:
    typedef std::function<real_number(real_number)> real_function;
    real_function potential;
    real_function dpotential_dr;
    real_function pot_kernel;
    virtual void set_potential(void) = 0;
    virtual void set_dpotential_dr(void) = 0;
    virtual void set_pot_kernel(void) = 0;
  public:
    inline real_function& get_potential( void ) { return potential; }
    inline const real_function& get_potential( void ) const { return potential; }
    inline real_function& get_dpotential_dr( void ) { return dpotential_dr; }
    inline const real_function& get_dpotential_dr( void ) const { return dpotential_dr; }
    inline real_function& get_pot_kernel( void ) { return pot_kernel; }
    inline const real_function& get_pot_kernel( void ) const { return pot_kernel; }
    virtual ~NondirectionalPairPotential() = default;
};

class SutherlandMie : public NondirectionalPairPotential
{
  private:
    real_number pot_well;
    real_number mol_radius;
    real_number exponent;
    virtual void set_potential(void) override;
    virtual void set_dpotential_dr(void) override;
    virtual void set_pot_kernel(void) override;
  public:
    SutherlandMie(real_number, real_number, real_number);
    virtual ~SutherlandMie() = default;
};

class SutherlandMorse : public NondirectionalPairPotential
{
  private:
    real_number pot_well;
    real_number mol_radius;
    real_number alpha;
    virtual void set_potential(void) override;
    virtual void set_dpotential_dr(void) override;
    virtual void set_pot_kernel(void) override;
  public:
    SutherlandMorse(real_number, real_number, real_number);
    virtual ~SutherlandMorse() = default;
};

#endif /* EV_POTENTIAL_HPP */
