#ifndef EV_SPECIES_HPP
#define EV_SPECIES_HPP

class Species
{

private:

// Fluid particles
real_number diam_fluid;       /*!< Fluid particle diameter */
real_number hdiam_fluid;      /*!< Fluid particle radius */

// Wall particles
real_number diam_solid;       /*!< Solid particle diameter */
real_number hdiam_solid;      /*!< Solid particle radius */

// Pair cross-sections
real_number diam_gw;      /*!< Solid/fluid cross section */
real_number hdiam_gw;     /*!< Solid/fluid half cross section */

// Masses
real_number mass_fluid;   /*!< Fluid mass */
real_number mass_solid;   /*!< Solid mass */

public:

// Constructor
Species(real_number dg, real_number dw, real_number mg, real_number mw):
diam_fluid(dg),
hdiam_fluid(dg/2.0),
diam_solid(dw),
hdiam_solid(dw/2.0),
diam_gw((dw+dg)/2.0),
hdiam_gw((dw+dg)/4.0),
mass_fluid(mg),
mass_solid(mw)
{
}

~Species() = default;

inline const real_number& get_diam_fluid(void) const { return diam_fluid; }

inline real_number get_hdiam_fluid(void) const { return hdiam_fluid; }

inline real_number get_diam_solid(void) const { return diam_solid; }

inline real_number get_hdiam_solid(void) const { return hdiam_solid; }

inline real_number get_diam_gw(void) const { return diam_gw; }

inline real_number get_hdiam_gw(void) const { return hdiam_gw; }

inline const real_number& get_mass_fluid(void) const { return mass_fluid; }

inline real_number get_mass_solid(void) const { return mass_solid; }

};

#endif /* EV_SPECIES_HPP */
