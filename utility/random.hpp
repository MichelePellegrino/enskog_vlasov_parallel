/*! \file ev_random.hpp
 *  \brief Module containing random number generators
 */

#ifndef EV_RANDOM
#define EV_RANDOM

#include "types.hpp"

#include <random>
#include <functional>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdint>

/*! \def DEFAULT_SEED
 *  \brief Default seed computed using __TIME__
 */
#ifndef DEFAULT_SEED
#define DEFAULT_SEED ( (__TIME__[7] - '0') * 1  + (__TIME__[6] - '0') * 10  + \
          (__TIME__[4] - '0') * 60   + (__TIME__[3] - '0') * 600 + \
          (__TIME__[1] - '0') * 3600 + (__TIME__[0] - '0') * 36000 ) + \
          (__LINE__ * 100000)
#endif

/*! \namespace ev_random
 *  \brief A namespace containing classes for RNG
 *
 *  This namespace contains 1) an abstract RNG interface, 2) a wrapper for STL RNG,
 *  3) template classes for custom RNG
 */
namespace ev_random
{

// ABSTRACT RNG INTERFACE

/*! \class RngAbstract
 *  \brief An abstract interface for RNG
 *
 *  This class contains the basic samplig function to be overridden as well as the
 *  definition of the algorithm for sampling a unit sphere and the Box_muller algorithm
 */
class RngAbstract
{
protected:
  int seed;
public:
  RngAbstract(int seed_ = DEFAULT_SEED):
    seed( seed_ ) { }
  inline int get_seed(void)
  {
    return seed;
  }
  inline virtual void reseed(int new_seed)
  {
    seed = new_seed;
  }
  /*! \fn sample_uniform
   *  \brief algorithm for sampling a uniform between [0,1]
   *
   *  This is the samplig algorithm to be overridden by derived RNG classes
   */
  virtual real_number sample_uniform (void) = 0;
  /*! \fn sample_unit_sphere
   *  \brief algorithm for sampling the vector k uniformly on a unit sphere
   */
  inline virtual void sample_unit_sphere (real_number& kx, real_number& ky, real_number& kz)
  {
    real_number phi;
    kx = 2.0 * sample_uniform() - 1.0;
    ky = sqrt(1.0-kx*kx);
    phi = ev_const::pi2 * sample_uniform();
    kz = ky * sin(phi);
    ky = ky * cos(phi);
    phi = sqrt(kx*kx + ky*ky + kz*kz);
    kx = kx/phi;
    ky = ky/phi;
    kz = kz/phi;
  }
  /*! \fn sample_box_muller
   *  \brief Box-Muller algorithm for sampling velocities from a Maxwellian
   */
  inline virtual void sample_box_muller (real_number mass, real_number ux, real_number uy, real_number t,
    real_number& vx, real_number& vy, real_number& vz)
  {
    real_number sigma, r , ro, teta;
    sigma = sqrt(t / mass);
    do { r = sample_uniform(); } while (r == 0);
    ro = sqrt( -2.0 * log(r) );
    r = sample_uniform();
    teta = ev_const::pi2 * r;
    vy = ro * cos(teta);
    vz = ro * sin(teta);
    do { r = sample_uniform(); } while (r == 0);
    ro = sqrt(-2.0 * log(r));
    r = sample_uniform();
    teta = ev_const::pi2 * r;
    vx = ro * cos(teta);
    vx = ux + vx * sigma;
    vy = uy + vy * sigma;
    vz = vz * sigma;
  }
  inline virtual ~RngAbstract() {}
};

// STL RNG WRAPPER

/*! \class StdRngObject
 *  \brief An wrapper for RNG form the Standard Template Library
 *
 *  This class contains the basic samplig function to be overridden as well as the
 *  definition of the algorithm for sampling a unit sphere and the Box_muller algorithm
 */
template <class engineType = std::default_random_engine>
class StdRngObject : public RngAbstract
{
private:
  engineType randon_engine;
  std::uniform_real_distribution<real_number> uniform_dist;
  /*! \var uniform_rv
   *  \brief Functional modelling a uniform r.v. defined through std::bind (efficient call)
   */
  std::function<real_number(void)> uniform_rv = std::bind(uniform_dist, randon_engine);
public:
  StdRngObject(int seed_ = DEFAULT_SEED):
    RngAbstract(seed_), randon_engine( seed_ ), uniform_dist(0.0, 1.0) {}
  inline virtual void reseed(int new_seed = DEFAULT_SEED) override
  {
    seed = new_seed;
    randon_engine.seed(seed);
  }
  inline virtual real_number sample_uniform(void) override
  {
    return uniform_rv();
  }
  inline virtual ~StdRngObject() {}
};

// CUSTOM RNG TEMPLATES

/*! \enum RngType
 *  \brief Enumeration of template arguments for custom RNG template class specification
 */
enum RngType
{
  Knuth,
  ParkMiller,
  Marsiglia,
  Splitmix,
  Xorshift64
};

/*! \class CustomRngObject
 *  \brief Prototype template function for custom RNG
 */
template <RngType dummy_rng_type>
class CustomRngObject : public RngAbstract
{
public:
  CustomRngObject(int seed_) {
    throw "Invalid spacification for CustomRngObject template";
  }
};

// SPECIFICATIONS:

/*! \class CustomRngObject<Knuth>
 *  \brief Implementation of Knuth random engine (see Numerical Recipes)
 */
template <>
class CustomRngObject<Knuth> : public RngAbstract
{
private:
  static constexpr int MBIG = 1000000000;
  static constexpr int MSEED = 161803398;
  static constexpr int MZ = 0;
  static constexpr real_number FAC = 1.0 / ( (real_number)MBIG );
  int inext;
  int inextp;
  int ma[56];
  void init(void)
  {
    int i, ii, k, mj, mk;
    mj = MSEED - abs(seed);
    mj= mj % MBIG;
    ma[55] = mj;
    mk = 1;
    for ( i = 1; i<55; ++i ) {
      ii = (i*21) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if ( mk < MZ ) mk += MBIG;
      mj = ma[ii];
    }
    for ( k=1; k<=4; ++k ) {
      for ( i=1; i<56; ++i ) {
        ma[i] -= ma[ 1 + (i+30) % 55 ];
        if ( ma[i] < MZ ) ma[i] += MBIG;
      }
    }
    inext = 0;
    inextp = 31;
  }
public:
  CustomRngObject<Knuth>(int seed_ = DEFAULT_SEED):
    RngAbstract(seed_) {
      init();
    }
  inline virtual void reseed(int new_seed = DEFAULT_SEED) override
  {
    seed = new_seed;
    init();
  }
  inline virtual real_number sample_uniform(void) override
  {
    int mj;
    if (++inext == 56)
      inext = 1;
    if (++inextp == 56)
      inextp = 1;
    mj = ma[inext-1] - ma[inextp-1];
    if(mj < MZ)
      mj += MBIG;
    ma[inext] = mj;
    return mj * FAC;
  }
  virtual ~CustomRngObject<Knuth>() { }
};

/*! \class CustomRngObject<ParkMiller>
 *  \brief Implementation of Park-Miller random engine (see SPARTA)
 */
template<>
class CustomRngObject<ParkMiller> : public RngAbstract
{
private:
  static constexpr int IA = 16807;
  static constexpr int IM = 2147483647;
  static constexpr real_number AM = (1.0/IM);
  static constexpr int IQ = 127773;
  static constexpr int IR = 2836;
  // int save;
  // int second;
  // WARMUP...?
  void init(void)
  {
    if (seed < 0)
      seed = -seed;
    if (seed == 0)
      seed = 1;
  }
public:
  CustomRngObject<ParkMiller>(int seed_ = DEFAULT_SEED):
    RngAbstract(seed_) {
      // save = 0;
      init();
    }
  inline virtual void reseed(int new_seed = DEFAULT_SEED) override
  {
    seed = new_seed;
    init();
  }
  inline virtual real_number sample_uniform(void) override
  {
    int k = seed/IQ;
    seed = IA*(seed-k*IQ) - IR*k;
    if (seed < 0) seed += IM;
    real_number ans = AM*seed;
    return ans;
  }
  virtual ~CustomRngObject<ParkMiller>() { }
};

/*! \class CustomRngObject<Marsiglia>
 *  \brief Implementation of Marsiglia random engine (see SPARTA)
 */
template<>
class CustomRngObject<Marsiglia> : public RngAbstract
{
private:
  int i97,j97;
  real_number c,cd,cm;
  real_number *u;
  void init(void)
  {
    int ij,kl,i,j,k,l,ii,jj,m;
    real_number s,t;
    while (seed > 900000000)
      seed -= 900000000;
    u = new real_number[97+1];
    ij = (seed-1)/30082;
    kl = (seed-1) - 30082*ij;
    i = (ij/177) % 177 + 2;
    j = ij %177 + 2;
    k = (kl/169) % 178 + 1;
    l = kl % 169;
    for (ii = 1; ii <= 97; ii++) {
      s = 0.0;
      t = 0.5;
      for (jj = 1; jj <= 24; jj++) {
        m = ((i*j) % 179)*k % 179;
        i = j;
        j = k;
        k = m;
        l = (53*l+1) % 169;
        if ((l*m) % 64 >= 32) s = s + t;
          t = 0.5*t;
      }
      u[ii] = s;
    }
    c = 362436.0 / 16777216.0;
    cd = 7654321.0 / 16777216.0;
    cm = 16777213.0 / 16777216.0;
    i97 = 97;
    j97 = 33;
  }
public:
  CustomRngObject<Marsiglia>(int seed_ = DEFAULT_SEED):
    RngAbstract(seed_) {
      init();
    }
  inline virtual void reseed(int new_seed = DEFAULT_SEED) override
  {
    seed = new_seed;
    init();
  }
  inline virtual real_number sample_uniform(void) override
  {
    real_number uni = u[i97] - u[j97];
    if (uni < 0.0)
      uni += 1.0;
    u[i97] = uni;
    i97--;
    if (i97 == 0)
      i97 = 97;
    j97--;
    if (j97 == 0)
      j97 = 97;
    c -= cd;
    if (c < 0.0)
      c += cm;
    uni -= c;
    if (uni < 0.0)
      uni += 1.0;
    return uni;
  }
  virtual ~CustomRngObject<Marsiglia>() { }
};

/*! \class CustomRngObject<Splitmix>
 *  \brief Implementation of Splitmix random engine
 */
/*
template<>
class CustomRngObject<Splitmix> : public RngAbstract
{
private:
  static constexpr uint64_t C1 = 0x9e3779b97f4a7c15;
  static constexpr uint64_t C2 = 0xbf58476d1ce4e5b9;
  static constexpr uint64_t C3 = 0x94d049bb133111eb;
  static constexpr uint64_t MX = 18446744073709551615;  // WARINING: INTEGER TOO LARGE TO BE REPRESENTED
  uint64_t state;
public:
  CustomRngObject<Splitmix>(int seed_ = DEFAULT_SEED):
    RngAbstract(seed_) {
      state = (uint64_t)seed;
    }
  inline virtual void reseed(int new_seed = DEFAULT_SEED) override
  {
    seed = new_seed;
    state = (uint64_t)seed;
  }
  inline virtual real_number sample_uniform(void) override
  {
    uint64_t z = ( state += C1 );
	  z = (z ^ (z >> 30)) * C2;
	  z = (z ^ (z >> 27)) * C3;
	  z = z ^ (z >> 31);
    return ( (real_number)z / MX );
  }
  virtual ~CustomRngObject<Splitmix>() { }
};
*/

/*! \class CustomRngObject<Xorshift64>
 *  \brief Implementation of Xorshift 64 random engine
 */
/*
template<>
class CustomRngObject<Xorshift64> : public RngAbstract
{
private:
  static constexpr uint64_t MX = 18446744073709551615;  // WARNING: INTEGER TOO LARGE TO BE REPRESENTED
  uint64_t state;
public:
  CustomRngObject<Xorshift64>(int seed_ = DEFAULT_SEED):
    RngAbstract(seed_) {
      state = (uint64_t)seed;
    }
  inline virtual void reseed(int new_seed = DEFAULT_SEED) override
  {
    seed = new_seed;
    state = (uint64_t)seed;
  }
  inline virtual real_number sample_uniform(void) override
  {
    uint64_t x = state;
	  x ^= (x << 13);
	  x ^= (x >> 7);
	  x ^= (x << 17);
	  state = x;
    return ( (real_number)x / MX );
  }
  virtual ~CustomRngObject<Xorshift64>() { }
};
*/


} /* namespace ev_random */

#endif /* EV_RANDOM */
