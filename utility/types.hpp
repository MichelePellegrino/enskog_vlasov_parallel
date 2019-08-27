#ifndef EV_TYPES_HPP
#define EV_TYPES_HPP

#include <limits>
#include <string>
#include <memory>
#include <mpi.h>

#ifndef DefaultPointer
#define DefaultPointer std::shared_ptr
#endif

#ifndef DefaultString
#define DefaultString std::string
#endif

#ifndef real_number
#define real_number double
#endif

#ifndef RNG
#define RNG ev_random::RngType::ParkMiller
#endif

#ifndef CORR
#define CORR ev_correlations::CorrelationType::CarnahanStarling
#endif

#ifndef TM
#define TM MarchingType::Standard
#endif

#ifndef MPI_DEFAULT_TAG
#define MPI_DEFAULT_TAG 0
#endif

#ifndef DefaultWatchPrecision
#define DefaultWatchPrecision std::chrono::milliseconds
#endif

enum MarchingType
{
  Standard,
  Verlet,
  Leapfrog,
  RK4,
};

namespace ev_mpi
{

  template<class data_type>
  inline MPI_Datatype MPI_Data_convert(void)
  {
    throw "Invalid type: undefined conversion from MPI_Datatype to C++ type";
  }
  template<> inline MPI_Datatype MPI_Data_convert<char>(void)                 { return MPI_CHAR; }
  template<> inline MPI_Datatype MPI_Data_convert<short int>(void)            { return MPI_SHORT; }
  template<> inline MPI_Datatype MPI_Data_convert<int>(void)                  { return MPI_INT; }
  template<> inline MPI_Datatype MPI_Data_convert<long int>(void)             { return MPI_LONG; }
  template<> inline MPI_Datatype MPI_Data_convert<unsigned char>(void)        { return MPI_UNSIGNED_CHAR; }
  template<> inline MPI_Datatype MPI_Data_convert<unsigned short int>(void)   { return MPI_UNSIGNED_SHORT; }
  template<> inline MPI_Datatype MPI_Data_convert<unsigned int>(void)         { return MPI_UNSIGNED; }
  template<> inline MPI_Datatype MPI_Data_convert<unsigned long int>(void)    { return MPI_UNSIGNED_LONG; }
  template<> inline MPI_Datatype MPI_Data_convert<float>(void)                { return MPI_FLOAT; }
  template<> inline MPI_Datatype MPI_Data_convert<double>(void)               { return MPI_DOUBLE; }
  template<> inline MPI_Datatype MPI_Data_convert<long double>(void)          { return MPI_LONG_DOUBLE; }

}

namespace ev_const
{

  constexpr double  pi      = 3.141592653589793238462643383279502884197;
  constexpr double  pio2    = 1.57079632679489661923132169163975144209858;
  constexpr double  pi2     = 6.283185307179586476925286766559005768394;
  constexpr double  sqrt2   = 1.41421356237309504880168872420969807856967;
  constexpr double  euler   = 2.71828182845904523536028747135266249776;
  constexpr float   pi_s    = 3.141592653589793238462643383279502884197;
  constexpr float   pio2_s  = 1.57079632679489661923132169163975144209858;
  constexpr float   pi2_s   = 6.283185307179586476925286766559005768394;

  constexpr real_number pinfty = std::numeric_limits<real_number>::max();
  constexpr real_number minfty = std::numeric_limits<real_number>::lowest();

  constexpr real_number pzero = std::numeric_limits<real_number>::min();
  constexpr real_number mzero = -std::numeric_limits<real_number>::min();

} /* namespace ev_const */

#endif /* EV_TYPES_HPP */
