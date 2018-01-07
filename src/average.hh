//  average.hh
//
// Routines for averaging observations

#ifndef AVERAGE_HEADER
#define AVERAGE_HEADER

#include "hkl_datatypes.hh"
#include "hkl_unmerge.hh"
#include "normalise.hh"

namespace scala
{
  // Average all observations of a reflection
  // Note if there are no valid observations, this will
  //  return 0,0
  IsigI average_Is(const reflection& this_refl);
  void average_I(const reflection& this_refl, Rtype& I, Rtype& sigI);
  IsigI average_Es(const reflection& this_refl,
		   const Normalise& NormRes);
  void average_E(const reflection& this_refl, Rtype& I, Rtype& sigI,
		 const Normalise& NormRes);
} // scala
#endif

