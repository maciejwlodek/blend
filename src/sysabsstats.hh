//
// sysabsstats.hh

#ifndef SYSABSSTATS_HEADER
#define SYSABSSTATS_HEADER

//#include "latsym.hh"
//#include "sysabszones.hh"
#include "zone.hh"
//#include "intensitystatistics.hh"
#include "Output.hh"

namespace scala
{
  //--------------------------------------------------------------
  void SysAbsStats(const hkl_unmerge_list& hkl_list,
		   std::vector<Zone>& SZones);
  // Accumulate statistics for systematic absence test
  //
  // Uses averaged intensities, <I>/sig(<I>)
  //--------------------------------------------------------------
}
#endif
