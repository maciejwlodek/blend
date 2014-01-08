// sysabszones.hh

#ifndef SYSABSZONES_HEADER
#define SYSABSZONES_HEADER

#include "zone.hh"
#include "pointgroup.hh"
using CCtbxSym::PointGroup;


namespace scala
{
  class HKLTest
  {
  public:
    // Characteristic indices for axes and zones
    // static constants
    static Hkl Qh00;
    static Hkl Q0k0;
    static Hkl Q00l;
    
    static Hkl Q0kl;
    static Hkl Qh0l;
    static Hkl Qhk0;
    
    static Hkl Qhhl;
    static Hkl Qhkk;
    static Hkl Qhkh;
  };

  std::vector<Zone> SysAbsZones(const PointGroup& PG, const Chirality& chiral);
  // Make list of zones for testing for systematic absences
  // A zone may be an axis corresponding to a screw axis
  // or (if chiral != CHIRAL) a plane corresponding to a glide plane
}

#endif
