// checklatticecentering.hh
//

#ifndef CHECKLATTICECENTERING_HEADER
#define CHECKLATTICECENTERING_HEADER

#include <string>
#include "hkl_unmerge.hh"

namespace scala {
  // Returns character indicating lattice type  P, A, B, C, I, F, R (hexagonal setting)
  char CheckLatticeCentering (hkl_unmerge_list& hkl_list);
  // returns the minimal space group symbol for a group with given lattice type
  // viz. P 1; A,B,C,I 2; F 2 2 2; R 3
  std::string MinimalLatticeGroup(const char& latType);
}

#endif
