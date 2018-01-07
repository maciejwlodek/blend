// completeness.hh
//
// Test completeness (usually of merged file)
//

#ifndef COMPLETENESS_HEADER
#define COMPLETENESS_HEADER

#include <vector>
#include "hkl_unmerge.hh"

namespace scala {
  std::vector<double> Completeness(const hkl_unmerge_list& hkl_list, const int& Nbins);
} // namespace scala
#endif
