//
// writesummedlist.hh

#ifndef WRITESUMMEDLIST_HEADER
#define WRITESUMMEDLIST_HEADER

#include <string>
#include "hkl_unmerge.hh"

namespace scala {
  void WriteSummedList(const hkl_unmerge_list& TestList,
		       const std::string& filename);

}
#endif
