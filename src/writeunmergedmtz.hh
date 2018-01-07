// writeunmergedmtz.hh

#ifndef WRITEUNMERGEDMTZ_HEADER
#define WRITEUNMERGEDMTZ_HEADER

#include <string>

#include "hkl_unmerge.hh"
#include "Output.hh"
#include "hkl_symmetry.hh"

namespace MtzIO
{
  //--------------------------------------------------------------
  bool WriteUnmergedMTZ(scala::hkl_unmerge_list& hkl_list,
			const scala::hkl_symmetry& NewSymm,
			const std::string& filename_out,
			const std::string& title);
}
#endif
