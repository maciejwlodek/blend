// copymergedmtz.hh
//

#ifndef COPYMERGEDMTZ_HEADER
#define COPYMERGEDMTZ_HEADER

// Clipper
#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"
using clipper::Message;
using clipper::Message_fatal;


// CCP4
#include "csymlib.h"    // CCP4 symmetry stuff
#include "ccp4_general.h"
#include "cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)

#include "hkl_datatypes.hh"
#include "Output.hh"

namespace MtzIO
{
  // Copy all MTZ file changing hkl according to reindex
  // returns true if OK
  //  reduce true to reduce to asymmetric unit, else don't
  bool CopyMergedMTZ(const std::string& filename_in,
		     const std::string& filename_out,
		     const std::string& SpaceGroup,
		     const scala::ReindexOp& reindex,
		     const bool& reduce,
		     const bool& verbose,
		     phaser_io::Output& output);
}
#endif
