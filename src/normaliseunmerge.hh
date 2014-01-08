// normaliseunmerge.hh


#ifndef NORMALISEUNMERGE_HEADER
#define NORMALISEUNMERGE_HEADER

#include "hkl_unmerge.hh"
#include "globalcontrols.hh"
#include "controls.hh"
#include "Output.hh"
#include "normalise.hh"

namespace scala {
//--------------------------------------------------------------
Normalise NormaliseUnmerge(hkl_unmerge_list& hkl_list,
                           const GlobalControls& GC,
			   const all_controls& controls, 
			   const bool& verbose,
                           phaser_io::Output& output);
}

#endif
