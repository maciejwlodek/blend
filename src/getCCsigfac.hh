// getCCsigfac.hh
//
// Make list of unrelated pairs, use these to get relationship between sd(CC) and N
//

#ifndef GETCCSIGFAC_HEADER
#define GETCCSIGFAC_HEADER

#include "hkl_unmerge.hh"
#include "range.hh"
#include "normalise.hh"
#include "Output.hh"

// ************************************************************
double GetCCsigFac(hkl_unmerge_list& hkl_list,
		   const ResoRange& ResRange,
		   const Normalise& NormRes,
		   double& ECC0,  // estimated E(CC) returned
		   const bool& print,
		   phaser_io::Output& output);

#endif
