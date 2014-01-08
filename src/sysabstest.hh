// sysabstest.hh

#ifndef SYSABSTEST_HEADER
#define SYSABSTEST_HEADER

#include "Output.hh"
#include "hkl_unmerge.hh"
#include "hkl_symmetry.hh"
#include "controls.hh"
#include "zone.hh"
#include "pgscore.hh"

#include <vector>

using namespace scala;

int CheckPossibleSpacegroups(std::vector<scala::PossibleSpaceGroup>& groups,
			     const double& Pcutoff);

// Get list of possible spacegroups belonging to this
// Laue group SGscore
// Returns also zone list in SZones
std::vector<scala::PossibleSpaceGroup> 
SysAbsTest(const PGscore& SGscore,
	   const Chirality& chiral,
	   std::vector<Zone>& SZones);

#endif
