// testlauegroup.hh

#ifndef TESTLAUEGROUP_HEADER
#define TESTLAUEGROUP_HEADER

#include "hkl_unmerge.hh"
#include "controls.hh"
#include "globalcontrols.hh"
#include "Output.hh"
#include "pgscore.hh"
#include "io_files.hh"

//--------------------------------------------------------------
ValidElementObserved TestLaueGroup(hkl_unmerge_list& hkl_list,
				   const GlobalControls& GC,
				   const all_controls& controls, 
				   const RefListType& HklinIsMerged, 
				   phaser_io::Output& output,
				   std::vector<PGscore>& SGscores);
//--------------------------------------------------------------
ValidElementObserved ElementObserved(const hkl_symmetry& Symmetry,
				     const std::vector<SetScores>& scores,
				     const std::vector<PGscore>& SubgroupScores);
// On entry:
//  Symmetry    lattice symmetry
//  scores      element scores
//  SubgroupScores subgroups with scores
//
// Returns info about whether sufficient elements were observed


#endif

