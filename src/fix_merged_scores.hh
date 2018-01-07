// fix_merge_scores.hh
//
// If a merged file is being tested for undermerging, the scores for symmetry elements present
// in the the merging pointgroup group need to be set to indicate that that element is
// definitely assumed to be present
//


#ifndef FIX_MERGE_SCORES_HEADER
#define FIX_MERGE_SCORES_HEADER

// Clipper
#include <clipper/clipper.h>

#include "setscores.hh"
#include "scsignificance.hh"
#include "hkl_symmetry.hh"
#include "Output.hh"
#include "io_files.hh"
using phaser_io::LOGFILE;


namespace scala {
  RefListType FixMergedSymmScores(std::vector<SetScores>& scores,
                        std::vector<SCsignificance>& CCsig,
			const RefListType& HklinIsMerged,
			const hkl_symmetry& fileSymm,
			const hkl_symmetry& latticeSymm,
		        const ReindexOp& latticereindex,
			phaser_io::Output& output);
  //
  // Returns  MERGED or UNMERGED
  // 
  // On entry:
  //  scores       for each symmetry element, CC & Rfactor
  //  CCsig        for each symmetry element, significance etc
  //  HklinIsMerged  MERGED or UNMERGED, or NONE if unknown
  //  fileSymm     symmetry in input file, ie symmetry used for previous merging
  //  latticeSymm  lattice symmetry used for symmetry testing
  //  latticereindex reindex operator from original file to current lattice
  //  output 
  //
  // On exit:
  //  scores       for symmetry elements in fileSymm, CC=1.0, R = 0.0
  //  CCsig        for symmetry elements in fileSymm, likelihood = 1.0

  //--------------------------------------------------------------
  // Return true if the islat'th symmetry element of latticeSymm is present in reindexed mergeGroup
  bool ElementIsInGroup(const hkl_symmetry& latticeSymm, const int& islat,
			const hkl_symmetry& mergeGroup, const ReindexOp& reindex);

} // namespace scala
#endif
