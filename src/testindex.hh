// testindex.hh

#ifndef TESTINDEX_HEADER
#define TESTINDEX_HEADER

#include "hkl_unmerge.hh"
#include "reindexscore.hh"
#include "hkl_symmetry.hh"
#include "hkl_datatypes.hh"
#include "Output.hh"
#include "controls.hh"
#include "globalcontrols.hh"


namespace scala {
  // Note that "ref_list"s cannot be const because normalisation sets
  // resolution limits
  //--------------------------------------------------------------
  std::vector<ReindexScore> TestIndexUnmerged
  (hkl_merged_list& RefList,
   hkl_unmerge_list& hkl_list,
   const std::vector<ReindexOp>& ReindexList,
   const GlobalControls& GC,
   const all_controls& controls, 
   phaser_io::Output& output);
  //--------------------------------------------------------------
  std::vector<ReindexScore> TestIndexMerged
  (hkl_merged_list& RefList,
   const hkl_merged_list& TestList,
   const std::vector<ReindexOp>& ReindexList,
   const GlobalControls& GC,
   phaser_io::Output& output);
  //--------------------------------------------------------------
  // Test unmerged data  against unmerged reference dataset
  // for alternative indexing schemes
  std::vector<ReindexScore> TestIndexBothUnmerged
      (hkl_unmerge_list& ref_list,
       hkl_unmerge_list& test_list,
       const std::vector<ReindexOp>& ReindexList,
       const GlobalControls& GC,
       const all_controls& controls, 
       const int& PrintLevel,    
       phaser_io::Output& output);
  //--------------------------------------------------------------
  // Update ReindexScoreList with significance score (likelihood)
  void GetTestUnmergedSignificance(hkl_unmerge_list& test_list,
				   const Normalise& NormRes,
				   std::vector<ReindexScore>& ReindexScoreList,
				   phaser_io::Output& output);
    
}
#endif
