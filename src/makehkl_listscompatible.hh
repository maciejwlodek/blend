//
// makehkl_listscompatible.hh
//
// try to make a test hkl_unmerge_list compatible with a reference one
//

#ifndef MAKELISTSCOMPATIBLE_HEADER
#define MAKELISTSCOMPATIBLE_HEADER


#include "hkl_unmerge.hh"
#include "reindexscore.hh"
#include "hkl_symmetry.hh"
#include "hkl_datatypes.hh"
#include "Output.hh"
#include "controls.hh"
#include "globalcontrols.hh"


namespace scala {

  //--------------------------------------------------------------
  // Returns list of reindex operators to match Test on to Ref
  std::vector<ReindexOp> AlternativeReindexList
     (const Scell& ref_cell, const Scell& test_cell, 
      const hkl_symmetry& TestSym,
      const double& LatticeTolerance, const int& AllowI2);
  //--------------------------------------------------------------
  // return status = 0 if OK
  class MakeHKL_listscompatible
  {
  public:
    MakeHKL_listscompatible() : status(0), altIndex(false) {}

    MakeHKL_listscompatible(hkl_unmerge_list& ref_list,
			    hkl_unmerge_list& test_list,
			    const GlobalControls& GC,
			    const all_controls& controls, 
			    phaser_io::Output& output);
    // Unmerged reference
    // CheckIndex   true to check for alternative indexing, otherwise just check symmetry is
    // compatible
    void init(hkl_unmerge_list& ref_list,
	      hkl_unmerge_list& test_list,
	      const GlobalControls& GC,
	      const bool& CheckIndex,
	      const all_controls& controls, 
	      phaser_io::Output& output);
    // Merged reference
    void init(hkl_merged_list& ref_mrgd_list,
	      hkl_unmerge_list& test_list,
	      const GlobalControls& GC,
	      const all_controls& controls, 
	      phaser_io::Output& output);


    void MakeBatchNumbersUnique(hkl_unmerge_list& ref_list,
				hkl_unmerge_list& test_list);

    int Status() const {return status;}
    // return best score object
    ReindexScore BestReindexScore() const;
    // return all scores
    std::vector<ReindexScore> BestReindexScores() const {return ReindexScoreList;}
    // return "confidence" = sqrt(p1(p1-p2)) for best and next best
    float Confidence() const;


    bool AltIndex() const {return altIndex;}
    bool IsOffset() const;  // true if any non-zero offsets
    std::vector<int> RunOffsets() const {return runOffsets;}

  private:
    void ReportDifferentCells(const Scell& RefCell,
			      const Scell& TestCell,
			      const double& Tolerance,
			      phaser_io::Output& output);

    std::vector<ReindexOp> CheckSymmetry(const hkl_symmetry& RefSymm,
					 const hkl_symmetry& TestSymm,
					 const Scell& RefCell,
					 const Scell& TestCell,
					 const double& Tolerance,
					 const int& AllowI2,
					 phaser_io::Output& output);

    bool SameLaueGroup;
    bool ChangeIndex;
    int status;
    std::vector<ReindexScore> ReindexScoreList;
    bool altIndex; // true if possibility of alternative indexing
		   // (even if identity is used)
    std::vector<int> runOffsets;
  };
}

#endif
