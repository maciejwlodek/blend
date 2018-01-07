//
// fix_merge_scores.cpp
//
// If a merged file is being tested for undermerging, the scores for symmetry elements present
// in the the merging pointgroup group need to be set to indicate that that element is
// definitely assumed to be present
//

#include "fix_merged_scores.hh"


namespace scala {
  RefListType FixMergedSymmScores(std::vector<SetScores>& scores,
                        std::vector<SCsignificance>& CCsig,
			const RefListType& hklinListType,
			const hkl_symmetry& fileSymm,
			const hkl_symmetry& latticeSymm,
		        const ReindexOp& latticereindex,
			phaser_io::Output& output)
  //
  // Returns   MERGED or UNMERGED
  // 
  // On entry:
  //  scores       for each symmetry element, CC & Rfactor
  //  CCsig        for each symmetry element, significance etc
  //  hklinListType  MERGED or UNMERGED or NONE
  //  fileSymm     symmetry in input file, ie symmetry used for previous merging
  //  latticeSymm  lattice symmetry used for symmetry testing
  //  latticereindex reindex operator from original file to current lattice
  //  output 
  //
  // On exit:
  //  scores       for symmetry elements in fileSymm, CC=1.0, R = 0.0
  //  CCsig        for symmetry elements in fileSymm, likelihood = 1.0
  {
    ASSERT (latticeSymm.Nelement() == int(scores.size()));
    ASSERT (scores.size() == CCsig.size());  // sanity checks
    RefListType IsMerged = hklinListType;
    //
    if (IsMerged == MERGED) {
      //^
      //      std::cout << "FixMergedSymmScores: file symmetry: " << fileSymmtoLat.symbol_xHM() <<"\n";
      //      fileSymmtoLat.print_elements();
      //      std::cout << "\nLattice symmetry: " << latticeSymm.symbol_xHM() <<"\n";
      //      latticeSymm.print_elements();
      //^-
      // We have now decided that the file is merged, reset scores for symmetry elements
      // previously merged, ie in fileSymm reindexed
      for (int islat=0;islat<latticeSymm.Nelement();++islat) {
	//	std::cout << "\nLattice element: "
	//		  << latticeSymm.format_element(islat) << "\n";  //^
	// Is this lattice symmetry element in the merged group
	if (ElementIsInGroup(latticeSymm, islat, fileSymm, latticereindex)) {
	  //	  std::cout << "   is in group\n"; //^
	  scores[islat].SetScoresTrue(); // set CC=1.0, R=0.0
	  CCsig[islat].SetScoresTrue();  // flag for likelihood = 1.0
	}
      }
    }
    return IsMerged;
  }
  //--------------------------------------------------------------
  bool ElementIsInGroup(const hkl_symmetry& latticeSymm, const int& islat,
			const hkl_symmetry& mergeGroup, const ReindexOp& reindex)
  // Return true if the islat'th symmetry element of latticeSymm is present in reindexed mergeGroup
  {
    // First element of lattice element
    clipper::Symop elementOps = latticeSymm.RotSymopsInElement(islat)[0];
    //^    std::cout << "element symop " << elementOps[0].format() << "\n";
    std::vector<clipper::Symop> mergeOps = mergeGroup.PrimRotSymops();
    const double tol = 0.001;
    for (size_t i=0;i<mergeOps.size();++i) {
      clipper::Symop symop = reindex.Symop(mergeOps[i]);  // reindexed symop
      if (symop.equals(elementOps, tol)) { // only need to test first element op
	//^	std::cout << "symop in group " << mergeOps[i].format() << "\n";
	return true;
      }
    }
    return false;
  }

} // namespace scala

