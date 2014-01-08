// testlauegroup.cpp

#include "testlauegroup.hh"
#include "latsym.hh"
#include "normaliseunmerge.hh"
#include "pointless.hh"
#include "printthings.hh"
#include "fix_merged_scores.hh"
#include "getCCsigfac.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;
using namespace scala;

//--------------------------------------------------------------
ValidElementObserved TestLaueGroup(hkl_unmerge_list& hkl_list,
				   const GlobalControls& GC,
				   const all_controls& controls,
				   const RefListType& HklinIsMerged, 
				   phaser_io::Output& output,
				   std::vector<PGscore>& SGscores)
// If no information, returns single Laue group == input group from file
//
// On entry:
//  HklinIsMerged      == MERGED if hklinfile is known to be pre-merged
//
// On exit:
//  Returns ValidElementObserved 
//  SGscores  updated
{
  // *************************************************************
  // Get average cell from all datasets
  Scell OriginalCell = hkl_list.cell(PxdName());
  //  OriginalCell.dump(); //^

  // Lattice type character (P, I, etc) from spacegroup in mtz file
  char LatticeType = hkl_list.symmetry().lattice_type();

  // Store file (original) symmetry
  hkl_symmetry OrigSymm = hkl_list.symmetry();
  std::string OrigName = OrigSymm.symbol_xHM();

  // For XML
  Scell LatCell = OriginalCell;
  std::string LatGroupName = OrigName;
  double LatDelta = 0.0;

  // Unless OriginalLattice option selected,
  // generate lattice symmetry from unit cell (using CCtbx routines)
  ReindexOp reindex_op;

  if (!GC.OriginalLattice()) {
    double max_delta = GC.LatticeTolerance();
    CCtbxSym::LatticeSymmetry lat(OriginalCell.UnitCell(), LatticeType, GC.AllowI2(), max_delta);
    
    // "Best" spacegroup symbol for lattice symmetry
    std::string BestSpaceGroupSymbol = lat.best_spacegroup_symbol();
    
    // Lattice symmetry from name
    hkl_symmetry LatSymm(BestSpaceGroupSymbol);
    // Compare Lattice Symmetry with current symmetry in hkl_list
    // Ignore translations
    // (ie test number of operators & rotation parts only)
    bool change_symm = true;
    
    if (LatSymm == hkl_list.symmetry()) {  // compare rotations only
      // Note that although lattice & MTZ symmetry are the same
      // the symmetry operators may be in different orders,
      // so be careful to use hkl_list symmetry
      change_symm = false;
      output.logTab(0,LOGFILE, "Lattice symmetry == HKLIN symmetry\n");
    }
        
    output.logTabPrintf(0,LOGFILE,
			"\nUnit cell (from HKLIN file) used to derive lattice symmetry");
    output.logTabPrintf(0,LOGFILE,
			" with tolerance %5.1f degrees\n", max_delta);
    for (int i=0;i<6;i++) output.logTabPrintf(0,LOGFILE,"%7.2f",OriginalCell.UnitCell()[i]);
    
    output.logTab(0,LOGFILE,
		  "\n\nTolerance (and delta) is the maximum deviation from the");
    output.logTab(0,LOGFILE,
		  " expected angle between two-fold axes in the lattice group\n");
    output.logTab(0,LOGFILE, "\nLattice point group: "
		  + BestSpaceGroupSymbol + "\n");
    LatGroupName = BestSpaceGroupSymbol;
    
    // Reindexing operator for "best" spacegroup
    reindex_op = MVutil::SetCMat33(lat.best_sg_reindex_op());
    bool reindex = false;
    if ( ! reindex_op.IsIdentity()) {
      reindex = true;
    }
        
    if (change_symm || reindex)	{
      output.logTab(0,LOGFILE, "Reindexing or changing symmetry\n");
      output.logTab(0,LOGFILE, "Reindex operator from input cell to lattice cell: "
		    + reindex_op.as_hkl() + "\n\n"
		    + reindex_op.as_matrix() + "\n");
      //^
      //      std::cout << "LatSymm\n" << LatSymm.GetSpaceGroup().formatAllSymops_as_hkl() << "\n";
      //^-
      
      // Transform hkl_list to lattice symmetry
      // reindex, sort & reorganise hkl list
      hkl_list.change_symmetry(LatSymm, reindex_op);
      // Get average cell from all datasets
      Scell OverallCell = hkl_list.cell(PxdName());
      output.logTabPrintf(0,LOGFILE,"\nLattice unit cell after reindexing:");
      output.logTabPrintf(0,LOGFILE," deviation%5.2f degrees\n", lat.Delta());
      for (int i=0;i<6;i++)
	output.logTabPrintf(0,LOGFILE,"%7.2f",OverallCell.UnitCell()[i]);
      output.logTabPrintf(0,LOGFILE,"\n\n");
      
      LatCell = OverallCell.UnitCell();
      LatDelta = lat.Delta();
    }
  }  // End generate lattice option
  
  output.logTab(0,LXML,"\n<LatticeSymmetry>\n");
  output.logTab(1,LXML,"<LatticegroupName>"+LatGroupName+"</LatticegroupName>");
  output.logTab(1,LXML,LatCell.xml());
  output.logTabPrintf(1,LXML,"<CellDelta>%8.2f</CellDelta>", LatDelta);
  output.logTab(0,LXML,"\n</LatticeSymmetry>\n");

  Normalise NormRes =  scala::NormaliseUnmerge(hkl_list, GC, controls, true, output);
  ResoRange ResRange = hkl_list.ResRange();

  // ************************************************************
  // Make list of unrelated pairs
  double ECC0;  // estimated E(CC)
  double CCsigFac = GetCCsigFac(hkl_list, ResRange, NormRes,
				ECC0, true, output);
  if (CCsigFac > 0.001) {
    output.logTabPrintf(0,LOGFILE, "    Estimated sd(CC) = %5.3f / Sqrt(N)\n\n", CCsigFac);
  } else {
    CCsigFac = 1.0;
    output.logTabPrintf(0,LOGFILE, "    Set sd(CC) = %5.3f / Sqrt(N)\n\n", CCsigFac);
  }
  // ************************************************************
  // Get scores for reflections related by lattice symmetry
  int maxMult;
  int reducedMult = GC.MaxMult();
  // get list of scores for each symmetry element
  std::vector<SetScores> scores =
    MergeSym(hkl_list, ResRange, NormRes, GC.GetChiral(),
	     maxMult, reducedMult);
  if (reducedMult < maxMult) {
    output.logTabPrintf(0,LOGFILE,
			"  NOTE: high multiplicity data. Maximum multiplicity%4d reduced to maximum %3d to save time\n",
			maxMult, reducedMult);
  }

  if (hkl_list.num_reflections_icering() > 0)  {
    output.logTabPrintf(0,LOGFILE,
			"  Number of reflections omitted from ice rings: %8d\n",
			hkl_list.num_reflections_icering());
  }
  // Store CCs in significance list
  ASSERT(int(scores.size())==hkl_list.symmetry().Nelement());
  // Make list of CC-significance for each symmetry element
  // from scores list
  std::vector<SCsignificance> CCsig = CCsum(scores);
  // Assess E(CC) based on E(CC) from unrelated scores and from score
  //  for identity operator (just the mean, at present)
  double ECC = ECC0;
  if (hkl_list.symmetry().IsElementIdent(0)) {
    ECC = EstimateECC(CCsigFac, ECC0, CCsig[0]);
    output.logTab(0,LOGFILE,"");
    output.logTabPrintf(1,LOGFILE,
			"\nEstimated E(CC) of true correlation coefficient from identity = %6.3f\n\n",
			ECC);
    
  }
  // Get scores for unrelated reflections to compare with
  // related ones: update CCsig (vector<SCsignificance>)
  ScoreSig<correl_coeff>(CCsigFac, ECC, CCsig);
  // Fix up undetermined scores for merged data
  RefListType isMerged =
    FixMergedSymmScores(scores, CCsig, HklinIsMerged,
			OrigSymm, hkl_list.symmetry(), hkl_list.TotalReindex(), output);
  PrintSymElementScores(hkl_list, scores, CCsig, isMerged, output);
  // ************************************************************
  // Make list of sub-groups of lattice group
  // Combine symmetry element scores into scores for each
  // sub-pointgroup
  SGscores =
    ScoreSubGroups(hkl_list.symmetry(), OrigSymm, HklinIsMerged,
		   OrigName, OriginalCell.UnitCell(),
		   reindex_op,
		   CCsigFac, ECC, scores, CCsig, GC.AllowI2());
  // Store "confidence" values for each group relative to next one
  SetConfidenceScores(SGscores);
  return ElementObserved(hkl_list.symmetry(), scores, SGscores);
}

//--------------------------------------------------------------
ValidElementObserved ElementObserved(const hkl_symmetry& Symmetry,
				     const std::vector<SetScores>& scores,
				     const std::vector<PGscore>& SubgroupScores)
// On entry:
//  Symmetry    lattice symmetry
//  scores      element scores
//  SubgroupScores subgroups with scores
//
// Returns info about whether sufficient elements were observed
{
  int Nelement = Symmetry.Nelement();
  ASSERT(int(scores.size()) == Nelement);
  int NumOrigElement = 0;
  int NumOrigElementPresent = 0;
  int NumElementPresent = 0;

  for (int i=0;i<Nelement;i++)
    {
      if (! Symmetry.IsElementIdent(i))
	if (scores[i].OverallCC().Number() > 0) NumElementPresent++;
    }

  // We have a list of "subgroups" constructed from
  // symmetry elements - loop round subgroups
  for (size_t k=0; k<SubgroupScores.size();k++)
    {
      // Is this the original Laue group?
      if (SubgroupScores[k].Original())
	{
	  // List of symmetry elements belonging to this subgroup
	  std::vector<int> elements = SubgroupScores[k].Elements();
	  NumOrigElement += elements.size() - 1; // exclude identity
	  
	  for (size_t i=0;i<elements.size();i++)  // loop elements
	    {
	      // Count valid elements for this subgroup,
	      // Exclude identity unless triclinic
	      if (SubgroupScores[k].crystal_system() == TRICLINIC  ||
		  ! Symmetry.IsElementIdent(elements[i]))
		{
		  if (scores[elements[i]].OverallCC().Number() > 0)
		    {
		      NumOrigElementPresent += 1;
		    }
		}
	    }
        }
    }
  // How well have we sampled possible subgroups & in particular
  // the original Laue group?
  return ValidElementObserved(Nelement, NumElementPresent,
			      NumOrigElement, NumOrigElementPresent);
}

