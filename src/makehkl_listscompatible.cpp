//
// makehkl_listscompatible.cpp
//
// try to make a test hkl_unmerge_list compatible with a reference one
//

#include "makehkl_listscompatible.hh"
#include "checkcompatiblesymmetry.hh"
#include "pointgroup.hh"
#include "scala_util.hh"
#include "average.hh"
#include "normaliseunmerge.hh"
#include "printthings.hh"
#include "testindex.hh"

namespace scala {
  //--------------------------------------------------------------
  std::vector<ReindexOp> AlternativeReindexList
     (const Scell& ref_cell, const Scell& test_cell, 
      const hkl_symmetry& TestSym,
      const double& LatticeTolerance,
      const int& AllowI2)
  // Returns list of reindex operators to match Test on to Ref
  {
    // Allowed cell discrepancy for alternative indexing
    // Very rough: convert angular tolerance to "length" using
    // average cell edge
    //	double edge = Max(Max(test_cell[0], test_cell[1]), test_cell[2]);
    double edge = 0.333 * (test_cell[0] + test_cell[1] + test_cell[2]);
    double max_diff = clipper::Util::d2rad(LatticeTolerance) * edge;

    // Reference Laue group, identity reindex op
    CCtbxSym::PointGroup TestLG(TestSym.symbol_xHM());
    TestLG.SetCell(test_cell.UnitCell(), ReindexOp(), AllowI2);
    return  CCtbxSym::AlternativeIndexing(TestLG,false,ref_cell,max_diff, AllowI2);
  }
  //--------------------------------------------------------------
MakeHKL_listscompatible::MakeHKL_listscompatible(hkl_unmerge_list& ref_list,
						 hkl_unmerge_list& test_list,
						 const GlobalControls& GC,
						 const all_controls& controls, 
						 phaser_io::Output& output)
  :  altIndex(false) 
{
  bool CheckIndex = !GC.AssumeSameIndexing();
  init(ref_list, test_list, GC, CheckIndex, controls, output);
}
//--------------------------------------------------------------
void MakeHKL_listscompatible::init(hkl_merged_list& ref_mrgd_list,
				   hkl_unmerge_list& test_list,
				   const GlobalControls& GC,
				   const all_controls& controls,
				   phaser_io::Output& output)
// Merged reference list
{
  status = 0;
  
  // Set altIndex if alternatives are possible
  std::vector<ReindexOp>  ReindexList =
    CheckSymmetry(ref_mrgd_list.symmetry(), test_list.symmetry(),
		  ref_mrgd_list.Cell(), test_list.Cell(),
		  GC.LatticeTolerance(), GC.AllowI2(), output);

  //***  altIndex = true;
  ReindexOp BestReindex;
  //  int PrintLevel = 0;
  ChangeIndex = altIndex;
  if (altIndex) {
    ReindexScoreList =
      TestIndexUnmerged(ref_mrgd_list, test_list, ReindexList, GC, controls,
			output);
    //^    PrintIndexScores(ReindexScoreList, false, output);  //^
    BestReindex = ReindexScoreList[0];
  }
    
  // Change indexing of test_list if needed, and change to ref_mrgd_list symmetry
  // Always test the data
  int Nfract_index = test_list.change_symmetry(ref_mrgd_list.symmetry(), BestReindex, true);
  ASSERT (Nfract_index == 0); // should not be any fractional indices
  // Check if cells are compatible after reindexing
  ReportDifferentCells(ref_mrgd_list.Cell(), test_list.Cell(), GC.LatticeTolerance(), output);
}
//--------------------------------------------------------------
void MakeHKL_listscompatible::init(hkl_unmerge_list& ref_list,
				   hkl_unmerge_list& test_list,
				   const GlobalControls& GC,
				   const bool& CheckIndex,
				   const all_controls& controls, 
				   phaser_io::Output& output)
{
  // Unmerged reference list
  // CheckIndex   true to check for alternative indexing, otherwise just check symmetry is
  // compatible
  status = 0;
  
  // Set altIndex if alternatives are possible
  std::vector<ReindexOp>  ReindexList =
    CheckSymmetry(ref_list.symmetry(), test_list.symmetry(),
		  ref_list.Cell(), test_list.Cell(), GC.LatticeTolerance(), GC.AllowI2(), output);

  ReindexOp BestReindex;
  int PrintLevel = 0;
  ChangeIndex = altIndex && CheckIndex;
  if (ChangeIndex) {
    ReindexScoreList =
      TestIndexBothUnmerged(ref_list, test_list, ReindexList, GC, controls,
      			    PrintLevel, output);
    //^PrintIndexScores(ReindexScoreList, false, output);  //^
    BestReindex = ReindexScoreList[0];
    if (BestReindex .IsIdentity()) {
      ChangeIndex = false;  // clear flag if best reindex is identity
    }
  }
    
  // Change indexing of test_list if needed, and change to ref_list symmetry
  if (ChangeIndex || !SameLaueGroup) {
    int Nfract_index = test_list.change_symmetry(ref_list.symmetry(), BestReindex, true);
    ASSERT (Nfract_index == 0); // should not be any fractional indices
  }
  // Check if cells are compatible after reindexing
  ReportDifferentCells(ref_list.Cell(), test_list.Cell(), GC.LatticeTolerance(), output);
  // Rest altIndex false if not checking index
  if (!CheckIndex) {altIndex = false;}
}
//--------------------------------------------------------------
std::vector<ReindexOp> MakeHKL_listscompatible::CheckSymmetry(const hkl_symmetry& RefSymm,
							      const hkl_symmetry& TestSymm,
							      const Scell& RefCell,
							      const Scell& TestCell,
							      const double& Tolerance,
							      const int& AllowI2,
							      phaser_io::Output& output)
// Returns ReindexList
// Sets SameLaueGroup = true or false
// Sets altindex
{
  // Check for equivalent symmetry (same crystal system even if
  // different Laue group)
  // Fails in here if not compatible
  SameLaueGroup = 
    CheckCompatibleSymmetry(TestSymm, RefSymm, false);
  // * * *
  // Check for alternative indexing if necessary
  //  double max_diff;
  std::vector<ReindexOp> ReindexList =
    AlternativeReindexList(RefCell, TestCell, RefSymm, Tolerance, AllowI2);
  //^  PrintAlternativeIndexing(ReindexList, Tolerance, output); //^

  altIndex = false;	
  if ((ReindexList.size() > 1) ||
      ((ReindexList.size() == 1) && 
       (!ReindexList[0].IsIdentity())))
    {altIndex = true;}  // true if we need to test for alternative indexing
  return ReindexList;
}
//--------------------------------------------------------------
void MakeHKL_listscompatible::ReportDifferentCells(const Scell& RefCell,
						  const Scell& TestCell,
						  const double& Tolerance,
						  phaser_io::Output& output)
{    
  // Check that cells are the same after reindexing
  if (! RefCell.equalsTol(TestCell, Tolerance)) {
    // Unit cells are too different
    std::string r = "";
    if (ChangeIndex) {r = " reindexed";}
    output.logTab(0,LOGFILE,
		  "\nUnit cell of"+r+" new HKLIN is too different from reference cell:-");
    output.logTabPrintf(0,LOGFILE,"   Reference:   ");
    for (int i=0;i<6;i++) output.logTabPrintf(0,LOGFILE,"%6.1f",RefCell[i]);
    output.logTabPrintf(0,LOGFILE,"\n        Test:   ");
    for (int i=0;i<6;i++) output.logTabPrintf(0,LOGFILE,"%6.1f",TestCell[i]);
    output.logTab(0,LOGFILE,"\n\n");
    status = +2;
  }
}
//--------------------------------------------------------------
void MakeHKL_listscompatible::MakeBatchNumbersUnique(hkl_unmerge_list& ref_list,
						     hkl_unmerge_list& test_list)
{
  // Make sure that test_list has unique batch numbers
  // Get list of offsets
  runOffsets =  CompareRunRanges(ref_list.RunList(), test_list.RunList());
  // Apply them
  test_list.OffsetBatchNumbers(runOffsets);
}
//--------------------------------------------------------------
bool MakeHKL_listscompatible::IsOffset() const  // true if any non-zero offsets
{
  bool offset = false;
  for (size_t i=0;i<runOffsets.size();i++) {if (runOffsets[i] != 0) offset = true;}
  return offset;
}
//--------------------------------------------------------------
ReindexScore MakeHKL_listscompatible::BestReindexScore() const
{
  if (ReindexScoreList.size() > 0) {
    return ReindexScoreList[0];
  }
  return ReindexScore();
}
//--------------------------------------------------------------
float MakeHKL_listscompatible::Confidence() const
// return "confidence" = sqrt(p1(p1-p2)) for best and next best
{
  if (ReindexScoreList.size() <= 1) {
    return 1.0; // complete confidence in one solution
  }
  float p1 = ReindexScoreList[0].Likelihood();
  float p2 = ReindexScoreList[1].Likelihood();
  if (p2 > p1) return 0.0; // shouldn't happen
  return sqrt(p1 * (p1 - p2));
}
//--------------------------------------------------------------
}  // namespace scala
