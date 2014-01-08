// outputunmergeddata.cpp
//

#include "outputunmergeddata.hh"
#include "writeunmergedmtz.hh"
#include "pointgroup.hh"

namespace scala{
//--------------------------------------------------------------
void OutputUnmergedData(hkl_unmerge_list& hkl_list,
                        const std::string& hklout_filename,
                        const hkl_symmetry& SPGsymm,
                        const ReindexOp& reindex_sg,
			const bool& reindexedTestData,
                        const std::string& Title,
                        const std::string& groupflag,
			phaser_io::Output& output)
//  Write out data in best pointgroup or space group
{
  output.logTabPrintf(0,LOGFILE,
		       "\nWriting unmerged data to file %s ",hklout_filename.c_str());
  output.logTab(0,LOGFILE,
		 " in "+groupflag+" group "+SPGsymm.symbol_xHM()+"\n\n");
  if (!reindexedTestData) {
    output.logTabPrintf(1,LOGFILE,
			"Reindexing operator         %s\n\n",
			reindex_sg.as_hkl().c_str());
    output.logTabPrintf(1,LOGFILE,
			"Real space transformation   %s\n\n",
			CCtbxSym::Reindex_as_xyz(reindex_sg).c_str());
  }

  if (!reindex_sg.Positive()) {
    output.logTab(0,LOGFILE,
		  "\nWARNING! reindex operator inverts hand WARNING!\n\n");
  }

  // Reindex list if necessary to match first Laue group
  // Note that "reindex" operator here is total reindexing from
  // original setting to spacegroup
  //
  //  h(hkl_list)T = h(orig)T H1       H1 = TotalReindex
  //  h(wanted)T   = h(orig)T H2       H2 = reindex_sg
  //  h(wanted)T = h(hkl_list)T H1^-1 H2
  ReindexOp reindex;
  if (!reindexedTestData) {
    reindex = hkl_list.TotalReindex().inverse() * reindex_sg;
  }
  //^
  //  std::cout << "OutputUnmergedData: reindex operator from TestList to solution: "
  //	    << reindex.as_hkl() << "\n\n";
  //  std::cout <<"Cell before reindex: " << hkl_list.Cell().formatPrint() << "\n";
  //^-

  if (hklout_filename != "") {
    // Only test rotational part of symmetry, but check for same order
    // don't bother to change symmetry just for translations
    //   but use correct space group in writing
    if (!SPGsymm.equals_r_order(hkl_list.symmetry()) ||
	! reindex.IsIdentity()) {
      int NfractIdx = hkl_list.change_symmetry(SPGsymm, reindex, true);
      if (NfractIdx > 0) {
	output.logTabPrintf(0,LOGFILE,
			    "\n%8d reflections with fractional indices discarded\n\n",
			    NfractIdx);
      }
    }
    // Clear Icerings
    hkl_list.SetIceRings(Rings());
    // Write out file with chosen space group
    bool status =
      MtzIO::WriteUnmergedMTZ(hkl_list, SPGsymm, hklout_filename, Title);
    status = status;
  }
}
}
