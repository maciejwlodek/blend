// choosesolution.cpp

#include "choosesolution.hh"
#include "pointgroup.hh"
#include "latsym.hh"
#include "lattice.hh"
#include "string_util.hh"

using phaser_io::LOGFILE;

namespace scala {
  //------------------------------------------------------------
  ChooseSolution::ChooseSolution(const int& solution,
				 const std::string& Lauegroup,
				 const std::string& Spacegroup,
				 phaser_io::Output& output)
  : solutionnumber(solution),
    lauegroup(Lauegroup), spacegroup(Spacegroup)
  {
    std::string outstring;
    if (lauegroup != "") {
      lauegroup = CheckInputRhombohedralSGname("In CHOOSE LAUEGROUP command: ",
					       lauegroup, outstring);
      if (outstring != "") output.logTab(0,LOGFILE, outstring);
      lauegroup = SpaceGroup(lauegroup).PattersonGroup().Symbol_hm();
    }
    if (spacegroup != "") {
      spacegroup = CheckInputRhombohedralSGname("In CHOOSE SPACEGROUP command: ",
					       spacegroup, outstring);
      if (outstring != "") output.logTab(0,LOGFILE, outstring);
      // If Space group set, also set equivalent Laue group
      SpaceGroup SG(spacegroup);
      std::string sgn = SG.Symbol_hm();
      if (sgn != "") spacegroup = sgn;
      lauegroup= SG.PattersonGroup().Symbol_hm();
      //      std::string sgn = CCtbxSym::SpaceGroupName(spacegroup);
      //      if (sgn != "") spacegroup = sgn;
      //      lauegroup = CCtbxSym::PointGroup
      //	(CCtbxSym::SpaceGroupName(spacegroup)).LGname();
    }
  }
  //------------------------------------------------------------
  bool ChooseSolution::ChooseLaueGroup() const
  // Return true if Laue group choice is specified, either as
  // solution number or as Laue group
  {
    return (solutionnumber > 0) || (lauegroup != "");
  }
  //------------------------------------------------------------
  int ChooseSolution::MarkAcceptedSubgroup(std::vector<PGscore>& SGscores)
  // Mark chosen group in Laue group solution list,
  // =  0 if not found as solution number
  // = -1 if not found as Laue group
  {
    int found = 0;
    if (lauegroup != "") {
      // choosing by group name
      found = -1;
      CCtbxSym::PointGroup LG(lauegroup);
      for (size_t k=0; k<SGscores.size();k++) {
	//^
	//	std::cout << "MASG " << lauegroup << " " 
	//		  <<  LG.LGname() << " :: " <<  SGscores[k].RefLGname() << "\n";
	if (found < 0 && LG.LGname() == SGscores[k].RefLGname()) {
	  SGscores[k].SetAccept(true);  // mark as accepted, 1st one only
	  found = k+1;
	} else {
	  SGscores[k].SetAccept(false);  // mark as not accepted
	}
      }
    } else if (solutionnumber > 0) {
      // choosing by number
      for (size_t k=0; k<SGscores.size();k++) {
	if (int(k+1) == solutionnumber) { 
	  SGscores[k].SetAccept(true);  // mark as accepted
	  found = k+1;
	} else {
	  SGscores[k].SetAccept(false);  // mark as not accepted
	}
      }
    }
    return found;
  }
  //------------------------------------------------------------
  int ChooseSolution::ChooseSpacegroup(const std::vector<scala::PossibleSpaceGroup>& AllGroups) const
  // Return chosen group number in space group solution list,
  // return solution index chosen, = -1 if none
  {
    int found = -1;
    if (spacegroup == "") return found;  // no check to do
    // choosing by group name: use H lattice to check reference group
    std::string sgname = StringUtil::Strip(SGnameHtoR(spacegroup, 'H')); 
    sgname = SpaceGroup(sgname).Symbol_hm();
    for (size_t k=0; k<AllGroups.size();k++) {
      if (StringUtil::Strip(sgname) == StringUtil::Strip(AllGroups[k].Refname())) {
	found = k;
	return found;
      }
    }
    return found;
  }
  //------------------------------------------------------------
  bool ChooseSolution::IsI2() const
  // true if either Laue group or space group is I-centred monoclinic
  {
    if (lauegroup != "") {
      char LT = lauegroup[0];
      CCtbxSym::PointGroup LG(lauegroup);
      if (LT == 'I' && LG.crystal_system() == MONOCLINIC) {
	return true;
      }
    } else if (spacegroup != "") {
      char LT = spacegroup[0];
      CCtbxSym::PointGroup SG(spacegroup);
      if (LT == 'I' && SG.crystal_system() == MONOCLINIC) {
	return true;
      }
    }
    return false;
  }
  //------------------------------------------------------------
  bool ChooseSolution::IsRhombohedralR() const
  // true if either Laue group or space group is rhombohedral R-lattice
  {
    if (lauegroup != "") {
      if (lauegroup[0] == 'R') return true;
    }
    if (spacegroup != "") {
      if (spacegroup[0] == 'R') return true;
    }
    return false;
  }
}
