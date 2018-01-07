//
// sysabstest.cpp

#include "sysabstest.hh"
#include "sysabszones.hh"
#include "sysabsstats.hh"

using namespace scala;


//--------------------------------------------------------------
int CheckPossibleSpacegroups(std::vector<scala::PossibleSpaceGroup>& groups,
			      const double& Pcutoff)
{			      
  // Maximum score
  double Pmax = groups[0].SysAbsProb();
  double Plimit = Pmax * Pcutoff;
  int naccept = 0;
  for (size_t i=0;i<groups.size();i++) {
    //^
    //^      std::cout << "GroupProb " << i << " "
    //^		<< groups[i].Name().c_str() << " "
    //^		<< groups[i].SysAbsProb() << "\n";
    //^-
    if (groups[i].SysAbsProb() > Plimit) {
      groups[i].Accept();
      naccept++;
    }
  }
  return naccept;
}
//--------------------------------------------------------------
std::vector<scala::PossibleSpaceGroup> SysAbsTest(const PGscore& SGscore,
						  const Chirality& chiral,
						  std::vector<Zone>& SZones)
// Get list of possible spacegroups belonging to this Laue group SGscore
// Returns also updated zone list in SZones
{
  //?  CCtbxSym::PointGroup PG(SGscore.RefPGname());

  // Get list of possible spacegroups with "probabilities"
  std::vector<scala::PossibleSpaceGroup> groups
    = SGscore.TestSpaceGroupList(SZones, chiral);

  // Store things from Laue group into spacegroup
  for (size_t i=0;i<groups.size();i++) {
    groups[i].StoreLaueGroupName(SGscore.RefLGname());
    groups[i].StorePointGroupName(SGscore.RefPGname());
    groups[i].StoreLaueGroupProb(SGscore.Likelihood(), SGscore.Confidence());
    // Original to Laue group frame
    groups[i].StoreOrigReindex(SGscore.RefSGreindex());
    //^
    //            std::cout << "SGposs # " << i << " " << groups[i].Name()
    //      		<< " reindex "
    //      		<< SGscore.RefSGreindex().as_hkl() << "\n"; //^-
    // Special kludge for SG 205, in case it's a non-cyclic permutation of
    // the standard setting
    if ((groups[i].Sgnumber() == 205) && 
	groups[i].Name(true) == "P a -3 (c,b,-a)") {
      // if it is the non-cyclic, then reset reindex & the name
      groups[i].StoreOrigReindex(ReindexOp("-l,k,h")*
				 SGscore.RefSGreindex());
      groups[i].StoreSGname(groups[i].Refname());
    }
  }
  
  return groups;
}

