// getCCsigfac.cpp
//
// Make list of unrelated pairs, use these to get relationship between sd(CC) and N
//

#include "getCCsigfac.hh"
#include "pointless.hh"

// ************************************************************
double GetCCsigFac(hkl_unmerge_list& hkl_list,
		   const ResoRange& ResRange,
		   const Normalise& NormRes,
		   double& ECC0,  // estimated E(CC) returned
		   const bool& print,
		   phaser_io::Output& output)  
{
  // Make list of unrelated pairs
  std::vector<IsPair> pair_list = NullStats(hkl_list, ResRange, NormRes,
					    ECC0, print, output);
  //  sigma(CC) = CCsigFac/sqrt(n)
  double CCsigFac =  FitScoreSig<correl_coeff>(pair_list);
  return CCsigFac;
}
