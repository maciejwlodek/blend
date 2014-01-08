// setscores.hh
//

#ifndef SCALA_SETSCORES
#define SCALA_SETSCORES

//#include "scala.hh"
//#include "util.hh" 
#include "score_datatypes.hh"
#include "spline.hh"
//#include "scala_util.hh"
#include "range.hh"
#include "controls.hh"
//#include "hash.hh"
#include "mtz_unmerge_io.hh"
#include "hkl_symmetry.hh"
#include "hkl_unmerge.hh"
#include "normalise.hh"

namespace scala
{
//--------------------------------------------------------------
class SetScores
// All scores, pairwise correlation coefficients, R-factors etc
// for observations with a set
{
private:
  correl_coeff OverallCC_;
  Rfactor OverallRfactor_;

public:
  SetScores() {}  // Constructor clears all sums
  // Statistics within unmerged data
  //    statistics relative to mean (Rfactors)
  void addIstats(const IsigI mnIsigI, const reflection& this_refl,
		 const PairSet& set);
  //    pairwise statistics
  void addPairStats(const reflection& this_refl,
		    const PairSet& set,
		    const Normalise& NormRes);
  void DumpPairs(FILE* file,
		 const reflection& this_refl,
		 const PairSet& set,
		 const Normalise& NormRes);

  // Statistics between datasets (pairwise)
  void addDsetIstats(const IsigI& Is1, const IsigI& Is2,
		     const Normalise& NormRes1,
		     const Normalise& NormRes2,
		     const float& sSqr);
  // Statistics between datasets (pairwise), already normalised
  void addDsetIstats(const IsigI& Is1, const IsigI& Is2);

  // Set CC = 1, R = 0.0
  void SetScoresTrue() {OverallCC_.SetCC(); OverallRfactor_.SetR();}

  void print() const;
  // ... and we need some more return functions
  correl_coeff OverallCC() const {return OverallCC_;}
  Rfactor OverallRfactor() const {return OverallRfactor_;}

  SetScores& operator+=(const SetScores& other);
  friend SetScores& operator+ (const SetScores& a, const SetScores& b);

};
}

#endif
