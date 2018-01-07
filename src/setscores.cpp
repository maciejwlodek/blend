// setscores.cpp

#include "scala.hh"

namespace scala
{
//--------------------------------------------------------------
void SetScores::addIstats(const IsigI mnIsigI, const reflection& this_refl,
			  const PairSet& set)
{
  std::vector<int> list = set.ObsIndices();

  int n = list.size();
  double an = n;
  // total weight = (multiplicity weight) * (pair weight)
  //  multiplicity weight = Sqrt(n/n-1)
  //  pair weight < 1 for centric reflections which may be related
  //   by multiple symmetry operators
  double wtmult = set.Weight() * sqrt(an/(an-1.0));

  for (int i = 0; i < n; i++)
    {
      // Note that a check on obs.valid_I() is not needed
      // as it has already been checked
      Rtype I = this_refl.get_observation(list[i]).I();
      OverallRfactor_.add((I-mnIsigI.I()), mnIsigI.I(), wtmult);
    }
}

//--------------------------------------------------------------
void SetScores::addPairStats(const reflection& this_refl,
			     const PairSet& set,
			     const Normalise& NormRes)
// Accumulate pairwise statistics (correlation coefficients etc)
// "correct" intensities by B-factor
{
  std::vector<IndexPair> pairs = set.Pairs();
  std::vector<RPair> fac = set.Fac();


  double an = Rtype(set.Nobs());
  int Npairs = pairs.size();
  // Weight for each pair = (Nobs / 2*Npair ) * (pair weight)
  //  pair weight < 1 for centric reflections which may be related
  //   by multiple symmetry operators
  double w = set.Weight() * an/(2.0*double(Npairs));
  observation obs1, obs2;

  for (int i = 0; i < Npairs; i++) {   
    //  unweighted version
    // Note that a check on obs.valid_I() is not needed
    // as it has already been checked
    obs1 = this_refl.get_observation(pairs[i].first);
    obs2 = this_refl.get_observation(pairs[i].second);
    IsigI Is1 = obs1.I_sigI().scaleIs(fac[i].first);
    IsigI Is2 = obs2.I_sigI().scaleIs(fac[i].second);
    if (NormRes.NotTooLarge(Is1) && NormRes.NotTooLarge(Is2)) {
      OverallCC_.add(Is1.I(), Is2.I(), w);
      //  weighted version (worse!)
      //	  OverallCC_.add(Is1.I(), Is2.I(), w, Is1.sigI(), Is2.sigI());
    }
  }
}

//--------------------------------------------------------------
void SetScores::DumpPairs(FILE* file,
			  const reflection& this_refl,
			  const PairSet& set,
			  const Normalise& NormRes)
// Dump pairs to file for testing
{
  std::vector<IndexPair> pairs = set.Pairs();

  //  double an = Rtype(set.Nobs());
  int Npairs = pairs.size();
  // Weight for each pair = Nobs / 2*Npair
  //  double w = an/(2.0*double(Npairs));

  //  double res = sqrt(1.0/this_refl.invresolsq());
  //  fprintf(file, "** %s   %5.2f %7.3f\n",
  //	  this_refl.hkl().format().c_str(), w, res);
  for (int i = 0; i < Npairs; i++)
    {
      IsigI I1 = this_refl.get_observation(pairs[i].first).I_sigI();
      IsigI I2 = this_refl.get_observation(pairs[i].second).I_sigI();
      //      I1 = NormRes.apply(I1, sSqr);
      //      I2 = NormRes.apply(I2, sSqr);
      //      fprintf(file, " %9.3f  %9.3f %9.3f  %9.3f %8.3f\n",
      //              I1.I(), I1.sigI(),
      //              I2.I(), I2.sigI(), w); 
      fprintf(file, " %10s  %10s%3d%6d %8.0f   %10s%3d%6d %8.0f\n",
	      this_refl.hkl().format().c_str(),
      	      this_refl.get_observation(pairs[i].first).hkl_original().format().c_str(),
      	      this_refl.get_observation(pairs[i].first).Isym(),
      	      this_refl.get_observation(pairs[i].first).Batch(),
	      I1.I(), 
      	      this_refl.get_observation(pairs[i].second).hkl_original().format().c_str(),
      	      this_refl.get_observation(pairs[i].second).Isym(),
      	      this_refl.get_observation(pairs[i].second).Batch(),
	      I2.I());
      //	      I2.I(), I2.sigI(), w); 
    }
}

//--------------------------------------------------------------
void SetScores::addDsetIstats(const IsigI& IsigI1, const IsigI& IsigI2,
		     const Normalise& NormRes1,
		     const Normalise& NormRes2,
		     const float& sSqr)
  // Statistics between datasets (pairwise)
  // Use average normalisation, just function of s^2
{
  IsigI Is1 = NormRes1.applyAvg(IsigI1, sSqr);
  IsigI Is2 = NormRes2.applyAvg(IsigI2, sSqr);
  if (NormRes1.NotTooLarge(Is1) && NormRes2.NotTooLarge(Is2))  {
    // R factors on normalised data
    double w = 1.0;
    //^
    //    std::cout << "Add: " << IsigI1.I() << "  " << IsigI2.I() << "\n";
    //    std::cout << "Add: " << Is1.I() << "  " << Is2.I() << "\n";
    //^-
    OverallRfactor_.add((Is1.I()-Is2.I()), (Is1.I()+Is2.I()), w);

    // Correlation coefficients
    OverallCC_.add(Is1.I(), Is2.I(), w);
  } 
}
//--------------------------------------------------------------
void SetScores::addDsetIstats(const IsigI& Is1, const IsigI& Is2)
  // Statistics between datasets (pairwise)
  // Already normalised
{
  // R factors on normalised data
  double w = 1.0;
  //^
  //    std::cout << "Add: " << IsigI1.I() << "  " << IsigI2.I() << "\n";
  //    std::cout << "Add: " << Is1.I() << "  " << Is2.I() << "\n";
  //^-
  OverallRfactor_.add((Is1.I()-Is2.I()), (Is1.I()+Is2.I()), w);
  
  // Correlation coefficients
  OverallCC_.add(Is1.I(), Is2.I(), w);
}
//--------------------------------------------------------------
void SetScores::print() const
{
  printf("\nOverall Correlation Coefficient: %8.3f    Number of contributions:%10d\n",
	 OverallCC_.result().val, OverallCC_.result().count);
  printf(  "Overall R-factor               : %8.3f    Number of contributions:%10d\n\n",
	   OverallRfactor_.result().val, OverallRfactor_.result().count);
}
//--------------------------------------------------------------
SetScores& SetScores::operator +=(const SetScores& other)
{
  OverallCC_ += other.OverallCC_;
  OverallRfactor_ += other.OverallRfactor_;
  return *this;
}
//--------------------------------------------------------------
SetScores& operator +
(const SetScores& a, const SetScores& b)
{
  SetScores c = a;
  return c += b;
}

}
