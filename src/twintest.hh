// twintest.hh
//
// Tests for twinning
// 
// (1)  Padilla & Yeates L test
// Acta Cryst. D59,1124-1130 (2003)
// 

#ifndef TWINTEST_HEADER
#define TWINTEST_HEADER

#include "hkl_unmerge.hh"
#include "Output.hh"
#include "scala_util.hh"
#include "range.hh"
#include "hkl_symmetry.hh"

namespace scala {
  class Twintest
  {
  public:
    Twintest() :numl(0){}

    Twintest(const hkl_unmerge_list& hkl_list,
	  const double& MinIsigRatio,
	  phaser_io::Output& output);

    void init(const hkl_unmerge_list& hkl_list,
	      const double& MinIsigRatio,
	      phaser_io::Output& output);

    //! return true if there is enough data for analysis
    bool Enoughdata(const int& minimumcount) const;

    //! returns cumulative distribution of |L|
    std::vector<double> NL() const;

    //! return <|L|>
    double MeanL() const;

    //! return <L^2>
    double MeanL2() const;

    //! "Best" (maximum) estimate of twin fraction
    double Alpha() const;

    //! estimate alpha from <|L|>
    double EstimateAlphaFromAbsL() const;

    //! estimate alpha from <L^2>
    double EstimateAlphaFromL2() const;

    //! estimate alpha from N(|L|) distribution
    MeanSD EstimateAlphaFromNL() const;

    //! returns true if might be twinned
    bool PrintTwinResult(phaser_io::Output& output) const;

    void PrintCumulativeNL(phaser_io::Output& output) const;

    //! write summary data
    void WriteResultSummary(phaser_io::Output& output) const;

  private:
    MeanSD absl;   // <|L|>
    MeanSD l2;     // <L^2>
    Range lbins;
    std::vector<int> lbcount;
    int numl; // total count
    ResoRange TwinResRange;   // maybe cut back

    bool GetNeighbour(const hkl_unmerge_list& hkl_list,
		      const Hkl& hkl2, const int& irun,
		      const hkl_symmetry& symm,
		      Rtype& Ic2) const;

    // interpolate between lowerval and upperval according to where
    // val lies between lowerbin and upperbin
    double Interpolate(const double& val,
		       const double& lowerbin, const double& upperbin,
		       const double& lowerval, const double& upperval) const;

    static const double ALPHATHRESHOLD;  // threshold for twin warning

  };
} // namespace scala
#endif
