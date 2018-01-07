//  intensitystatistics.hh
//

#ifndef INTENSITYSTATISTICS_HEADER
#define INTENSITYSTATISTICS_HEADER

#include "hkl_unmerge.hh"
#include "range.hh"
#include "normalise.hh"

namespace scala
{
  //=================================================================
  class IntensityStatistics
  {
  public:
    IntensityStatistics(){}
    IntensityStatistics(const hkl_unmerge_list& hkl_list, const double& MinIsigRatio,
			phaser_io::Output& output);
    
    // Access
    double MeanE2minus1() const {return TotE2minus1;}
    double SD_MeanE2minus1() const {return SD_E2m1;}
    double MeanE() const {return Emean;}
    double MeanE2() const {return E2mean;}
    double ProbAcen() const {return pacentric;}
    
    void Print(phaser_io::Output& output) const; 
    
  private:
    Normalise NormRes;
    ResoRange Rrange;
    std::vector<double> E2minus1;
    std::vector<int> num;

    double TotE2minus1;  // <|E^2-1|>
    int n;
    double SD_E2m1;
    double Emean;   // <E>
    double E2mean;  // <E^2>

    double pacentric;  // probability of non-centrosymmetric
		       // p(centrosymmetric) = 1 - pacentric

  };
}
#endif
