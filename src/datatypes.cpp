// datatypes.cpp

#include <math.h>
#include <cmath>


#include "datatypes.hh"

namespace scala {
  // minimum fraction of input variance for correction
  const double SDcorrection::MINVARINFRAC = 0.1;
  const double SDcorrection::MINSDFRAC = 0.2;     // minimum SDfac

  //--------------------------------------------------------------
  void SDcorrection::Correct(observation& Observation, const float& Iav) const
  // in-place correction, using Iav as intensity
  {
    double var = Observation.sigI()*Observation.sigI();
    double gi = Observation.Gscale() * Iav;
    double si = sdadd * gi;
    double si2 = si*si;
    if (sdadd < 0.0) si2 = -si2;
    double ssc = Max(MINVARINFRAC*var, (var + sdb*gi + si2));
    double sd = Max(sdfac, MINSDFRAC) * sqrt(ssc);
    double corr = sd/Observation.sigI();
    mincorr = Min(corr, mincorr);
    maxcorr = Max(corr, maxcorr);
    Rtype fsd = sd;
    if (fsd <= 0.0) {
      std::cout << "sdfac " << sdfac << " " << fsd <<"\n";;
    }  //^
    Observation.set_sigI(fsd);
  }
  //--------------------------------------------------------------
  IsigI SDcorrection::Corrected(const IsigI& Is, const float& Iav) const
  {
    double var = Is.sigI()*Is.sigI();
    double si = sdadd * Iav;
    double si2 = si*si;
    if (sdadd < 0.0) si2 = -si2;
    double ssc = Max(MINVARINFRAC*var, (var + sdb*Iav + si2));
    double sd = Max(sdfac, MINSDFRAC) * sqrt(ssc);
    double corr = sd/Is.sigI();
    mincorr = Min(corr, mincorr);
    maxcorr = Max(corr, maxcorr);
    return IsigI(Is.I(), sd);
  }
  //--------------------------------------------------------------
  void SDcorrection::Correct(observation& Observation) const
  // in-place correction
  {
    double var = Observation.sigI()*Observation.sigI();
    double si = sdadd * Observation.I();
    double si2 = si*si;
    if (sdadd < 0.0) si2 = -si2;
    double ssc = Max(MINVARINFRAC*var, (var + sdb*Observation.I() + si2));
    double sd = Max(sdfac, MINSDFRAC) * sqrt(ssc);
    double corr = sd/Observation.sigI();
    mincorr = Min(corr, mincorr);
    maxcorr = Max(corr, maxcorr);
    Rtype fsd = sd;
    Observation.set_sigI(fsd);
  }
  //--------------------------------------------------------------
  IsigI SDcorrection::Corrected(const IsigI& Is) const
  {
    double var = Is.sigI()*Is.sigI();
    double si = sdadd * Is.I();
    double si2 = si*si;
    if (sdadd < 0.0) si2 = -si2;
    double ssc = Max(MINVARINFRAC*var, (var + sdb*Is.I() + si2));
    double sd = Max(sdfac, MINSDFRAC) * sqrt(ssc);
    double corr = sd/Is.sigI();
    mincorr = Min(corr, mincorr);
    maxcorr = Max(corr, maxcorr);
    return IsigI(Is.I(), sd);
  }
  //--------------------------------------------------------------
  std::string SDcorrection::format() const
  {
    return
      "SdFac = "+clipper::String(Max(sdfac, MINSDFRAC),6,3)+
      ",  SdB = "+clipper::String(sdb,6,3)+
      ",  SdAdd = "+clipper::String(sdadd,6,3);
  }
  //--------------------------------------------------------------
  std::vector<double> SDcorrection::GetParameters() const
  {
    std::vector<double> params;
    params.push_back(sdfac);
    if (!fixsdb) params.push_back(sdb);
    params.push_back(sdadd);
    return params;
  }
  //--------------------------------------------------------------
  std::vector<double> SDcorrection::GetShifts(const double& scale) const
  // Initial shifts for each parameter type, scaled by "scale"
  {
    const double SDFAC_SHIFT = 0.5;
    const double SDB_SHIFT =  5;
    const double SDADD_SHIFT = 0.05;
    std::vector<double> shifts;
    shifts.push_back(scale * SDFAC_SHIFT);
    if (!fixsdb) shifts.push_back(scale * SDB_SHIFT);
    shifts.push_back(scale * SDADD_SHIFT);
    return shifts;
  }
  //--------------------------------------------------------------
  // Set all parameters from vector (2 or 3)
  void SDcorrection::SetParameters(const std::vector<double>& params)
  {
    sdfac = Max(params[0], MINSDFRAC);
    int k=1;
    if (!fixsdb) sdb = params[k++];	 
    sdadd = params[k];
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
}
