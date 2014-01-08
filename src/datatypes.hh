// datatypes.hh

#ifndef DATATYPES_HEADER
#define DATATYPES_HEADER

#include "hkl_datatypes.hh"
#include "hkl_unmerge.hh"

namespace scala
{
  class SDcorrection
  // sd' = SDfac * Sqrt(sd^2 + SdB * I + (SDadd * I)^2)
  {
  public:
    SDcorrection() : sdfac(1.0), sdadd(0.0) {}
    SDcorrection(const double& SDfac, const double& SDb, const double& SDadd,
		 const bool& fixSDb=false)
      : sdfac(SDfac), sdb(SDb), sdadd(SDadd), fixsdb(fixSDb) {}

    void FixSDb() {fixsdb = true;}
    void UnFixSDb() {fixsdb = false;}

    // Reset minimum & maximum
    void ResetRange() {mincorr = +1.0e10; maxcorr = -mincorr;}

    // Multiply SDfac by update
    void UpdateFactor(const double& update) {sdfac *= update;}
    // Multiply SDadd by update
    void UpdateAdd(const double& update) {sdadd *= update;}

    void Set(const double& SDfac, const double& SDb, const double& SDadd)
    {sdfac = SDfac; sdb = SDb; sdadd = SDadd;}
    double SDfac() const {return sdfac;} // get SDfac
    double SDb() const {return sdb;} // get SDb
    double SDadd() const {return sdadd;} // get SDadd
    double& SDfac() {return sdfac;}  // set SDfac
    double& SDb() {return sdb;}  // set SDb
    double& SDadd() {return sdadd;}  // set SDadd

    // Apply appropriate SD correction 
    void Correct(observation& Observation) const;
    // in-place correction, using Iav as intensity (scaled by Gscale)
    void Correct(observation& Observation, const float& Iav) const;

    // return corrected sd as IsigI
    IsigI Corrected(const IsigI& Is) const;
    //  NB Iav should be pre-scaled by g(hl)
    IsigI Corrected(const IsigI& Is, const float& Iav) const;

    std::string format() const;
    // Range of values applied
    std::pair<double,double> MinMax() const
    {return std::pair<double,double>(mincorr, maxcorr);}

    // Get vector of parameters (2 or 3)
    std::vector<double> GetParameters() const;
    // Initial shifts for each parameter type, scaled by "scale"
    std::vector<double> GetShifts(const double& scale) const;

    // Set all parameters from vector
    void SetParameters(const std::vector<double>& params);

    int Nparams() const {return (fixsdb) ? 2 : 3;}
    
  private:
    double sdfac;  // NB number of parameters in Npar() above
    double sdb;
    double sdadd;
    bool fixsdb;   // true if SdB is fixed at 0.0
    mutable double mincorr, maxcorr;  // minimum & maximum values 
    static const double MINVARINFRAC;  // minimum fraction of input variance for correction
    static const double MINSDFRAC;     // minimum SDfac
  };
}

#endif
