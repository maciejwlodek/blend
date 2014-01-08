//  scsignificance.hh
//
#ifndef SCSIGNIFICANCE
#define SCSIGNIFICANCE

#include "probfunctions.hh"


namespace scala
{
  //--------------------------------------------------------------
  class SCsignificance
  // Significance of score (correlation coefficient, mean
  // square difference etc)
  {
  public:
    SCsignificance(): SDUnrelSC(-1.0), SDrelSC(-1.0),
		      truncated(false), dummy(true) {}
    SCsignificance(const float& SCin, const int& Nsample_in);

    // Store sigma value for unrelated pairs
    void SetUnrSC(const float& SD_SC_unrel);
    // Store sigma value for related pairs
    void SetRelSC(const float& SD_SC_rel);
    // SetUnrSC & SetRelSC both set the other SD value equal
    // if it has not already been set

    // Store E(score|true)
    void SetTrueScore(const double& TrueSc) {truescore = TrueSc;}
    // Store range of score for possible truncation
    void SetScoreRange(const double& SCmin, const double& SCmax)
    {scmin = SCmin; scmax = SCmax; truncated = true;}

    // Access
    float SC() const {return SC_;}
    int Nsample() const {return NsampleSC;}
    float UnrelSd() const {return SDUnrelSC;}
    float RelSd() const {return SDrelSC;}
    bool Valid() const {return (SDUnrelSC>0.0);}
    double TrueScore() const {return truescore;}

    // Likelihood
    //  P(S|score)/[P(S|score) + P(!S|score)]
    //  E(score|S) = TrueScore  (eg 1.0 for correlation coefficient)
    //  E(score|!S) is in range 0.0 to TrueScore
    //  P(S|E(S)) = Gaussian possibly truncated at SCmin and SCmax
    float Likelihood() const;
    // Z-score = SC/SD(SC)
    float Z() const;
    //! set exponent
    static void SetExponent(const double& exp) {exponent = exp;}

    // Set dummy scores to "true", ie score = 1.0 etc, set dummy flag
    void SetScoresTrue();
    // True if dummy value has been set
    bool Dummy() const {return dummy;}

  private:
    float SC_;           // sample score
    int NsampleSC;      // number of pairs in sample score
    float SDUnrelSC;    // SD of score for unrelated pairs
    float SDrelSC;    // SD of score for related pairs
    float truescore;    // E(score) if condition is true
    float scmin, scmax;
    bool truncated;
    bool dummy;
    static double exponent;    // exponent for P(m) model, default = 2
    static float SDmin;  // minimum value for unrelated sigma
			 //   (set in scsignificance.cpp)
  };
}
#endif
