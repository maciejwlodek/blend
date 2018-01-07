//  scsignificance.cpp
//

#include "scsignificance.hh"
#include "probfunctions.hh"
#include "pgscore.hh"
#include "util.hh"


namespace scala
{
//--------------------------------------------------------------
  // minimum value for unrelated sigma
  float SCsignificance::SDmin = 0.10;  
  double SCsignificance::exponent = 2.0; // exponent in P(m;!S) model
//--------------------------------------------------------------
  SCsignificance::SCsignificance(const float& SCin, const int& Nsample_in)
    : SC_(SCin), NsampleSC(Nsample_in), SDUnrelSC(-1.0), SDrelSC(-1.0),
      truescore(1.0), truncated(false), dummy(false)
      {}
  //--------------------------------------------------------------
  void SCsignificance::SetUnrSC(const float& SD_SC_unrel)
  {
    SDUnrelSC = Max(SD_SC_unrel, SDmin);
    if (SDrelSC <= 0.0) SDrelSC = SDUnrelSC;
  }
  //--------------------------------------------------------------
  void SCsignificance::SetRelSC(const float& SD_SC_rel)
  {
    // No minimum
    SDrelSC = SD_SC_rel;
    if (SDUnrelSC <= 0.0) SDUnrelSC = SDrelSC;
  }
  //--------------------------------------------------------------
  float SCsignificance::Likelihood() const
  // Likelihood
  //  P(S|score)/[P(S|score) + P(!S|score)]
  //  E(score|S) = TrueScore  (eg 1.0 for correlation coefficient)
  //  E(score|!S) is in range 0.0 to TrueScore
  //
  //  Returns = 0.0 if SDUnrelSC = 0.0
  //          = 1.0 if SDUnrelSC < 0
  {
    if (SDUnrelSC < -1.9) return 1.0;  // special flag
    if (SDUnrelSC <= 0.0) return 0.0;

    IntgrtProb IP;

    if (Close<double>(exponent, 2.0)) {
      DM_2sqrt  DM(0.0,truescore); // Model for CC- if true, most likely 0.0
      //    DM_2sqrt  DM(0.0,1.0); // Model for CC- if true, most likely 0.0
      //    DM_lin  DM(0.0,1.0);
      //    DM_cubic  DM(0.0,1.0);
      //    DM_cubicsu  DM(0.0,1.0);
      //    DM_1minusmSq  DM(0.0,1.0);
      //    DM_1minusmCu  DM(0.0,1.0);
      IP.init(DM);
    } else {
      DM_power DM(0.0,truescore);
      DM.SetPower(exponent);
      IP.init(DM);
    }

    double PS, PnotS;
    // P(S|score) Lorentzian  (usual case, see class ScoreSig)
    if (truncated) {  
      PS = TruncatedLorentzianProb(SC_, truescore, SDrelSC, scmin, scmax);
      // P(!S|score) allowing for pseudosymmetry
      PnotS = IP.LorentzProb(SC_,SDUnrelSC, scmin, scmax);
    } else {
      PS = GaussProb(SC_, truescore, SDrelSC);
      PnotS = IP.Prob(SC_,SDUnrelSC);       // P(!S|score) allowing for pseudosymmetry
    }
    //^
    //    std::cout << "Likelihood (score, sd, TrueS, PS, P!S, P) " << SC_ << " "
    //    	      << SDrelSC << " " << truescore << " " << PS << " " << PnotS << " "
    //    	      << PS/(PS+PnotS) << "\n";
    //^-

    return PS/(PS+PnotS);  // Normalise
  }
  //--------------------------------------------------------------
  float SCsignificance::Z() const
  {
    if (SDUnrelSC <= 0.0) return 0.0;
    //^
    //    std::cout << "ZZ " << SC_ << " " <<  MeanUnrelSC << " " 
    //	      << SDUnrelSC << " " << (SC_ - MeanUnrelSC)/SDUnrelSC << "\n";
    return SC_/SDUnrelSC;
  }
  //--------------------------------------------------------------
  void SCsignificance::SetScoresTrue()
  // Set dummy scores to "true", ie score = 1.0 etc
  {
    SDUnrelSC = -2.0;  // flag
    SC_ = 1.0;
    NsampleSC = 4;
    dummy = true;
  }
}
