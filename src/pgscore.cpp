// pgscore.cpp

#include "pgscore.hh"
#include "scala_util.hh"
#include "probfunctions.hh"
#include "spacegroupreindex.hh"


namespace scala
{
  //--------------------------------------------------------------
  double PGscore::exponent = 2.0;
  //--------------------------------------------------------------
  PGscore::PGscore(const CCtbxSym::PointGroup& PGin)
    : CCtbxSym::PointGroup(PGin),
      CCTrueScore(1.0),
      likelihood(-1.0),
      AveragedZCCfor(0.0),
      AveragedZCCagainst(0.0),
      AveragedZCCnet(0.0),
      NAveragedZCCfor(0),
      NAveragedZCCagainst(0),
      AverCalculated(false),
      accepted(false),
      original(false)
  {}
  //--------------------------------------------------------------
  void PGscore::SetScore(const SetScores& SGscore) {score = SGscore;}
  //--------------------------------------------------------------
  void PGscore::SetCC(const SCsignificance& CCfor,
		      const SCsignificance& CCagainst)
  {
    CCsig_for = CCfor;
    CCsig_against = CCagainst;
    CC_NetZ = CCsig_for.Z() - CCsig_against.Z();
    CCTrueScore = CCsig_for.TrueScore();
  }
  //--------------------------------------------------------------
  void PGscore::SetMSD(const SCsignificance& MSDfor,
		      const SCsignificance& MSDagainst)
  {
    MSDsig_for = MSDfor;
    MSDsig_against = MSDagainst;
    MSD_NetZ = MSDsig_for.Z() - MSDsig_against.Z();
  }
  //--------------------------------------------------------------
  void PGscore::SetRfactor(const Rfactor& rfor, const Rfactor& ragainst)
  {
    Rfor = rfor;
    Ragainst = ragainst;
  }
  //--------------------------------------------------------------
  void PGscore::StoreCCfor(const SCsignificance& sig)
  {
    if (sig.Nsample() > 2)
      {CCsig_for_list.push_back(sig);}
  }
  //--------------------------------------------------------------
  void PGscore::StoreCCagainst(const SCsignificance& sig)
  {
    if (sig.Nsample() > 2)
      {CCsig_against_list.push_back(sig);}
  }
  //--------------------------------------------------------------
  // Averaged (RMS) scores from just adding Z-CC for elements
  void PGscore::AddAveragedZCCfor(const double& CZfor)
  {
    if (CZfor >= 0.0)
      AveragedZCCfor += CZfor*CZfor;
    else
      AveragedZCCfor -= CZfor*CZfor;
    NAveragedZCCfor++;
    AverCalculated = false;
  }
  //--------------------------------------------------------------
  // Averaged (RMS) scores from just adding Z-CC for elements
  void PGscore::AddAveragedZCCagainst(const double& CZagainst)
  {
    if (CZagainst >= 0.0)
      AveragedZCCagainst += CZagainst*CZagainst;
    else
      AveragedZCCagainst -= CZagainst*CZagainst;
    NAveragedZCCagainst++;
    AverCalculated = false;
  }
  //--------------------------------------------------------------
  void PGscore::CalcAver() const
    // Function to calculate RMS values from totals 
 {
   if (!AverCalculated)
     {
       if (NAveragedZCCfor > 0)
	 AveragedZCCfor = AveragedZCCfor/double(NAveragedZCCfor);
       else
	 AveragedZCCfor = 0.0;

       if (NAveragedZCCagainst > 0)
	 AveragedZCCagainst = AveragedZCCagainst/double(NAveragedZCCagainst);
       else
	 AveragedZCCagainst = 0.0;

       AveragedZCCfor = SafeSqrt(AveragedZCCfor);
       AveragedZCCagainst = SafeSqrt(AveragedZCCagainst);
       AveragedZCCnet = AveragedZCCfor - AveragedZCCagainst;
     }
   AverCalculated = true;
  }
  //--------------------------------------------------------------
  void PGscore::CalcLikelihood()
  {
    // return if no data
    if (NAveragedZCCfor == 0 &&  NAveragedZCCagainst == 0) return;

    IntgrtProb IP1;
    if (Close<double>(exponent, 2.0)) {
      DM_2sqrt DM(0.0,CCTrueScore); // Model for CC- if true, most likely 0.0
      //DM_2sqrt DM(0.0,1.0); // Model for CC- if true, most likely 0.0
      //    DM_lin  DM(0.0,CCTrueScore);
      //    DM_cubic  DM(0.0,CCTrueScore);
      //    DM_cubicsu  DM(0.0,CCTrueScore);
      //    DM_1minusmSq  DM(0.0,CCTrueScore);
      //    DM_1minusmCu  DM(0.0,CCTrueScore);
      IP1.init(DM);
    } else {
      DM_power DM(0.0,CCTrueScore);
      DM.SetPower(exponent);
      IP1.init(DM);
    }

    const bool DEBUG = false;

    double CC;
    double SdCC;
    double PL;
    double PLfor = 0.0;
    double PLagainst = 0.0;
    double CCmin = -1.0;   // range of possible CC
    double CCmax = +1.0;

    if (DEBUG) std::cout << "\n" << "TrueScore " << CCTrueScore << "\n";

    // If Laue group is not true, then CC(for) should be CCTrueScore, but
    // CC(against) might be > 0.0
    // Loop elements in group
    for (size_t i=0;i<CCsig_for_list.size();i++) {
      if (CCsig_for_list[i].Nsample() > 2) {
	CC = CCsig_for_list[i].SC();
	SdCC = CCsig_for_list[i].RelSd();
	if (SdCC < 0.0) {
	  PL = 1.0;
	} else {
	  PL = TruncatedLorentzianProb(CC, CCTrueScore, SdCC, CCmin, CCmax);
	}
	PLfor += log(PL);
	//^
	if (DEBUG)
	  {printf("    In group: %8.4f %5.2f %6.3f\n",
		  CC, SdCC, PL);}
      }
    }
    // Loop elements not in group
    for (size_t i=0;i<CCsig_against_list.size();i++) {
      if (CCsig_against_list[i].Nsample() > 2) {
	CC = CCsig_against_list[i].SC();
	SdCC = CCsig_against_list[i].UnrelSd();
	if (SdCC < 0.0) {
	  PL = 1.0;
	} else {
	  PL = IP1.LorentzProb(CC, SdCC, CCmin, CCmax);
	}
	PLagainst += log(PL);
	//^
	if (DEBUG)
	  {printf("Not in group: %8.4f %5.2f %6.3f\n",
		  CC, SdCC, PL);}
      }
    }

    likelihood = exp(PLfor+PLagainst);
    //^
    if (DEBUG) {
      printf(
	     "Name: %10s CC+ %5.2f CC- %5.2f\n",
	     CCtbxSym::PointGroup::RefPGname().c_str(), CCsig_for.SC(), CCsig_against.SC());
      printf(
	     "       PL+ %5.2f  PL- %5.2f   Lk %6.3f\n",
	     exp(PLfor),  exp(PLagainst), likelihood);
    }
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  Solution::Solution(const int& Index,
		     const std::string& Groupname,
		     const bool& SGorPG,
		     const ReindexOp& Reindex,
		     const double& TotalProb,
		     const double& Confidence
		     ) : index(Index), name(Groupname), spaceorpointgroup(SGorPG),
			 reindex(Reindex), totprob(TotalProb), confidence(Confidence)
  {}
  //--------------------------------------------------------------

}
