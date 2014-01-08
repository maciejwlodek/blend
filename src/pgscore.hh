// pgscore.hh
//

#ifndef PGSCORE
#define PGSCORE

#include "pointgroup.hh"
#include "setscores.hh"
#include "scsignificance.hh"

namespace scala
{
  //--------------------------------------------------------------
  class PGscore : public  CCtbxSym::PointGroup
  // A pointgroup and all its scores etc
  {
  public:
    PGscore(){};
    PGscore(const CCtbxSym::PointGroup& PGin);

    void SetScore(const SetScores& SGscore);

    void SetCC(const SCsignificance& sfor, const SCsignificance& sagainst);
    void SetMSD(const SCsignificance& sfor, const SCsignificance& sagainst);
    void SetRfactor(const Rfactor& rfor, const Rfactor& ragainst);

    // store CC significance for each symmetry element separately, for & against
    void StoreCCfor(const SCsignificance& sig);
    void StoreCCagainst(const SCsignificance& sig);
    void SetCCTrueScore(const double& CCTrue) {CCTrueScore = CCTrue;}

    // Averaged (RMS) scores from just adding Z-CC for elements
    void AddAveragedZCCfor(const double& CZfor);
    void AddAveragedZCCagainst(const double& CZagainst);
    void CalcAver() const;  // Calculate RMS scores

    float CC_Zfor() const {return CCsig_for.Z();}
    float CC_Zagainst() const {return CCsig_against.Z();}
    float CC_Zgain() const {return CC_NetZ;}

    float MSD_Zfor() const {return MSDsig_for.Z();}
    float MSD_Zagainst() const {return MSDsig_against.Z();}
    float MSD_Zgain() const {return MSD_NetZ;}

    float CC_Zafor() const {CalcAver(); return AveragedZCCfor;}
    float CC_Zaagainst() const {CalcAver(); return AveragedZCCagainst;}
    float CC_Zagain() const {CalcAver(); return AveragedZCCnet;}

    SCsignificance CC_against() const {return CCsig_against;}
    Rfactor Rfactor_against() const {return Ragainst;}

    // Just store a likelihood value
    void SetLikelihood(const double& Likelihood) {likelihood = Likelihood;}
    // Return calculated or stored value
    double Likelihood() const {return likelihood;}
    // Calculate likelihood
    void CalcLikelihood();
    //! set exponent
    static void SetExponent(const double& exp) {exponent = exp;}

    void SetConfidence(const double& ConfidenceValue) {confidence = ConfidenceValue;}
    double Confidence() const {return confidence;}

    correl_coeff OverallCC() const {return score.OverallCC();}
    Rfactor OverallRfactor() const {return score.OverallRfactor();}

    void SetAccept(bool accept=true) {accepted=accept;}
    bool Accepted() const {return accepted;}

    void SetOrig(const bool& origflag) {original = origflag;}
    bool Original() const {return original;}

    static double CCSDmin() {return SDmin;}

    friend bool operator < (const PGscore& a,const PGscore& b)
    // Sort on likelihood
    {return (a.likelihood < b.likelihood);}

      // for sorting by rank on AveragedZCCnet
      //    {return (a.AveragedZCCnet < b.AveragedZCCnet);}
      // for sorting by rank on CC_NetZ
    //    {return (a.CC_NetZ < b.CC_NetZ);}
    //!      // for sorting by rank on MSD_NetZ
    //!    {return (a.MSD_NetZ < b.MSD_NetZ);}


  private:
    static double SDmin;
    double CCTrueScore;    // expected CC if symmetry true, allowing for errors
    double likelihood;   // Total probability (see Prob())
    SetScores score;
    // correlation coefficient & significance
    SCsignificance CCsig_for;      // for
    SCsignificance CCsig_against;  // against
    std::vector<SCsignificance> CCsig_for_list;  // each symmetry element for
    std::vector<SCsignificance> CCsig_against_list;  // each symmetry element against
    // mean square deviation & significance
    SCsignificance MSDsig_for;      // for
    SCsignificance MSDsig_against;  // against
    float CC_NetZ, MSD_NetZ;                   // 
    // Averaged (RMS) scores from just summing Z-scores for elements
    mutable double AveragedZCCfor;
    mutable double AveragedZCCagainst;
    mutable double AveragedZCCnet;
    int NAveragedZCCfor;
    int NAveragedZCCagainst;
    mutable bool AverCalculated; // true after RMS values calculated
    // R-factors
    Rfactor Rfor;
    Rfactor Ragainst;
    bool accepted;
    bool original; // true is this corresponds to the original pointgroup    
    double confidence;
    static double exponent;
  };
  //--------------------------------------------------------------
  // For sorting on pointgroup "order"
  class  PGcompareOrder
  {
  public:
    bool operator()(const CCtbxSym::PointGroup pg1,
		    const CCtbxSym::PointGroup pg2)
      //!    {
      //!      return (pg1.OrderValue() > pg2.OrderValue());
      //!    }
      //!    bool operator()(const PGscore pg1, const PGscore pg2)
    {
      bool Order = true;  // pg1 > pg2
      if (pg1.OrderValue() == pg2.OrderValue())
	{
	  // Compare deviation from ideal cell
	  Order = (pg1.Delta() < pg2.Delta());
	}
      else
	{
	  Order = (pg1.OrderValue() > pg2.OrderValue());
	}
      return Order;
    }
  };
  //--------------------------------------------------------------
  class Solution
  // A solution, space group or point group, with some scores
  {
  public:
    Solution(){}
    Solution(const int& Index,
	     const std::string& Groupname,
	     const bool& SGorPG,
	     const ReindexOp& Reindex,
	     const double& TotalProb,
	     const double& Confidence
	     );

    // Access
    int Index() const {return index;}
    std::string Name() const {return name;}
    bool SpaceOrPointGroup() const {return spaceorpointgroup;}
    ReindexOp Reindex() const {return reindex;}
    double TotalProb() const {return totprob;}
    double Confidence() const {return confidence;}

  private:
    int index;
    std::string name;   // group name
    bool spaceorpointgroup;    // true if space group, else point group
    ReindexOp reindex; // reindex
    double totprob;
    double confidence;
  };
}
#endif
