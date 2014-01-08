// pointless.hh

// Routines in pointless.cpp needed elsewhere (mainly testlauegroup)

#ifndef POINTLESS_HEADER
#define POINTLESS_HEADER

#include <vector>
#include "io_files.hh"
#include "normalise.hh"
#include "hkl_unmerge.hh"
#include "score_datatypes.hh"
#include "range.hh"
#include "Output.hh"
#include "pgscore.hh"

using namespace scala;
//--------------------------------------------------------------
enum PointlessTask {NOTASK, COPY, LAUE, REFERENCE, CENTRE};

class PointlessStatus {
public:
  // Records what Pointless has done & success or failure
  PointlessStatus () : success(false), NumHKLIN(0), 
		       IsHKLOUT(false) {
    task = NOTASK;
    HKLrefFileType = NONE;
  }
  PointlessStatus (const bool& ssuccess, const int& numhklin,
		   const scala::RefListType& reftype, const bool& ishklout,
		   const PointlessTask& Task)
    : success(ssuccess), NumHKLIN(numhklin), HKLrefFileType(reftype),
      IsHKLOUT(ishklout), task(Task)  {}

  bool Success() const { return success;}

private:
  bool success;
  int NumHKLIN;       // Number of HKLIN files
  scala::RefListType HKLrefFileType;    //  NONE, MERGED, or UNMERGED
  bool IsHKLOUT;                         //  true if output file written
  PointlessTask task;
};

//--------------------------------------------------------------
std::vector<IsPair> NullStats(hkl_unmerge_list& ref_list,
                              const ResoRange& ResRange,
                              const Normalise& NormRes,
			      double& ECC0,
			      const bool& print,
                              phaser_io::Output& output);

//--------------------------------------------------------------
template<class T> std::vector<T> Score_Pairs(const std::vector<IsPair>& pair_list,
                                             const int& Ngroup)
// Get scores, eg correlation
// in groups of Ngroup pairs
// On entry:
//  pair_list    list of unrelated pairs, normalised
//  Ngroup       number per group
//
// On exit, returns:
//  sc_group     vector of scores for groups of size Ngroup
{
  std::vector<T> sc_group;
  T sc;
  //  const float w = 1.0;

  int k = -1;
  for (size_t i=0;i<pair_list.size();i++) {
    if (++k%Ngroup == 0)  {
      if (i > 0) {
	sc_group.push_back(sc);
      }
      sc.zero();
      k = 0;
    }
    sc.add(pair_list[i].Is1(), pair_list[i].Is2());
  }
  return sc_group;
}
//--------------------------------------------------------------
template<class T> RPair MeanScore(const std::vector<T>& sc_list)
// Return mean & SD of list of scores eg correlation coefficients
// Class T must provide .result().val & .result().count members
{
  int n = sc_list.size();
  double sum_sc = 0.0;
  double sum_sc2 = 0.0;
  double sc;
  double sd;
  int count = 0;

  if (n<1) return RPair(0.0,0.0);
  if (n==1) {
    sc = sc_list[0].result().val;
    count = sc_list[0].result().count;
    sd = 0.0;
  } else {
    for (int i = 0;i<n;i++) {
      sc = sc_list[i].result().val;
      sum_sc  += sc;
      sum_sc2 += sc*sc;
      count += sc_list[i].result().count;
    }
    // Mean
    sc = sum_sc/double(n);
    // SD
    sd = sqrt((sum_sc2 - double(n)*sc*sc)/double(n-1));
  }

  return RPair(sc, sd);
}
//--------------------------------------------------------------
template <class T> double FitScoreSig(const std::vector<IsPair>& pair_list)
// Get scores between unrelated pairs within same resolution bin
// as a control to compare with symmetry elements or point groups
//
// Divide pair list into groups of ascending size Ngroup and fit
// to function   sigma(score) = SigFac/Sqrt(Ngroup)
// Returns SigFac
//
// On entry:
//  pair_list    list of unrelated pairs, normalised
//
// 
{
  std::vector<T> SC_group;
  int Npairs = pair_list.size();

  // Correlate them in pairs, within groups of ascending size
  int MinNumGroup = 10;  // Minimum number of groups
  int MaxNgroup = Min(Npairs/MinNumGroup, 200); // Maximum number in group
  int MinNgroup = Min(5, MaxNgroup);

  if(MaxNgroup-MinNgroup < 4) {return 0.0;}
  
  LinearFit SigFunction;
  const float w = 1.0;
  
  for (int Ngroup=MinNgroup;Ngroup<=MaxNgroup;Ngroup++) {
    SC_group = Score_Pairs<T>(pair_list, Ngroup);
    // calculate mean & sd of SCs
    RPair SCmnsd = MeanScore<T>(SC_group);
    // x = Ngroup^-0.5, y = sd(score)
    SigFunction.add(1.0f/sqrt(float(Ngroup)), SCmnsd.second, w);
  }
  return SigFunction.slope();
}
//--------------------------------------------------------------
std::vector<SetScores> MergeSym(hkl_unmerge_list& ref_list,
                                const ResoRange& ResRange,
                                const Normalise& NormRes,
				const Chirality& ChiralFlag,
				int& maxMult, const int& reducedMult);
//--------------------------------------------------------------
std::vector<SCsignificance> CCsum(const std::vector<SetScores>& scores);
//--------------------------------------------------------------
// Assess E(CC) based on E(CC) from unrelated scores and from score
//  for identity operator
double EstimateECC(const double& SigFac, const double& ECC0,
		   const SCsignificance& SCsig1);
//--------------------------------------------------------------
template <class T> void ScoreSig(const double& SigFac, const double& ECC,
                                 std::vector<SCsignificance>& SCsig)
//
// On entry:
//   SigFac  factor used to determine sd(CC) based on sample number N
//          sd(CC) = Sigfac/Sqrt(N)
//   ECC    expectation value for CC E(CC) if symmetry is present
//   SCsig  contains for each symmetry item tested, score and
//          the number of pairs in the sample
//
// On exit:
//   SCsig  contains for each symmetry item tested,
//          the estimated Score & SD(Score) for this number of pairs
//          in the sample
{
  const int MaxGroup = 200;  // flatten off SD estimate at this group size
			     // to stop SD getting saller & smaller for large N
  if (SCsig.size() == 0) return;

  for (size_t i=0;i<SCsig.size();i++) {   // Loop symmetry elements
    int Ngroup = Min(MaxGroup, SCsig[i].Nsample());  // size of sample
    if (Ngroup > 2) {
      float SD = SigFac/sqrt(float(Ngroup));  
      SCsig[i].SetUnrSC(SD);  // unrelated SD
      // set related sd = unrelated SD
      float SDrel = SD;
      SCsig[i].SetRelSC(SDrel);
      SCsig[i].SetTrueScore(ECC);  // E(CC) if element present
      SCsig[i].SetScoreRange(-1.0, +1.0);   // set range of possible CC for truncation
    }
  }
}
//--------------------------------------------------------------
std::vector<PGscore> ScoreSubGroups(const hkl_symmetry& Symmetry,
				    const hkl_symmetry& OrigSymm,
				    const RefListType& HklinIsMerged,
				    const std::string& OrigName,
				    const Scell& OriginalCell,
				    const ReindexOp& reindex_op,
				    const double& CCsigFac,
				    const double& ECC,
				    const std::vector<SetScores>& scores,
				    const std::vector<SCsignificance>& CCsig,
				    const int& AllowI2);
//--------------------------------------------------------------
void  SetConfidenceScores(std::vector<PGscore>& SGscores);
// Store "confidence" values for each group relative to next one

#endif
