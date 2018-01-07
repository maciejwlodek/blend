// pointless.cpp
//
//==================================================================
//
// Pointless
//
// Calculate statistics for merging observations related by each
// symmetry element belonging to the highest lattice symmetry
// compatible with the cell dimensions.
//
// Use this to determine the "best" pointgroup & spacegroup
//
// Phil Evans 2003-2006 etc
// MRC Laboratory of Molecular Biology, Cambridge, UK
//
//==================================================================


#include "scala.hh"
#include "pointless.hh"
#include "latsym.hh"
#include "pointgroup.hh"
#include "version.hh"
#include "matvec_utils.hh"
#include "sysabszones.hh"
#include "sysabsstats.hh"
#include "printthings.hh"
#include "icering.hh"
#include "observationflags.hh"
#include "util.hh"
#include "io_files.hh"
#include "outputunmergeddata.hh"
#include "xds_unmerge.hh"
#include "sca_unmerge.hh"
#include "makehkl_listscompatible.hh"
#include "normaliseunmerge.hh"
#include "readallhklinfiles.hh"
#include "testindex.hh"
#include "testlauegroup.hh"
#include "choosesolution.hh"
#include "spacegroupreindex.hh"
#include "writesummedlist.hh"
#include "checkcompatiblesymmetry.hh"
#include "string_util.hh"
#include "filetype.hh"
#include "completeness.hh"
#include "timer.hh"
#include "twintest.hh"

#if defined _MSC_VER
#include <io.h>
#define isatty _isatty
#define STDIN_FILENO _fileno(stdin)
#endif

using namespace scala;
using phaser_io::LOGFILE;
using phaser_io::LXML;
using phaser_io::RESULT;

int Verbose()
{
  // Returns = 0 silent (default) > 0 verbosity level
  static int verbose = -1;
  char* env;

  if (verbose < 0) {
    env = getenv("SCALA_VERBOSE");
    if (env == NULL) {
      verbose = 0;
    } else {
      verbose = +1;
    }
  }
  return verbose;
}
//--------------------------------------------------------------
bool IsOnline()
// Returns true if stdin is not connected to a file
{
  if (isatty(STDIN_FILENO)== 0)
    return false;
  else
    return true;  
}

//--------------------------------------------------------------
IsigI average_Iset(const reflection& this_refl, const PairSet& set)
// Average I for set of observations
// Partials must have been added first
// Variance weighted & no scaling or selections for now
{
  int n = 0;
  Rtype w;

  Rtype I = 0.0;
  Rtype sigI = 0.0;

  std::vector<int> list = set.ObsIndices();

  for (size_t i = 0; i < list.size(); i++) {
    observation this_obs = this_refl.get_observation(list[i]);
      if (this_obs.IsAccepted()) {
	n += 1;
	w = 1./(this_obs.sigI()*this_obs.sigI());
	I += w * this_obs.I();
	sigI += w;
      }
  }
  if (n > 0) {
    I = I/sigI;
    sigI = sqrt(1.0/sigI);
  }
  return IsigI(I,sigI);
}
//--------------------------------------------------------------
std::vector<SetScores> MergeSym(hkl_unmerge_list& ref_list,
                                const ResoRange& ResRange,
                                const Normalise& NormRes,
				const Chirality& ChiralFlag,
				int& maxMult, const int& reducedMult)
// Merge observations in pairs related by symmetry elements,
//   accumulate statistics
// Correlation-coefficient based on "corrected" ("normalised")
//  intensities
//
// On entry:
//  ref_list    the data
//  ResRange    resolution range & binning information
//  NormRes     normalisation object
//  ChiralFlag      chirality flag
//      If ChiralFlag != CHIRAL then exclude zonal reflections
//  reducedMult maximum multiplicity
//
// On exit:
//  maxMult       maximum reflection multiplicity ...
//  reducedMult   ... reduced to this
{
  observation this_obs1, this_obs2;
  int ksymelmt, Nsets;
  bool found;
  IsigI IsigImean;
  maxMult = 0;
  int irun2;
  float time2;
  float fac2;

  // Number of symmetry elements
  hkl_symmetry Symmetry = ref_list.symmetry();
  int Nelement = Symmetry.Nelement();
  std::vector<int> listsymelmt(Symmetry.NsymP());
  int nm;
  // A set of scores for each symmetry element
  std::vector<SetScores> AllScores(Nelement);
  // **** Open files for dumping pairs for debug analysis ****
  //      One file for each symmetry element
  std::vector<FILE*> file(Nelement);
  if (Verbose()>0)  {
    for (int i=0;i<Nelement;i++) {
      int w = 1;         // number of digits in element number
      if (i+1 >=10) w=2; // no more than 99 elements
      std::string name = "related_pairs_"+String(i+1,w)+".pair";
      file[i] = fopen(name.c_str(), "w");
      if (file[i] == NULL)
	Message::message(Message_fatal("Can't open file "+name));
    }
  }
  // ****
  // list of pairs 
  std::vector<PairSet> sets;

  reflection this_refl;

  double weight = 1.0;
  ref_list.rewind();

  int Maxmult = reducedMult + 1;

  // loop all reflections
  while (ref_list.next_reflection(this_refl) >= 0)   {
    int Nobs = this_refl.num_observations();
    // Some statistics (only if more than one observation)
    if (Nobs > 1) {
      sets.clear();
      if (ChiralFlag != CHIRAL)
	if (!this_refl.hkl().IsGeneral()) continue;
      bool IsCentric = Symmetry.is_centric(this_refl.hkl());
      float sSqr = this_refl.invresolsq();
      // For very large multiplicity (> Maxmult), select a subset of pairs
      int linc = Nobs/Maxmult + 1;
      maxMult = Max(maxMult, Nobs);
      // Check all pairs, accumulate into sets by symmetry element
      for (int lobs1 = 0; lobs1 < Nobs-1; lobs1+=linc) {  //**
	this_obs1 = this_refl.get_observation(lobs1);
	if (this_obs1.IsAccepted()) {
	  int irun1 = this_obs1.run();
	  float time1 = this_obs1.time();
	  // multiplying normalisation factor
	  float fac1 = NormRes.Corr(sSqr, irun1, time1);
	  for (int lobs2 = lobs1+1; lobs2 < Nobs; lobs2+=linc) { //**
	    this_obs2 = this_refl.get_observation(lobs2);
	    if (this_obs2.IsAccepted()) {
	      irun2 = this_obs2.run();
	      time2 = this_obs2.time();
	      fac2 = NormRes.Corr(sSqr, irun2, time2);
	      // Identify symmetry element relating this pair
	      // For centric reflections more than one symmetry operation
	      // relates a pair, so use all of them, but weight down by
	      // the number found
	      if (IsCentric) {
		listsymelmt = Symmetry.get_symelmt
		  (this_obs1.hkl_original(), this_obs2.hkl_original());
		nm = listsymelmt.size();
	      }
	      else {
		// acentric
		listsymelmt[0] = Symmetry.get_symelmt(this_obs1.Isym(),
						 this_obs2.Isym());
		nm = 1;
	      }
	      weight = 1./double(nm);  // weight
	      for (int k1=0;k1<nm;k1++) {
		ksymelmt = listsymelmt[k1];
		if (ksymelmt >= 0) {
		  // add to list of sets
		  if (sets.size() == 0) {
		    sets.push_back(PairSet
				   (lobs1,lobs2,ksymelmt,weight,fac1,fac2));
		  } else {
		    Nsets = sets.size();
		    found = false;
		    for (int k2=0;k2<Nsets;k2++) {
		      if (sets[k2].AddPair(lobs1,lobs2,ksymelmt,fac1,fac2)) {
			found = true;
			break;
		      }
		    }
		    if (!found) {
		      sets.push_back(PairSet
				     (lobs1,lobs2,ksymelmt,weight,fac1,fac2));
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      // Loop sets
      Nsets = sets.size();
      for (int k=0;k<Nsets;k++)
	{
	  // Get Imean
	  IsigImean = average_Iset(this_refl, sets[k]);
	  // Symmetry element number
	  ksymelmt = sets[k].SymElement();
	  // Rfactors
	  AllScores[ksymelmt].addIstats(IsigImean, this_refl, sets[k]);
	  // Correlation coefficients & mean square difference
	  AllScores[ksymelmt].addPairStats(this_refl, sets[k], NormRes);
	  // **** Dump pairs for testing ****
	  if (Verbose()>0) 
	    AllScores[ksymelmt].DumpPairs(file[ksymelmt], this_refl,sets[k], NormRes);
	}  // end loop sets
    }
  }  // end reflection loop
  
  return AllScores;
}
//--------------------------------------------------------------
void DumpPairList(const std::string& name,
                  const std::vector<IsPair>& pair_list)
// Dump pair list to file
{
  FILE* file = fopen(name.c_str(), "w");
  if (file == NULL)
    Message::message(Message_fatal("Can't open file "+name));

  float w = 1.0;
  int nl = pair_list.size();

  for (int i=0;i<nl;i++)
    {
      fprintf(file, " %9.3f  %9.3f %9.3f  %9.3f %8.3f\n",
              pair_list[i].Is1().I(), pair_list[i].Is1().sigI(),
              pair_list[i].Is2().I(), pair_list[i].Is2().sigI(), w); 
    }
  fclose(file);
}


//--------------------------------------------------------------
bool NextPair(std::vector<IKode>& IKlist, int& l, int& m)
// Get next pair from list, with different hkl-codes
// On entry:
//  IKlist   IKode list
//           initially all HKL-codes are positive numbers
//   l  1st member from last time here
//   m  2nd member from last time here
// On exit:
//  IKlist   items used flagged with negative hkl-code
//    l,m pair found
//  Returns false if list exhausted
{
  int Nlist = IKlist.size();
  if (l >= Nlist) return false;

  while (l < Nlist && IKlist[l].kode < 0) l++;
  if (l >= Nlist) return false;

  int lkode = IKlist[l].kode;

  m = l+1;
  while (m < Nlist && ((IKlist[m].kode < 0) || 
		       (IKlist[m].kode == lkode))) m++;
  if (m >= Nlist) return false;
  IKlist[l].kode = -lkode;
  IKlist[m].kode = -IKlist[m].kode;
  return true;
}
//--------------------------------------------------------------
int NormalisePair(std::vector<IsPair>& pair_list)
// Normalise pair list, omitting pairs where normalised
// intensity E^2 is too large for either pair
// Returns number of rejected pairs
{
  int k=0;
  int N = pair_list.size();

  float E2max = 20.0;

  for (int i=0;i<N;i++)
    {
      IsigI Is1 = pair_list[i].Is1norm();
      IsigI Is2 = pair_list[i].Is2norm();
      //      if (NormRes.NotTooLarge(Is1) && NormRes.NotTooLarge(Is2))
      if (Is1.I() < E2max && Is2.I() < E2max)
        { 
	  // Reset normalisation factors to 1.0
          pair_list[k].Set(Is1, 1.0, Is2, 1.0);
          k++;
        }
    }
  pair_list.resize(k);
  return (N-k);
}
//--------------------------------------------------------------
std::vector<IsPair> NullStats(hkl_unmerge_list& ref_list,
                              const ResoRange& ResRange,
                              const Normalise& NormRes,
			      double& ECC0,
			      const bool& print,
                              phaser_io::Output& output)
// Make list of unrelated pairs within same resolution bin
// as a control to compare with symmetry elements
// No print if false
// 
// Returns:
//  pair_list  contains vector of pairs, normalised by resolution
//              (using NormRes)
//              Length of list is truncated to MAXPAIRLIST to save memory
//              (defined here)
//  ECC0        estimated E(CC) allowing for errors (<1.0)
{
  const int MAXPAIRLIST = 20000;
  std::vector<IsPair> pair_list;

  reflection this_refl;
  observation this_obs;
  IsigI IsigImean, Is;
  int kode;  // hkl code, packed reduced hkl
  // IKode is intensity + hkl packed code + InvResolution^2
  std::vector<IKode> IKlist;
  // Objects for sigma(E^2) and sigma(sd(E^2))
  MeanSD SigmaE2;
  MeanSD MeanSqSdE2;
  MeanSD SigmaSdE2;
  // Correct SDs by "typical" SD factors
  SDcorrection SDC(2.0, 0.0, 0.03);
  ref_list.rewind();

  // loop all reflections
  while (ref_list.next_reflection(this_refl) >= 0) {
      kode = this_refl.hkl().code();
      // Resolution 
      float sSqr = this_refl.invresolsq();      
      // Next valid observations
      while (this_refl.next_observation(this_obs) >= 0)
        {
          // Store intensity, packed version of reduced hkl,
          // & resolution (sSqr)
	  float time = this_obs.time();
          IKlist.push_back(IKode(this_obs.I_sigI(), kode, sSqr,
				 this_obs.run(), time));
	  // Sums for sigmas, using normalised E with corrected sigmas
	  Is = NormRes.apply(SDC.Corrected(this_obs.I_sigI()), sSqr,
			     this_obs.run(), time);
	  SigmaE2.Add(Is.I());
	  SigmaSdE2.Add(Is.sigI());
	  MeanSqSdE2.Add(Is.sigI()*Is.sigI());
	}
    }

  // Sort list on resolution
  std::sort(IKlist.begin(), IKlist.end());

  int l = 0;
  int m = 0;

  std::vector<float> sSqrList;

  // Make list of intensity pairs roughly matched in resolution
  // (from list sorted on resolution) but with different indices
  double sum_dsSqr = 0.0;
  double max_dsSqr = 0.0;
  int n_dsSqr = 0;

  pair_list.clear();

  while (NextPair(IKlist,l,m)) {
    // Average resolution
    float sSqr = 0.5*(IKlist[l].sSqr+IKlist[m].sSqr);
    sSqrList.push_back(sSqr);
    pair_list.push_back
      (IsPair(IKlist[l].Is, NormRes.Corr(sSqr, IKlist[l].irun, IKlist[l].time), 
	      IKlist[m].Is, NormRes.Corr(sSqr, IKlist[m].irun, IKlist[m].time)));
    // resolution difference
    sum_dsSqr += std::abs(IKlist[l].sSqr-IKlist[m].sSqr);
    max_dsSqr = Max(max_dsSqr, std::abs(IKlist[l].sSqr-IKlist[m].sSqr));
    n_dsSqr++;
  }
  // Normalise (roughly) intensities based on resolution
  //  Apply same normalisation factor to both members of pair
  int Nrej = NormalisePair(pair_list);
  if (print && Nrej > 0)
    output.logTabPrintf(0,LOGFILE,"\n%8d pairs rejected for E^2 too large\n",Nrej);
  if (pair_list.size() <= 0) {
    Message::message(Message_fatal("All reflection pairs rejected"));
  }

  // **** Write pairlist to file for testing ****
  if (Verbose()>0) DumpPairList("NullPair.pair", pair_list);

  // Scramble intensity pairs to mix them up at different resolutions
  std::random_shuffle(pair_list.begin(), pair_list.end());

  // Truncate length of list If unnecessarily long
  if (int(pair_list.size()) > MAXPAIRLIST) pair_list.resize(MAXPAIRLIST);


  // Overall correlation coefficient
  int Npairs = pair_list.size();
  correl_coeff cc_overall;
  const float w = 1.0;
  for (int i=0;i<Npairs;i++) {
    cc_overall.add(pair_list[i].Is1().I(), pair_list[i].Is2().I(), w);
  }
  if (print) 
    output.logTabPrintf(1,LOGFILE,
			"\n\nOverall CC for %5d unrelated pairs: %7.3f  N= %d\n\n",
                       Npairs,
		       cc_overall.result().val, cc_overall.result().count);
  float   sigE2 = SigmaE2.SD();
  float VarDelE2 = MeanSqSdE2.Mean();
  //^
  //  std::cout 
  //    << "SigmaE2 " << SigmaE2.format()
  //    << "\nSigmaSdE2 " << SigmaSdE2.format()
  //    << "\nSigmaSqSd " << MeanSqSdE2.format() << "\n";
  // var(E2)/(var(E2) + var(delE2))
  ECC0 = (sigE2*sigE2)/(sigE2*sigE2+VarDelE2);
  //^
  //  output.logTabPrintf(1,LOGFILE,
  //   " sigma(E^2) %6.3f  Variance(E^2) %6.3f  VarianceDelta(E^2) %6.3f\n",
  //	      sigE2, sigE2*sigE2, VarDelE2);
  if (print) 
    output.logTabPrintf(1,LOGFILE,
			" Estimated expectation value of true correlation coefficient E(CC) = %6.3f\n\n",
			ECC0);

  return pair_list;
}

//--------------------------------------------------------------
// Assess E(CC) based on E(CC) from unrelated scores and from score
//  for identity operator
double EstimateECC(const double& SigFac, const double& ECC0,
		   const SCsignificance& SCsig1)
// Assess E(CC) based on E(CC) from unrelated scores and from score
//  for identity operator (as long as there are a significant
//  number of identity observations
//
// On entry:
//   SigFac  factor used to determine sd(CC) based on sample number N
//          sd(CC) = Sigfac/Sqrt(N)
//   ECC0    expectation value for CC E(CC) if symmetry is present
//   SCsig1 for identity operator, contains score and
//          the number of pairs in the sample
// 
// Returns updated ECC
{
  const int MIN_SAMPLE = 10;
  const double MIN_SD = 0.05;   // minimum SD of CC for weighting
  // Number of observations to estimate SDunrelated E(CC)
  const double NsdUnrelated = 200;
  // Guess at SD of ECC from unrelated pairs
  double SD0 = Max(MIN_SD, SigFac/sqrt(NsdUnrelated));
  double w0 = 0.0;  // weight
  if (SD0 > 0.0001) {
    w0 = 1./(SD0*SD0);  // weight
  }
  // Sigma CC(identity)
  double w1 = 0.0;
  if (SCsig1.Nsample() > MIN_SAMPLE) {
    double SD1 = Max(MIN_SD, SigFac/sqrt(double(SCsig1.Nsample()))); 
    w1 = 1./(SD1*SD1);  // weight
  }
  double ECCestimate =(w0*ECC0 + w1*SCsig1.SC())/(w0+w1); 
  return  ECCestimate;
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
				    const int& AllowI2)
// Create list of all possible sub-pointgroups of lattice group
// Combine symmetry element scores into scores for each sub-pointgroup
// Rank subgroups
//
// On entry:
//  Symmetry    lattice symmetry
//  OrigSymm    original symmetry
//  HklinIsMerged  NONE, MERGED, UNMERGED  file merge type
//  OrigName    original spacegroup name
//  OriginalCell original cell before any symmetry constraints
//               were applied
//  reindex_op  reindex operator applied from original cell to
//                    lattice cell
//  CCsigfac    sigma factor for CC
//  ECC          estimated CC for symmetry "true"
//  scores      scores for each symmetry element
//  AllowI2     =  0 don't allow I2 setting of centred monoclinic
//              = +1 allow it
//              = -1 force it
{
  const bool DEBUG = false;
  int Nelement = Symmetry.Nelement();
  ASSERT(int(scores.size()) == Nelement);
  CCtbxSym::PointGroup OrigLG(OrigName);  // Original Laue group

  std::vector<CCtbxSym::PointGroup> subgroups = scala::GetSubGroups(Symmetry);

  if (HklinIsMerged == MERGED) {
    // file is merged, remove all subgroups which are subgroups of original group
    //  since these are implicit
    hkl_symmetry mergeSymm = OrigSymm;
    RemoveImplicitSubgroups(subgroups, mergeSymm, reindex_op);
  }
  //^
  //  std::cout << "=== Subgroup size " <<subgroups.size()  <<"\n";
  //^-
  std::vector<PGscore> SubgroupScores;

  // We have a list of "subgroups" constructed from
  // symmetry elements - loop round subgroups
  for (size_t k=0; k<subgroups.size();k++)  {
    // List of symmetry elements belonging to this subgroup
    std::vector<int> elements = subgroups[k].Elements();
    
    SetScores SGscore;
    for (size_t i=0;i<elements.size();i++) {  // loop elements
      // Add together scores belonging to all elements
      // for this subgroup,
      // Now always identity even if triclinic
      //^ testing: always include identity
      //^	  if (subgroups[k].crystal_system() == TRICLINIC  ||
      //^	      ! Symmetry.IsElementIdent(elements[i]))
      SGscore += scores[elements[i]];
    }
    // Store combined score
    PGscore s(subgroups[k]);
    s.SetScore(SGscore);
    SubgroupScores.push_back(s);
  }
  ASSERT (SubgroupScores.size() == subgroups.size());  
  // Significance scores for & against for each subgroup
  // Correlation coefficients
  std::vector<SCsignificance> CCfor(SubgroupScores.size());
  std::vector<SCsignificance> CCagainst(SubgroupScores.size());
  // Mean squared differences
  //  std::vector<SCsignificance> MSDfor(SubgroupScores.size());
  //  std::vector<SCsignificance> MSDagainst(SubgroupScores.size());

  // "Original" Laue group is either identical or failing that one with the same
  // reference Laue group
  int sameSGindex = -1;
  bool OrigLGfound = false;

  for (size_t k=0; k<SubgroupScores.size();k++) {  // loop subgroups
    // Test initial cell against symmetry constraints for
    // this pointgroup
    // store change-of-basis from original indexing
    // to lattice symmetry
    // Cell is in original indexing convention (basis)
    //  so must call this after changing basis
    //^
    //    if (DEBUG) {
    //          std::cout << "\n-----\nOriginal cell " << OriginalCell.format() << " "
    //		    << reindex_op.as_hkl() << "\n";
    //	  SubgroupScores[0].dump();
    //	  std::cout << "\nDumped\n";
    //}
    //^-

    double delta = SubgroupScores[k].SetCell(OriginalCell.UnitCell(),reindex_op, AllowI2);
    delta = delta;
    // Is this the original Laue group?
    if (OrigLG.Equals(SubgroupScores[k])) {
      // Yes
      SubgroupScores[k].SetOrig(true);
      OrigLGfound = true;
      if (DEBUG) {
	std::cout << " **** Orig group " << k << "\n";
      }
    } else {
      // Not yet
	SubgroupScores[k].SetOrig(false);
	if (OrigLG.EqualsRef(SubgroupScores[k])) {
	  //	  std::cout << "Same group as orig " << k << " "
	  //		    << SubgroupScores[k].RefSGreindex().as_hkl() << "\n";
	  sameSGindex = k;  // same reference group, fallback
	}
    }
    //^
    if (DEBUG) {
      std::cout << "\nSubGroupScores: Original Laue group " << OrigLGfound
		<< " " << OrigLG.RefLGname()
		<< " subgroup " << k << " " << SubgroupScores[k].RefLGname() 
		<< "\n";
      //^ SubgroupScores[k].dump();
    } //^-
    
    // Sum scores for each symmetry element in this
    // subgroup, and for those not in this subgroup
    std::vector<int> symel = SubgroupScores[k].Elements();
    SetScores score_for;
    SetScores score_against;
    
    if (DEBUG) {
      printf("\n*** Subgroup: %s\n",
	     subgroups[k].RefPGname().c_str());
    }
    
    for (int i=0;i<Nelement;i++) {   // loop all elements, including identity
      // We must include the identity, otherwise triclinic gets too high a score
      // is this element in this subgroup?
      if (std::find(symel.begin(),symel.end(),i) != symel.end()) {
	// yes it is
	score_for += scores[i];
	SubgroupScores[k].AddAveragedZCCfor(CCsig[i].Z());
	SubgroupScores[k].StoreCCfor(CCsig[i]); // add CC+ to list
	if (DEBUG) {
	  printf("    Element in group, %3d: CC = %7.3f N %6d R = %7.3f N %6d\n",
		 i+1,
		 scores[i].OverallCC().result().val,
		 scores[i].OverallCC().result().count,
		 scores[i].OverallRfactor().result().val,
		 scores[i].OverallRfactor().result().count);
	}
      } else {
	// oh no it isn't
	score_against += scores[i];
	SubgroupScores[k].AddAveragedZCCagainst(CCsig[i].Z());
	SubgroupScores[k].StoreCCagainst(CCsig[i]); // add CC- to list
	if (DEBUG) {
	  printf("Element not in group, %3d: CC = %7.3f N %6d R = %7.3f N %6d\n",
		 i+1, 
		 scores[i].OverallCC().result().val,
		 scores[i].OverallCC().result().count,
		 scores[i].OverallRfactor().result().val,
		 scores[i].OverallRfactor().result().count);
	}
      }
    }
    // Store CC for k'th subgroup
    CCfor[k] = SCsignificance(score_for.OverallCC().result().val,
			      score_for.OverallCC().result().count);
    CCagainst[k] = SCsignificance(score_against.OverallCC().result().val,
				  score_against.OverallCC().result().count);
    
    if (DEBUG) {
      printf("\nCC+ %5.2f %6d CC- %5.2f %6d  R+ %5.2f %6d R- %5.2f %6d\n",
	     score_for.OverallCC().result().val,
	     score_for.OverallCC().result().count,
	     score_against.OverallCC().result().val,
	     score_against.OverallCC().result().count,
	     score_for.OverallRfactor().result().val,
	     score_for.OverallRfactor().result().count,
	     score_against.OverallRfactor().result().val,
	     score_against.OverallRfactor().result().count);
    }
    
    SubgroupScores[k].CalcAver();  // Calculate RMS scores for sorting
    
    // & Rfactors
    SubgroupScores[k].SetRfactor(score_for.OverallRfactor(), score_against.OverallRfactor());
  }  // end loop subgroups
  // Try to set original LG if not done
  if (!OrigLGfound && sameSGindex >= 0) {
    SubgroupScores[sameSGindex].SetOrig(true);
  }
  // Calculate significance from unrelated pairs
  //   updates CCfor & CCagainst
  ScoreSig<correl_coeff>(CCsigFac, ECC, CCfor);
  ScoreSig<correl_coeff>(CCsigFac, ECC, CCagainst);
  //
  double total_likelihood = 0.0;
  // Store significance values for each subgroup
  for (size_t k=0; k<SubgroupScores.size();k++) {
    SubgroupScores[k].SetCC(CCfor[k], CCagainst[k]); // store CCfor & CCagainst
    //      SubgroupScores[k].SetMSD(MSDfor[k], MSDagainst[k]);
    SubgroupScores[k].CalcLikelihood();   // calculate likelihood
    total_likelihood += SubgroupScores[k].Likelihood();
  }
  // Normalise likelihoods to add up to 1.0
  for (size_t k=0; k<SubgroupScores.size();k++) {
    SubgroupScores[k].SetLikelihood(SubgroupScores[k].Likelihood()/total_likelihood);
  }

  // Sort them by rank on Likelihood
  // & reverse so that best is first
  std::sort(SubgroupScores.begin(),SubgroupScores.end());
  std::reverse(SubgroupScores.begin(),SubgroupScores.end());

  return SubgroupScores;
}
//--------------------------------------------------------------
std::vector<SCsignificance> CCsum(const std::vector<SetScores>& scores)
//  Make vector of CCs for significance test
{
  std::vector<SCsignificance> CCsig(scores.size());
  for (size_t k=0;k<scores.size();k++) {
    CCsig[k] = SCsignificance(scores[k].OverallCC().result().val,
			      scores[k].OverallCC().result().count);
  }
  return CCsig;
}
//--------------------------------------------------------------
void  SetConfidenceScores(std::vector<PGscore>& SGscores)
// Store "confidence" values for each group relative to next one
{
  //  double LGprob_first = SGscores[0].Likelihood();
  double LG_score = -1.0;
  double nextLG_score = -1.0;
  double LGconfidence;

  if (SGscores.size() > 1)  {
    for (size_t i=0;i<SGscores.size()-1;i++) {
      LG_score = SGscores[i].Likelihood();
      nextLG_score = SGscores[i+1].Likelihood();
      if (LG_score > 0.0 && nextLG_score > 0.0) {
	double lgc = LG_score*(LG_score - nextLG_score);
	LGconfidence = sqrt(std::abs(lgc));
	if (lgc < 0.0) LGconfidence = -LGconfidence;
      } else {LGconfidence = 0.0;}
      SGscores[i].SetConfidence(LGconfidence);
    }
    LGconfidence = 0.0;
  } else {LGconfidence = 1.0;}  // Only one Laue group

  SGscores[SGscores.size()-1].SetConfidence(LGconfidence); // last one
}
//--------------------------------------------------------------
void PrintSummary(const Solution& solution,
		  const std::string& groupflag,
		  const Scell& cell,
		  const ResoRange& ResRangeLG, const ResoRange& ResRangeFile,
		  const int& Nbatches,
		  const Twintest& TwinTest,
		  phaser_io::Output& output)
//  ResRangeLG      resolution range used for Laue group test, perhaps limited
//                  unset if explcit resolution range given
//  ResRangeFile    resolution range used for systematic absence check, full range
{
  output.logTab(0,LXML,"\n<BestSolution Type=\""+groupflag+"group\">"); 
  output.logTab(1,LXML,"<GroupName>"+solution.Name()+"</GroupName>\n");
  output.logTab(1,LXML, solution.Reindex().as_hkl_XML()+"\n");
  output.logTab(1,LXML, solution.Reindex().as_XML());
  output.logTabPrintf(1,LXML,
		       "<Confidence>%9.3f</Confidence>\n", solution.Confidence());
  output.logTabPrintf(1,LXML,
		       "<TotalProb>%7.3f</TotalProb>\n", solution.TotalProb());
  output.logTab(0,LXML,"</BestSolution>\n"); 

  // Print summary
  //r  output.logTab(0,LOGFILE,"\n<!--SUMMARY_BEGIN--> $TEXT:Result: $$ $$\n\n");
  output.logTab(0,RESULT,"Best Solution     "+groupflag+" group "+solution.Name()+"\n\n"); 
  output.logTab(1,RESULT,"Reindex operator: "+
		StringUtil::CentreString(solution.Reindex().as_hkl(),42)+"\n");
  output.logTabPrintf(1,RESULT, "Laue group probability:          %8.3f\n",
		       solution.TotalProb());
  output.logTabPrintf(1,RESULT, "Confidence:                      %8.3f\n",
		       solution.Confidence());
  output.logTabPrintf(0,RESULT, "\n   Unit cell: %s\n",
		      cell.change_basis(solution.Reindex()).format().c_str());
  if (ResRangeLG.isSet()) {
    output.logTabPrintf(0,RESULT, "\n%8.2f to %6.2f   - Resolution range used for Laue group search\n",
			ResRangeLG.ResLow(), ResRangeLG.ResHigh());
    output.logTabPrintf(0,RESULT, "\n%8.2f to %6.2f   - Resolution range in file\n",
			ResRangeFile.ResLow(), ResRangeFile.ResHigh());
  } else {
    output.logTabPrintf(0,RESULT, "\n%8.2f to %6.2f   - Resolution range selected on input\n",
			ResRangeFile.ResLow(), ResRangeFile.ResHigh());
  }
  output.logTabPrintf(0,RESULT, "\n   Number of batches in file: %6d\n", Nbatches);

  // Write twin assessment
  TwinTest.WriteResultSummary(output);

  output.WriteResult();
  //r  output.logTab(0,LOGFILE,"\n$$ <!--SUMMARY_END-->\n\n");
}
//--------------------------------------------------------------
bool ConventionalSetting(const scala::PossibleSpaceGroup& group, const int& settingOption)
{
  // Return true if the 'conventional setting' is to be used
  //  Conventional setting is:
  //   1) for primitive orthorhombic, a<b<c
  //   2) for centred monoclinic, I2 if "better" than C2
  //
  //  settingOption     = 0 CELL-BASED = +1 SYMMETRY-BASED = -1 C2
  //                    <= 0 to keep primitive orthorhombic spacegroups in the Lauegroup setting
  //                              ie a < b < c
  //                    == 0 to allow C2 as possibly I2
  bool LGsetting = true;  // conventional setting (a<b<c) and allow I2
  if (group.IsPIForthorhombic() && settingOption > 0)
    {LGsetting = false;} // oP && SYMMETRY-BASED
  if (group.IsI2() && settingOption != 0)
    {LGsetting = false;}  // mC/I && C2
  return LGsetting;
}
//--------------------------------------------------------------
Solution ChooseBestGroup(const ChooseSolution& SolutionChoice,
			 std::vector<scala::PossibleSpaceGroup>& AllGroups,
			 const int& settingOption, const Scell& cell,
			 const char& origLatType,
			 const ResoRange& ResRangeLG, const ResoRange& ResRangeFile,
			 const int& Nbatches,
			 const Twintest& TwinTest,
			 const std::string& hklinSpacegroup,
			 phaser_io::Output& output)
// Select best group, space group or point group
// If there are multiple groups with the same score, use the Laue group
//
// On entry:
//  SolutionChoice  solution choice settings from input
//  AllGroups       possible space groups with scores
//  settingOption     = 0 CELL-BASED = +1 SYMMETRY-BASED = -1 C2
//                    <= 0 to keep primitive orthorhombic spacegroups in the 
//                         Lauegroup setting  ie a < b < c
//                    == 0 to allow C2 as possibly I2
//  cell            original cell
//  origLatType     original lattice type character
//  ResRangeLG        resolution range used for Laue group test, perhaps limited
//                   unset if specified explicitly
//  ResRangeFile    resolution range used for systematic absence check, full range
//                    
// Returns best solution 
{
  // If a solution has been explicitly chosen on input (CHOOSE), find it
  //  = -1 if no selection or not found
  int ChosenSolution = SolutionChoice.ChooseSpacegroup(AllGroups);
  if (SolutionChoice.ChooseSpaceGroup() && ChosenSolution < 0) {
    Message::message(Message_fatal("Cannot find the chosen spacegroup "+
				   SolutionChoice.Spacegroup()));
  }

  output.logTabPrintf(0,LOGFILE,
       "\n\n---------------------------------------------------------------\n\n");
  double prob_first = AllGroups[0].Prob();
  double LGprob_first = AllGroups[0].LaueGroupProb();
  int ntop = 1;
  int nextSG = -1;
  double nextLG_score = -1.0;
  int ig = -1;
  std::string groupflag = "space";
  std::string Sname = AllGroups[0].Refname();
  
  double confidence = 1.0;
  double LGconfidence = AllGroups[0].LaueGroupConfidence();
  bool LGsetting = false;

  if (ChosenSolution >= 0) {  //  ===  Solution specified
    // Solution chosen on input
    ig = ChosenSolution;
    // set flag true for keeping Laue group setting for P, I, F orthorhombic groups, or I2
    LGsetting = ConventionalSetting(AllGroups[ig], settingOption);

    Sname = SolutionChoice.Spacegroup();
    confidence = 1.0;
    LGconfidence = 1.0;
    if (SolutionChoice.IsRhombohedralR() && origLatType != 'R') {
      // We are choosing a rhombohedral R setting, and the original lattice was not
      // rhombohedral, work out reindexing
      AllGroups[ig].ReindexRhombohedral(false);
    }
    output.logTabPrintf(0,LOGFILE,
	"\nChoosing specified space group %s with reindex operator %s\n",
			Sname.c_str(),
			AllGroups[ig].Reindex(LGsetting).as_hkl().c_str());
  } else {  //  ===  No solution specified, so choose one
    // Choose solution
    for (size_t i=0;i<AllGroups.size();i++) {
      //^      std::cout << "Name " << i << " " << AllGroups[i].Refname() << "\n";; //^
      if (i > 0) {
	if (nextSG < 0) {nextSG = i;}
	if (std::abs(prob_first - AllGroups[i].Prob()) <= 0.001)
	  {ntop++;} // Count solutions with same score
	else {
	  if (nextLG_score < 0.0 && AllGroups[i].LaueGroupName() != AllGroups[0].LaueGroupName())
	    {nextLG_score = AllGroups[i].LaueGroupProb();}
	}
      }
    }
    // ntop (at least 1) have the same score
    bool sameLG = true;
    bool sameSG = true;
    bool enantiomorphic = false;
    if (AllGroups.size() > 1) {
      if (nextLG_score > 0.0) {
	double lgc = LGprob_first*(LGprob_first-nextLG_score);
	LGconfidence = sqrt(std::abs(lgc));
	if (lgc < 0.0) {LGconfidence = -LGconfidence;}
      }
    }
    
    if (ntop > 1) {
      output.logTab(0,LOGFILE,"\nChoosing between possible best groups:\n");
      output.logTab(0,LOGFILE,"\n    Space group     Point group            Reindex\n\n");
      // Do any have the identity reindex?
      // Do they all have the same Lauegroup & same spacegroup?
      for (int i=0;i<ntop;i++)	{
	LGsetting = ConventionalSetting(AllGroups[i], settingOption);
	if (AllGroups[i].Reindex(LGsetting).IsIdentity()) {
	  ig = i;
	}
	if (i>0 && AllGroups[i].LaueGroupName() != AllGroups[0].LaueGroupName()) {
	  sameLG = false;}
	if (i>0 && AllGroups[i].Sgnumber() != AllGroups[0].Sgnumber()) {
	  sameSG = false;
	}
	Sname = AllGroups[i].Name(LGsetting);
	output.logTabPrintf(0,LOGFILE,"%15s %15s %20s\n",
			    Sname.c_str(),
			    CCtbxSym::PointGroupName(Sname).c_str(),
			    AllGroups[i].Reindex(LGsetting).as_hkl().c_str());
      }
      bool matchingInputSpacegroup = false;
      if (sameLG && ntop == 2) {
	// If there are just two possibilities, check if they enantiomers
	// as this is acceptable
	if (CCtbxSym::EnantiomorphicGroups(AllGroups[0].Refname(), AllGroups[1].Refname())) {
	  enantiomorphic = true;
	  ig = 0;
	  //  Pick second one if same as input group
	  if (hklinSpacegroup == AllGroups[0].Refname()) {
	    ig = 0;
	    matchingInputSpacegroup = true;
	  }
	  if (hklinSpacegroup == AllGroups[1].Refname()) {
	    ig = 1;
	    matchingInputSpacegroup = true;
	  }
	  LGsetting = ConventionalSetting(AllGroups[ig], settingOption);
	  Sname = AllGroups[ig].Name(LGsetting);
	  if (AllGroups.size() > 2)
	    {nextSG = 2;}
	  else
	    {nextSG = -1;}
	}
      }
      if (sameLG) {
	// ig >= 0 found one with identity operator
	if (ig < 0) ig = 0;  // Use first operator if none of them are the identity
	if (enantiomorphic) {
	  groupflag = "space";
	  std::string inputmatch = "";
	  if (matchingInputSpacegroup) {
	    inputmatch = " and it matches the input file";
	  }
	  output.logTabPrintf(0,LOGFILE,
			      "\nSelecting space group %s as solutions are enantiomorphic%s\n",
			      Sname.c_str(), inputmatch.c_str());
	} else if (sameSG) {
	  // If there are two indistinguishable solutions with the same space group,
	  // we cannot tell the point group either, so use the Laue group
	  LGsetting = ConventionalSetting(AllGroups[ig], settingOption);
	  Sname = CCtbxSym::LaueGroupName(AllGroups[ig].Refname());
	  groupflag = "Laue";
	  output.logTabPrintf(0,LOGFILE,
      "\nSelecting Laue group %s as multiple solutions have indistinguishable settings\n    of the same space group\n",
			      Sname.c_str());
	  if (AllGroups.size() > 2)
	    {nextSG = 2;}
	  else
	    {nextSG = -1;}
	} else {
	  Sname = CCtbxSym::PointGroupName(AllGroups[ig].Refname());
	  groupflag = "point";
	  output.logTabPrintf(0,LOGFILE,
			      "\nSelecting point group %s as multiple space groups have the same score\n",
			      Sname.c_str());
	}
      }	else {
	// Different Laue groups, probably can't happen but trap anyway
	Message::message(Message_fatal("ERROR: cannot decide on which Laue group to select"));
      }
    } else {
      ig = 0;
      LGsetting = ConventionalSetting(AllGroups[ig], settingOption);
      Sname = AllGroups[ig].Name(LGsetting);
      if (SolutionChoice.ChooseLaueGroup()) {
	output.logTabPrintf(0,LOGFILE,
      "\nSelecting space group %s as there is a single space group with the highest score\n   in the chosen Laue group\n",
			    Sname.c_str());
      } else {
	output.logTabPrintf(0,LOGFILE,
		 "\nSelecting space group %s as there is a single space group with the highest score\n",
			    Sname.c_str());
      }
    }
    
    if (nextSG >= 0)
      {confidence = sqrt(prob_first*(prob_first-AllGroups[nextSG].Prob()));}
    else
      {confidence = LGconfidence;}
    
    output.logTabPrintf(0,LOGFILE,"\nSpace group confidence (= Sqrt(Score * (Score - NextBestScore))) = %8.2f\n",
			confidence);
    output.logTabPrintf(0,LOGFILE,"\nLaue group confidence  (= Sqrt(Score * (Score - NextBestScore))) = %8.2f\n",
			LGconfidence);
  } // === end choose solution

  double prob = AllGroups[ig].Prob();

  ReindexOp reindex = AllGroups[ig].Reindex(LGsetting);
  output.logTab(0,LXML,"\n<BestSolution Type=\""+groupflag+"group\">"); 
  output.logTab(1,LXML,"<GroupName>"+Sname+"</GroupName>\n");
  output.logTab(1,LXML, reindex.as_hkl_XML()+"\n");
  output.logTab(1,LXML, reindex.as_XML());
  output.logTabPrintf(1,LXML,
		       "<Confidence>%9.3f</Confidence>\n", confidence);
  output.logTabPrintf(1,LXML,
		       "<TotalProb>%7.3f</TotalProb>\n", prob);
  output.logTab(0,LXML,"</BestSolution>\n"); 

  // Print summary
  //r  output.logTab(0,LOGFILE,"\n<!--SUMMARY_BEGIN--> $TEXT:Result:$$ $$\n\n");
  if (ChosenSolution >= 0) {
    output.logTab(0,RESULT,"Chosen Solution:    "+groupflag+" group "+Sname+"\n\n"); 
  } else {
    output.logTab(0,RESULT,"Best Solution:    "+groupflag+" group "+Sname+"\n\n"); 
  }
  output.logTab(1,RESULT,"Reindex operator: "+
		StringUtil::CentreString(reindex.as_hkl(),42)+"\n");
  output.logTabPrintf(1,RESULT, "Laue group probability:          %8.3f\n", AllGroups[ig].LaueGroupProb());
  output.logTabPrintf(1,RESULT, "Systematic absence probability:  %8.3f\n", AllGroups[ig].SysAbsProb());
  if (AllGroups[0].Chiral() != CHIRAL)
    output.logTabPrintf(1,RESULT, "Centrosymmetric probability:     %8.3f\n", AllGroups[ig].CentroProb());
  output.logTabPrintf(1,RESULT, "Total probability:               %8.3f\n", prob);
  output.logTabPrintf(1,RESULT, "Space group confidence:          %8.3f\n", confidence);
  output.logTabPrintf(1,RESULT, "Laue group confidence            %8.3f\n", LGconfidence);
  output.logTabPrintf(0,RESULT, "\n   Unit cell: %s\n", cell.change_basis(reindex).format().c_str());
  if (ResRangeLG.isSet()) {
    output.logTabPrintf(0,RESULT, "\n%8.2f to %6.2f   - Resolution range used for Laue group search\n",
			ResRangeLG.ResLow(), ResRangeLG.ResHigh());
    output.logTabPrintf(0,RESULT, "\n%8.2f to %6.2f   - Resolution range in file, used for systematic absence check\n",
			ResRangeFile.ResLow(), ResRangeFile.ResHigh());
  } else {
    output.logTabPrintf(0,RESULT, "\n%8.2f to %6.2f   - Resolution range selected on input\n",
			ResRangeFile.ResLow(), ResRangeFile.ResHigh());
  }

  output.logTabPrintf(0,RESULT, "\n   Number of batches in file: %6d\n", Nbatches);
  // Write twin assessment
  TwinTest.WriteResultSummary(output);

  output.WriteResult();
  //r  output.logTab(0,LOGFILE,"\n$$ <!--SUMMARY_END-->\n\n");


  bool SGorPG = true; // space group
  if (groupflag == "point") SGorPG = false;
  return Solution(ig, Sname, SGorPG, reindex, prob, confidence);
}
//--------------------------------------------------------------
void WriteOutputFile(const std::string hklout_filename,
		     const Solution& BestSolution,
		     hkl_unmerge_list& hkl_list,
		     phaser_io::Output& output)
{
  std::string Sname = BestSolution.Name();
  std::string groupflag;

  if (BestSolution.SpaceOrPointGroup())
    {groupflag = "space";}
  else
    {groupflag = "point";}
  // Get point group name, use to make symmetry
  hkl_symmetry PGsymm(Sname);
  OutputUnmergedData(hkl_list, hklout_filename,
		     PGsymm, BestSolution.Reindex(), false,
		     hkl_list.Title(), groupflag, output);
}
//--------------------------------------------------------------
void MarkAcceptedSubgroup(std::vector<PGscore>& SGscores,
			  const ScoreAccept& Accept)
// Mark acceptable groups in SGscores, including original
{
  if (Accept.IfSet())
    {
      for (size_t k=0; k<SGscores.size();k++)
	{
	  // Don't now accept original Laue group
	  //	  if (SGscores[k].Original()) 
	  //	    {
	  //	      SGscores[k].SetAccept(true);  // mark as accepted
	  //	    }
	  //alternative// if (Accept.Accept(SGscores[k].CC_Zagain(), SGscores[k].CC_Zgain()))
	  if (Accept.Accept(SGscores[k].Likelihood()))
	    {
	      SGscores[k].SetAccept(true);  // mark as accepted
	    }
	}
    }
}
//--------------------------------------------------------------
void SetChiralityFromHKLIN(const hkl_unmerge_list& hkl_list, GlobalControls& GC)
// Set chirality in GC from hkl_list space group
{
  if (hkl_list.symmetry().IsCentro())
    {GC.set_Chiral(CENTROSYMMETRIC);}
  else if (hkl_list.symmetry().IsChiral())
    {GC.set_Chiral(CHIRAL);}
  else
    {GC.set_Chiral(NONCHIRAL);}
}
//--------------------------------------------------------------
bool ResetInputLaueGroup(std::vector<PGscore>& SGscores,
			 const std::string& hklinSpacegroup,
			 const bool& symSet,
			 phaser_io::Output& output)
// Insufficient data to determine Laue group, reset to HKLIN one if symSet true
// Return false if fails
{
  if (symSet) {
    output.logTab(0,LOGFILE,
		  "\n****** WARNING! WARNING! WARNING! ******\n\n");
    output.logTab(0,LOGFILE,
		  "Too few of the symmetry elements have been observed to distinguish the correct Laue group.\n");
    output.logTabPrintf(0,LOGFILE,
			"The Laue group from the HKLIN file %s will be used.\n",
			hklinSpacegroup.c_str());
    output.logTab(0,LOGFILE,
		  "Use the command LAUEGROUP ALL to override this\n");
    output.logTab(0,LOGFILE,
		  "\n****** WARNING! WARNING! WARNING! ******\n\n");
    // Find original one in list
    int k = -1;
    for (size_t i=0;i<SGscores.size();i++) {
      SGscores[i].SetAccept(false); // mark all as not accepted
      if (SGscores[i].Original())  {
	k = i;
      }
    }
    if (k < 0) {
      output.logTab(0,LOGFILE,
	    " Can't find original Laue group in subgroup list, revert to original\n");
      return false;
    }
    SGscores[k].SetAccept(true); // mark selected one as accepted
    SGscores[k].SetLikelihood(1.0); // store likelihood = 1.0
  } else {
    int k = -1;	
    // Find the group with the maximumm number of symmetry elements, highest symmetry
    int nsym = -1;
    int ntop = 0;
    double minDiff = 0.001;
    int ktopIdent = -1;  // top score with identity operator
    for (size_t i=0;i<SGscores.size();i++) {
      SGscores[i].SetAccept(false); // mark all as not accepted
      if (i>0 &&
	  ((SGscores[0].Likelihood() - SGscores[i].Likelihood()) < minDiff)) {
	// score not significantly different from first (top) one
	ntop++;
	if (ktopIdent < 0 && SGscores[i].RefSGreindex().IsIdentity()) {
	  ktopIdent = i;
	}
      }
      if (SGscores[i].NElements() > nsym) {
	nsym = SGscores[i].NElements();
	k = i; // maximum number of elements
      }
    }
    // Accept the top one with no reindex [h,k,l] as long as it has no worse score than the best
    if (ktopIdent >= 0) {
      k = ktopIdent;
    } else if (k > 0) {
      // Accept highest symmetry (k) as long as it has no worse score than the best
      if (SGscores[k].Likelihood() <  SGscores[0].Likelihood()) {
	k = 0;  // it is worse so accept first
      }
    } else if (k < 0) {
      k = 0;
    }
    SGscores[k].SetAccept(true); // mark highest symmetry one as accepted
    SGscores[k].SetLikelihood(1.0); // store likelihood = 1.0
    output.logTab(0,LOGFILE,
		  "\n****** WARNING! WARNING! WARNING! ******\n\n");
    output.logTab(0,LOGFILE,
		  "No symmetry information available in HKLIN file: Laue group "+
		  SGscores[k].RefLGname()+" accepted even though \n"+
	  "too few of the symmetry elements have been observed to distinguish the correct Laue group\n\n");
  }
  return true;
}
//--------------------------------------------------------------
std::vector<Zone>  MakeUniqueZones(const std::vector<PGscore>& SGscores,
				   const Chirality& Chiral)
{
  std::vector<Zone> AllUniqueZones;
  // Loop all possible (accepted) Laue groups (SGscores)
  //    Laue groups were constructed in lattice frame
  for (size_t k=0; k<SGscores.size();k++)  {
    if (SGscores[k].Accepted())	{
      //^
      //      std::cout << "\nAccepted Laue group " << SGscores[k].LGname()
      //		    << " reindex " << SGscores[k].RefSGreindex().as_hkl() << "\n";
      //^-
      // Reindex operator from constructor (lattice) to reference
      ReindexOp reindex = SGscores[k].SGreindex();
      // List of possible zones for this Laue group
      CCtbxSym::PointGroup PG(SGscores[k]);
      std::vector<Zone> SZones = SysAbsZones(PG, Chiral);
      if (SZones.size() > 0) {
	// Accumulate list of unique zones
	for (size_t iz=0;iz<SZones.size();iz++) { // loop new zones
	  // reindex from lattice to sub-Laue-group reference frame
	  SZones[iz].StoreReindex(reindex);
	  AllUniqueZones.push_back(SZones[iz]);
	  //^
	  //		std::cout << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n"
	  //			  << "Zone " << iz << " reindex " << reindex.as_hkl() 
	  //			  << "\n";
	  //		SZones[iz].dump();
	  //^-
	}
      }
    }
  }
  // Sort to put axes first
  std::sort(AllUniqueZones.begin(),AllUniqueZones.end());

  return AllUniqueZones;
}
//--------------------------------------------------------------
PointlessStatus Pointless(int argc, char* argv[], phaser_io::Output& output)
// Main program
{
  GlobalControls GC;
  std::string hklin_filename = "";
  ReflectionFileType  hklin_filetype;
  hklin_filetype.SetType("MTZ");
  PointlessTask pointlessTask = LAUE;  // default task
  int N_hklin = 0;
  scala::RefListType reflisttype = NONE;
  bool IsHKLOUT = false;
  bool success = true;
  ResoRange ResRange0;   // resolution range in file or from RESOLUTION command
  ResoRange ResRangeCut; // resolution range cut by normalisation
  char origLatType = ' ';      // original lattice type from file, for R v H testing
  Timer cputime;

  try {
    // Input from command line: optional HKLIN filename etc
    phaser_io::InterpretCommandLine CL(argc, argv);
    // Read input & store if not "online"
    phaser_io::InputAll input(IsOnline(), output);
    input.Analyse();

    // Store filenames
    scala::IO_files AllFiles;

    // Get XML file names if relevant first, because we need
    // this before printing header
    AllFiles.StoreFilename("XMLOUT", CL.getXMLOUT());
    AllFiles.StoreFilename("XMLOUT", input.getXMLOUT());
    bool IsXMLOUT = AllFiles.IsFile("XMLOUT");
    if (IsXMLOUT)
      {output.setXmlout(AllFiles.Filename("XMLOUT"));}

    output.logHeader(LOGFILE);
    PrintTitle(output);
    std::string outputstring;

    //  N_hklin is number of input files from input commands
    //  nf is number on the command line
    N_hklin = input.HKLINnumber() + input.XDSINnumber() + input.SCAINnumber();
    //  hklin, xdsin or scain
    // Note that all HKLIN files on command line will be considered as fileSeries 1
    // Put them into the input object, for compatibility with command input
    int nf = input.AddHKLINfilename(CL.getHKLIN());
    nf += input.AddXDSINfilename(CL.getXDSIN());
    nf += input.AddSCAINfilename(CL.getSCAIN());

    N_hklin += nf;    // total number of HKLIN filenames
    if (!input.AllowOutofSequenceFiles()) {
      // Optionally [by default] check file series for files out of time
      // sequence (MTZ files only)
      std::vector<std::string> rejectedFiles =
	input.RejectOutofSequenceHKLINFiles();
      if (rejectedFiles.size() > 0) {
	output.logTab(0,LOGFILE,
		      std::string("\n$TEXT:Warning:$$ $$\nWARNING:")+
		      " some files have been rejected "+
		      "because they are older than a previous file in the series.\n"+
		      "You can override this using the command \"ALLOW OUTOFSEQUENCEFILES\"\n"+
		      " or use the system command \"touch\" to update the file times\n"+
		      "$$\nRejected files:\n");
	for (size_t i=0;i<rejectedFiles.size();++i) {
	  output.logTab(1,LOGFILE, rejectedFiles[i]);
	}
	output.logTabPrintf(0,LOGFILE,"\n");
	N_hklin -= rejectedFiles.size();
      }
    }
    ReflectionFileType fltype;
    fltype.SetType("MTZ");
    int nfm = AllFiles.AddHKLINfilename(input.getHKLINnames(),
					input.getHKLINpxdNames(), fltype,
					input.getHKLINseriesNames(),
					input.getHKLINfileseries());
    nfm = nfm;
    // other types of XDS files will be distinguished later
    fltype.SetType("XDS");
    int nfx = AllFiles.AddHKLINfilename(input.getXDSINnames(),
					input.getXDSINpxdNames(), fltype);
    nfx = nfx;
    // other file types will be distinguished later
    fltype.SetType("SCA");
    int nfc = AllFiles.AddHKLINfilename(input.getSCAINnames(),
					input.getSCAINpxdNames(), fltype);
    nfc = nfc;

    // If there are multiple files defined on the command line,
    // or on commands with wild-cards,
    // assume that they are consistently indexed, but not necessarily between series
    int ASI = input.ASSUMESAMEINDEXING();
    if (ASI < 0) {
      // not set on input, set default
      if (input.HKLIN_wild()) GC.set_AssumeSameIndexing(true);
    } else {
      GC.set_AssumeSameIndexing(ASI);
    }

    ASSERT (N_hklin == AllFiles.NumberOfHKLINnames());
    if (N_hklin <= 0) {
      Message::message(Message_fatal("No HKLIN filename given"));}
    // Filenames from command line
    AllFiles.StoreFilename("HKLREF", CL.getHKLREF());
    AllFiles.StoreFilename("XYZIN", CL.getXYZIN());
    AllFiles.StoreFilename("HKLOUT", CL.getHKLOUT());
    // ... but may be overwritten by command input
    // [which should take precedence? Here it is the command input]
    AllFiles.StoreFilename("HKLREF", input.getHKLREF());
    AllFiles.StoreFilename("XYZIN", input.getXYZIN());
    AllFiles.StoreFilename("HKLOUT", input.getHKLOUT());

    bool IsHKLREF = AllFiles.IsFile("HKLREF");
    bool IsXYZIN =  AllFiles.IsFile("XYZIN");
    IsHKLOUT = AllFiles.IsFile("HKLOUT");

    GC.set_Chiral(input.getCHIRAL());
    GC.set_MinIsig(input.getISIG());
    GC.set_OriginalLattice(input.getOriginalLatFlag());
    GC.set_LatticeTolerance(input.getTOLERANCE());
    GC.set_SysAbsCheck(input.getSYSABSCHECK());
    GC.set_TestFirstFile(input.TestFirstFile());
    GC.set_MaxMult(input.getScoreMaxMultiplicity());
    if (input.getNeighbourFraction() >= 0.0) {
      Zone::SetNeighbourFraction(input.getNeighbourFraction());
    }
    ChooseSolution SolutionChoice(input.getSolutionNumber(),
				  input.getChosenLaueGroup(),
				  input.getChosenSpaceGroup(), output);

    // Setting options from input
    // 0 CELL-BASED +1 SYMMETRY-BASED -1 C2
    int settingOption = input.getSetting();
    bool KeepOrthorhombicSetting = true;  // Default to ABC order of primitive orthorhombic
    int AllowI2 = +1;  // Use I2 setting if appropriate
    if (settingOption > 0) {
      KeepOrthorhombicSetting = false;  // SYMMETRY-BASED oP
    }
    if (settingOption != 0) {
      AllowI2 = 0; // don't allow I2 (SYMMETRY-BASED or C2)
    }
    if (SolutionChoice.IsI2()) {
      AllowI2 = -1; // force I2 setting
    }
    // If "choose spacegroup C2", or  "choose lauegroup C2" etc,
    // don't allow I2 solution
    if ((SolutionChoice.Spacegroup() == "C 1 2 1") ||
	(SolutionChoice.Lauegroup() == "C 1 2/m 1")) {
      AllowI2 = 0;
    }
    //^    std::cout << "Solutionchoice:"<<SolutionChoice.Spacegroup()
    //	      << "  SolutionchoiceLaue:"<<SolutionChoice.Lauegroup()
    //	      <<" settingOption "<<settingOption
    //^-      <<" AllowI2 " <<AllowI2<<"\n";

    GC.set_AllowI2(AllowI2);

    std::string outstring;
    std::string LaueGroup = input.getLAUEGROUP(); // blank unless LAUEGROUP command given
    LaueGroup = CheckInputRhombohedralSGname("In LAUEGROUP command: ",
					     LaueGroup, outstring);
    output.logTab(0,LOGFILE, outstring);

    GC.set_Lauegroup(input.getLAUEGROUP());
    std::string SpaceGroup = input.getSPACEGROUP();
    SpaceGroup = CheckInputRhombohedralSGname("In SPACEGROUP command: ",
					      SpaceGroup, outstring);
    output.logTab(0,LOGFILE, outstring);
    GC.set_Spacegroup(SpaceGroup);

    // Reindex operator to apply in output (ie not previously applied to test data)
    ReindexOp ReindexOut;
    // true if test data list has been reindexed
    bool reindexedTestData = false;

    // Only set GC.reindex if set
    if (input.REINDEXisSet())
      {
	GC.set_Reindex(input.getREINDEX());
	ReindexOut = GC.Reindex();
      }

    // Conditions for COPY option
    bool CopyFlag = false;
    // just write output file (ie copy) if
    //  1) SPACEGROUP or REINDEX given
    //  2) REINDEX given unless LAUEGROUP given too
    //  2) "-c[opy]" on command line
    //  3) COPY command given
    if (SpaceGroup != "" ||
	(GC.IsReindexSet() && GC.Lauegroup() == "") ||
	CL.CopyFlag() || input.getCOPY()) {
      CopyFlag = true;
      pointlessTask = COPY;
    }
    
    // Set Bitflag control, to reject any flagged observations
    // apart from overloads
    ObservationFlagControl ObsFlagControl; 
    ObsFlagControl.SetAcceptOverload();

    if (SpaceGroup != "" && IsHKLREF)
      {Message::message(Message_fatal("Cannot have both HKLREF filename given and SPACEGROUP"));}
    if (SpaceGroup != "" && !IsHKLOUT && !IsXYZIN)
      {Message::message(Message_fatal("If SPACEGROUP given you must define an HKLOUT file"));}
    if (GC.IsReindexSet() && IsHKLREF)
      {Message::message(Message_fatal("Cannot have both HKLREF filename given and REINDEX"));}
    if (IsXYZIN && IsHKLREF)
      {Message::message(Message_fatal("Cannot have both HKLREF and XYZIN filenames given"));}
    if (GC.IsReindexSet() && !IsHKLOUT)
      {Message::message(Message_fatal("If REINDEX given you must define an HKLOUT file"));}
    for (int i=0;i<AllFiles.NumberOfHKLINnames();i++) {
      if ((AllFiles.HKLINfilename(i) == AllFiles.Filename("HKLREF"))
	  && IsHKLOUT)
	{Message::message(Message_fatal
			  ("Cannot write output file if HKLIN & HKLREF files are the same"));}
    }
    double Pcutoff = 0.01;  // FIXME

    // *************************************************************
    //  Setup up controls for reflection & column selection etc
    //    these should be set from command input (or default)
    
    // Set Profile-fitted [default] or integrated intensity
    col_controls column_selection; 
    
    // File selection flags (resolution, datasets, batches etc)
    file_select file_sel(input, AllFiles.NumberOfFileSeries());
    // scale factor to apply on input
    file_sel.InputScale() = input.getMULTIPLY();
    // If resolution range given explicitly, reset minI/sigI to -1 to suppress
    // resolution cutoff in normalisation
    if (file_sel.reslimits().isSet()) {
      GC.set_MinIsig(-1);
    }

    // Scala control classes (default settings)
    //  run controls
    //  partials controls
    all_controls controls;
    
    controls.partials = partial_controls(input.getFracLimMin(),
					 input.getFracLimMax(),
					 input.getSclMinLim(),
					 input.getCheck(),
					 input.getMaxGap());

    if (input.IsPolarisationSet()) {
      // Note default constructor flags polarisation as unset, ie default
      controls.polarisationcontrol.SetFactor(input.Polarisation());
    }

    Scell input_cell = input.getCELL();     // mainly for SCALEPACK
    double wavelength = input.Wavelength();  // mainly for SAINT

    // Make list of required columns in unmerged HKLIN file
    MtzIO::column_labels column_list_unmrg = MtzIO::setup_columns();
    // label for merged HKLIN for I or F column (may be set by LABIN command)
    MtzIO::column_labels column_list_mrg;
    column_list_mrg.addLabin(input.getLABI_I(), input.getLABI_sigI());

    //<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
    // Read reference list, if required
    //
    hkl_merged_list RefList;  // reference list may be merged
    hkl_unmerge_list RefUnmrgList;     //   ... or unmerged
    hkl_symmetry RefSymm;
    Scell RefCell;        // Read file into list

    // Reference list from HKLREF or XYZIN
    if (IsXYZIN) {
      // Input coordinate file, calculate structure factors to make reference reflection list
      pointlessTask = REFERENCE;
      // should take this from HKLIN but more complicated!
      output.logTab(0,LOGFILE,
		    "\n---------------------------------------------------------------\n");
      output.logTab(0,LOGFILE,
		    "\nReading reference coordinate list from file "+
		    AllFiles.Filename("XYZIN")+"\n");
      double FCresolution = GetResolutionFromHKLIN(AllFiles);
      RefList.CreateFromAtoms(AllFiles.Filename("XYZIN"),SpaceGroup,input_cell,
			      FCresolution, true, output);
      PrintFileInfoToXML("XYZIN",AllFiles.Filename("XYZIN"),
			 RefList.Cell(), RefList.SpaceGroupSymbol(),
			 output);
      reflisttype = MERGED;
      RefSymm = RefList.symmetry();
      RefCell = RefList.Cell();
      SpaceGroup = RefSymm.symbol_xHM();  // output file should have same spacegroup as reference
    } else if (IsHKLREF) {
      // Create reference reflection list
      pointlessTask = REFERENCE;
      output.logTab(0,LOGFILE,
		    "\n---------------------------------------------------------------\n");
      output.logTab(0,LOGFILE,
		    "\nReading reference data set from file "+
		    AllFiles.Filename("HKLREF")+"\n");
      
      RefList.init(AllFiles.Filename("HKLREF"),false,output);
      // Resolution limits from input (if given)
      if (RefList.Unmerged()) {
	int verbose = +1;
	ReflectionFileType hklref_filetype;
	hklref_filetype.SetType("MTZ");  // for now anyway
	file_select ref_file_select; // dummy
	double Tol = 0.0; // dummy
	AddHKLIN(0,AllFiles.Filename("HKLREF"),
		 hklref_filetype, true, ref_file_select,
		 column_selection, column_list_unmrg,
		 controls, PxdName(), Scell(), Tol, wavelength,
		 output, verbose, RefUnmrgList);
	verbose = +1;
	PrintUnmergedHeaderStuff(RefUnmrgList, output, verbose);
	output.logTabPrintf(0,LOGFILE,
			    "Space group from HKLREF file : %s\n",
			    RefUnmrgList.symmetry().symbol_xHM().c_str());
	PrintFileInfoToXML("HKLREF",AllFiles.Filename("HKLREF"),
			   RefUnmrgList.Cell(),
			   RefUnmrgList.symmetry().symbol_xHM(),
			   output);
	reflisttype = UNMERGED;
	RefSymm = RefUnmrgList.symmetry();
	RefCell = RefUnmrgList.Cell();
      } else {
	MtzIO::column_labels ref_column_list;
	// label for merged HKLIN for I or F column (may be set by LABREF command)
	ref_column_list.addLabin(input.getLABR_I(), input.getLABR_sigI());
	RefList.read(ref_column_list, file_sel.reslimits().ResHigh(), true, output);
	RefList.PrintHeaderStuff(output);
	PrintFileInfoToXML("HKLREF",AllFiles.Filename("HKLREF"),
			   RefList.Cell(), RefList.SpaceGroupSymbol(),
			   output);
	reflisttype = MERGED;
	RefSymm = RefList.symmetry();
	RefCell = RefList.Cell();
      }
      SpaceGroup = RefSymm.symbol_xHM();  // output file should have same spacegroup as reference
      if (RefSymm.lattice_type() == 'R') {
	{Message::message(Message_fatal
	  ("At present, HKLREF file cannot be in a rhombohedral R lattice setting"));}
      }
    }
    //<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
    // Construct test reflection list by reading one or more MTZ (or XDS) files
    // Test reflection file may be merged or unmerged, but
    // a merged test file only allowed if reference dataset is given
    hkl_unmerge_list TestList;              // the test if unmerged
    scala::hkl_merged_list TestListMerged;  // the test list if merged

    // Merged HKLIN files can either
    //    1) be treated as merged, and used purely for comparison with reference,
    //     change of space group, or reindexing, ie not Laue group determination
    // or 2) be treated as unmerged, in which case a limited Laue group determination
    //     may be possible

    bool TestDataAsMerged = false;  // true if used NOT for Laue group determination
    RefListType HklinIsMerged = NONE;     // unknown
    RefListType mergeStatus = NONE;
    bool symSet;  // false if no symmetry info was available from the input files

    // check file type for all hklin files, reassign them if necessary
    AllFiles.CheckHklinFileType();

    // Check all input files to see if they are merged or not
    // Must be all merged or all unmerged, not mixed
    for (int ifl=0;ifl<AllFiles.NumberOfHKLINnames();++ifl) {
      hklin_filename = AllFiles.HKLINfilename(ifl, hklin_filetype);
      //    output.logTab(0,LOGFILE, "\nReading test dataset from HKLIN filename: " +
      //		  hklin_filename + "\n");
      bool fileIsUnmerged = true; // true if unmerged
      bool fileValid;
      if (hklin_filetype.IsTypeMTZ()) {
	MtzIO::MtzMrgFile mtzin;
	fileValid = mtzin.open_read(hklin_filename);
	fileIsUnmerged = !mtzin.Merged();  // invert to true if unmerged
      } else if (hklin_filetype.IsTypeSCA()) {
	SCAIO::SCAunmergeFile scain(hklin_filename, outputstring);
	output.logTab(0,LOGFILE,outputstring);
	if (scain.FileType() == SCAIO::SCA_MERGED) {
	  fileIsUnmerged = false;
	} else if (scain.FileType() == SCAIO::SCA_UNMERGED) {
	  fileIsUnmerged = true;
	} else if (scain.FileType() == SCAIO::SCA_SHELX) {
	  fileIsUnmerged = true;  // assume ShelX file is unmerged ...
	  // ... unless we are told otherwise
	  if (input.getHklinMergedGroup() != "") {
	    fileIsUnmerged = false;
	    if (input.getHklinMergedGroup() != "UNKNOWN") {
	      // Store merging space(point)group
	      AllFiles.StoreHKLINsymmetry(ifl, hkl_symmetry(input.getHklinMergedGroup()));
	    }
	  }
	}
      } else if (hklin_filetype.IsTypeXDS()) {
	int verbose = 0;
	XDSIO::XDSunmergeFile(hklin_filename, verbose, outputstring);
	output.logTab(0,LOGFILE,outputstring);
	fileIsUnmerged = true; //  XDS files must be unmerged (checked in XDSunmergeFile)
      }

      if ((mergeStatus == UNMERGED && !fileIsUnmerged) ||
	  (mergeStatus == MERGED && fileIsUnmerged)) {
	// Mixed file types
	Message::message
	  (Message_fatal("HKLIN files must be all merged or all unmerged"));
      }
      mergeStatus = fileIsUnmerged ? UNMERGED: MERGED;
      AllFiles.SetMergedFlag(ifl, fileIsUnmerged);  // flag file as merged or unmerged
    } // end loop input files
    ASSERT (mergeStatus != NONE);
    if (mergeStatus == UNMERGED) {
      // HKLIN file(s) are unmerged
      TestDataAsMerged = false;
      HklinIsMerged = UNMERGED;
    } else {
      if (AllFiles.NumberOfHKLINnames() > 1) {
	Message::message
	  (Message_fatal("You can only have one merged HKLIN file"));
      }
      HklinIsMerged = MERGED;
      // HKLIN file(s) are merged
      if (IsHKLREF || IsXYZIN || CopyFlag) {
	// Merged file will be used just for copying with or without references file
	// ie not Laue group determination
	TestDataAsMerged = true;
      } else {
	// Merged file will be treated as if it were unmerged
	TestDataAsMerged = false;
      }
    }

    if (TestDataAsMerged) {
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      // Test data is merged
      if (!(IsHKLREF || IsXYZIN) && !CopyFlag) {
	Message::message
	  (Message_fatal("HKLIN is merged and no HKLREF file is defined, SPACEGROUP or REINDEX"));}
      if (AllFiles.NumberOfHKLINnames() > 1) {
	Message::message(Message_fatal("Only one merged HKLIN file is allowed"));}

      // Merged HKLIN file, read data
      // Resolution limits from input (if given)
      // Print on
      TestListMerged.init(hklin_filename, true, output);
      TestListMerged.read(column_list_mrg, file_sel.reslimits().ResHigh(),true,output);
      PrintFileInfoToXML("HKLIN",hklin_filename,
			 TestListMerged.Cell(), TestListMerged.SpaceGroupSymbol(),
			 output);
      TestListMerged.PrintHeaderStuff(output);
      scala::hkl_symmetry TestSym = TestListMerged.symmetry();
      scala::Scell test_cell = TestListMerged.Cell();

      if (IsHKLREF || IsXYZIN) {
	// HKLREF defined: test merged data for alternative
	// indexing schemes against reference dataset
	
	// Fails in here if symmetry incompatible
	CheckCompatibleSymmetry(RefSymm, TestSym, TestDataAsMerged);
	
	std::vector<scala::ReindexOp> ReindexList =
	  AlternativeReindexList(RefCell, test_cell,
				 TestSym,GC.LatticeTolerance(), AllowI2);
	
	// Single alternative is only allowed if it is not the identity
	bool altIndex = true;
	if ((ReindexList.size() > 1) ||
	    ((ReindexList.size() == 1) && 
	     (!ReindexList[0].IsIdentity())))
	  {
	    // Yes we have some possible alternative indexing schemes
	    PrintAlternativeIndexing(ReindexList, GC.LatticeTolerance(), output);
	  }
	else
	  {
	    output.logTab(0,LOGFILE,"\nNo alternative indexing possible\n\n");
	    output.logTabPrintf(0,LXML,
				"<IndexScores>\n  <ScoreCount> 0</ScoreCount>\n</IndexScores>\n");
	    altIndex = false;
	  }
	
	std::vector<ReindexScore> Scores;
	
	if (altIndex) {
	  Scores = TestIndexMerged(RefList, TestListMerged, ReindexList,
				   GC, output);
	  output.logTab(0,LOGFILE,"\nAlternative indexing relative to reference file "+
			AllFiles.Filename("HKLREF"));
	  PrintIndexScores(Scores, true, output);  // Print to RESULT block
	  ReindexOut =  Scores[0];  // just get ReindexOp part of derived class
	} else {
	  ReindexOut = ReindexOp();
	}
	
	// Merged Test list has not been reindexed
	reindexedTestData = false;
	
	if (IsHKLOUT) {
	  CopyFlag = true;
	}   // set flag for output
	else {
	  // If no output file all done, stop here 
	  // Testing merged test data against reference, no output
	  return PointlessStatus(true, N_hklin, reflisttype, IsHKLOUT, pointlessTask);
	}
      } else {
	// No reference: if a space group is input (SPACEGROUP) but no REINDEX,
	// generate reindex operator ReindexOut, provided that symmetries belong
	// to same lattice group
	// Note that this is probably really only useful (or indeed valid) for
	// C2 <-> I2 & R3 <-> H3
	// If REINDEX command is given, ReindexOut is already set & is not changed here
	// Returns true if ReindexOut is set 
	// fails if the symmetries do not belong to same lattice group
	bool reindexSet = SpacegroupReindex(GC, TestSym, test_cell, ReindexOut, output);
	reindexSet = reindexSet; // just to remove warning of unused variable
      }
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    else {
      // Unmerged HKLIN file(s), read them into list, putting them on the
      // same indexing system etc if necessary
      //  Also comes here for merged files treated as unmerged (with a different column_list)
      int istat = ReadAllHKLINfiles(AllFiles, file_sel, column_selection,
	    (HklinIsMerged==UNMERGED) ? column_list_unmrg : column_list_mrg,
				    controls, GC, input_cell, wavelength,
				    RefList, RefUnmrgList, reflisttype,
				    HklinIsMerged,
				    output, TestList);
      //^      std::cout << TestList.symmetry().symbol_xHM() << "\n"; //^
      origLatType = TestList.symmetry().lattice_type();
      if (istat == -1) {
	symSet = false; // no symmetry in input files
	istat = 0;
      } else {
	symSet = true;
      }
      // store resolution range from file, modified by input RESOLUTION limits,
      // for later use
      ResRange0 = TestList.ResRange();

      reindexedTestData = true;  // Test data is already reindexed, if necessary
      if (istat != 0) {
	if (istat == 2) {
	  output.logTab(0,LOGFILE,
			"\n***** Stopping because cell discrepancy between files is too large ****");
	  output.logTab(0,LOGFILE,
			"*****          Increase TOLERANCE to accept all cells              ****");
	}
	// Error exit
	return PointlessStatus(false, N_hklin, reflisttype, IsHKLOUT, pointlessTask);
      }

      output.logTab(0,LOGFILE,
		    "\n===============================================================");
      //      output.logTab(0,LOGFILE,
      //		    "\n<!--SUMMARY_BEGIN-->\n");
      output.logTab(0,LOGFILE,"\n>*> Summary of test data read in:");
      int verbose = +3;
      PrintUnmergedHeaderStuff(TestList, output,verbose);
      // reset Obs flag controls & print counts
      TestList.ResetObsAccept(ObsFlagControl);
      output.logTab(0,LOGFILE,ObsFlagControl.PrintCounts());
      output.logTab(0,LOGFILE,
		    "\n===============================================================");
      // Print reindexing information
      if (GC.IsReindexSet()) {
	output.logTab(0,LOGFILE,
		      "\nData reindexed by input operator "+GC.Reindex().as_hkl());
      } else if (reflisttype != NONE || N_hklin > 1) {
	if (GC.AssumeSameIndexing()) {
	  if (AllFiles.AreThereFileSeries()) {
	    output.logTab(0,LOGFILE,
			  "\nAssuming all file series have the same indexing");
	  } else if (!AllFiles.IsHKLINreindex()) {
	    output.logTab(0,LOGFILE,"\nNo alternative indexing");
	  } else {
	    output.logTab(0,LOGFILE,
			  "\nAssuming all files have the same indexing");
	  }
	} else {
	  if (AllFiles.IsHKLINreindex()) {
	    PrintReindexSummary(reflisttype, AllFiles, output);
	  } else {
	    output.logTab(0,LOGFILE,"\nNo alternative indexing to test");
	  }
	}
      }
      // 
      if (reflisttype != NONE) {
	if (IsHKLOUT) {
	  CopyFlag = true;
	}   // set flag for output
	else {
	  // if reference and no output file all done, stop here 
	  return PointlessStatus(true, N_hklin, reflisttype, IsHKLOUT, pointlessTask);
	}
      }
      // Set Chiral flag if necessary
      if (GC.GetChiral() == UNKNOWN)
	{SetChiralityFromHKLIN(TestList, GC);}
    }
    // At this point we have a test data set. It is either:
    // 1. unmerged, in TestList, or 
    // 2. merged, in TestListMerged (if TestDataAsMerged == true)
    //<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
    // Option to output summed list
    if (input.getSUMMEDLIST() != "") {
      WriteSummedList(TestList, input.getSUMMEDLIST());
    }

    if (CopyFlag) {
      // just write output file (ie copy) if
      //  1) SPACEGROUP or REINDEX given
      //  2) REINDEX given unless LAUEGROUP given too
      //  2) "-c[opy]" on command line
      //  3) COPY command given
      // Merged HKLIN
      // Always reduce to asymmetric unit
      if (TestDataAsMerged) {
	bool reduce = true;
	bool verbose = true;
	MtzIO::CopyMergedMTZ(hklin_filename,
			     AllFiles.Filename("HKLOUT"),
			     SpaceGroup, ReindexOut, reduce,
			     verbose, output);
      }
      else {    // Unmerged
	hkl_symmetry symm;
	if (SpaceGroup == "HKLIN" || SpaceGroup == "") {
	  symm = TestList.symmetry();  // Symmetry from HKLIN (merged or unmerged)
	} else {
	  symm = hkl_symmetry(SpaceGroup);
	}
	OutputUnmergedData(TestList, AllFiles.Filename("HKLOUT"),
			   symm, ReindexOut, reindexedTestData,
			   TestList.Title(), "space", output);
      }
      // Just copying
      return PointlessStatus(true, N_hklin, reflisttype, IsHKLOUT, pointlessTask);
    }
    // *************************************************************
    if (input.Centre()) {
      // CENTRE search for centre of lattice
      FindCentre(TestList, GC.Lauegroup(), GC.Reindex(), GC.GetMinIsig(),
		 input.CentreGridMax(), AllFiles.Filename("HKLOUT"), output);
      pointlessTask = CENTRE;
      return PointlessStatus(true, N_hklin, reflisttype, IsHKLOUT, pointlessTask);
    }
    // *************************************************************
    // *************************************************************
    // **** MODE POINTGROUP ****
    // Score all possible Laue groups
    std::vector<PGscore> SGscores;
    int ChosenSolution = 0;

    ScoreAccept Accept;
    // Save hklin space group (this will be the symmetry from the 1st input file)
    std::string hklinSpacegroup = TestList.symmetry().symbol_xHM();
    ValidElementObserved  ElementsObserved;
    
    output.logTabPrintf(0,LOGFILE,"\nTime for reading file(s): %8.3f secs\n", cputime.Stop());    
    // The calculation of likelihood from CC depends among other things on the weighting of
    // the P(CC;!S) term, based on modelling the expectation of CC;!S (ie with symmetry absent)
    // This is done (at present) with a model E(CC;!S) = m; P(m;!S) = (1 - m^k)^(1/k)
    // The value of k defaults to PMEXPONENTMIN, but is increased depending on the estimated twin factor,
    // up to a maximum of PMEXPONENTMAX
    const double PMEXPONENTMAX = 5.0;
    const double PMEXPONENTMIN = 2.0;
    double PMexponent = PMEXPONENTMIN; // default value
    // This is stored as a static variable in classes SCsignificance and PGscore
    SCsignificance::SetExponent(PMexponent);
    PGscore::SetExponent(PMexponent);

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    Twintest TwinTest;
    if (GC.GetChiral() == CHIRAL) {
      // twinning test only for chiral spacegroups
      // (I don't have the distributions for centric data)
      output.logTab(0,LOGFILE,
		    "\nChecking for possible twinning");
      cputime.Start();
      double MinIsigIstatsTwin = 8.0;  // minimum bin <I>/<sigma>
      TestList.ResetResoLimits();  // reset to file limits
      TwinTest.init(TestList, MinIsigIstatsTwin, output);
      TestList.ResetResoLimits();  // reset to file limits
      std::vector<double> nl = TwinTest.NL();
      TwinTest.PrintTwinResult(output);
      // Work out a suitable value for PMexponent
      double alpha = TwinTest.Alpha();  // estimated twin fraction
      double f = alpha*(4-alpha*(6.0-4.0*alpha));
      PMexponent = PMEXPONENTMIN + f*(PMEXPONENTMAX-PMEXPONENTMIN);
      SCsignificance::SetExponent(PMexponent);
      PGscore::SetExponent(PMexponent);
      output.logTabPrintf(0,LOGFILE,
 "Model for expectation(CC) = E(m) if symmetry is absent P(m;!S) = (1-m^k)^(1/k) with k = %4.1f\n\n",
			  PMexponent);

      output.logTabPrintf(0,LOGFILE,
			  "Time for twinning test %8.3f secs\n\n",cputime.Stop());
      output.logTab(0,LOGFILE,
       "======================================================================\n\n");
    }
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    if (GC.Lauegroup() != "" && GC.Lauegroup() != "ALL") {
      // LAUEGROUP given on input, set with dummy scores
      std::string LGname = GC.Lauegroup();
      // option to keep original Laue group
      if (GC.Lauegroup() == "HKLIN") {
	LGname = TestList.symmetry().symbol_xHM();
      }
      CCtbxSym::PointGroup LG(LGname);
      LG.SetCell(TestList.cell().UnitCell(), GC.Reindex(), GC.AllowI2());
      SGscores.push_back(PGscore(LG));
      // dummy probability score
      SGscores.back().SetLikelihood(1.0); 
      SGscores.back().SetAccept(true); // mark as accepted

      output.logTab(0,LOGFILE,
		    "Reindexing or changing symmetry to selected Laue group "+
		    LG.RefLGname());
      output.logTab(0,LOGFILE, "Reindex operator from input cell to given Lauegroup: "
		    + GC.Reindex().as_hkl() + "\n\n"
		    + GC.Reindex().as_matrix() + "\n");
      
      // reindex, sort & reorganise hkl list
      int Junk = TestList.change_symmetry(hkl_symmetry(LG.RefLGname()), GC.Reindex());
      Junk = Junk;
      // Get average cell from all datasets
      Scell OverallCell = TestList.cell(PxdName());
      output.logTabPrintf(0,LOGFILE,"\nUnit cell after reindexing:");
      for (int i=0;i<6;i++) output.logTabPrintf(0,LOGFILE,"%7.2f",OverallCell.UnitCell()[i]);
      output.logTabPrintf(0,LOGFILE,"\n\n");
      Normalise NormRes =  scala::NormaliseUnmerge(TestList, GC, controls,
						   true, output);
    } else { // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      // Check merged files for gross incompleteness
      double overallcompleteness = 0.0;
      if (HklinIsMerged == MERGED) {
	std::vector<double> complete =
	  scala::Completeness(TestList, 1);  // all in one bin
	overallcompleteness = complete[0];
      }

      // Test all possible Lauegroups
      ElementsObserved =
	TestLaueGroup(TestList, GC, controls, HklinIsMerged, output, SGscores);
      output.logTabPrintf(0,LOGFILE,"\nTime to determine pointgroup: %8.3f secs\n", cputime.Stop());

      if (GC.GetMinIsig() > 0.0) {
	// the resolution range may have been cut back in the normalisation, but not if
	// explicitly set on input (flagged by MinIsig < 0)
	ResRangeCut = TestList.ResLimRange();
      }
      if (SolutionChoice.ChooseLaueGroup()) {
	// CHOOSE solution or Lauegroup option given
	ChosenSolution = SolutionChoice.MarkAcceptedSubgroup(SGscores);
      } else {
	Accept = ScoreAccept(input.getACCEPTANCELIMIT(),
			     SGscores[0].Likelihood());
	// Test for acceptance for all candidates
	if (Accept.IfSet()) MarkAcceptedSubgroup(SGscores, Accept);

	if (!ElementsObserved.Valid() && !(HklinIsMerged == MERGED && symSet))  {
	  // Insufficient symmetry elements observed to determine Laue group:
	  // If merged and input symmetry given, then use the results anyway
	  //   - accept original Laue group from input file unless	  
	  //     there was no symmetry in the input file (eg SHELX),
	  //     in which case use the best group with no reindexing
	  if (!ResetInputLaueGroup(SGscores, hklinSpacegroup, symSet, output))  {
	    // Can't find it in list (because of reindexing)
	    // Just reset it to original
	    std::string LGname = TestList.symmetry().symbol_xHM();
	    CCtbxSym::PointGroup LG(LGname);
	    LG.SetCell(TestList.cell().UnitCell(), ReindexOp(), false);
	    SGscores.clear();
	    SGscores.push_back(PGscore(LG));
	    // dummy probability score
	    SGscores.back().SetLikelihood(1.0); 
	    SGscores.back().SetAccept(true); // mark as accepted
	  }
	} else {
	  if (Accept.IfSet()) {
	    if (symSet) {
	      output.logTabPrintf(0,LOGFILE,
				  "\nAcceptable Laue groups have scores above %5.2f\n",
				  Accept.Threshold());
	    }
	  } else if (ChosenSolution > 0) {
	    output.logTabPrintf(0,LOGFILE,
				"Selected solution number %3d chosen, Laue group %s\n",
				ChosenSolution,
				SGscores[ChosenSolution-1].RefLGname().c_str());
	  } else {
	    output.logTab(0,LOGFILE, "\n\nAll Laue groups accepted\n");
	  }
	}
      }

      PrintSubGroupScores(SGscores, HklinIsMerged, overallcompleteness, output);
      if (SolutionChoice.ChooseLaueGroup()) {
	if (ChosenSolution < 0) {
	  Message::message(Message_fatal
			   ("Cannot find chosen Laue group solution "+
			    SolutionChoice.Lauegroup()));
	} else if (ChosenSolution == 0) {
	  Message::message(Message_fatal
			   ("Cannot find chosen Laue group solution number "+
			    clipper::String(SolutionChoice.SolutionNumber(),3)));
	}
      }
    }
    //............................................................
    // If not chiral, score intensity statistics for possible
    // centrosymmetric spacegroup, even if defined as centrosymmetric
    double ProbAcentric = 1.0;
    if (GC.GetChiral() != CHIRAL)
      {
	output.logTabPrintf(0,LOGFILE,
			    "\n********************************************************\n");
	double MinIsigIstats = 6.0;  // minimum bin <I>/<sigma>
	TestList.ResetResoLimits();  // reset to file limits
	IntensityStatistics Istats(TestList, MinIsigIstats, output);
	Istats.Print(output);
	ProbAcentric = Istats.ProbAcen();
      }
    //............................................................
    // We now have a ranked list in SGscores of candidate
    // Laue groups
    //  Create list of possible space groups
    // Option to skip systematic absence checks
    Solution BestSolution;
    if ( !GC.SysAbsCheck()) {
      int k = 0;  // choose top score
      // Point group name
      std::string Sname =
	CCtbxSym::PointGroupName(SGscores[k].RefPGname());
      BestSolution = Solution(k, Sname,
			      false, SGscores[k].RefSGreindex(),
			      SGscores[k].Likelihood(),1.0);  
      Scell cellOrig =  TestList.cell().change_basis(TestList.TotalReindex().inverse());
      PrintSummary(BestSolution, "point", cellOrig, ResRangeCut, ResRange0, TestList.num_batches(),
		   TwinTest, output);
    } else {
      // Check if we have sufficient potentially symmetry-related
      // observations to choose the Laue group. If not, just use
      // the input group
      output.logTabPrintf(0,LOGFILE,
			  "\n********************************************************");
      output.logTabPrintf(0,LOGFILE,
			  "\n\nTesting Lauegroups for systematic absences\n");
      output.logTabPrintf(0,LOGFILE,
			  "------------------------------------------\n\n");
      
      std::vector<scala::PossibleSpaceGroup> AllGroups;
      //    and of zone lists for each Lauegroup
      // Unique zones
      std::vector<Zone>  AllUniqueZones = MakeUniqueZones(SGscores, GC.GetChiral());
      
      // Test for systematic absences in lattice group, return list of space groups
      // & zones
      // Get statistics
      if ( AllUniqueZones.size() > 0) {
	// Set list to accept profile-fitted overloads
	ObsFlagControl.SetAcceptOverload();
	TestList.ResetObsAccept(ObsFlagControl);
	ObsFlagControl.Clear();
	// Reset resolution limits to accept all data for systematic absence test
	TestList.ResetResoLimits();  // reset to file limits
	SysAbsStats(TestList, AllUniqueZones);
	int Noverload = ObsFlagControl.NumAccOverload();
	if (Noverload > 0) {
	  output.logTabPrintf(0,LOGFILE,
	      "\n\nNumber of accepted overloads in systematic absence analysis = %8d\n\n",
			      Noverload);
	}
      }
	      
      // test against space group possibilities
      for (size_t k=0;k<SGscores.size();k++) {
	if (SGscores[k].Accepted()) {
	  std::vector<scala::PossibleSpaceGroup> groups =
	    SysAbsTest(SGscores[k], GC.GetChiral(), AllUniqueZones);
	  //   and accumulate
	  AllGroups.insert(AllGroups.end(), groups.begin(), groups.end());
	}
      }

      if (GC.GetChiral() != CHIRAL) {
	for (size_t i=0;i<AllGroups.size();i++) {
	  AllGroups[i].StoreAcenProb(ProbAcentric);
	}
      }
      
      // Sort on total score
      std::sort(AllGroups.begin(),AllGroups.end());
      
      // Axis zone data for loggraph
      if (AllUniqueZones.size() > 0) {
	output.logTabPrintf(0,LOGFILE,
			    "I' is intensity adjusted by subtraction of a small fraction (%4.2f, NEIGHBOUR)\n of the neighbouring intensities, to allow for possible overlap\n",
			    Zone::GetNeighbourFraction());
      }
      for (size_t iz=0;iz<AllUniqueZones.size();iz++)  {
	// Skip zones we've done already
	bool done = false;
	if (iz > 0) {
	  for (size_t iz2=0;iz2<iz;iz2++)  {
	    if (AllUniqueZones[iz].Compare(AllUniqueZones[iz2])) {done = true;}
	  }
	}
	if (!done) {
	  OutputZoneTable(AllUniqueZones[iz], output);
	}
      }
      PrintZoneScores(AllUniqueZones, output);
      
      // restore rhombohedral lattice if in original file or explicitly required
      if (origLatType == 'R') {
	if (!SolutionChoice.ChooseSpaceGroup() ||
	    SolutionChoice.IsRhombohedralR()) {
	  // check all subgroups, reindex H lattices to R
	  ResetRhombohedralHexLatticestoR(AllGroups);
	}}

      // Flag accepted groups
      int naccept = CheckPossibleSpacegroups(AllGroups, Pcutoff);
      naccept = naccept;

      output.logTabPrintf(0,LOGFILE,"\nTime for systematic absence tests: %8.3f secs\n",
			  cputime.Stop());
      
      PrintPossibleSpaceGroups
	(AllGroups, Pcutoff, GC.GetChiral(), settingOption, output);
      
      Scell cellOrig =  TestList.cell().change_basis(TestList.TotalReindex().inverse());
      
      BestSolution = ChooseBestGroup(SolutionChoice, AllGroups, settingOption,
				     cellOrig, origLatType, ResRangeCut, ResRange0,
				     TestList.num_batches(), TwinTest, hklinSpacegroup, output);
    }
    
    output.logTab(0,LOGFILE,"\n\n");
    
    // Check if all input files had same spacegroup etc
    bool SameCrystalType = true;
    bool SameSymm = true;
    // CrystalType is crystal system + lattice centre
    std::vector<CrystalType> CrysTypes;
    CrystalType Ctyp0;
    if (!AllFiles.HKLINsymmetry(0).IsNull()) {
      Ctyp0 = CrystalType(AllFiles.HKLINsymmetry(0).CrysSys(),
			  AllFiles.HKLINsymmetry(0).lattice_type());
      CrysTypes.push_back(Ctyp0);
    }

    if (N_hklin > 1) {
      for (int i=1;i<AllFiles.NumberOfHKLINnames();i++) {  // loop from 2nd onwards
	if (!AllFiles.HKLINsymmetry(i).IsNull()) {
	  // test for equality of rotational operators, ignoring translations
	  // ie point group + any lattice centring
	  if (AllFiles.HKLINsymmetry(i) != AllFiles.HKLINsymmetry(0)) {
	    SameSymm = false;
	  }
	  // ... and test just for crystal type
	  CrystalType Ctyp(AllFiles.HKLINsymmetry(i).CrysSys(),
			   AllFiles.HKLINsymmetry(i).lattice_type());
	  // For test of equality, trigonal == hexagonal since they have the 
	  // same cell constraints
	  if (Ctyp != Ctyp0) {
	    SameCrystalType = false;
	  }
	  CrysTypes.push_back(Ctyp);  // store all of them even if the same
	}
      }
    }
    if (N_hklin == 1 || SameSymm) {
      if (CrysTypes.size() > 0) {
	output.logTab(0,LOGFILE,
	      "HKLIN spacegroup: "+hklinSpacegroup+"  "+CrysTypes[0].format(false)+"\n");
      } else {
	output.logTab(0,LOGFILE,"No HKLIN spacegroup specified");
      }
    } else {
      output.logTab(0,LOGFILE,"HKLIN spacegroups:");
      for (int i=0;i<AllFiles.NumberOfHKLINnames();i++) {
	if (!AllFiles.HKLINsymmetry(i).IsNull()) {
	  output.logTabPrintf(1,LOGFILE,"File%4d: %15s  %s\n",i+1,
			      AllFiles.HKLINsymmetry(i).symbol_xHM().c_str(),
			      CrysTypes[i].format(false).c_str());
	} else {
	  output.logTabPrintf(1,LOGFILE,"File%4d: unspecified\n",i+1);
	}
      }
    }

    CrystalType OutputCrystalType(hkl_symmetry(BestSolution.Name()));
    // Various warnings if the lattice type has change from the input 
    if (!SameCrystalType) {
      output.logTab(0,LOGFILE,
       "\n$TEXT:Warning:$$ $$\nWARNING:");
      output.logTab(0,LOGFILE,
		    "Different input files were integrated in different crystal systems.");
      output.logTab(0,LOGFILE,
       "You should rerun the integration in the chosen crystal system because the cell constraints differ.");
      output.logTabPrintf(0,LOGFILE,"Input crystal systems are:");
      CrystalType Ctyp = Ctyp0;
      for (int i=0;i<N_hklin;i++) {
	if (i==0 || CrysTypes[i] != Ctyp) {
	  output.logTabPrintf(0,LOGFILE," %s", CrysTypes[i].format(false).c_str());
	  if (i != N_hklin-1) {output.logTabPrintf(0,LOGFILE,",");}
	}
	Ctyp =  CrysTypes[i];
      }
      output.logTabPrintf(0,LOGFILE,".\n");
      output.logTab(0,LOGFILE,
		    "The crystal system chosen for output is "+OutputCrystalType.format(false)+".");
      output.logTab(0,LOGFILE,"$$\n");
    } else if (Ctyp0.isSet() && OutputCrystalType != Ctyp0) {
      if (EquivalentCentredMonoclinic(OutputCrystalType, Ctyp0)) {
	output.logTab(0,LOGFILE,
		      "\n$TEXT:Warning:$$ $$\nWARNING:");
	output.logTab(0,LOGFILE,
		      "Input and output crystal systems are different but equivalent centred monoclinic systems.");
      } else if (EquivalentRhombohedral(OutputCrystalType, Ctyp0)) {
	output.logTab(0,LOGFILE,
		      "\n$TEXT:Warning:$$ $$\nWARNING:");
	output.logTab(0,LOGFILE,
		      "Input and output crystal systems are different but equivalent rhombohedral systems.");
      }	else {
	output.logTab(0,LOGFILE,
		      "\n$TEXT:Warning:$$ $$");
	if (HklinIsMerged == MERGED) {
	  output.logTab(0,LOGFILE, std::string("\nWARNING: undermerged data\n")+
	    "The chosen output crystal system is different from that given in the input merged file.\n"
	   +"You should remerge the data in the correct chosen crystal system.");
	} else {
	  output.logTab(0,LOGFILE, std::string("\nWARNING:\n")+
	    "The chosen output crystal system is different from that used for integration of the input file(s).\n"
	    +"You should rerun the integration in the chosen crystal system because the cell constraints differ.");
	}
      }
      output.logTab(0,LOGFILE,
		    "\nThe input crystal system is "+Ctyp0.format(false));
      output.logTab(0,LOGFILE," (Cell: "+AllFiles.HKLINcell(0).format()+")");
      output.logTab(0,LOGFILE,
		    "The crystal system chosen for output is "+OutputCrystalType.format(false));      
      Scell cellOrig =  TestList.cell().change_basis(TestList.TotalReindex().inverse());
      output.logTab(0,LOGFILE," (Cell: "+cellOrig.change_basis(BestSolution.Reindex()).format()+")");
      output.logTab(0,LOGFILE,"$$\n");
    }

    if (N_hklin == 1) {
      std::string fname = TestList.Filename();
      std::cout << fname << "\n";
      output.logTab(0,LOGFILE,"Filename: "+TestList.Filename()+"\n");
    }
    output.logTab(0,LOGFILE,
		  "\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n\n");

    // Optionally write out data in best space group
    if (IsHKLOUT)
      {
	WriteOutputFile(AllFiles.Filename("HKLOUT"), BestSolution,
			TestList, output);
      }
  } // end try

  catch (phaser_io::PreprocessorError& capErr) {
    output.logWarning(LOGFILE, capErr.partialEcho()+"PREPROCESSOR ERROR: " + capErr.Message()
                      + "\n");
    success = false;
  }
  catch (phaser_io::SyntaxError& ccp4Err) {
    output.logWarning(LOGFILE, ccp4Err.Echo()+"\nSYNTAX ERROR: " + ccp4Err.Message()
                      + "\n");
    success = false;
  }
  catch (phaser_io::InputError& inpErr) {
    output.logWarning(LOGFILE, inpErr.Echo()+"\nINPUT ERROR: " + inpErr.Message()
                      + "\n");
    success = false;
  }

  catch (Message_fatal& message) {
    output.logWarning(LOGFILE, "\nFATAL ERROR message: \n"
		      + message.text() + "\n");
    success = false;
  }

  catch (std::exception const& err) {
    output.logWarning(LOGFILE, "\nUNHANDLED EXCEPTION: " + std::string(err.what())+"\n");
    success = false;
  }
  catch (...) {
    output.logWarning(LOGFILE, "\nUNKNOWN EXCEPTION TYPE\n");
    success = false;
  }

  return PointlessStatus(success, N_hklin, reflisttype, IsHKLOUT, pointlessTask);
}

//--------------------------------------------------------------
int main_two(int argc, char* argv[])
// Close down files (XML) & return ExitStatus to shell
{
  // Initialise CCP4 command line parser
  CCP4::ccp4fyp(argc, argv);

  CCP4::ccp4ProgramName (PROGRAM_NAME.c_str());
  std::string rcsdate = "$Date: "+std::string(PROGRAM_DATE2)+"$";
  CCP4::ccp4RCSDate     (rcsdate.c_str());
  CCP4::ccp4_prog_vers(PROGRAM_VERSION);
  CCP4::ccp4_banner();

  // Initialise output object, CCP4 mode, write header
  phaser_io::Output output;
  output.setPackageCCP4();
  output.SetMaxLineWidth(100);

  // Call main program
  PointlessStatus status = Pointless(argc, argv, output);
  // Close XML file if necessary
  if (output.doXmlout()) output.logTab(0,LXML,"</POINTLESS>");
  output.logTab(0,LOGFILE,
		"$TEXT:Reference: $$ Please reference $$");
  output.logTab(0,LOGFILE,
	std::string("P.R.Evans, 'Scaling and assessment  of data quality'")+
		" Acta Cryst. D62, 72-82  (2006).");
  output.logTab(0,LOGFILE,
		std::string("P.R.Evans, 'An introduction to data reduction: space-group determination, scaling and intensity statistics'")+
			    " Acta Cryst. D67, 282-292 (2011)");

  output.logTab(0,LOGFILE,"$$");
  int ExitStatus = 0;
  if (!status.Success()) {ExitStatus = +1;}
  return ExitStatus;
}
