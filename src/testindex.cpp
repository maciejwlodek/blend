//  testindex.cpp

#include "pointless.hh"
#include "testindex.hh"
#include "normalise.hh"
#include "normaliseunmerge.hh"
#include "average.hh"
#include "getCCsigfac.hh"

namespace scala {
  //--------------------------------------------------------------
  std::vector<ReindexScore> TestIndexUnmerged
  (hkl_merged_list& RefList,
   hkl_unmerge_list& hkl_list,
   const std::vector<ReindexOp>& ReindexList,
   const GlobalControls& GC,
   const all_controls& controls, 
   phaser_io::Output& output)

  // Test unmerged data for alternative
  // indexing schemes against reference dataset (HKLREF)
  // MODE INDEX
  {
    std::vector<ReindexScore>
      ReindexScoreList(ReindexList.begin(),ReindexList.end());
    int Nalt = ReindexScoreList.size();
    // *    std::vector<IsigI> IsigRef(Nalt);
    IsigI IsigRef;

//    int PrintLevel = 0;
    ResoRange RRange = RefList.ResRange();
    output.logTab(0,LOGFILE, "\n>>>> Normalising reference dataset");
    Normalise NormResRef = SetNormaliseMerged(RefList, -1.0, RRange);
    output.logTabPrintf(0,LOGFILE,
			"\nLog(<I>) fit for intensity normalisation: B (slope) %10.2f\n",
			NormResRef.Bcorr());
  
    output.logTab(0,LOGFILE, "\n>>>> Normalising test dataset");
    Normalise NormResTest =  scala::NormaliseUnmerge(hkl_list, GC, controls, true, output);
    output.logTabPrintf(0,LOGFILE,
			"\nLog(<I>) fit for intensity normalisation: B (slope) %10.2f\n",
			NormResTest.Bcorr());

    // Set mininum resolution range for test list
    ResoRange minrange = RefList.ResRange().MinRange(hkl_list.ResLimRange());
    hkl_list.SetResoLimits(minrange.ResLow(), minrange.ResHigh());

    reflection this_refl;
    observation this_obs;
    IsigI Isig;
    int isym;
    hkl_symmetry RefSymm = RefList.symmetry();
    hkl_list.rewind();

    // loop all reflections from test list
    while (hkl_list.next_reflection(this_refl) >= 0) {
      /*  --- Using average test data: not appropriate if symmetry is wrong!
	  Hkl hkl = this_refl.hkl();
	  // fetch equivalent IsigI for possible hkl's from reference list
	  for (int i=0;i<Nalt;i++) {
	  IsigRef[i] = RefList.Isig(RefSymm.put_in_asu
	  (hkl.change_basis(ReindexScoreList[i]),isym));
	  }
	  // Average unmerged I's
	  average_I(this_refl, Isig.I(), Isig.sigI());
	  float sSqr = this_refl.invresolsq();      
	  for (int i=0;i<Nalt;i++)
	  // add pairs of test, reference into scores
	  {
	  if (IsigRef[i].sigI() > 0.001) {
	  ReindexScoreList[i].addDsetIstats(Isig, IsigRef[i],
	  NormResTest, NormResRef,
	  sSqr);
	  }
	  }*/
      ///*---- Unaveraged alternative
      float sSqr = this_refl.invresolsq();      
      int Nobs = this_refl.num_observations();
	for (int lobs = 0; lobs < Nobs; lobs++) {
	  this_obs = this_refl.get_observation(lobs);
	  if (this_obs.IsAccepted()) {
	    Hkl hkl = this_obs.hkl_original();
	    IsigI Isig = this_obs.I_sigI();
	    // fetch equivalent IsigI for possible hkl's from reference list
	    for (int i=0;i<Nalt;i++) {
	      IsigRef = RefList.Isig(RefSymm.put_in_asu
				     (hkl.change_basis(ReindexScoreList[i]),isym));
	      if (IsigRef.sigI() > 0.001) {
		// add pairs of test, reference into scores
		ReindexScoreList[i].addDsetIstats(Isig, IsigRef,
						  NormResTest, NormResRef,
						  sSqr);
	      }
	    }
	  }
	}
	// ----*/
    }
    // update ReindexScoreList to add likelihood
    GetTestUnmergedSignificance(hkl_list, NormResTest, ReindexScoreList, output);
    // Restore resolution limits which were reset by normalisation
    hkl_list.ResetResoLimits();
    std::sort(ReindexScoreList.begin(), ReindexScoreList.end());
    return ReindexScoreList;
  }
  //--------------------------------------------------------------
  std::vector<ReindexScore> TestIndexMerged
  (hkl_merged_list& RefList,
   const hkl_merged_list& TestList,
   const std::vector<ReindexOp>& ReindexList,
   const GlobalControls& GC,
   phaser_io::Output& output)
  // Test merged data for alternative
  // indexing schemes against reference dataset (HKLREF)
  // MODE INDEX
  {
    std::vector<ReindexScore>
      ReindexScoreList(ReindexList.begin(),ReindexList.end());
    int Nalt = ReindexScoreList.size();
    std::vector<IsigI> IsigRef(Nalt);

    //    int PrintLevel = 0;
    ResoRange RRange = RefList.ResRange();
    output.logTab(0,LOGFILE, "\n>>>> Normalising reference dataset");
    Normalise NormResRef = SetNormaliseMerged(RefList, -1.0, RRange);
    output.logTabPrintf(0,LOGFILE,
			"\nLog(<I>) fit for intensity normalisation: B (slope) %10.2f\n",
			NormResRef.Bcorr());
  
    output.logTab(0,LOGFILE, "\n>>>> Normalising test dataset");
    Normalise NormResTest =  SetNormaliseMerged(TestList, GC.GetMinIsig(), RRange);
    output.logTabPrintf(0,LOGFILE,
			"\nLog(<I>) fit for intensity normalisation: B (slope) %10.2f\n",
			NormResTest.Bcorr());

    // loop all reflections in TestList
    TestList.start();
    IsigI Isig, Is;
    int isym;
    hkl_symmetry RefSymm = RefList.symmetry();
    bool OK;

    while (TestList.next(Isig)) {
      if ( ! Isig.missing()) {
	Hkl hkl = TestList.hkl();
	float sSqr = TestList.invresolsq();      
	IsigI Is1 = NormResTest.applyAvg(Isig, sSqr);
	if (NormResTest.NotTooLarge(Is1)) {
	  // fetch equivalent IsigI for possible hkl's from reference list
	  OK = true;
	  for (int i=0;i<Nalt;i++) {
	    Is = RefList.Isig(RefSymm.put_in_asu
			      (hkl.change_basis(ReindexScoreList[i]),isym));
	    if (Is.sigI() < 0.001) {
	      OK = false;
	    } else {
	      IsigRef[i] = NormResRef.applyAvg(Is, sSqr);
	      if (!NormResRef.NotTooLarge(IsigRef[i])) {
		OK = false;
	      }
	    }
	  }
	  if (OK) {
	    for (int i=0;i<Nalt;i++) {
	      ReindexScoreList[i].addDsetIstats(Is1, IsigRef[i]);
	    }
	  }
	}
      }
    }
    std::sort(ReindexScoreList.begin(), ReindexScoreList.end());
    return ReindexScoreList;
  }
  //--------------------------------------------------------------
  std::vector<ReindexScore> TestIndexBothUnmerged
  (hkl_unmerge_list& ref_list,
   hkl_unmerge_list& test_list,
   const std::vector<ReindexOp>& ReindexList,
   const GlobalControls& GC,
   const all_controls& controls, 
   const int& PrintLevel,    
   phaser_io::Output& output)
  // Test unmerged data  against unmerged reference dataset
  // for alternative indexing schemes
  {
    std::vector<ReindexScore>
      ReindexScoreList(ReindexList.begin(),ReindexList.end());
    int Nalt = ReindexScoreList.size();
    // *    std::vector<IsigI> IsigRef(Nalt);
    IsigI IsigRef;

    ResoRange RRange = ref_list.ResRange();
    bool verbose = false;
    if (PrintLevel > 0) verbose = true;
    // Note that after the following normalise operations, both reflection lists
    // have resolution limits set
    // These are reset at end of this procedure (ResetResoLimits)
    Normalise NormResRef =  NormaliseUnmerge(ref_list, GC, controls, verbose, output);
    Normalise NormResTest =  NormaliseUnmerge(test_list, GC, controls, verbose, output);
    // Set mininum resolution range for test list
    ResoRange minrange = ref_list.ResLimRange().MinRange(test_list.ResLimRange());
    test_list.SetResoLimits(minrange.ResLow(), minrange.ResHigh());

    reflection this_refl, ref_refl;
    observation this_obs;
    IsigI Isig;
    IsigI Isig0(0.0,0.0);
    int idx_refl;
    int isym;
    hkl_symmetry RefSymm = ref_list.symmetry();
    test_list.rewind();

    // loop all reflections from test list
    // Note that we assume that the symmetry of the reference list is correct,
    // but we cannot assume that the test list symmetry is, so don't merge test list Is
    while (test_list.next_reflection(this_refl) >= 0) {
      int Nobs = this_refl.num_observations();
      float sSqr = this_refl.invresolsq();      
      
      for (int lobs = 0; lobs < Nobs; lobs++)  {
	this_obs = this_refl.get_observation(lobs);
	if (this_obs.IsAccepted()) {
	  Hkl hkl = this_obs.hkl_original();
	  IsigI Isig = this_obs.I_sigI();
	  // fetch equivalent IsigI for possible hkl's from reference list
	  for (int i=0;i<Nalt;i++) {
	    idx_refl =
	      ref_list.get_reflection(ref_refl,
				      RefSymm.put_in_asu
				      (hkl.change_basis(ReindexScoreList[i]), isym));
	    
	    if (idx_refl >= 0) {
	      IsigRef = average_Is(ref_refl);
	    } else {
	      IsigRef = Isig0;
	    }
	    // add pairs of test, reference into scores
	    if (IsigRef.sigI() > 0.001) {
	      ReindexScoreList[i].addDsetIstats(Isig, IsigRef,
						NormResTest, NormResRef,
						sSqr);
	    }
	  }
	}
      }
    }
    // update ReindexScoreList to add likelihood
    GetTestUnmergedSignificance(test_list, NormResTest, ReindexScoreList, output);
    // Restore resolution limits which were reset by normalisation
    ref_list.ResetResoLimits();
    test_list.ResetResoLimits();
    std::sort(ReindexScoreList.begin(), ReindexScoreList.end());
    return ReindexScoreList;
  }
  //--------------------------------------------------------------
  void GetTestUnmergedSignificance(hkl_unmerge_list& test_list,
				   const Normalise& NormRes,
				   std::vector<ReindexScore>& ReindexScoreList,
				   phaser_io::Output& output)
  // Update ReindexScoreList with significance score (likelihood)
  {
    // ************************************************************
    // ***  Get significance of scores (cf testlauegroup)
    //  Get factor for sd(CC)
    ResoRange ResRange = test_list.ResRange();
    
    // Make list of unrelated pairs
    double ECC0;  // estimated E(CC)
    double CCsigFac = GetCCsigFac(test_list, ResRange, NormRes,
				  ECC0, false, output);
    //^
    //    std::cout << "ECC0, CCsigFac " << ECC0 <<" "<< CCsigFac << "\n"; //^-
    // Make list of CC-significance for each reindex op from scores list
    std::vector<scala::SetScores> scorelist;
    for (size_t i=0;i<ReindexScoreList.size();++i) {
      scorelist.push_back(ReindexScoreList[i]);
      //^
      //      std::cout << "ReindexScoreList " << i << " " << ReindexScoreList[i].OverallCC().result().val <<"\n";
      //^-
    }
    std::vector<SCsignificance> CCsig = CCsum(scorelist);
    // update CCsig with sd(CC) etc
    ScoreSig<correl_coeff>(CCsigFac, ECC0, CCsig);
    // Store probabilities back into ReindexScoreList
    ASSERT (ReindexScoreList.size() == CCsig.size());
    float totprob = 0.0; // total likelihood = 1.0
    for (size_t i=0;i<ReindexScoreList.size();++i) {
      //      std::cout << "Lkl " << CCsig[i].Likelihood() <<"\n"; //^
      totprob += CCsig[i].Likelihood();
    }
    for (size_t i=0;i<ReindexScoreList.size();++i) {
      ReindexScoreList[i].SetLikelihood(CCsig[i].Likelihood()/totprob);
    }
  }
  //--------------------------------------------------------------
} // namespace scala
