// twintest.cpp
//
// Padilla & Yeates L test
// Acta Cryst. D59,1124-1130 (2003)
// 

#include "twintest.hh"
#include "normalise.hh"
#include "printthings.hh"
#include "tablegraph.hh"

using phaser_io::LOGFILE;
using phaser_io::RESULT;

namespace scala {
  //--------------------------------------------------------------
  const double Twintest::ALPHATHRESHOLD = 0.09;  // threshold for twin warning
  //--------------------------------------------------------------
  Twintest::Twintest(const hkl_unmerge_list& hkl_list,
	       const double& MinIsigRatio,
	       phaser_io::Output& output)
  {
    init(hkl_list, MinIsigRatio, output);
  }
  //--------------------------------------------------------------
  void Twintest::init(const hkl_unmerge_list& hkl_list,
		   const double& MinIsigRatio,
		   phaser_io::Output& output)
  {
    // Intensity statistics only work for strong data, set I/sigma limit high
    // & cut resolution range
    TwinResRange = hkl_list.ResRange();
    TwinResRange.SetRange(hkl_list.num_observations());
    Rings DummyIceRing;
    Normalise NormRes = SetNormalise(hkl_list, MinIsigRatio, false, TwinResRange,
				     DummyIceRing, 0);
    //    PrintNormalisationResult(NormRes, TwinResRange, MinIsigRatio, output);

    // bins of |L| for cumulative distribution
    // |L| ranges from 0 to 1
    const int NLBINS = 20;
    lbins.SetRange(0.0, 1.0, true, NLBINS);
    lbcount.assign(NLBINS+1, 0);  // counts for cumulative distribution
    numl = 0;

    reflection this_refl;
    observation this_obs;
    Rtype Ic, Ic2;
    int irun;
    hkl_symmetry symm = hkl_list.symmetry();

    // Set up hkl increment list for looping over neighbours
    // Increment by 2 to avoid lattice absences and at least some
    // non-crystallographic translations
    std::vector<Hkl> hklinc;
    int maxinc = 2;
    int inc = 2;     //  +-2 by 2
    for (int h=-maxinc;h<=+maxinc;h+=inc) {
      for (int k=-maxinc;k<=+maxinc;k+=inc) {
	for (int l=-maxinc;l<=+maxinc;l+=inc) {
	  Hkl hkl(h,k,l);
	  if (hkl != Hkl(0,0,0)) {
	    hklinc.push_back(hkl);
	  }
	}}}

    hkl_list.rewind();
    // loop all reflections
    for (int index = 0;index<hkl_list.num_reflections();++index) {
      this_refl = hkl_list.get_reflection(index);
      // check for valid reflection
      // use only general reflections h!=k!=l!=0 to avoid centrics
      if ((this_refl.Status() == 0) && (this_refl.NvalidObservations() > 0) &&
	  (this_refl.hkl().IsGeneral())) {
	float sSqr = this_refl.invresolsq();
	int rbin = TwinResRange.tbin(sSqr);
	if (rbin >= 0) {  // check (restricted) resolution range
	  while (this_refl.next_observation(this_obs) >= 0) {
	    Hkl hkl = this_obs.hkl_original();
	    Ic = this_obs.kI();
	    irun = this_obs.run();
	    if (Ic != 0.0) {
	      // Loop neighbours
	      for (size_t i=0;i<hklinc.size();++i) {
		Hkl hkl2 = hkl + hklinc[i];
		if (hkl2.IsGeneral()) {
		  // find observation for index hkl2, in same run
		  // return true if found, Ic2 is intensity
		  bool found = GetNeighbour(hkl_list, hkl2, irun, symm, Ic2);
		  if (found) {
		    // Add into sums
		    Rtype L = Ic + Ic2;
		    if (L != 0.0) {
		      L = std::abs((Ic - Ic2)/L);  // |L|
		      if (L <= 1.0) {
			absl.Add(L);  // <|L|>
			l2.Add(L*L);            // <L^2>
			numl++;
			// Cumulative distribution of |L|
			int jb = lbins.bin(L);
			for (int i=jb;i<=NLBINS;++i) {
			  lbcount[i]++;
			}
		      }
		    }
		  }  // found
		} // general
	      } // neighbour loop
	    }
	  } // observation loop
	}   // resolution
      }
      //      std::cout << "Reflection " << this_refl.hkl().format() << " " <<numl
      //		<<" " << index <<"\n";
    } // reflection loop
    //^    std::cout <<"NumL" << numl <<"\n";
  }
  //--------------------------------------------------------------
  // find observation for index hkl2, in same run
  // return true if found, Ic2 is intensity
  bool Twintest::GetNeighbour(const hkl_unmerge_list& hkl_list,
			   const Hkl& hkl2, const int& irun,
			   const hkl_symmetry& symm,
			   Rtype& Ic2) const
  // Although there are multiple observations in each reflection, it is likely
  // that only one will match, and we can safely ignore any others, so we don't
  // need to buffer the reflection
  {
    int isym;
    reflection current_refl;  // current reflection
    observation this_obs;

    int index = hkl_list.get_reflection(current_refl, symm.put_in_asu(hkl2, isym));
    if (index < 0) {
      // no such reflection in list
      return false;
    }
    if ((current_refl.Status() == 0) && (current_refl.NvalidObservations() > 0)) {
      while (current_refl.next_observation(this_obs) >= 0) {
	if (this_obs.run() == irun) {  // check run
	  if (this_obs.hkl_original() == hkl2) {
	    // matching observation found
	    Ic2 = this_obs.kI();
	    if (Ic2 == 0.0) return false;
	    return true;
	  }}
      }
    }
    return false;  // not found
  }
  //--------------------------------------------------------------
  std::vector<double> Twintest::NL() const
  //! returns cumulative distribution of |L|
  {
    int nbins = lbins.Nbins();
    std::vector<double> nl(nbins);
    if (numl <= 0) return nl;
    for (int i=0;i<nbins;++i) {
      nl[i] = double(lbcount[i])/double(numl);
      //^      std::cout <<"i,N(|L|) " << i+1 <<" " <<nl[i]<<"\n"; //^
    }
    return nl;
  }
  //--------------------------------------------------------------
  //! return <|L|>
  double Twintest::MeanL() const
  {
    return absl.Mean();
  }
  //--------------------------------------------------------------
  //! return <L^2>
  double Twintest::MeanL2() const
  {
    return l2.Mean();
  }
  //--------------------------------------------------------------
  double Twintest::Alpha() const
  //! "Best" (maximum) estimate of twin fraction
  {
    return Max(EstimateAlphaFromAbsL(),
	       Max(EstimateAlphaFromAbsL(), EstimateAlphaFromL2()));
 
  }
  //--------------------------------------------------------------
  double Twintest::EstimateAlphaFromAbsL() const
  //! estimate alpha from <|L|>
  {
    const int NALPHAVALS = 101; // length of table on alpha 0.0 -> 0.5
    std::vector<double> MnLalpha(NALPHAVALS); // <|L|> for each value of alpha
    double alphainc = 0.5/double(NALPHAVALS-1); // 
    double alpha;
    double L = MeanL();
    double twinned = 0.375;  // N(|L|) for perfect twin

    if (L < twinned) {
      alpha = 0.5;
    } else {

      MnLalpha[0] = 0.5;
      //^
      //    alpha = 0.0;
      //    std::cout << "<|L|>(alpha) " << alpha <<" "<<MnLalpha[0]<<"\n";
      //^-
      for (int i=1;i<NALPHAVALS;++i) {
	alpha = Min(double(i) * alphainc, 0.5);  // 0 -> 0.5
	double den = (1.0 - 2.0*alpha);
	if (den <= 0.000001) {
	  MnLalpha[i] = 0.375;  // limiting value
	} else {
	  den = 2.0*den*den*den*den;
	  // P & Y equation 31
	  MnLalpha[i] = ((1.0 - 2.0*alpha)*(1.0 - 2.0*alpha)*
			 (1.0 - 6.*alpha +6.*alpha*alpha) -
			 8.0*(1. - alpha)*(1. - alpha)*alpha*alpha*
			 log(4.*alpha*(1. - alpha)))/den;
	}
	//^
	//^      std::cout << "<|L|>(alpha) " << alpha <<" "<<MnLalpha[i]<<"\n";
      }
      
      alpha = 0.0;
      for (int i=1;i<NALPHAVALS;++i) { // from 1
	if (L > MnLalpha[i]) {
	  //^	std::cout << i <<" " << L << " " << MnLalpha[i] <<"\n";
	  if (i == 1) {  // reset to 0 twin if in 1st bin
	    alpha = 0.0;
	  } else {
	    // Linear interpolation
	    alpha = Interpolate(L, MnLalpha[i], MnLalpha[i-1],
				double(i) * alphainc, double(i-1) * alphainc);
	    alpha = Min(alpha, 0.5);
	  }
	  break;
	}
      }
    }
    return alpha;
  }
  //--------------------------------------------------------------
  double Twintest::EstimateAlphaFromL2() const
  //! estimate alpha from <L^2>
  {
    const int NALPHAVALS = 101; // length of table on alpha 0.0 -> 0.5
    std::vector<double> L2alpha(NALPHAVALS); // <|L|> for each value of alpha
    double alphainc = 0.5/double(NALPHAVALS-1); // 
    double alpha = 0.0;    
    double L2 = MeanL2();
    double twinned = 0.2;  // N(L^2) for perfect twin

    if (L2 < twinned) {
      alpha = 0.5;
    } else {

      L2alpha[0] = 1.0/3.0;  // for untwinned
      //^
      //    std::cout << "<L^2>(alpha) " << alpha <<" "<<L2alpha[0]<<"\n";
      //^-
      for (int i=1;i<NALPHAVALS;++i) {
	alpha = Min(double(i) * alphainc, 0.5);  // 0 -> 0.5
	double den = (1.0 - 2.0*alpha);
	if (den <= 0.000001) {
	  L2alpha[i] = 0.2;  // limiting value at alpha = 0.5
	} else {
	  den = 3.0*den*den*den*den*den;
	  // P & Y equation 32
	  L2alpha[i] = ((1.0 - 2.0*alpha)*
			(1. + alpha*(-12. + alpha*(-4 +alpha*(+32. - 16.*alpha)))) +
			24.*(1. - alpha)*(1. - alpha)*alpha*alpha*
			log((1. - alpha)/alpha))/den;
	}
	//^
	//      std::cout << "<|L|>(alpha) " << alpha <<" "<<L2alpha[i]<<"\n";
      }
      
      alpha = 0.0;
      for (int i=1;i<NALPHAVALS;++i) { // from 1
	if (L2 > L2alpha[i]) {
	  if (i == 1) {  // reset to 0 twin if in 1st bin
	    alpha = 0.0;
	  } else {
	    alpha = Interpolate(L2, L2alpha[i], L2alpha[i-1],
				double(i) * alphainc, double(i-1) * alphainc);
	    alpha = Min(alpha, 0.5);
	  }
	  break;
	}
      }
    }
    return alpha;
  }
  //--------------------------------------------------------------
  MeanSD Twintest::EstimateAlphaFromNL() const
  //! estimate alpha from N(|L|) distribution
  //  mean & SD alpha for N(L) bins bewteen L 0.2 & 0.8
  {
    const int NALPHAVALS = 101; // length of table on alpha 0.0 -> 0.5
    std::vector<double> NLalpha(NALPHAVALS); // <|L|> for each value of alpha
    double alphainc = 0.5/double(NALPHAVALS-1); // 
    double alpha = 0.0;    
    double maxalpha = 0.0;
    // don't use the extreme |L| bins for estimating alpha
    const double MINL = 0.1999;
    const double MAXL = 0.7999;

    std::vector<double> nldistrib = NL();  // cumulative distribution N(|L|)

    int nbins = lbins.Nbins();
    std::vector<double> alphaestimate(nbins,0.0);
    MeanSD meanAlpha;

    for (int j=0;j<nbins;++j) { // loop over L bins
      double L = lbins.bounds(j).second;   // |L|
      double nl = nldistrib[j];            //N(|L|)
      //      double twinned = L * (3.0 - L*L)/2.0;  // N(|L| for perfect twin
      //      NLalpha[0] = L;  // for untwinned
      
      //^
      //      std::cout << "\nL = " << L << " N(L) " << nl <<"\n";
      //<< "\nN(L)(alpha) " << alpha <<" "<<NLalpha[0]<<"\n";
      //^-
      for (int i=0;i<NALPHAVALS;++i) { // loop alpha values for this |L|
	alpha = Min(double(i) * alphainc, 0.5);  // 0 -> 0.5
	
	double q  = (1.0 - alpha);
	double r = (1.0 - 2.0*alpha);
	if (r <= 0.000001) {
	  NLalpha[i] = 1.0;  // limiting value at alpha = 0.5
	} else {
	  // 2 * Integral of P(L), from P & Y equation 30
	  NLalpha[i] = (2.0/(r*r))*
	    (0.5*L*(alpha*alpha+q*q)
	     + alpha*q*q*(-1.0/(1.0+r*L) + 1.0/(r*(1.0+r*L)))
	     - alpha*alpha*q*(1.0/(1.0-r*L) + 1.0/(r*(1.0-r*L))));
	}
	//^
	//	std::cout << "N(L)(alpha) " << alpha <<" "<<NLalpha[i]<<"\n";
      } // alpha values for this L
      
      alpha = 0.0;
      for (int i=0;i<NALPHAVALS;++i) { // from 1
	if (nl < NLalpha[i]) {
	  if (i == 0) {
	    alpha = 0.0;
	  } else {
	    // Linear interpolation
	    alpha = Interpolate(nl, NLalpha[i], NLalpha[i-1],
				double(i) * alphainc, double(i-1) * alphainc);
	    alpha = Min(alpha, 0.5);
	  }
	  //	  std::cout << i <<" "<<alpha<<" "<<nl<<" "<<NLalpha[i]<<" iann\n";
	  break;
	}
      }
      alphaestimate[j] = alpha;
      if (L > MINL && L < MAXL) {
	// don't use the extreme |L| bins for estimating alpha
	maxalpha = Max(alpha, maxalpha);
	meanAlpha.Add(alpha);
      }
      //      std::cout << "Alpha " << alpha << " for L " << L <<"\n"; //^
    }  // loop |L|
    //^
    //    std::cout << "Mean & max alpha in range "
    //	      << MINL <<" - "<< MAXL
    //	      <<": " << meanAlpha.Mean() << " +/- " << meanAlpha.SD()
    //	      <<" " << maxalpha << "\n";
    //^-
    return meanAlpha;
  }
  //--------------------------------------------------------------
  bool Twintest::Enoughdata(const int& minimumcount) const
  //! return true if there is enough data for analysis
  {
    for (int i=1;i<lbins.Nbins();++i) {  // not 1st bin
      if (lbcount[i] < minimumcount) {
	return false;
      }
    }
    return true;
  }
  //--------------------------------------------------------------
  bool Twintest::PrintTwinResult(phaser_io::Output& output) const
  // returns true if might be twinned
  {
    if (numl == 0) return false;  // no data
    output.logTabPrintf(0,LOGFILE,
			"L-test for twinning (acentrics only) to maximum resolution %8.3f\n ",
			TwinResRange.ResHigh());

    // Check we have enough data
    const int MINIMUMCOUNT = 100;
    if (!Enoughdata(MINIMUMCOUNT)) {
      output.logTab(0,LOGFILE,
		    "Too few data to analyse twinning\n\n");
      return false;
    }

    bool twinned = false;

    PrintCumulativeNL(output);

    output.logTabPrintf(2,LOGFILE,
      "Estimated twin fraction alpha from cumulative N(|L|) plot%6.3f (%5.3f)\n",
			EstimateAlphaFromNL().Mean(),
			EstimateAlphaFromNL().SD());

    output.logTabPrintf(1,LOGFILE,
			"<|L|>: %8.3f (0.5 untwinned, 0.375 perfect twin)\n",
			MeanL());
    output.logTabPrintf(2,LOGFILE,
			"Estimated twin fraction alpha from <|L|> %8.3f\n",
			EstimateAlphaFromAbsL());
    output.logTabPrintf(1,LOGFILE,
			"<L^2>: %8.3f (0.333 untwinned, 0.2 perfect twin)\n",
			MeanL2());
    output.logTabPrintf(2,LOGFILE,
			"Estimated twin fraction alpha from <L^2> %8.3f\n",
			EstimateAlphaFromL2());

    if (EstimateAlphaFromAbsL() > ALPHATHRESHOLD || EstimateAlphaFromL2() > ALPHATHRESHOLD) {
      output.logTab(0,LOGFILE,
		    "\nWARNING: the L-test suggests that the data may be twinned,\nso the indicated Laue symmetry may be too high");
    } else {
      output.logTab(0,LOGFILE,
		    "\nThe L-test suggests that the data is not twinned");
    }
    output.logTab(0,LOGFILE,
		  "Note that the estimate of the twin fraction from the L-test is not very accurate,");
    output.logTab(0,LOGFILE,
		  "  particularly for high twin fractions. Better estimates from other test need knowledge of");
    output.logTab(0,LOGFILE,
		  "  the point group and the twin operator, which are not available here");
    output.logTab(0,LOGFILE,
		  "Also these statistics come from unscaled (and unmerged), so may be inaccurate for that reason\n\n");

    return twinned;
  }
  //--------------------------------------------------------------
  void Twintest::PrintCumulativeNL(phaser_io::Output& output) const
  //! print cumulative distribution of N(|L|)
  //  as TableGraph
  {
    int nbins = lbins.Nbins();

    std::string title = "L-test for twinning, twin fraction "+
      StringUtil::ftos(EstimateAlphaFromAbsL(),5,3);
    TableGraph table(title);
    output.logTab(0,LOGFILE,table.formatTitle());
    int c[] = {1,2,3,4};
    std::vector<int> cln(c,c+4);
    output.logTab(0,LOGFILE,
		table.Graph("Cumulative distribution function for |L|","N",cln));
    std::vector<std::string> collabels;
    collabels.push_back("|L|");          // 1
    collabels.push_back("N(|L|)");       // 2
    collabels.push_back("Untwinned");    // 3
    collabels.push_back("Twinned");      // 4
    collabels.push_back("Number");       // 5
    std::vector<bool> Zero(5, false);
    output.logTab(0,LOGFILE,
		  table.ColumnFields(collabels, Zero,
				     " %10.4f %10.4f %10.4f %10.4f %8d\n")); 
    int nc = collabels.size();
    output.logTab(0,LOGFILE,
		  table.Line(nc, 0.0, 0.0, 0.0, 0.0, 0));
    for (int i=0;i<nbins;++i) {
      double L = lbins.bounds(i).second;
      double twinned = L * (3.0 - L*L)/2.0;  // N(|L| for perfect twin
      double nl = double(lbcount[i])/double(numl);
      output.logTab(0,LOGFILE,
		    table.Line(nc, L, nl, L,twinned, lbcount[i]));
    }
    output.logTab(0,LOGFILE,
		  table.CloseTable());
  }
  //--------------------------------------------------------------
  void Twintest::WriteResultSummary(phaser_io::Output& output) const
  //! write summary data
  {
    if (numl == 0) return;  // no data
    // Check we have enough data
    const int MINIMUMCOUNT = 100;
    if (!Enoughdata(MINIMUMCOUNT)) {
      output.logTab(0,RESULT,
		    "\nToo few data to analyse twinning\n");
      return;
    }

    if (EstimateAlphaFromAbsL() > ALPHATHRESHOLD ||
	EstimateAlphaFromL2() > ALPHATHRESHOLD ||
	EstimateAlphaFromNL().Mean() > ALPHATHRESHOLD) {
      output.logTab(0,RESULT,
		    "\nWARNING: the L-test suggests that the data may be twinned,\n    so the indicated Laue symmetry may be too high");
      output.logTabPrintf(1,RESULT,
	"Rough estimated twin fraction alpha from cumulative N(|L|) plot%6.3f (%5.3f)\n",
			EstimateAlphaFromNL().Mean(),
			EstimateAlphaFromNL().SD());
      output.logTabPrintf(1,RESULT,
			  "Rough estimated twin fraction alpha from <|L|> %8.3f\n",
			  EstimateAlphaFromAbsL());
      output.logTabPrintf(1,RESULT,
			  "Rough estimated twin fraction alpha from <L^2> %8.3f\n",
			  EstimateAlphaFromL2());
    } else {
      output.logTab(0,RESULT,
		    "\nThe data do not appear to be twinned, from the L-test");
    }
  }
  //--------------------------------------------------------------
  double Twintest::Interpolate(const double& val,
			       const double& lowerbin, const double& upperbin,
			       const double& lowerval, const double& upperval) const
  // interpolate between lowerval and upperval according to where
  // val lies between lowerbin and upperbin
  {
    double frac = (val - lowerbin)/(upperbin - lowerbin);
    return frac*(upperval - lowerval) + lowerval;
  }
  //--------------------------------------------------------------
} // namespace scala
