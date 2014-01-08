//
// sysabsstats.cpp


#include "scala.hh"
#include "hkl_merged_list.hh"

#include "latsym.hh"
#include "sysabszones.hh"
#include "zone.hh"
#include "intensitystatistics.hh"

//#include "printthings.hh"

namespace scala
{
  //--------------------------------------------------------------
  double SigmaMod(const double& sigma, const double& sigmin)
  // Return modified sigma inflating small values, using sigmin
  {
    if (sigma <= 0.0) return sigmin;
    return (sigma*sigma + sigmin)/sigma;
  }
  //--------------------------------------------------------------
  void SysAbsStats(const hkl_unmerge_list& hkl_list,
		   std::vector<Zone>& SZones)
  // Accumulate statistics for systematic absence test
  //
  // Uses unaveraged intensities, I/sig(I)
  {
    if (SZones.size() == 0) return;
    // For the  control set of non-zonal data, there is no point
    // in using all of a huge file
    // Get approximately MAXGENLIST reflections by selected every
    // nchoose'th one (gives < 2 * MAXGENLIST)
    const int MAXGENLIST = 10000;
    int nchoose = Max(hkl_list.num_reflections()/MAXGENLIST, 1);
    // Maximum number of reflections to be considered in
    // each control sample
    const int MAXSAMPLE = 200;
    // Minimum number of samples
    const int MINNUMSAMPLE = 8;
    // Minimum number in sample
    const int MINNUMBER = 8;

    double sigmin = hkl_list.MinSigma();
    sigmin = 2.0*sigmin+1.0;

    reflection this_refl;
    IsigI Isig;
    bool axisinzone;
    bool inzone;

    int izr = 0;
    int kode;  // hkl code, packed reduced hkl
    // IKode is intensity + hkl packed code + InvResolution^2
    std::vector<IKode> IKlist;
    hkl_list.rewind();

    // Set up zone list by the Laue group they belong to
    std::vector<int> ZoneIndexLG;
    ZoneIndexLG.push_back(0);  // 1st group
    hkl_symmetry symm = SZones[0].LGsymm();
    for (size_t iz=1;iz<SZones.size();iz++) { // loop rest of zones
      if (SZones[iz].LGsymm() != symm) {
	// Changed Laue group
	ZoneIndexLG.push_back(iz);
      }
    }

    while (hkl_list.next_reflection(this_refl) >= 0) {
      float sSqr = this_refl.invresolsq();
      int Nobs = this_refl.num_observations();
      for (int l = 0; l < Nobs; l++) {
	observation this_obs = this_refl.get_observation(l);
	if (this_obs.IsAccepted()) {
	  axisinzone = false;
	  inzone = false;
	  Isig = this_obs.I_sigI();
	  Hkl hkl = this_obs.hkl_original();
	  int lg=0;
	  // Loop zones for systematic absence test
	  for (size_t iz=0;iz<SZones.size();iz++) {
	    // Is this for a new Laue group?
	    if (lg < int(ZoneIndexLG.size())) {
	      if (int(iz) == ZoneIndexLG[lg]) {
		// Yes, reset axisFound flag
		lg++;
		axisinzone = false;
	      }
	    }
	    if (!axisinzone || !SZones[iz].Axis()) {
	      // For each Laue group, axes come before glides,
	      // so an axial reflection
	      // will not be included in glide zone
	      // Reflection may belong to more than one glide
	      
	      // test original hkl in lattice frame
	      if (SZones[iz].InZone(hkl)) {
		if (SZones[iz].Axis()) axisinzone = true;		  
		inzone = true;
		if (Isig.sigI() > 0.0) {
		  // Modify sigma to avoid very small values
		  Isig.sigI() = SigmaMod(Isig.sigI(), sigmin);
		  SZones[iz].AddRef(hkl, Isig, sSqr);
		  //^
		  //		  std::cout << "In Zone " << iz << SZones[iz].formatLGFrame()
		  //			    << " " << hkl.format() << " "
		  //			    << Isig.I() << " " <<  Isig.sigI() << "\n";
		  //^-
		}
	      }
	    }
	  }
	  if (! inzone && izr%nchoose == 0) {
	    // Not in zone & selected by every nchoose'th one
	    // Store non-zone observations as a control
	    if (Isig.sigI() > 0.0) {
	      kode = hkl.code();
	      IKlist.push_back(IKode(Isig, kode, sSqr));
	    }
	  }
	  izr++;
	}
      }
    }
    
    // Sort control list on resolution
    std::sort(IKlist.begin(), IKlist.end());
  
    //*  std::cout << IKlist.size() << "  reflections in control list\n";
    // Loop Zones  
    for (size_t iz=0;iz<SZones.size();iz++) {
      // Get a sample of null scores for the same number of reflections
      // as were included in the scores
      if (SZones[iz].Nobs() > 0) {
	std::vector<int> IndxList = SZones[iz].Indices();
	int llength = IndxList.size();
	if (llength > 1) {
	  // Number of reflections to be used in each sample
	  int Nr = Min(SZones[iz].Nobs(), MAXSAMPLE);
	  Range InvResRange = SZones[iz].InvResoRange();
	  //			  std::cout << "\nZone " << iz
	  //			  << "  Number of observations  " << SZones[iz].Nobs()
	  //			  << "   resolution range " << InvResRange.min()
	  //			  << " to " << InvResRange.max() << "\n";
	  // How many reflections in control list between resolution
	  // limits of this zone sample?
	  int I1 = -1;
	  int I2 = -1;
	  for (size_t i=0;i<IKlist.size();i++) {
	    if (I1 < 0) {
	      if (IKlist[i].sSqr > InvResRange.min()) I1 = i;
	    } else {
	      if (I2 < 0) {
		if (IKlist[i].sSqr > InvResRange.max()) {
		  I2 = i;
		  break;
		}
	      }
	    }
	  }
	  if (I1 < 0) I1 = 0;
	  if (I2 < 0) I2 = IKlist.size();
	  int NumControl = I2-I1;
	  
	  //		std::cout << "Number of control reflections in resolution range : "
	  //			  <<  NumControl << " " << I1 << " " << I2 << "\n";
	  int minControl = MINNUMSAMPLE * MINNUMBER;
	  if (NumControl < minControl) {
	    // Too few control reflections in resolution range,
	    // try to expand range to find some more
	    if (int(IKlist.size()) > minControl) {
	      I1 = Max(0, I1 - minControl/2);
	      I2 = I1 + minControl;
	      NumControl = I2-I1;
	    } else {
	      NumControl = 0; // can't do it
	    }
	  }
	  
	  double cSD = 0.4; // default value in case of insufficient information
	  double cMean  = 0.0;
	  
	  int NgPoints = SZones[iz].NgridPoints();
				      
	  std::vector<double> controlSD(NgPoints, cSD);;
	  std::vector<double> controlMean(NgPoints, cMean);
	  
	  
	  if (NumControl > 0) {
	    int NSample = NumControl/Nr;  // Number of samples
	    if (NSample < MINNUMSAMPLE)  {
	      NSample = MINNUMSAMPLE;
	      Nr = NumControl/NSample;
	    }
	    // Now divide up the control (non-zone) reflections into
	    // NSample groups of size Nr, each running over the equivalent
	    // resolution range to the data, and get the Fourier scores
	    
	    //		  std::cout << "Number & size of control groups "
	    //			    << NSample << " " << Nr << "\n";
	    
	    // If very few control reflections, just use default SD
	    
	    std::vector<std::vector<double> > sampleFourierlist;
	    // Loop samples
	    for (int is=0;is<NSample;is++) {
	      int ir = I1+is;
	      OneDFourier F(SZones[iz].Ngrid());  // initialise 1-D Fourier
	      // Loop reflections in sample	  
	      for (int i=0;i<Nr;i++) {
		int j = i%llength;
		// Modify sigma to avoid very small values
		IKlist.at(ir).Is.sigI() = SigmaMod(IKlist.at(ir).Is.sigI(), sigmin);
		// Add "score value" into Fourier series
		// negative observations are truncated to zero
		F.AddRef(IndxList[j], IscoreVal(IKlist.at(ir).Is));
		ir += NSample;
		//^
		//		      std::cout << iz << " " << i << " : " << IndxList[j]
		//				<< " " << IKlist.at(ir).Is.I() << " "
		//				<< IscoreVal(IKlist.at(ir).Is) << "\n";
		//^-
	      }
	      std::vector<double> Fval = F.FourierVal();
	      if (Fval[0] > 0.0) {
		sampleFourierlist.push_back(Fval);  // store only if origin is non-zero
		//^
		//		      std::cout << "Zone: " << iz;
		//		      for (int ii=0;ii<Fval.size();++ii) {
		//			std::cout << " " << Fval[ii];
		//		      }
		//		      std::cout << "\n";
		//^-
	      }
	    }
	    // Mean & sd of Fourier scores at each Fourier spacing
	    // Mean should be around 0
	    int NsamplePoints = sampleFourierlist.size();
	    
	    for (int i=1;i<NgPoints;i++) {
	      std::vector<double> val(NsamplePoints);
	      for (int j=0;j<NsamplePoints;j++) {
		val[j] = sampleFourierlist[j][i];
	      }
	      controlSD[i] = MeanSD(val).SD();
	      controlMean[i] = MeanSD(val).Mean();
	      //^
	      //	      std::cout << SZones[iz].formatRefFrame(i)
	      //			<< " nobs= " << SZones[iz].Nobs()
	      //			<< "\n";
	      //	      printf("  %4d %8.3f control sd %8.3f  Nsample %6d size %6d\n\n",
	      //		     SZones[iz].Ngrid()[i], SZones[iz].FourierVal()[i],
	      //		     controlSD[i], NSample, Nr);
	      //^-	      
	    }
	  }		
	  SZones[iz].StoreMeanSD(controlMean, controlSD);
	}
      }  // if no observations in zone
    }  // end loop zones
  }
} // namespace scala
