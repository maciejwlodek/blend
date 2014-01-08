//  intensitystatistics.cpp
//

#include "intensitystatistics.hh"
#include "scala_util.hh"
#include "probfunctions.hh"
#include "icering.hh"
#include "printthings.hh"

namespace scala
{
  //=================================================================
  IntensityStatistics::IntensityStatistics(const hkl_unmerge_list& hkl_list,
					   const double& MinIsigRatio,
					   phaser_io::Output& output)
  // get mean |E^2-1| for general reflections, in resolution bins
  // MinIsigRatio   minimum I/sigI ratio on unaveraged data, default 6
  {
    // Intensity statistics only work for strong data, set I/sigma limit high
    // & cut resolution range
    Rrange = hkl_list.ResRange();
    Rrange.SetRange(hkl_list.num_observations());
    Rings DummyIceRing;
    NormRes = SetNormalise(hkl_list, MinIsigRatio, false, Rrange, DummyIceRing, 0);
    PrintNormalisationResult(NormRes, Rrange, MinIsigRatio, output);
    //    printf("\nResolution range reset to %8.2f to %8.2f\n",
    //	   Rrange.ResLow(), Rrange.ResHigh());

    reflection this_refl;
    observation this_obs;
    IsigI Isig;

    // Divide randomly into NRbins
    int NRbins = 20;

    E2minus1 = std::vector<double>(NRbins, 0.0);
    num = std::vector<int>(NRbins, 0);
    E2mean = 0.0;
    int nm = 0;
    Emean = 0.0;

    hkl_list.rewind();
    // loop all reflections
    while (hkl_list.next_reflection(this_refl) >= 0) {
      // use only general reflections h!=k!=l!=0 to avoid
      // problems with unknown epsilon
      if (this_refl.hkl().IsGeneral()) {
	float sSqr = this_refl.invresolsq();
	int rbin = Rrange.tbin(sSqr);
	if (rbin >= 0) {
	  while (this_refl.next_observation(this_obs) >= 0) {
	    IsigI E2sig = NormRes.applyAvg(this_obs.I_sigI(), sSqr);
	    //^
	    //^		    printf("%5d%5d%5d%10.3f%10.3f\n",
	    //^		    this_refl.hkl().h(),this_refl.hkl().k(),this_refl.hkl().l(),
	    //^			   sqrt(std::abs(E2sig.I())), std::abs(E2sig.I() - 1.0));
	    
	    if (NormRes.NotTooLarge(E2sig)) {
	      unsigned int rbin = IRandom(NRbins);
	      E2minus1[rbin] += std::abs(E2sig.I() - 1.0);
	      num[rbin]++;
	      E2mean += E2sig.I();
	      Emean += sqrt(std::abs(E2sig.I()));
	      nm++;
	    }
	  }
	}
      }
    }

    TotE2minus1 = 0.0;
    n=0;

    for (size_t i=0;i<E2minus1.size();i++) {
      if (num[i] > 0) {
	TotE2minus1 += E2minus1[i];
	n += num[i];
	E2minus1[i] /= num[i];
      }
    }
    TotE2minus1 /= n;
    // Mean & SD of samples
    RPair meansd_E2m1 = MnSd(E2minus1);
    SD_E2m1 = meansd_E2m1.second;
    if (nm > 0)  {
      E2mean /= nm;
      Emean /= nm;
    }
    // Probability
    double mcent = 0.968;
    double macen = 0.736;
    DM_2sqrt DMcentric(macen, mcent);
    double pcentro = ProbBiassed(TotE2minus1, SD_E2m1, mcent, DMcentric);


    //    double dc = Max(mcent - TotE2minus1, 0.0)/SD_E2m1;
    //    double da = Max(TotE2minus1 - macen, 0.0)/SD_E2m1;
    //    double pcent = exp(-0.5*dc*dc);
    //    pacentric = exp(-0.5*da*da);
    //    double sn = pcent+pacentric;
    pacentric = 1.0 - pcentro;
  }
  //--------------------------------------------------------------
  void IntensityStatistics::Print(phaser_io::Output& output) const
  {
    output.logTab(0,LOGFILE,
      "\n\nIntensity statistics from general reflections for centrosymmetric test\n");
    output.logTab(0,LOGFILE,
		   "\n            Observed     Theoretical: acentric  centric\n");

    output.logTabPrintf(0,LOGFILE,
			 "\n<|E^2-1|> %8.3f sd %5.3f            0.736     0.968 \n", TotE2minus1, SD_E2m1);
    output.logTabPrintf(0,LOGFILE,
		  	   "<|E|>     %8.3f                     0.886     0.798 \n", Emean);
    output.logTabPrintf(0,LOGFILE,
			   "<|E^2|>   %8.3f                     1.0       1.0\n ", E2mean);
    output.logTabPrintf(0,LOGFILE,
	 "\np(non-centrosymmetric) %8.3f  p(centrosymmetric) %8.3f\n",
			 pacentric,1.0-pacentric);
    if (pacentric > 0.5)
      output.logTabPrintf(0,LOGFILE,"                          =====\n");
    else
      output.logTabPrintf(0,LOGFILE,"                                                       =====\n");
  }
}
