// normaliseunmerge.cpp

#include "normaliseunmerge.hh"
#include "printthings.hh"

namespace scala {
//--------------------------------------------------------------
Normalise NormaliseUnmerge(hkl_unmerge_list& hkl_list,
                           const GlobalControls& GC,
			   const all_controls& controls, 
			   const bool& verbose,
                           phaser_io::Output& output)
// Organise hkl_list & return normalisation object
{
  //  output.logTab(0,LOGFILE, "\n>>>> Normalising test dataset");
  // sort & organise reflection list as required
  int Nrefl = hkl_list.prepare();
  // Sum partials, return number of partials
  int Npart = hkl_list.sum_partials();
  Npart = Npart;
  int Nobs = hkl_list.num_observations();

  //  if (Npart > 0)
  //    controls.partials.Print(output);
  float ratio = float(Nobs)/float(Nrefl);

  if (verbose) {
    output.logTabPrintf(0,LOGFILE,"\nNumber of reflections  =        %10d\n",
			hkl_list.num_reflections_valid());
    output.logTabPrintf(0,LOGFILE,    "Number of observations =        %10d\n",Nobs);
    if (hkl_list.num_observations_scaled() > 0)  {
      output.logTabPrintf(0,LOGFILE,    "Number of scaled observations = %10d\n",
			  hkl_list.num_observations_scaled());
    }
    output.logTabPrintf(0,LOGFILE,
			"Average multiplicity =            %8.1f\n",ratio);
    output.logTabPrintf(0,LOGFILE,"\nResolution range in list:  %9.2f ->%7.2f\n",
			hkl_list.ResRange().ResLow(), hkl_list.ResRange().ResHigh());
  }

  // ************************************************************
  // Set up resolution bins
  ResoRange ResRange(hkl_list.RRange().min(), hkl_list.RRange().max(),Nobs);
  // ************************************************************
  // Go through data to get normalisation correction
  int PrintLevel = 0;   //^
  // Ice rings
  Rings Icerings;
  Icerings.DefaultIceRings();

  // Resolution limits may be reset if Min(I/sigI) > 0
  Normalise NormRes = SetNormalise(hkl_list, GC.GetMinIsig(), false,
                                   ResRange, Icerings, PrintLevel);
  hkl_list.SetIceRings(Icerings);

  if (verbose) {
    PrintNormalisationResult(NormRes, ResRange, GC.GetMinIsig(), output);
  }
  // Set resolution limits for reflection list
  hkl_list.SetResoLimits(ResRange.ResLow(), ResRange.ResHigh());
  return NormRes;
}
}
