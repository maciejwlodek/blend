// centre.cpp
//
// Find centre of diffraction pattern
// This Friedel centre should be at h=(0,0,0) if the intensities
// were properly indexed
// Use Rmeas as score function (CC also calculated, but commented out)


#include "centre.hh"
#include "pointgroup.hh"
#include "outputunmergeddata.hh"
#include "scala_util.hh"

#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;

using phaser_io::LOGFILE;
using phaser_io::LXML;

namespace scala
{
  void FindCentre(hkl_unmerge_list& hkl_list,
		  const std::string& LaueGroupName,
		  const ReindexOp& reindex,
		  const double& MinIsigRatio,
		  const clipper::Vec3<int>& GridMax,
		  const std::string& hklout_filename,
		  phaser_io::Output& output)
  // On entry:
  //  hkl_list            data list
  //  LaueGroupName       Laue group name
  //  reindex             reindex operator
  //  MinIsigratio        for resolution cut-off
  //  GridMax             maximum indices for grid search
  //  hklout_filename     filename for optional output file
  //  output              output object
  // 
  {
    std::string LGname = LaueGroupName;
    if (LGname != "" && LGname != "ALL")
      {
	// LAUEGROUP given on input, set with dummy scores
	// option to keep original Laue group
	if (LGname == "HKLIN")
	  {
	    LGname = hkl_list.symmetry().symbol_xHM();
	  }
	CCtbxSym::PointGroup LG(LGname);
	output.logTab(0,LOGFILE, "Reindexing or changing symmetry\n");
	output.logTab(0,LOGFILE, "Reindex operator from input cell to given Lauegroup: "
		      + reindex.as_hkl() + "\n\n"
		      + reindex.as_matrix() + "\n");
	
	// reindex, sort & reorganise hkl list
	int Junk = hkl_list.change_symmetry(hkl_symmetry(LG.RefLGname()), reindex);
	Junk = Junk;
      }
    // Search for centre
    ReindexOp reindex_cen = Centre(hkl_list, MinIsigRatio, GridMax, output);
    if (hklout_filename != "")
      {
	// Write out unmerged reindexed file with best reindexing,
	output.logTabPrintf(0,LOGFILE,
			    "\n\n---------------------------------------------------------------\n\n");
	output.logTab(0,LOGFILE,
		      "Writing unmerged test data to file "+hklout_filename);
	output.logTab(0,LOGFILE,
		      " with reindexing operator "+
		      reindex_cen.as_hkl()+"\n");
	output.logTab(0,LOGFILE,
		      " in space group "+hkl_list.symmetry().symbol_xHM());
	output.logTab(0,LOGFILE,"\n");
	OutputUnmergedData(hkl_list, hklout_filename,
			   hkl_list.symmetry(), reindex_cen, true,
			   hkl_list.Title(),"space",output);
      }
  }
  // ------------------------------------------------------------------------
  
ReindexOp Centre(hkl_unmerge_list& hkl_list,
		   const double& MinIsigRatio,
		   const clipper::Vec3<int>& GridMax,
		   phaser_io::Output& output)
  // Returns reindexing operation to fix up data
  //
  // On entry:
  //  hkl_list            data list
  //  MinIsigratio        for resolution cut-off
  //  GridMax             maximum indices for grid search
  //  output              output object
  {
    output.logTab(0,LOGFILE, "\nSearching for centre of diffraction pattern\n\n");
    output.logTabPrintf(0,LOGFILE,
			 "\nGrid search limits %4d %4d %4d\n",
			 GridMax[0], GridMax[1], GridMax[2]);
    double Zmin = 1.5;  // Significance level for Z-score
    // sort & organise reflection list as required
    int Nrefl = hkl_list.prepare();
    Nrefl = Nrefl;
    // Sum partials, return number of valid observations
    int Nobs = hkl_list.sum_partials();
    Nobs = Nobs;
    ResoRange ResRange = hkl_list.ResRange();
    Rings DummyRings;
    // Resolution limits may be reset
    Normalise NormRes = SetNormalise(hkl_list, MinIsigRatio, false,
				     ResRange, DummyRings, 0);
    // Set resolution limits for reflection list
    hkl_list.SetResoLimits(ResRange.ResLow(), ResRange.ResHigh());
    output.logTabPrintf(0,LOGFILE, "\nResolution range reset to %8.2f to %8.2f\n",
			 ResRange.ResLow(), ResRange.ResHigh());

    int hmax=0;
    int kmax=0;
    int lmax=0;
    reflection this_refl;
    observation this_obs;
    Hkl hkl;
    clipper::Coord_grid cg;
    Scell cell = hkl_list.cell();
    //    float smax = hkl_list.ResRange().SResHigh();
    //    const float unity = 1.0;

    int isymt;
    hkl_symmetry symm = hkl_list.symmetry();
    bool centro = symm.IsCentro();  // centrosymmetric
    // number of ISYM flags 
    //   = 2* Nsymop for non-centrosymmetric space group
    //   =    Nsymop for centrosymmetric space group
    int nisym = symm.NsymP();
    if (!centro) nisym *=2;

    hkl_list.rewind();  // point to beginning of list

    // Maximum hkl (original) indices in list
    while (hkl_list.next_reflection(this_refl) >= 0)
      {
	while (this_refl.next_observation(this_obs) >= 0)
	  {
	    hkl = this_obs.hkl_original();
	    hmax = Max(hmax, std::abs(hkl.h()));
	    kmax = Max(kmax, std::abs(hkl.k()));
	    lmax = Max(lmax, std::abs(hkl.l()));
	  }
      }

    Array3D<float> Data(2*hmax+1,2*kmax+1,2*lmax+1,-hmax,-kmax,-lmax);
    Data.Zero();

    output.logTabPrintf(0,LOGFILE,
			 "Storing data for maximum h,k,l %5d %5d %5d\n",
			 hmax, kmax, lmax);

    // Loop observations & store
    while (hkl_list.next_reflection(this_refl) >= 0)
      {
	while (this_refl.next_observation(this_obs) >= 0)
	  {
	    hkl = this_obs.hkl_original();
	    Data(hkl) = this_obs.I();
	  }
      }


    // Search around origin
    int hgrid = GridMax[0];
    int kgrid = GridMax[1];
    int lgrid = GridMax[2];

    int nh = 2*hgrid+1;
    int nk = 2*kgrid+1;
    int nl = 2*lgrid+1;
    Array3D<Rfactor> Rfac(nh,nk,nl, -hgrid, -kgrid, -lgrid);
    //    double w=1.0;

    // Loop search grid
    for (int ih=-hgrid;ih<=hgrid;ih++) {
      for (int ik=-kgrid;ik<=kgrid;ik++) {
	for (int il=-lgrid;il<=lgrid;il++) {
	  //  test for lattice absence in grid
	  if (symm.LatticePresent(Hkl(ih,ik,il))) {
	    // reset origin of data array
	    Data.Offset(-hmax+ih,-kmax+ik,-lmax+il);
	    
	    // Loop data
	    for (int kh=-hmax;kh<=hmax;kh++) {
	      for (int kk=-kmax;kk<=kmax;kk++) {
		for (int kl=-lmax;kl<=lmax;kl++) {
		  Hkl hkl(kh,kk,kl);
		  // matching pairs of data
		  if (Data.Test(hkl) && 
		      Data(hkl) > 0.0) {
		    Hkl hklt = symm.put_in_asu(hkl, isymt);
		    // Find symmetry-related pairs
		    std::vector<Hkl> symh;
		    std::vector<int> isyms;
		    symh.push_back(hkl);
		    isyms.push_back(isymt);
		    for (int is=1;is<=nisym;is++) {  // is 1:nisym
		      Hkl hkls = symm.get_from_asu(hklt, is);
		      if (hkls != hkl &&
			  Data.Test(hkls) &&
			  Data(hkls) > 0.0) {
			bool newIsym = true;
			for (size_t js=0;js<symh.size();js++) {
			  if (is == isyms[js]) {newIsym = false; break;}
			}
			if (newIsym) {
			  symh.push_back(hkls);
			  isyms.push_back(is);
			}
		      }
		    }
		    // Calculate Rfactor (using unweighted mean)
		    // weight for Rfactor (Rmeas) = sqrt(n/n-1)
		    int ns = symh.size();
		    if (ns > 1) {
		      double an = ns;
		      double w = 1.0/(an-1.0);
		      double wtm = sqrt(an*w);
		      double sumi = 0.0;
		      for (int i=0;i<ns;i++)
			{sumi += Data(symh[i]);}
		      double Imean = sumi/double(ns);  // unweighted mean
		      for (int i=0;i<ns;i++) {
			Rfac(ih,ik,il).add((Data(symh[i])-Imean), Imean, wtm);
		      }
		    }
		  }
		}
	      }
	    }  // data loop
	  }
	}
      }
    }  // grid loop

    // Print region around origin
    //cc    double ccmax=0.0;
    double rfacmin =  1000.;
    Hkl hkltop(0,0,0);
    const int NRmin = 4;  // minimum number of of observations to count
    std::vector<float> Rlist;

    // Loop grid to find best point
    for (int il=-lgrid;il<=lgrid;il++) {
      for (int ik=-kgrid;ik<=kgrid;ik++) {
	double R;
	int N;
	for (int ih=-hgrid;ih<=hgrid;ih++)  {
	  R = Rfac(ih,ik,il).result().val;
	  N = Rfac(ih,ik,il).result().count;
	  if (N > NRmin) {
	    Rlist.push_back(R);
	    if (R < rfacmin)  {
	      rfacmin = R;
	      hkltop = Hkl(ih,ik,il);
	    }
	  }
	}
      }
    }

    MeanSD Rmean(Rlist);
    double Rfacmean = Rmean.Mean();
    double sdR = Rmean.SD();
    float Zscore = 0.0;
    float Zcentre = 0.0;
    if (sdR > 0.00001) 
      {
	// Z relative to mean
	Zscore = std::abs(rfacmin - Rfacmean)/sdR;
	// Z relative to centre
	if (hkltop != Hkl(0,0,0))
	  {
	    if (Rfac(0,0,0).result().count > NRmin)
	      {Zcentre = std::abs(rfacmin - Rfac(0,0,0).result().val)/sdR;}
	  }
      }
    // Print
    output.logTab(0,LOGFILE,"\nSearch grid on hkl around 0,0,0\n\n");
    output.logTabPrintf(0,LOGFILE,
			 "                   R-factors                          Number of observations");

    for (int il=-lgrid;il<=lgrid;il++)
      {
	output.logTabPrintf(0,LOGFILE,"\nl = %3d\n    h  ",il);
	for (int ih=-hgrid;ih<=hgrid;ih++)
	  {output.logTabPrintf(0,LOGFILE," %5d ", ih);}
	output.logTabPrintf(0,LOGFILE,"\n");
	for (int ik=-kgrid;ik<=kgrid;ik++)
	  {
	    bool cenkl = (ik==0 && il==0);
	    bool bestkl = (ik== hkltop.k() && il==hkltop.l());
	    output.logTabPrintf(0,LOGFILE," k %3d  ",ik);
	    std::vector<double> Rh(nh);
	    std::vector<int> Nh(nh);
	    for (int ih=-hgrid;ih<=hgrid;ih++)
	      {
		Rh[ih+hgrid] = Rfac(ih,ik,il).result().val;
		Nh[ih+hgrid] = Rfac(ih,ik,il).result().count;
	      }
	    for (int j=0;j<nh;j++)
	      {
		//cc		output.logTabPrintf(0,LOGFILE," %6.2f",CCh[j]);
		if (bestkl && j == hgrid+hkltop.h())
		  { // best
		    output.logTabPrintf(0,LOGFILE,"*%5.3f*",Rh[j]);}
		else if (cenkl && j == hgrid)
		  { // central point, 0,0,0
		    output.logTabPrintf(0,LOGFILE,">%5.3f<",Rh[j]);}
		else
		  {output.logTabPrintf(0,LOGFILE," %5.3f ",Rh[j]);}
	      }
	    output.logTabPrintf(0,LOGFILE,"   ");
	    for (int j=0;j<nh;j++)
	      {
		output.logTabPrintf(0,LOGFILE," %6d",Nh[j]);
	      }
	    output.logTabPrintf(0,LOGFILE,"\n");
	  }
      }

    output.logTab(0,LOGFILE,"\n<!--SUMMARY_BEGIN-->\n\n");
    //cc    output.logTabPrintf(0,LOGFILE, "\nMaximum CC %6.2f at %s\n",ccmax,hkltop.format().c_str());
    output.logTabPrintf(0,LOGFILE,
			 "\nMinimum Rfactor %6.3f at lattice shift %3d %3d %3d\n",
			 rfacmin, hkltop.h(), hkltop.k(), hkltop.l());
    output.logTabPrintf(1,LOGFILE,
			     "Z-score relative to mean   = %5.1f\n", Zscore);
    if (hkltop != Hkl(0,0,0))
      {
	output.logTabPrintf(1,LOGFILE,
			     "Z-score relative to centre = %5.1f\n", Zcentre);
      }

    clipper::Vec3<double> v(hkltop.h(), hkltop.k(), hkltop.l());
    ReindexOp reindex_cen(clipper::Mat33<double>::identity(), v);
    output.logTab(0,LXML, "<Centre>");
    output.logTabPrintf(1,LXML,
			 "<LatticeShift> %3d %3d %3d </LatticeShift>\n",
			 hkltop.h(), hkltop.k(), hkltop.l());
    std::string signif = "no";
    if (Zcentre > Zmin) {signif = "yes";}
    output.logTabPrintf(1,LXML,
			 "<Zscore significant=%s> %6.2f </Zscore>\n", signif.c_str(), Zcentre);

    output.logTab(1,LXML,reindex_cen.as_hkl_XML()+"\n");


    if (hkltop == Hkl(0,0,0))
      {
	output.logTabPrintf(0,LOGFILE,
			     "\n>>>> Correct centre of diffraction pattern confirmed <<<<\n\n");
	output.logTab(1,LXML, "<CorrectCentre> true </CorrectCentre>");
      }
    else
      {
	output.logTab(1,LXML, "<CorrectCentre> false </CorrectCentre>");
	output.logTabPrintf(0,LOGFILE,
			     "\n **** WARNING, WARNING WARNING ****\n\n");
	if (Zcentre > Zmin)
	  {
	    output.logTabPrintf(0,LOGFILE,
				 "**** Centre of diffraction pattern is WRONG ****\n\n");
	    output.logTabPrintf(0,LOGFILE,
				 "The Z-score %5.1f for miscentering seems significant\n\n",
				 Zcentre);
	  }
	else
	  {
	    output.logTabPrintf(0,LOGFILE,
				 "**** Centre of diffraction pattern may be WRONG ****\n\n");
	    output.logTabPrintf(1,LOGFILE,
				 " but the Z-score %5.1f for miscentering may not be significant\n",
				 Zcentre);
	    output.logTabPrintf(0,LOGFILE,
				 "         Check program output carefully\n\n");
	  }
	output.logTabPrintf(0,LOGFILE,
			     "**** Required reindexing operator is %s\n\n\n",
			     reindex_cen.as_hkl().c_str());
      }
    output.logTab(0,LOGFILE,"\n<!--SUMMARY_END-->\n\n");
    output.logTab(0,LXML, "</Centre>");
    return reindex_cen;
  }
}
