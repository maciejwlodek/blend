// checknullbatches.cpp
//
// Check for and eliminate batches with essentially zero intensity
// (ie from blank images)
//

#include "checknullbatches.hh"
#include "score_datatypes.hh"
#include "scala_util.hh"
using phaser_io::LOGFILE;

namespace scala {
  //--------------------------------------------------------------
  int CheckNullBatches(hkl_unmerge_list& hkl_list,
		       const file_select& file_sel,
		       phaser_io::Output& output)
  {
    if (file_sel.NullResFrac() <= 0.0) return -1;  // no check to be made
    //    int nruns = hkl_list.num_runs();
    int nbatches = hkl_list.num_batches();
    // By batch
    std::vector<double> sumIbatch(nbatches);
    std::vector<double> sumIovSbatch(nbatches);
    std::vector<int> numobs(nbatches);
    std::vector<int> numZero(nbatches);

    std::vector<int> rejectedBatches; // list of rejected batches

    // High resolution limit 1/d^2
    float highRes = hkl_list.Srange().max();
    // Use low resolution only
    float reslimit = file_sel.NullResFrac() * highRes;

    float I, sigI;
    int ib;

    observation this_obs;
    reflection this_refl;
    hkl_list.rewind();  // Just in case

    while (hkl_list.next_reflection(this_refl) >= 0) {
      if (this_refl.invresolsq() < reslimit) {  // low resolution only
	while (this_refl.next_observation(this_obs) >= 0) {
	  I = this_obs.I();
	  sigI = this_obs.ksigI();
	  ib = hkl_list.batch_serial(this_obs.Batch());
	  sumIbatch[ib] += I;
	  if (sigI > 0.0) {
	    sumIovSbatch[ib] += I/sigI;
	  }
	  numobs[ib]++;
	  if (I <= 0.0) numZero[ib]++;
	}
      }
    } // end loop reflections

    int nrej = 0;
    const int MIN_NUM_OBS = 20; // must have at least this many observations
    double maxZf = 0.0;
    int NbatZero = 0;  // batches with no observations
    
    // Split batches into up to 4 parts, each part with at least 10 batches
    int nb = Max(10,nbatches/4);
    int nbpart = Max(1,nbatches/nb); // number of parts
    nb = nbatches/nbpart;     // number of batches/part
    std::vector<MeanSD> meanZf(nbpart);
    std::vector<float> batZf(nbatches,0.0);

    for (int ib=0;ib<nbatches;++ib) {
      if (numobs[ib] <= 0) {
	NbatZero++;;
      } else {
	int mb = Min(ib/nb, nbpart-1);
	double Zf = double(numZero[ib])/double(numobs[ib]);
	sumIbatch[ib] /=double(numobs[ib]);
	sumIovSbatch[ib] /=double(numobs[ib]);
	meanZf.at(mb).Add(Zf);
	batZf[ib] = Zf;

	//^
	//	std::cout << "Batch, <I> <I/sig>  N N0 Zf"
	//			  << hkl_list.batch(ib).num()
	//			  << " " << sumIbatch[ib] << " " << sumIovSbatch[ib]
	//			  << " " << numobs[ib] << " " << numZero[ib]
	//			  << " " << Zf << "\n";
	if (numobs[ib] > MIN_NUM_OBS) {
	  maxZf = Max(maxZf, Zf); // ignore if too few observations
	}
      }
    }
    // Pick out ~3/4 the parts with the smallest SD(Zf)
    MeanSD meanZfm;
    // nbpart  1  2  3  4
    // mpart   1  1  2  3
    int mpart = Max(Min((nbpart*3)/4, nbpart-1),1);
    std::sort(meanZf.begin(), meanZf.end(), MeanSD::MeanSDsmallerSD);
    for (int i=0;i<nbpart;++i) {
      //            std::cout << "  Mean Zf, SD(Zf) "
      //      		<< meanZf[i].Mean() << " " << meanZf[i].SD()
      //      		<< "\n"; //^
      if (i < mpart) {
	meanZfm += meanZf[i];
      }
    }
    //        std::cout << "Selected Zf, SD " << meanZfm.Mean() << " " << meanZfm.SD() <<"\n"; //^
    // Accept batches with fraction(below 0) Zf < accept
    // NullNegativeReject default = ~0.3
    float accept = Max(file_sel.NullNegativeReject(), meanZfm.Mean() + 4.*meanZfm.SD());
    std::vector<bool> acceptBatch(nbatches, false);
    
    if (NbatZero > 0 || maxZf > accept) {
      if (maxZf > accept) {
	// Reject batches with Zf > accept
	output.logTabPrintf(0,LOGFILE,
			    "\nBatches with too many observations with intensities below zero (blank batches)\n  are rejected if Fraction<0 > %6.2f (average %6.2f)\n",
			    accept, meanZfm.Mean());
	output.logTabPrintf(1,LOGFILE,
			    "(BLANK  RESOFRACTION %7.2f NEGATIVE %7.2f)\n\n",
			    file_sel.NullResFrac(), accept);

	for (int ib=0;ib<nbatches;++ib) {
	  bool OK = true;
	  if (numobs[ib] <= 0) {
	    output.logTabPrintf(2,LOGFILE,
				"Batch %5d eliminated as it has no observations\n",
				hkl_list.batch(ib).num());
	    OK = false;
	  } else if (numobs[ib] > MIN_NUM_OBS && batZf[ib] > accept) {
	    output.logTabPrintf(2,LOGFILE,
				"Batch %5d eliminated, Fraction<0 = %6.2f\n",
				hkl_list.batch(ib).num(), batZf[ib]);
	    OK = false;
	  }
	  if (!OK) {
	    rejectedBatches.push_back(hkl_list.batch(ib).num());
	    nrej++;
	  } else {
	    acceptBatch[ib] = true;
	  }
	}
      }
      if (nrej > 0) {
	// Yes we have done some rejections
	// Calculate mean(I) & SD(I) for accepted batches
	MeanSD mnI;
	for (int ib=0;ib<nbatches;++ib) {
	  if (acceptBatch[ib]) {mnI.Add(sumIbatch[ib]);}
	}
	double Imean =  mnI.Mean();
	double acceptI =  mnI.Mean() - 4.*mnI.SD();
	//^      std::cout << "Recalculated mean, SD, threshold "
	//		<< Imean << " " << mnI.SD() << " " << acceptI << "\n";
	bool extraRejected = false;

	do {
	  bool rejected = false;
	  extraRejected = false;
	  // Check last batch before a rejected one to see if that should be rejected too
	  for (int ib=nbatches-1;ib>=0;--ib) {
	    if (acceptBatch[ib]) {
	      // this batch accepted
	      if (rejected) {
		// Last accepted batch before a rejected one, test it
		if (sumIbatch[ib] < acceptI) {
		  output.logTabPrintf(2,LOGFILE,
				      "Batch %5d eliminated, last batch before rejected batches\n          too weak <I> =%6.0f threshold %6.0f overallAverage<I> %7.0f\n",
				      hkl_list.batch(ib).num(), sumIbatch[ib], acceptI, Imean);
		  rejectedBatches.push_back(hkl_list.batch(ib).num());
		  nrej++;
		  extraRejected = true;
		  acceptBatch[ib] = false;
		}
	      }
	      rejected = false;
	    } else {
	      // This batch rejected
	      rejected = true;
	    }
	  }
	} while (extraRejected);

	int nrefrej = hkl_list.EliminateBatches(rejectedBatches);
	nrefrej = nrefrej;
	hkl_list.AutoSetRun();
	hkl_list.prepare();
	hkl_list.sum_partials();
      }
      return nrej;
    }
    return nrej;
  } // CheckNullBatches
}

