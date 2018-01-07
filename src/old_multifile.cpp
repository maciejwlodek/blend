// Multifile.cpp
//
// Outline program for comparing multiple MTZ files for
// alternative indexing

// For now, just reads two files & compares them

// Pointlees "library" things (namespace scala)
#include "hkl_unmerge.hh"
#include "mtz_unmerge_io.hh"
#include "Output.hh"
#include "controls.hh"
#include "hkl_controls.hh"
#include "openinputfile.hh"
#include "normalise.hh"
#include "pointgroup.hh"
#include "average.hh"

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;

using namespace scala;
using phaser_io::LOGFILE;


//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
//=================================================================



int main(int argc, char* argv[])
{
  // JAF - START
  std::cout << "This is the very beginning!" << std::endl;
  // JAF -   END

  // Initialise CCP4 command line parser
  CCP4::ccp4fyp(argc, argv);                                                       // In ccp4_general.h

  // JAF - START
  std::cout << "CCP4 has been initialized" << std::endl;
  // JAF -   END

  // Initialise output object, CCP4 mode, write header
  phaser_io::Output output;                                                        // In Output.hh
  output.setPackageCCP4();

  // JAF - START
  std::cout << "PHASER stuff has been initialized" << std::endl;
  // JAF -   END

  // Read two filenames
  const int Nfiles = 2;

  std::vector<std::string> files(Nfiles);
  std::cin >> files[0] >> files[1];
  std::cout << "Files: " << files[0] << " " << files[1] << "\n";

  // JAF - START
  std::cout << "The two files have been read" << std::endl;
  // JAF -   END

  // Everything in a try ... catch block
  try {
    //  Unmerged reflection data structures
    std::vector<hkl_unmerge_list> hkl_lists(Nfiles);
    std::vector<Normalise> NormRes(Nfiles);

    // *************************************************************
    //  Setup up controls for reflection & column selection etc
    //    these should be set from command input (or default)
    // Set Profile-fitted [default] or integrated intensity
    col_controls column_selection; 
    
    // File selection flags (resolution, datasets, batches etc)
    file_select file_sel;
    
    // Scala control classes (default settings)
    //  run controls
    //  partials controls
    all_controls controls;
    
    // Make list of required columns in unmerged HKLIN file
    MtzIO::column_labels column_list = MtzIO::setup_columns();
    
    // Define ice rings: these might be omitted from analysis
    Rings Icerings;
    Icerings.AddRing(3.90, 0.0025);   // resolution in A, width of ring
    Icerings.AddRing(3.67, 0.0025);
    Icerings.AddRing(3.44, 0.0025);
    
    hkl_symmetry Symmetry;  // Space group from mtz files
    int AllowI2 = +1;

    // *************************************************************
    //   Get highest low resolution limit & lowest high resolution
    double maxHiRes = 0.0;
    double minLoRes = 100000.0;
    // loop to read in unmerged data files
    for (int kdatafile=0;kdatafile<Nfiles;kdatafile++)
      {
	const int verbose = +1;
	MtzIO::MtzUnmrgFile mtzin;
	mtzin.AddHklList(kdatafile, files[kdatafile], file_sel, 
			  column_selection, column_list, controls,
			 PxdName(), Scell(), double(2.0),
			  output, verbose, hkl_lists[kdatafile]);
	// Symmetry check
	if (kdatafile == 0)
	  {
	    // first file, just store
	    Symmetry = hkl_lists[kdatafile].symmetry();
	  }
	else
	  {
	    if (Symmetry != hkl_lists[kdatafile].symmetry())
	      {
		// Files must have same symmetry (at least rotation part)
		Message::message(Message_fatal
				 ( "Inconsistent symmetry between files\n")); 
	      }
	  }
	// sort & organise reflection list as required
	int Nrefl = hkl_lists[kdatafile].prepare();
	// Sum partials, return number of valid observations
	int Nobs = hkl_lists[kdatafile].sum_partials();
	
	output.logTabPrintf(0,LOGFILE,"\nNumber of reflections  = %10d\n",Nrefl);
	output.logTabPrintf(0,LOGFILE,    "Number of observations = %10d\n",Nobs);
	float ratio = float(Nobs)/float(Nrefl);
	output.logTabPrintf(0,LOGFILE,    "  Average multiplicity =   %8.1f\n",ratio);
	output.logTabPrintf(0,LOGFILE,"\nResolution range in list:  %9.2f ->%7.2f\n",
			    hkl_lists[kdatafile].ResRange().ResLow(), hkl_lists[kdatafile].ResRange().ResHigh());
	
	// Get normalisation obkect to convert I -> E^2
	//   and [optionally] get resolution range cutoff at <I>/<sigma> = 10
	const double MinIsigRatio = 0.; // No resolution cutoff for now
	ResoRange Rrange = hkl_lists[kdatafile].ResRange();  
	Rrange.SetRange(hkl_lists[kdatafile].num_observations());
	NormRes[kdatafile] = SetNormalise(hkl_lists[kdatafile], MinIsigRatio, false,
					  Rrange,Icerings,0);
	output.logTabPrintf(0,LOGFILE,"\nResolution range selected %8.2f to %8.2f\n",
			    Rrange.ResLow(), Rrange.ResHigh());

	minLoRes = Min(minLoRes, Rrange.ResLow());
	maxHiRes = Max(maxHiRes, Rrange.ResHigh());
      }
    // *************************************************************
    // Set overall resolution limits to worst for any file
    output.logTabPrintf(0,LOGFILE,"\nOverall resolution range reset to %8.2f to %8.2f\n\n",
			    minLoRes, maxHiRes);
    ResoRange ResRange(minLoRes, maxHiRes);
    //  set limits for all file datasets 
    for (int kdatafile=0;kdatafile<Nfiles;kdatafile++)
      {
	hkl_lists[kdatafile].SetResoLimits(minLoRes, maxHiRes);
      }
    // *************************************************************
    // Possible alternative indexing schemes, within cell tolerance

    // Allowed cell discrepancy for alternative indexing
    // Very rough: convert angular tolerance to "length" using
    // average cell edge
    double angular_tolerance = 2.0;   // degrees
    Scell cell = hkl_lists[0].cell();
    double edge = 0.333 * (cell[0] + cell[1] + cell[2]);
    double max_diff = clipper::Util::d2rad(angular_tolerance) * edge;
    
    CCtbxSym::PointGroup PG(Symmetry.symbol_xHM());
    
    PG.SetCell(cell.UnitCell(), ReindexOp(), AllowI2);
    std::vector<scala::ReindexOp> ReindexList =
      AlternativeIndexing(PG,false, cell, max_diff, AllowI2);

    int Nalt = ReindexList.size();  // number of alternative indexing schemes
    if (Nalt <= 1)
      {
	std::cout << "No alternative indexing in this space group\n";
	return 0;  // Stop
      }

    // *************************************************************
    // Loop all reflections in file 1 and fetch matching one in file 2
    reflection this_ref1, this_ref2;
    std::vector<correl_coeff> CC(Nalt);  // CC for each reindexing
    int iref;    

    while (hkl_lists[0].next_reflection(this_ref1) >= 0)
      {
	float sSqr = this_ref1.invresolsq();
	
	for (int idx=0;idx<Nalt;idx++)
	  {
	    // Reindexed reflection, put back into asymmetric unit
	    int isym;
	    Hkl hkl2 = Symmetry.put_in_asu(this_ref1.hkl().change_basis(ReindexList[idx]), isym);
	    int jref2 = hkl_lists[1].get_reflection(this_ref2, hkl2);
	    if (jref2 >= 0)   // = -1 if missing
 	      {
		// Average I, sigI
		//		IsigI Is1 = average_Is(this_ref1);
		//		IsigI Is2 = average_Is(this_ref2);
		// Average E
		IsigI E2sig1 = average_Es(this_ref1, NormRes[0]);
		IsigI E2sig2 = average_Es(this_ref2, NormRes[1]);
		// Correlation on E
		CC[idx].add(E2sig1, E2sig2);
	      }
	  }
      }
    
    for (int idx=0;idx<Nalt;idx++)
      {
	output.logTab(0,LOGFILE,
		      "\n\nReindex operator: "+ReindexList[idx].as_hkl()+"\n");
	output.logTabPrintf(0,LOGFILE,
			"Correlation coefficient for %6d pairs = %10.3f\n",
			    CC[idx].Number(), CC[idx].result().val);
      }   
    std::cout << "\n\n********\n********\n********\n";

  }
  
  catch (Message_fatal& message)
    {
      output.logWarning(LOGFILE, "\nFATAL ERROR message: \n"
                        + message.text() + "\n");
    }
  catch (std::exception const& err) {
    output.logWarning(LOGFILE, "\nUNHANDLED EXCEPTION: " + std::string(err.what())+"\n");
  }
  catch (...) {
    output.logWarning(LOGFILE, "\nUNKNOWN EXCEPTION TYPE\n");
  }
}

