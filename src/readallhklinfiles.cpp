// readallhklinfiles.cpp

#include "readallhklinfiles.hh"
#include "scala_util.hh"
#include "printthings.hh"
#include "makehkl_listscompatible.hh"
#include "xds_unmerge.hh"
#include "sca_unmerge.hh"
#include "testlauegroup.hh"
#include "spacegroupreindex.hh"
#include "mtz_merge_io.hh"
#include "checklatticecentering.hh"
#include "checknullbatches.hh"
#include "string_util.hh"

using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

namespace scala {
  //--------------------------------------------------------------
  FileRead AddHKLIN(const int& fileSeries,
		    const std::string& hklin_filename,
		    const ReflectionFileType& hklin_filetype,
		    const bool& Unmerged,
		    file_select& file_sel,
		    col_controls& column_selection,
		    MtzIO::column_labels& column_list,
		    const all_controls& controls,
		    const scala::PxdName& InputPxdName,
		    const Scell& cell,
		    const double& cellTolerance,
		    const double& wavelength,
		    phaser_io::Output& output,
		    const int& verbose,
		    hkl_unmerge_list& hkl_list)
  // Add this file into hkl_list, if compatible with whatever is already there
  // Returns:
  //    FileRead.Opened() true if file was successfully opened
  //    FileRead.Read() true if file was read
  //
  // fileSeries is index number for file or file-series (from 1) for batch
  //  exclusion in file_sel (MTZ files only)
  // non-null cell overrides file cell, or adds it to SCALEPACK files
  {
    std::string outputstring;
    //    int Nrej_batch = 0;
    if (hklin_filetype.IsTypeMTZ()) {
      if (Unmerged) {
	MtzIO::MtzUnmrgFile mtzin;
	FileRead fr = mtzin.AddHklList(fileSeries, hklin_filename, file_sel, 
				       column_selection, column_list, controls,
				       InputPxdName, cell, cellTolerance, outputstring,
				       verbose, hkl_list);
	output.logTab(0,LOGFILE,outputstring);
	return fr;
      } else {
	// Merged file
	MtzIO::MtzMrgFile mtzmrgin;
	mtzmrgin.open_read(hklin_filename);
	FileRead fr = mtzmrgin.MakeHklList(hklin_filename, file_sel, column_selection,
					   column_list, InputPxdName, cell, 
					   outputstring, verbose, hkl_list);
	output.logTab(0,LOGFILE,outputstring);
	return fr;
      }
    }
    else if (hklin_filetype.IsTypeXDS()) {
      XDSIO::XDSunmergeFile xdsin(hklin_filename, verbose, outputstring);
      output.logTab(0,LOGFILE,outputstring);
      outputstring = "";
      xdsin.AddHklList(file_sel, controls,
		       InputPxdName, cell, fileSeries, outputstring, verbose, hkl_list);
      output.logTab(0,LOGFILE,outputstring);
      return FileRead(true, true, true, 0);
    }
    else if (hklin_filetype.IsTypeSSS()) {
      SCAIO::SCAunmergeFile scain(hklin_filename, outputstring);
      output.logTab(0,LOGFILE,outputstring);
      outputstring = "";
      FileRead fr = scain.AddHklList(file_sel, controls,
		  InputPxdName, cell, wavelength, outputstring, hkl_list);
      output.logTab(0,LOGFILE,outputstring);
      return fr;
    }
    return FileRead(false, false, false, 0);
  }
  //--------------------------------------------------------------
  int ExcludeBatches(hkl_unmerge_list& hkl_list,
		      file_select& file_sel)
  // Remove excluded batches from hkl_list
  // Excluded batches/ranges are specified in object file_sel
  // Here we are only concerned with those exclusions in file_sel
  // which refer to the possibly renumbered batches, ie those which are
  // flagged in file_sel with fileSeriesList == 0
  //
  // Returns number of observation parts rejected
  {
    int fileSeries = 0;
    BatchSelection BatSel = file_sel.BatchExclude();
    int Nbat = hkl_list.num_batches();
    // Loop all batches in data object
    for (int ib=0;ib<Nbat;ib++) {
      Batch batch = hkl_list.batch(ib);
      // Is this batch is the exclusion selection?
      if (BatSel.InSelection(batch.num(), fileSeries)) {
	// yes, flag as rejected
	hkl_list.RejectBatchSerial(ib);
	//^	std::cout << "Batch Rejected: " << batch.num() << "\n"; //^
      }
    }
    //
    int nrefrej = hkl_list.PurgeRejectedBatches();
    return nrefrej;
  }
  //--------------------------------------------------------------
  std::vector<int> MakeRunOffsets(hkl_unmerge_list& hkl_list1,
				  hkl_unmerge_list& hkl_list_new,
				  MakeHKL_listscompatible& makecompatible,
				  const bool& verbose,
				  phaser_io::Output& output)

  // Make batch numbers unique relative to first file (unmerged)
  // and append hkl_list_new to hkl_list1
  {
    //^
    //^    std::cout << "**** MakeRunOffsets\n";
    std::vector<int> runOffsets;
    makecompatible.MakeBatchNumbersUnique(hkl_list1, hkl_list_new);
    if (makecompatible.IsOffset()) {
      runOffsets = makecompatible.RunOffsets();
      if (runOffsets.size() > 1) {
	output.logTabPrintf(1,LOGFILE,
			    "Batch numbers for each run incremented by ");
	for (size_t i=0;i<runOffsets.size();i++) {
	  output.logTabPrintf(0,LOGFILE, " %6d", runOffsets[i]);
	}
	output.logTabPrintf(0,LOGFILE,"\n");
      }
      else {
	output.logTabPrintf(1,LOGFILE,
			    "Batch numbers incremented by %6d\n", runOffsets[0]);
      }
    }
    //	hkl_list_next.prepare();
    //	hkl_list_next.sum_partials();
    int iverbose = (verbose ? 1 : 0);
    PrintUnmergedHeaderStuff(hkl_list_new, output, iverbose);
    // append new list to old
    hkl_list1.append(hkl_list_new);
    return runOffsets;
  }
  //--------------------------------------------------------------
  void PrintIndexingResult(const std::vector<ReindexScore>& BestReindexScores,
			   const bool& AltTest,
			   const int& fileSeries,
			   const std::string& sseries,
			   const GlobalControls& GC,
			   const ReindexOp& reindex,
			   const bool& newSG, const std::string& NewSGname,
			   const RefListType& hklrefListType,		      
			   phaser_io::Output& output)
  // Print indexing result  
  {
    std::string sf = "files";
    if (sseries != "") sf = "file series";
    if (AltTest) {
      bool printscore = true;
      if (hklrefListType == NONE) {
	if (fileSeries == 1) {
	  printscore = false;
	  if (!reindex.IsIdentity()) {
	    output.logTabPrintf(0,LOGFILE,
				"\nData reindexed by input operator %s\n",
				reindex.as_hkl().c_str());
	  }		
	  if (newSG) {
	    output.logTab(0,LOGFILE,
			  "Space group changed to "+NewSGname);
	  }
	} else if (fileSeries == 2) {
	  output.logTab(0,LOGFILE,
			"\nAlternative index test relative to first file"+sseries);
	} else {
	  output.logTab(0,LOGFILE,
			"\nAlternative index test relative to "+sf+" so far");
	}
      } else {
	output.logTab(0,LOGFILE,
		      "\nAlternative index test relative to reference file");
      }
      if (printscore)
	PrintIndexScores(BestReindexScores, false, output);
    } else {
      // No alternative test (AltTest false)
      if (fileSeries > 1 && GC.AssumeSameIndexing()) {
	output.logTab(0,LOGFILE,
		      "\nAssuming all "+sf+" have the same indexing");
      } else if ((hklrefListType != NONE) ||
		 (fileSeries > 1 && !GC.AssumeSameIndexing())) {
	output.logTab(0,LOGFILE,
		      "\nNo possible alternative indexing");
      }
    }
  }
  //--------------------------------------------------------------
  MakeHKL_listscompatible CheckAgainstReference(const all_controls& controls,
						const GlobalControls& GC,
						hkl_merged_list& RefList,
						hkl_unmerge_list& RefUnmrgList,
						const RefListType& hklrefListType,
						phaser_io::Output& output,
						hkl_unmerge_list& hkl_list)
  {
    MakeHKL_listscompatible makecompatible;
    if (hklrefListType == MERGED) {
      // Reference dataset present and merged, check for consistent indexing
      makecompatible.init(RefList, hkl_list, GC,
			  controls, output);
    }
    else if (hklrefListType == UNMERGED) {
      // Reference dataset present and unmerged, check for consistent indexing
      bool CheckIndex = true;
      makecompatible.init(RefUnmrgList, hkl_list, GC, CheckIndex,
			  controls, output);
    }
    return makecompatible;
  }
  //--------------------------------------------------------------
  bool ReindexList(const GlobalControls& GC,
		   phaser_io::Output& output,
		   hkl_unmerge_list& hkl_list,
		   ReindexOp& reindex,
		   bool& newSG,
		   std::string& NewSGname)
  // Return true is list has been reindexed
  {	
    // No reference and first file, do we have an input reindex command?
    if (GC.IsReindexSet() ||
	(GC.Spacegroup() != "" && GC.Spacegroup() != "HKLIN")) {
      // REINDEX or SPACEGROUP given on input
      if (GC.IsReindexSet()) {
	reindex = GC.Reindex();
      }
      // If SPACEGROUP is given but REINDEX is not, generate reindex operator
      // within a lattice group
      // Note that this is probably really only useful (or indeed valid) for
      // C2 <-> I2 & R3 <-> H3
      bool reindexSet = SpacegroupReindex(GC,
					  hkl_list.symmetry(),
					  hkl_list.cell(), reindex, output);
      reindexSet = reindexSet;
      hkl_symmetry NewSymm;
      if (GC.Spacegroup() != "" && GC.Spacegroup() != "HKLIN") {
	NewSymm = hkl_symmetry(GC.Spacegroup());
	NewSGname = NewSymm.symbol_xHM();
	newSG = true;
      } else {
	NewSymm = hkl_list.symmetry();
	if (GC.IsReindexSet()) {
	  NewSymm.ChangeBasis(reindex);
	  NewSGname = NewSymm.symbol_xHM();
	}
      }
      hkl_list.change_symmetry(NewSymm, reindex);
      return true;
    }
    return false;
  }
  //--------------------------------------------------------------
  void ResetLatticeCentering(const bool& symSet, hkl_unmerge_list& hkl_list,
  			     phaser_io::Output& output)
  // Check if the lattice centering in the list matches the presumed symmetry
  // If not, reset symmetry to simplest Bravais lattice group
  {
    char latType = CheckLatticeCentering(hkl_list);

    // Allow primitive R setting
    if (symSet && (hkl_list.symmetry().lattice_type() == 'R' && latType == 'P')) {
      output.logTab(0,LOGFILE,"\nRhombohedral R lattice");
      latType = hkl_list.symmetry().lattice_type();
    } else if (latType != hkl_list.symmetry().lattice_type()) {
      // actual lattice type is not the same as that indicated in the symmetry
      // This is OK if the hklin files were ShelX or Saint files for which
      // the space group was set as P1
      if (symSet) {
	// Symmetry has been set from at least one input hklin file
	// Print warning
	Message::message
	  (Message_warn(std::string("\nWARNING: Lattice type ")+latType+
			" deduced from reflection data is different from that in the input symmetry "+
			hkl_list.symmetry().lattice_type()+"\n"));
      }
      // Reset symmetry to <latType> xxx
      hkl_symmetry lSym;
      std::string group;
      if (symSet) {
	scala::SpaceGroup SG = hkl_list.symmetry().GetSpaceGroup().NewLatticePointGroup(latType);
	if (SG.LatType() == latType) {
	  // New symmetry including point group
	  lSym = hkl_symmetry(SG);
	  group = lSym.symbol_xHM();
	} 
      }
      if (group == "") {
	// failed to change or no symSet, use minimal group
	group = MinimalLatticeGroup(latType);
	lSym = hkl_symmetry(group);
      }
      hkl_list.change_symmetry(lSym, ReindexOp());
      output.logTab(0,LOGFILE,
		    std::string("\nLattice centering type changed to ")+
		    latType+", space group "+group);
      hkl_list.prepare();
      hkl_list.sum_partials();
    }
  }
  //--------------------------------------------------------------
  void TestFirstFileSymmetry(const all_controls& controls,
			     const GlobalControls& GC,
			     const std::string& sseries,
			     const RefListType& hklinListType,
			     phaser_io::Output& output,
			     hkl_unmerge_list& hkl_list)
  // Test for symmetry in the first file in the usual way
  {
    output.logTabPrintf(0,LOGFILE,
			"\n********************************************************\n");
    output.logTabPrintf(0,LOGFILE,
			"\n    Determining Laue group for first file%s\n",
			sseries.c_str());
    std::vector<PGscore> SGscores;
    ValidElementObserved ElementsObserved =
      TestLaueGroup(hkl_list, GC, controls, hklinListType, output, SGscores);
    PrintSubGroupScores(SGscores, hklinListType, 0.0, output);
    // Choose best group if enough information & apply it
    if (ElementsObserved.Valid()) {
      ReindexOp reindex =
	hkl_list.TotalReindex().inverse() * SGscores[0].RefSGreindex();
      // Symmetry
      hkl_symmetry Symm(CCtbxSym::PointGroupName(SGscores[0].RefPGname()));
      output.logTabPrintf(0,LOGFILE,
			  "\nReindexing test data from first file%s by %s\n  in point group %s\n\n",
			  sseries.c_str(),
			  SGscores[0].RefSGreindex().as_hkl().c_str(),
			  Symm.symbol_xHM().c_str());
      

      hkl_list.change_symmetry(Symm, reindex, false);
    } else {
      output.logTab(0,LOGFILE,
		    "\nInsufficient data to determine Laue group of first input file\n");
    }
  }
  //--------------------------------------------------------------
  void RecordReindexing(const MakeHKL_listscompatible& makecompatible,
			IO_files& AllFiles,
			const int& fSeries, 
			const std::string& sseries,
			const bool& newSG,
			const std::string& NewSGname,
			const RefListType& hklrefListType,		      
			const GlobalControls& GC,
			const bool& verbose,
			phaser_io::Output& output)
  {
    ReindexOp reindex;
    bool AltTest = makecompatible.AltIndex();
    if (AltTest) {
      // Store reindex for all files in series
      reindex = makecompatible.BestReindexScore();
      AllFiles.StoreHKLINreindex(fSeries, reindex,
				 makecompatible.BestReindexScore().OverallCC(),
				 makecompatible.BestReindexScore().Likelihood(),
				 makecompatible.Confidence());
    }
    if (verbose > 0 ) {
      PrintIndexingResult(makecompatible.BestReindexScores(),
			  AltTest, fSeries, sseries, GC, reindex,
			  newSG, NewSGname,hklrefListType, output);
    }
  }
  //--------------------------------------------------------------
  int ReadAllHKLINfiles(IO_files& AllFiles,
			file_select& file_sel, 
			col_controls& column_selection,
			MtzIO::column_labels& column_list,
			const all_controls& controls,
			const GlobalControls& GC,
			const Scell& cell,
			const double& wavelength,
			hkl_merged_list& RefList,
			hkl_unmerge_list& RefUnmrgList,
			const RefListType& hklrefListType,		      
			const RefListType& hklinListType,		      
			phaser_io::Output& output,
			hkl_unmerge_list& hkl_list)
  // Read all HKL files into hkl_list
  // Returns status:
  //   =  0  OK
  //   >  0  error = +2 cell discrepancy
  //   = -1  no symmetry set from file (dummy value set in hkl_list)
  {
    int verbose = +3;   // verbosity of output for each file
    //    int verbose = +2;   // verbosity of output for each file
			// +1 minimal if file series defined
			// +2 medium if multiple files
			// +3 more verbose if only one file
    if (AllFiles.NumberOfHKLINnames() == 1) {verbose = +3;}
    if (AllFiles.AreThereFileSeries()) {
      verbose = +1;
      if (AllFiles.NumberOfFileSeries() == 1) {
	int nfiles = AllFiles.NumberOfFilesInSeries(1);
	output.logTabPrintf(0,LOGFILE,
			    "Reading %4d files from file series: %s\n",
			    nfiles, AllFiles.SeriesNames()[0].c_str());
      } else {
	output.logTab(0,LOGFILE,"Reading file series:");
	for (size_t i=0;i<AllFiles.SeriesNames().size();++i) {
	  int nfiles = AllFiles.NumberOfFilesInSeries(1);
	  output.logTabPrintf(1,LOGFILE,
			      "%4d files from series %s\n",
			      nfiles, AllFiles.SeriesNames()[i].c_str());
	}
      }
    }
    int status = 0;
    PxdName InputPxdName;
    int fileSeries;
    int fSeries;
    int Nfiles = AllFiles.NumberOfHKLINnames();
    bool AltTest = false;
    bool newSG = false;
    std::string NewSGname;
    bool firstSeries = true;
    ReindexOp reindex;
    int Nrej_batch = 0;
    bool CheckIndex = !GC.AssumeSameIndexing();
    MakeHKL_listscompatible makecompatible;
    std::vector<int> runOffsets;  // Empty if no non-zero offsets
    double cellTolerance = GC.LatticeTolerance();

    hkl_unmerge_list hkl_list_work;  // temporary store for new (2nd+) series
    hkl_unmerge_list* HklList = &hkl_list;  // input list, == _work for 2nd+ series
    ReflectionFileType hklin_filetype;

    // for printing, depending on whether there are any HKLIN file series or not
    std::string sseries = "";
    if (AllFiles.AreThereFileSeries()) {sseries = " series";}
    //    int Nseries = AllFiles.NumberOfFileSeries();
    FileRead fileRead(true, true, true,0);  // to force reading of 1st file
    int kfile = -1;
    std::string hklin_filename;
    std::vector<std::string> fileNames;
    bool Unmerged;
    bool symSet = false;

    // + + + + + + + + + + + + + + + + + + + + + + + + + + + + 
    // Loop reading files
    while (++kfile < Nfiles){
      AltTest = false;
      if (fileRead.Read()) {
	// File name, PXD name, file type & file series number for file number kfile
	// (NB fileSeries is numbered from 1)
	hklin_filename =
	  AllFiles.HKLINfilename(kfile, hklin_filetype, InputPxdName, fileSeries);
	Unmerged = AllFiles.MergedFlag(kfile);
	fSeries = fileSeries;
	if (Nfiles == 1) {fSeries = 0;}
      }
      // Try to read or re-read file
      fileRead = AddHKLIN(fSeries, hklin_filename, hklin_filetype, Unmerged,
			  file_sel, column_selection, column_list,
			  controls, InputPxdName, cell, cellTolerance,
			  wavelength, output, verbose, (*HklList));
      symSet = symSet || fileRead.SymmetrySet();  // true if symmetry set from any input file

      if (!fileRead.Opened()) {
	// Can't open file, try to keep going if in series
	if (fSeries > 0 && AllFiles.AreThereFileSeries()) {
	  output.logTab(0,LOGFILE,
			"**** Failed to open MTZ file "+hklin_filename);
	  // Remove this file from the list
	  AllFiles.RemoveHKLINfile(kfile);
	  Nfiles--;
	  kfile--;
	  fileRead = FileRead(true, true, true, 0);
	  continue;  // try next one
	} else {
	  Message::message
	    (Message_fatal("Failed to read file "+hklin_filename));
	}
      }

      if (fileRead.Read()) {
	// Store symmetry from input files
	if (fileRead.SymmetrySet()) {
	  AllFiles.StoreHKLINsymmetry(kfile, (*HklList).symmetry());
	}
	AllFiles.StoreHKLINcell(kfile, (*HklList).cell());
	Nrej_batch += fileRead.Nexcluded();
	std::string hklstream = "HKLIN";
	if (AllFiles.AreThereFileSeries()) {
	  hklstream = StringUtil::Strip("HKLIN"+clipper::String(fileSeries));
	}
	PrintFileInfoToXML(hklstream, hklin_filename,
			   (*HklList).cell(),
			   (*HklList).symmetry().symbol_xHM(),
			   output);
	fileNames.push_back(hklin_filename); // accumulate filenames
	// Current file has been added to list
	// Should we try to add next file?
	//  yes, if AssumeSameIndexing & not on last file
	if (!CheckIndex && (kfile < Nfiles)) {
	  continue;  // read next file
	} 
      } else {
	// File opened but not read (ie incompatible), re-read into next list
	kfile--;
	fileRead = FileRead(true, true, true, 0); // force next read
      }
      // Don't add next file to current list, switch lists
      if (HklList == &hkl_list) {
	// Main list full, check for
	//  (1) Reference file
	//  (2) Explicit reindex operation
	//  (3) first file & TestFirstFile
	// if so, deal with main list before going further
	if (hklrefListType != NONE) {
	  makecompatible = CheckAgainstReference(controls, GC, RefList, RefUnmrgList,
						 hklrefListType, output, (*HklList));
	  AltTest = makecompatible.AltIndex();
	  RecordReindexing(makecompatible, AllFiles, fileSeries,  sseries,
			   newSG,  NewSGname, hklrefListType, GC, verbose, output);
	  status += makecompatible.Status();
	  if (status > 0) {
	    // Incompatible unit cells, bail out
	    return status;
	  }
	} else if (firstSeries) {
	  // No reference and first file, do we have an input reindex command?
	  //  if so, reindex
	  AltTest = ReindexList(GC, output, (*HklList), reindex,
				newSG, NewSGname);
	  if (!AltTest) {
	    if (AllFiles.NumberOfFileSeries() > 1 &&
		(!symSet || GC.TestFirstFile())) {
	      // No reindex
	      // Optionally determine Laue group for first fileseries if more than one,
	      // since this will be used as a reference
	      // Test if either no symmetry in first file or asked for on input
	      if (!symSet) {  // reset lattice if needed
		hkl_list.prepare();
		hkl_list.sum_partials();
		ResetLatticeCentering(symSet, hkl_list, output);
	      }
	      TestFirstFileSymmetry(controls, GC, sseries, hklinListType,
				    output, (*HklList));
	    }
	  }
	}
	// finished with 1st (main) list
	// We have been filling the main list, switch to 2nd (work) list
	HklList = &hkl_list_work;
	// Read next file or previous incompatible file
	continue;
      } else {  // work list full
	// We have been filling the work list, now we have finished it
	// Is there anything in it?
	if (HklList->num_parts() > 0) {
	  // Merge work list into main list
	  //  First make it compatible with main list or reference
	  if (hklrefListType != NONE) {
	    makecompatible = CheckAgainstReference(controls, GC, RefList, RefUnmrgList,
						   hklrefListType, output, (*HklList));
	    AltTest = makecompatible.AltIndex();
	  } else {
	    // No reference and not the first file series
	    //    check for consistent indexing if CheckIndex set
	    makecompatible.init(hkl_list, (*HklList), GC, CheckIndex,
				controls, output);
	    if (CheckIndex && makecompatible.AltIndex()) AltTest = true;
	  }
	  
	  status += makecompatible.Status();
	  if (status > 0) {
	    // Incompatible unit cells, bail out
	    return status;
	  }
	  // Apply any necessary batch offsets, and append work list to main list
	  runOffsets = MakeRunOffsets(hkl_list, *(HklList),
				      makecompatible, verbose, output);
	  RecordReindexing(makecompatible, AllFiles, fileSeries,  sseries,
			   newSG,  NewSGname, hklrefListType, GC, verbose, output);
	}
	// Clear work list
	hkl_list_work.clear();
      }
    }  // End file read loop

    if (verbose == +1) {
      // Now list all file names read
      output.logTab(0,LOGFILE,"\nList of files read:");
      for (size_t i=0;i<fileNames.size();++i) {
	output.logTab(1,LOGFILE,fileNames[i]);
      }
      output.logTabPrintf(0,LOGFILE,"\n");
    }

    if (Nfiles > 1) {
      // Remove excluded batches from hkl_list
      // if they haven't already been removed on reading
      // ie those given on the EXCLUDE command with any FILE qualifier
      Nrej_batch += ExcludeBatches(hkl_list, file_sel);
    }
    if (Nrej_batch > 0) {
      output.logTabPrintf(1,LOGFILE,
			  "\nNumber of observation parts excluded: %9d\n",
			  Nrej_batch);
    }
    hkl_list.prepare();
    hkl_list.sum_partials();
    CheckNullBatches(hkl_list, file_sel, output);

    ResetLatticeCentering(symSet, hkl_list, output);

    if (!symSet) {status = -1;}
    return status;
  }
  //-------------------------------------------------------------- 
  double GetResolutionFromHKLIN(const IO_files& AllFiles)
  // Extract maximum resolution from HKLIN files
  // Only works for MTZ input: in other cases returns a default value
  {
    const double MAXRESDEFAULT = 2.5;  // default value if not available from file
    double maxres = -1.0;
    ReflectionFileType hklin_filetype;
    std::string hklin_filename;
    // Loop all input files
    for (int ifl=0;ifl<AllFiles.NumberOfHKLINnames();++ifl) {
      hklin_filename = AllFiles.HKLINfilename(ifl, hklin_filetype);
      if (hklin_filetype.IsTypeMTZ()) {
	MtzIO::MtzMrgFile mtzin;
	bool fileValid = mtzin.open_read(hklin_filename);
	if (fileValid) {
	  maxres = Max(maxres, double(mtzin.MtzResolution().limit()));
	}
      }
    }
    if (maxres < 0.0) {maxres = MAXRESDEFAULT;}
    return maxres;
  }
  //--------------------------------------------------------------
}

