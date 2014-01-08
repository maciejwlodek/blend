// sca_unmerge.cpp
//
// Read SCALEPACK data into hkl_unmerge object

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include "sca_unmerge.hh"
#include "controls.hh"
#include "hkl_controls.hh"
#include "observationflags.hh"
///#include "score_datatypes.hh"
#include "util.hh"
#include "scala_util.hh"
#include "string_util.hh"
///#include "latsym.hh"
///#include "jiffy.hh"
#define ASSERT assert

using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

namespace SCAIO {
  //--------------------------------------------------------------
  // Open file & check file type
  //  scaname            name of file (or logical name)
  SCAunmergeFile::SCAunmergeFile(const std::string& scaName,
				 std::string& output)
  {
    output = "";
    sca_type = SCA_NONE;
    std::vector<double> vcell(6);
    spaceGroupName = "P 1"; // default

    // Open input stream from file
    scaname = scaName;
    scain.open(scaname.c_str());
    if (!scain) {
      Message::message(Message_fatal
		       ("SCAIO: cannot open file "+scaname));
    }
    
    std::getline(scain, line);  // first line
    line = StringUtil::BuftoLine(line);  // remove any trailing lf or cr characters
    lineFull = false;  // stays false unless a reflection was read as header
    // Merged files have first line == "    1"  (i5)
    if (line.size() == 5) {
      if (line.i() == 1) {
	sca_type = SCA_MERGED;  // merged file
	// Skip next line, then read cell & spacegroup
	std::getline(scain, line);
	std::getline(scain, line);
	for (int i=0;i<6;i++) {
	  vcell[i] = clipper::String(line.substr(i*10, 10)).f();
	}
	cell = scala::Scell(vcell);
	spaceGroupName = line.substr(60,line.length()-1);
      }
    } else if (isalpha(line[6])) {
      // SCA unmerged   7th character is a letter
      // unmerged file starts with number of symops & spacegroup
      int Nsymop = line.i();
      spaceGroupName = line.substr(5, line.length()-1);
      // Skip 2 * Nsymop lines
      for (int i=0;i<2*Nsymop;i++) {std::getline(scain, line);}
      sca_type = SCA_UNMERGED;  // unmerged file
    } else if (line.size() >= 250) {
      // Probably Saint format
	sca_type = SCA_SAINT;
	lineFull = true;        // keep this reflection record for later
    } else if (line.size() >= 28) {
      // Probably ShelX, 3i4,2f8.2
	sca_type = SCA_SHELX;
	lineFull = true;        // keep this reflection record for later
    }

    if (sca_type == SCA_NONE) {
      Message::message(Message_fatal
		       ("Unrecognised file type for file "+scaname));
    }
  }
//--------------------------------------------------------------
  bool SCAunmergeFile::CheckCells(const scala::Scell& input_cell, const scala::Scell& file_cell,
				  const std::string& label, std::string& output)
  // Returns false if there is a discrepancy, else true
  {
    if (input_cell.null() || file_cell.null()) return true; // null cells
    const double TOLERANCE = 2.0;   // angular difference limit for warning
    bool OK = true;
    if (!file_cell.equalsTol(input_cell, TOLERANCE)) {
      output += FormatOutput::logTab(0, 
		    std::string("**** WARNING: input CELL is significantly different from cell from ")+
		    label+" file");
      OK = false;
    } else {
      output += FormatOutput::logTab(0, 
		    std::string("Cell from ")+label+" file replaced by input cell");
    }
    output += FormatOutput::logTabPrintf(1,"Cell from file: ");
    for (int i=0;i<6;i++) output += FormatOutput::logTabPrintf(0,"%7.2f",file_cell[i]);
    output += FormatOutput::logTab(0,"\n");
    output += FormatOutput::logTabPrintf(1,"Input cell:     ");
    for (int i=0;i<6;i++) output += FormatOutput::logTabPrintf(0,"%7.2f",input_cell[i]);
    output += FormatOutput::logTab(0,"\n\n");
    return OK;
  }
//--------------------------------------------------------------
  scala::FileRead SCAunmergeFile::FillHklList(std::string& output,
					     scala::hkl_unmerge_list& hkl_list)
  {
    scala::Scell dcell;
    return FillHklList(output, dcell, hkl_list);
  }
  //--------------------------------------------------------------
  scala::FileRead SCAunmergeFile::FillHklList(std::string& output,
					      const scala::Scell& input_cell,
					      scala::hkl_unmerge_list& hkl_list)
  {
    // File selection flags (resolution, datasets, batches etc)
    scala::file_select file_sel;
    // Scala control classes (default settings)
    //  run controls
    //  partials controls
    scala::all_controls controls;
    scala::PxdName InputPxdName;
    double wavelength = -1.0;
    return AddHklList(file_sel, controls, InputPxdName, input_cell, wavelength, output, hkl_list);
  }
//--------------------------------------------------------------
  scala::FileRead SCAunmergeFile::AddHklList(const scala::file_select& file_sel, 
					     const scala::all_controls& controls,
					     const scala::PxdName& InputPxdName,
					     const scala::Scell& input_cell,
					     const double& input_wavelength,
					     std::string& output,
					     scala::hkl_unmerge_list& hkl_list)
  // On entry:
  //  scaname            name of SCA file (or logical name)
  //  file_sel           flags for general selection
  //                     - dataset selection
  //                     - batch selection
  //                     - resolution limits
  //                     - detector coordinate rejection ranges
  //                     - input scale factor (MULTIPLY)
  //  controls           run controls, partial controls
  //  InputPxdName       PXD name (none in SCA file)
  //  output             output object for printing
  // 
  // On exit:
  //  hkl_list  has been filled, but not organised:
  //  this needs a call to "prepare" or "change_symmetry"
  {
    output = "";
    inputcell = input_cell;
    inputwavelength = input_wavelength;
    if (sca_type == SCA_UNMERGED || sca_type == SCA_SHELX) {
      if (input_cell.null()) {
	Message::message(Message_fatal
			 ("You must provide the CELL for unmerged SCALEPACK or SHELX data"));
      }
      cell = input_cell;
      wavelength = input_wavelength; // = -1 if not set
    } else if (sca_type == SCA_MERGED) {
      // Scalepack merged files have a cell, but you can override it here
      // print warning if they are too different
      if (!CheckCells(input_cell, cell, "Scalepack Merged", output)) {
	cell = input_cell;
      }
    wavelength = input_wavelength; // = -1 if not set
    } else if (sca_type == SCA_SAINT) {
    }

    output += FormatOutput::logTab(0,
		  "\nReading ascii file from file "+scaname+"\n\n");

    // Initialise hkl list
    // Cook up some sort of title
    std::string Title;
    if (sca_type == SCA_MERGED) {
      Title = "From merged SCALEPACK file "+scaname;
      output += FormatOutput::logTab(1,
		    "File type is Merged SCALEPACK\n");
    } else if (sca_type == SCA_UNMERGED) {
      Title = "From unmerged SCALEPACK file "+scaname;
      output += FormatOutput::logTab(1,
		    "File type is Unmerged SCALEPACK\n");
    } else if (sca_type == SCA_SHELX) {
      Title = "From ShelX file "+scaname;
      output += FormatOutput::logTab(1,
		    "File type is SHELX\n");
    } else if (sca_type == SCA_SAINT) {
      Title = "From SAINT file "+scaname;
      output += FormatOutput::logTab(1,
		    "File type is SAINT\n");      
    }

    // reserve space for some plausible number of observations, not critical
    int NreflReserve = 10000;
    
    // ----  Space group
    // Is this a rhombohedral spacegroup?
    //  if so, is this a rhombohedral setting?
    char RLatTyp = 'H';
    // Might be rhombohedral, false unless angles are 90,90,120
    if (scala::RhombohedralAxes(cell.UnitCell())) {RLatTyp = 'R';}
    // get extended Hermann-Mauguin symbol
    // RLatTyp does not matter except for rhombohedral cells
    spaceGroupName = StringUtil::ToUpper(StringUtil::Strip(spaceGroupName));
    output += FormatOutput::logTab(0, "Spacegroup from file: "+spaceGroupName);

    hkl_list.init(Title, NreflReserve,
		  scala::hkl_symmetry(spaceGroupName), controls);
    
    int NobsRej = 0;
    if (sca_type == SCA_SAINT) {
      // SAINT input
      ReadSaint(scain, file_sel, controls,
		output, hkl_list);
      // Make dataset information (just one dataset)
      int setid = MakeDataset(InputPxdName);
      setid = setid;
    } else {
      // Scalepack or ShelX
      NobsRej = ReadObservations(scain, file_sel, controls, output, hkl_list);
      // Store list of columns present in file, just the compulsory ones
      hkl_list.StoreDataFlags(scala::data_flags());
      // Make dataset information (just one dataset)
      int setid = MakeDataset(InputPxdName);
      setid = setid;
      // Make batch information, only for batches containing data
      int dataset_index = 0;  // only one dataset
      int nbatches_used = MakeBatches(hkl_list.symmetry(), dataset_index);
      nbatches_used= nbatches_used;
    }
    // All types
    // Sort batches into order on batch number
    std::sort(batches.begin(), batches.end());
    // Put dataset & batch information into hkl_list
    hkl_list.StoreDatasetBatch(datasets, batches);
    hkl_list.AppendFileName(scaname);
    bool symSet = true;
    if (sca_type == SCA_SHELX || sca_type == SCA_SAINT) {
      symSet = false;} // no symmetry from ShelX file
    return scala::FileRead(true, true, symSet, NobsRej);
  }
  //--------------------------------------------------------------
  int SCAunmergeFile::ReadObservations(std::ifstream& scain,
				       const scala::file_select& file_sel, 
				       const scala::all_controls& controls,
				       std::string& output,
				       scala::hkl_unmerge_list& hkl_list)
  // Read all observations into list
  // sca_type  file type

  // Returns number of observations
  {
    // Stores for observation part
    scala::Hkl hkl;
    int isym;
    int batch=1;
    Rtype I;  Rtype sigI;
    Rtype Im = 0.0;  Rtype sigIm = 0.0;
    Rtype Ipr;  Rtype sigIpr;
    Rtype Xdet=0.0;  Rtype Ydet=0.0;
    Rtype phi=0.0;  Rtype time=0.0;
    Rtype fraction_calc=0.0;  Rtype width=0.0;
    Rtype LP=0.0; 
    int Npart=1;  int Ipart=1;
    scala::ObservationFlag ObsFlag;

    scala::ResoRange rrange;
    rrange.clear();

    const int TABLE_SIZE = 5000;
    lookup.set_size(TABLE_SIZE);
    int table_size = lookup.Size();  // smallest prime >= set size
    Nbatches = 0;
    int nreso_rej = 0;

    // read and parse input lines
    while (true) {
      if (!lineFull) {
	std::getline(scain, line );
	if (scain.eof()) break;  // end of file
      }
      lineFull = false;  // true only for 1st line if it was read as header
      line = StringUtil::BuftoLine(line);  // remove any trailing lf or cr characters
      // Indices
      int h = std::atoi(line.substr(0, 4).c_str());
      int k = std::atoi(line.substr(4, 4).c_str());
      int l = std::atoi(line.substr(8, 4).c_str());
      if (h == 0 && k == 0 && l == 0) continue;
      // Reduce to asymmetric unit
      hkl = hkl_list.symmetry().put_in_asu(scala::Hkl(h,k,l), isym);
      Dtype sSqr = hkl.invresolsq(cell);
      if (!file_sel.in_reslimits(sSqr)) {
	nreso_rej++;
	continue;
      }
      rrange.update(sSqr);
      bool IsIm = false;

      if (sca_type == SCA_MERGED || sca_type == SCA_SHELX) {
	// merged or ShelX
	// I & sigmaI
	I    = std::atof(line.substr(13,8).c_str());
	sigI = std::atof(line.substr(21,8).c_str());
	// I- if present
	if (line.length() > 28 && sca_type == SCA_MERGED) {
	  Im    = std::atof(line.substr(29,8).c_str());
	  sigIm = std::atof(line.substr(37,8).c_str());
	  IsIm = true;
	}
      } else {
	// unmerged
	// I & sigmaI
	I    = std::atof(line.substr(37,8).c_str());
	sigI = std::atof(line.substr(45,8).c_str());
	// Use image number as batch number, numbered from 1
	batch =  std::atoi(line.substr(24, 6).c_str());
      }
      // do we have this batch already?
      int ib = lookup.lookup(batch);
      if (ib < 0) {
	// New batch, add to list
	Nbatches++;
	if (Nbatches > 0.8*table_size) {
	  Message::message(Message_fatal
			   ("SCAIO: increase TABLE_SIZE"));
	}
	lookup.add(batch, Nbatches-1);
      }
      // Apply input scale (MULTIPLY)
      I *= file_sel.InputScale();
      sigI *= file_sel.InputScale();
      //   store as IPR as well
      Ipr = I;
      sigIpr = sigI;
      // Fake phi = time = batch
      phi = batch;
      time = phi;
      
      // width = 0, Npart = Ipart = 1, ObsFlag = 0
      hkl_list.store_part(hkl, isym, batch, I, sigI, Ipr, sigIpr,
			  Xdet, Ydet, phi, time,
			  fraction_calc, width, LP, 
			  Npart, Ipart, ObsFlag);
      
      if (IsIm) {
	// I- if present
	hkl = hkl_list.symmetry().put_in_asu(-hkl, isym);
	I = Im;
	sigI = sigIm;
	Ipr = I;
	sigIpr = sigI;
	hkl_list.store_part(hkl, isym, batch, I, sigI, Ipr, sigIpr,
			    Xdet, Ydet, phi, time,
			    fraction_calc, width, LP, 
			    Npart, Ipart, ObsFlag);
      }
    }  // end loop reflection read
    rrange.Set(); // resolution range is set
    int nobs = hkl_list.close_part_list(rrange, false);
    
    output += FormatOutput::logTabPrintf(0,
			"\n%8d observations accepted\n",
			nobs);
    output += FormatOutput::logTabPrintf(0,
			"         Resolution range %8.3f %8.3f\n",
			rrange.ResLow(), rrange.ResHigh());
    if (file_sel.reslimits().isSet())  {
      output += FormatOutput::logTabPrintf(0,
			  "%8d observations rejected because outside resolution limits (%5.1f to %4.1fA)\n",
			  nreso_rej, file_sel.reslimits().ResLow(),
			  file_sel.reslimits().ResHigh());
    }
    return nreso_rej;
  }
//--------------------------------------------------------------
  int SCAunmergeFile::MakeDataset(const scala::PxdName& InputPxdName)
  // Enter the single dataset into list
  // Returns dataset ID number
  {
    // For now, make dummy Project, Crystal, Dataset names
    // These should be taken from user input
    int setid = 1; // just one dataset, numbered 1
    // Default names
    scala::PxdName pxd("SCAproject", "SCAcrystal", "SCAdataset");
    //  incorporate any non-blank elements
    pxd.update(InputPxdName);
    datasets.clear();
    datasets.push_back(scala::Xdataset(pxd, cell, wavelength, setid));
    return setid;
  }
  //--------------------------------------------------------------
  int SCAunmergeFile::MakeBatches(const scala::hkl_symmetry& Symm,
				  const int& Idataset)
  // Make list of batch data for all batches with observations
  // Returns number actually used
  {
    CMtz::MTZBAT BHeader;   // header structure
    // Store all common data, the bits that are the same for all batches
    // Items marked "//+" here are individual to each batch, but are listed
    // here for completeness
    //+     BHeader.int num;        /**< batch number */
    strcpy(BHeader.title, "SCALEPACK data");   /**< batch title */	      
    strcpy(BHeader.gonlab[0], "        ");    /**< names of the three axes */
    strcpy(BHeader.gonlab[1], "        ");    // 8-characters only! 
    strcpy(BHeader.gonlab[2], "        ");
    BHeader.iortyp = 0;                  /**< type of orientation block (for 
					    possible future use, now = 0) */
    /**< refinement flags for cell (constraint flags)*/
    for (int i=0;i<6;i++) {BHeader.lbcell[i] = Symm.CellConstraint()[i];}
    BHeader.misflg = 0;
    BHeader.jumpax = 0;
    BHeader.ncryst=1;                    /**< crystal number */
    BHeader.lcrflg = 0;
    BHeader.ldtype = 0;                  /**< type of data: 2D (1), 3D (2), or 
					    Laue (3) */
    BHeader.jsaxs = 0;
    BHeader.nbscal = 0;
    BHeader.ngonax = 0;
    BHeader.lbmflg = 0;
    BHeader.ndet = 1;                    /**< number of detectors (current maximum
					    2) */
    BHeader.nbsetid = 1;                 /**< dataset id - should be pointer? */
    for (int i=0;i<6;i++) {BHeader.cell[i] = cell[i];}   /**< cell dimensions */
    for (int i=0;i<9;i++) {BHeader.umat[i] = 0.0;}
    for (int i=0;i<2;i++) {
      for (int j=0;j<3;j++) {BHeader.phixyz[i][j] = 0.0;}}
    for (int i=0;i<12;i++) {BHeader.crydat[i] = 0.0;}
    for (int i=0;i<3;i++) {BHeader.datum[i] = 0.0;}
    BHeader.phistt = 0.0;
    BHeader.phiend = 0.0;
    for (int i=0;i<3;i++) {BHeader.scanax[i] = 0.0;}
    BHeader.time1 = 0.0;                /**< start time */
    BHeader.time2 = 0.0;                /**< stop time */
    BHeader.bscale = 0.0;               /**< batch scale */
    BHeader.bbfac = 0.0;                /**< batch temperature factor */
    BHeader.sdbscale = 0.0;             /**< sd bscale */
    BHeader.sdbfac = 0.0;               /**< sd bbfac */
    BHeader.phirange = 0.0;             /**< phi range */
    for (int i=0;i<3;i++) {BHeader.e1[i] = 0.0;}
    for (int i=0;i<3;i++) {BHeader.e2[i] = 0.0;}
    for (int i=0;i<3;i++) {BHeader.e3[i] = 0.0;}
    for (int i=0;i<3;i++) {BHeader.source[i] = 0.0;}
    for (int i=0;i<3;i++) {BHeader.so[i] = 0.0;}
    BHeader.alambd = wavelength;   /**< wavelength (A) */
    BHeader.delamb = 0.0;         /**< dispersion (deltalambda / lambda) */
    BHeader.delcor = 0.0;         /**< correlated component */
    BHeader.divhd = 0.0;          /**< horizontal beam divergence */
    BHeader.divvd = 0.0;          /**< vertical beam divergence */
    for (int i=0;i<2;i++) {
      BHeader.dx[i] = 0.0;
      BHeader.theta[i] = 0.0;
      for (int j=0;j<2;j++) {
	for (int k=0;k<2;k++) {
	  {BHeader.detlm[i][j][k] = 0.0;}}}
    }
    // End of MTZBAT

    int nbatches_used = 0;
    for (int ib=0;ib<Nbatches;ib++)
      {
	BHeader.num = lookup.number(ib);   /**< batch number */
	ASSERT (BHeader.num > 0);
	// Create and store batch header
	batches.push_back(scala::Batch(BHeader, true, Idataset));
	nbatches_used++;
	//^
	//^	    MtzPrintBatchHeader(&BHeader);
	
      }
    return nbatches_used;
  }
//--------------------------------------------------------------
  void SCAunmergeFile::ZeroRecVector(std::vector<clipper::Coord_orth>& vc)
  {
    vc.clear();
  // Add in origin vector to both lists (a few times!)
    const int NTIMES=20;
    for (int i=0;i<NTIMES;i++)  {
      vc.push_back(clipper::Coord_orth(0.,0.,0.));
    }
  }
//--------------------------------------------------------------
  int SCAunmergeFile::ReadSaint(std::ifstream& scain,
				const scala::file_select& file_sel, 
				const scala::all_controls& controls,
				std::string& output,
				scala::hkl_unmerge_list& hkl_list)
  // Read SAINT format
  //
  {
    bool first = true;
    const double TRNTOL = 0.0001; // Tolerance on translation length^2

    // Stores for observation part    
    scala::Hkl hkl;
    int isym;
    int batch;
    int sbatch=1;
    int batchOffset = 1;
    Rtype I;  Rtype sigI;
    Rtype Ipr;  Rtype sigIpr;
    Rtype Xdet=0.0;  Rtype Ydet=0.0;
    Rtype rot=0.0;  Rtype time=0.0;
    Rtype fraction_calc=0.0;  Rtype width=0.0;
    Rtype LP=0.0; 
    int Npart=1;  int Ipart=1;
    scala::ObservationFlag ObsFlag;

    // Storage for each reflection read
    //    int h,k,l;
    int iaxis, istl, pksum, bgsum, jcryst;
    double zfr, xdp, ydp, zfp, corr, accExp, swing;
    double rotrel, corL, xg, yg, chi, other;
    DVect3 s0r;   // unit s0 in reciprocal axis frame from file
    DVect3 s2r;   // unit s2 in reciprocal axis frame

    std::vector<SaintRun> SR;  // usually only one
    int iSrun = 0;            // index to current SR element
    scala::ResoRange rrange;
    rrange.clear();
    int nreso_rej = 0;

    clipper::String line;
    while (true) {
      std::getline(scain, line );
      if (scain.eof()) break;  // end of file
      line = StringUtil::BuftoLine(line);  // remove any trailing lf or cr characters
      // Indices
      int h = std::atoi(line.substr(0, 4).c_str());
      int k = std::atoi(line.substr(4, 4).c_str());
      int l = std::atoi(line.substr(8, 4).c_str());
      // Reduce to asymmetric unit
      hkl = hkl_list.symmetry().put_in_asu(scala::Hkl(h,k,l), isym);

      I    = std::atof(line.substr(12,8).c_str());   // intensity
      sigI = std::atof(line.substr(20,8).c_str());   // SD(intensity)
      sbatch = std::atoi(line.substr(28, 4).c_str()); // SAINT batch number 

      s0r[0] = std::atof(line.substr(32,8).c_str()); // incident beam
      s2r[0] = std::atof(line.substr(40,8).c_str()); // secondary beam
      s0r[1] = std::atof(line.substr(48,8).c_str()); // y
      s2r[1] = std::atof(line.substr(56,8).c_str());
      s0r[2] = std::atof(line.substr(64,8).c_str()); // z
      s2r[2] = std::atof(line.substr(72,8).c_str());
      
      Xdet = std::atof(line.substr(83,7).c_str()); // detector pixels
      Ydet = std::atof(line.substr(90,7).c_str());
      zfr = std::atof(line.substr(97,8).c_str()); // Observed frame number
    
      xdp = std::atof(line.substr(105,7).c_str()); // Predicted x,y
      ydp = std::atof(line.substr(112,7).c_str()); //   detector pixels
      zfp = std::atof(line.substr(119,8).c_str()); // Predicted frame
      
      LP   = std::atof(line.substr(127,6).c_str()); // Applied LP etc
      corr = std::atof(line.substr(133,5).c_str()); // profile correlation
    
      accExp = std::atof(line.substr(138,7).c_str()); // accumlated hours exposure
      swing = std::atof(line.substr(145,7).c_str());  // detector swing degrees
      rot = std::atof(line.substr(152,7).c_str()); // scan axis value (omega or phi) degrees
      iaxis = std::atoi(line.substr(159,2).c_str());  // scan axis 2=omega, 3=phi
      istl = std::atoi(line.substr(161,5).c_str());  // sin theta/lambda * 10000
      pksum = std::atoi(line.substr(166,9).c_str());  // peak sum
      bgsum = std::atoi(line.substr(175,7).c_str());  // background sum
      rotrel = std::atof(line.substr(182,7).c_str()); // relative rotation degrees
      jcryst = std::atoi(line.substr(189,4).c_str());  // crystal number
      corL = std::atof(line.substr(213,6).c_str());  // Lorentz factor
      xg = std::atof(line.substr(219,8).c_str());   // X pixel corrected
      yg = std::atof(line.substr(227,8).c_str());   // Y pixel corrected
      chi = std::atof(line.substr(235,8).c_str());   // Chi degrees
      other = std::atof(line.substr(243,8).c_str());   // Phi or Omega (non-scan) degrees

      Dtype sSqr = double(istl)/10000.; // sin theta/lambda
      sSqr = 4. * sSqr * sSqr;        // 1/d^2
      if (!file_sel.in_reslimits(sSqr)) {
	nreso_rej++;
	continue;
      }
      rrange.update(sSqr);

      batch = int(zfr) + batchOffset;

      // Check for discontinuities in goniostat angles, crystal number etc
      if (first) {
	// First
	SR.push_back(SaintRun(sbatch, iaxis, chi, other, jcryst, swing));
	first = false;
	iSrun = 0;
      } else {
	if (!SR[iSrun].SameRun(sbatch, iaxis, chi, other, jcryst, swing)) {
	  // Not in current run, is it in one we have already?
	  bool found = false;
	  for (size_t ir=0;ir<SR.size();++ir) {
	    if (SR[ir].SameRun(sbatch, iaxis, chi, other, jcryst, swing)) {
	      iSrun = ir;
	      found = true;
	      break;
	    }
	  }
	  if (!found) { 
	    // make a new one
	    SR.push_back(SaintRun(sbatch, iaxis, chi, other, jcryst, swing));
	    iSrun = SR.size()-1;
	    // Batch offset, increment by multiple of 1000
	    int maxBatchLast = SR[iSrun-1].MaxBatch();
	    while (batchOffset <= maxBatchLast) {
	      batchOffset += 1000;
	    }
	    batch = int(zfr) + batchOffset;
	  }
	}
      }
      // Store for orientation stuff etc
      // actual batch, not "sbatch" from file
      SR[iSrun].StoreReflection(IVect3(h,k,l), batch, s0r, s2r, rot, istl);

      //   store as IPR as well
      Ipr = I;
      sigIpr = sigI;

      hkl_list.store_part(hkl, isym, batch, I, sigI, Ipr, sigIpr,
			  Xdet, Ydet, rot, time,
			  fraction_calc, width, LP, 
			  Npart, Ipart, ObsFlag);
    } // End loop read reflections
    // ---
    rrange.Set();
    int nobs = hkl_list.close_part_list(rrange, false);

    // Store list of columns present in file
    hkl_list.StoreDataFlags(SetSaintDataFlags());

    int dataset_index = 0;  // only one dataset

    std::vector<Dtype> SumCell(6,0.0);
    scala::MeanSD meanLambda;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (iSrun=0;iSrun<int(SR.size());++iSrun) {  // loop Saint runs
      //  get cell from d* data
      cell = SR[iSrun].CalcCell();
      output += FormatOutput::logTabPrintf(0,"\nCell inferred from file: ");
      for (int i=0;i<6;i++) output += FormatOutput::logTabPrintf(0,"%7.2f",cell[i]);
      output += FormatOutput::logTab(0,"\n");
      bool cellOK =
	CheckCells(inputcell, cell, "SAINT", output); // check agreement with any input cell
      cellOK = cellOK;
      if (inputcell[0] > 0.01) {
	cell = inputcell;
	output += FormatOutput::logTab(0,
		      "Cell set from input\n");
	SR[iSrun].SetCell(cell);
      }

      Nbatches = SR[iSrun].Nbatches();
      if (SR.size() > 1) {
	output += FormatOutput::logTabPrintf(0,
			    "\nSAINT run %4d\n", iSrun+1);
      }
      // angle (degrees) between rotation axis (omega or phi) and beam
      double rotTos0 = SR[iSrun].RotAxistoBeamAngle();
      if (SR[iSrun].ScanAxis() == "OMEGA") {
	output += FormatOutput::logTabPrintf(0,
			    "\nScan axis OMEGA: chi = %8.3f phi = %8.3f\n",
			    SR[iSrun].Chi(), SR[iSrun].Other());
	output += FormatOutput::logTabPrintf(0,
			    "Angle between OMEGA and beam = %8.3f\n", rotTos0);
      } else {
	output += FormatOutput::logTabPrintf(0,
			    "\nScan axis PHI: chi = %8.3f omega = %8.3f\n",
			    SR[iSrun].Chi(), SR[iSrun].Other());
	output += FormatOutput::logTabPrintf(0,
			    "Angle between PHI and beam = %8.3f\n", rotTos0);
      }
      output += FormatOutput::logTabPrintf(0,
			  "Detector swing angle = %8.3f\n", SR[iSrun].Swing());
      output += FormatOutput::logTabPrintf(0,
			  "Number of batches:     %6d\n", Nbatches);

      DVect3 Utrn;
      DMat33 U = SR[iSrun].Umatrix(Utrn);
      output += FormatOutput::logTabPrintf(0,
       "\nOrientation matrix reconstructed from SAINT file:\n%s\n    Determinant: %5.3f\n",
			  U.format().c_str(),
			  U.det());

      if ((Utrn*Utrn) > TRNTOL) {
      output += FormatOutput::logTabPrintf(0,
			  "WARNING: Translation component: (%6.3 f%6.3f %6.3f) (should be 0,0,0)\n\n",
			  Utrn[0],Utrn[1],Utrn[2]);
      }

      wavelength = SR[iSrun].Wavelength();
      double sdwavelength = SR[iSrun].SdWavelength();
      output += FormatOutput::logTabPrintf(0,
			  "\nWavelength inferred from SAINT file: %8.4f SD %8.4f\n",
			  wavelength, SR[iSrun].SdWavelength());
      // Replace cell & wavelength with input values if given
      double WAVELENGTH_TOL = 3.*sdwavelength;
      if (inputwavelength > 0.0) {
	//  wavelength was explicitly given
	if (std::abs(inputwavelength-wavelength) > WAVELENGTH_TOL) {
	  output += FormatOutput::logTabPrintf(0,
    "\n**** WARNING: input WAVELENGTH is significantly different from wavelength inferred from SAINT file\n");
	}
	output += FormatOutput::logTabPrintf(0,
			    "Wavelength reset from input:         %8.4f\n",
			    inputwavelength);
	wavelength = inputwavelength;
	SR[iSrun].SetWavelength(wavelength);
      } else {
	// Check to see if it's close to CuKa
	double CuKa = 1.54182;
	if (std::abs(CuKa-wavelength) > WAVELENGTH_TOL) {
	  output += FormatOutput::logTab(0,
	      "\n**** WARNING: inferred wavelength is significantly different from CuKa");
	  output += FormatOutput::logTab(0,
			"You can override this with the WAVELENGTH command");
	}
      }


      // Make batch information
      std::vector<CMtz::MTZBAT> BH = SR[iSrun].MakeBatch(wavelength);
      // Store batch headers
      for (int i=0;i<Nbatches;++i) {
	batches.push_back(scala::Batch(BH[i], true, dataset_index));
      }

      for (int i=0; i<6; i++) {SumCell[i] += batches[0].cell().UnitCell()[i];}

    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if (SR.size() > 1) {
      for (int i=0; i<6; i++) {SumCell[i] /= SR.size();}
      cell = scala::Scell(SumCell);
    
      output += FormatOutput::logTab(0,
		    std::string("\nAverage Unit cell: ")+cell.format()+"\n");
    }

    Nbatches = batches.size();

    output += FormatOutput::logTabPrintf(0,
			"\n%8d observations accepted\n",
			nobs);
    output += FormatOutput::logTabPrintf(0,
			"         Resolution range %8.3f %8.3f\n",
			rrange.ResLow(), rrange.ResHigh());
    if (file_sel.reslimits().isSet()) {
      output += FormatOutput::logTabPrintf(0,
			  "%8d observations rejected because outside resolution limits (%5.1f to %4.1fA)\n",
			  nreso_rej, file_sel.reslimits().ResLow(),
			  file_sel.reslimits().ResHigh());
    }
    return nreso_rej;
  }
//-------------------------------------------------------------
  scala::data_flags SCAunmergeFile::SetSaintDataFlags()
  // flags for which columns are st from SAINT file
  {
    scala::data_flags flags;  // Constructor sets compulsory columns true, other false
    // Optional
    flags.is_Xdet = true;
    flags.is_Ydet = true;
    flags.is_Rot = true;
    flags.is_LP = true;
    return flags;
  }
  //--------------------------------------------------------------
  bool  SaintRun::CompareMtzBat::operator()(const CMtz::MTZBAT& a, const CMtz::MTZBAT& b)
    {return (a.num < b.num);}
  //--------------------------------------------------------------
 DVect3 SaintRun::e1(0., 0., 1.0);  // Omega along +z
 DVect3  SaintRun::e2(-1.0, 0.0, 0.);  // Chi along -x
 DVect3  SaintRun::e3(0., 0., -1.0);  // Phi along -z
 DVect3  SaintRun::s0(-1.0,0.0,0.0);  // s0 along x
 int SaintRun::nxpix = 512;
 int SaintRun::nypix = 512;
//--------------------------------------------------------------
  SaintRun::SaintRun(const int& Ibatch, const int& Iaxis,
		     const double& Chi, const double& Other,
		     const int& Icryst, const double& Swing)
{
  init(Ibatch, Iaxis, Chi, Other, Icryst, Swing);
}
//--------------------------------------------------------------
  void SaintRun::init(const int& Ibatch, const int& Iaxis, const double& Chi, const double& Other,
		    const int& Icryst, const double& Swing)
{
  ibatch = Ibatch;
  iaxis = Iaxis;
  chi = BasicAngleD(Chi);
  other = BasicAngleD(Other);
  icryst = Icryst;
  swing = BasicAngleD(Swing);

  double phi = 0.0;
  double omega = 0.0;
  if (iaxis == 2) {
    // Omega scan
    phi = other;
  } else {
    // Phi scan
    omega = other;
  }
  omegar = clipper::Util::d2rad(omega);
  chir = clipper::Util::d2rad(chi);
  phir = clipper::Util::d2rad(phi);

  const int TABLE_SIZE = 5000;
  lookup.set_size(TABLE_SIZE);
  table_size = lookup.Size();  // smallest prime >= set size
  nbatches = 0;
  batchoffset = 0;
  cellDone = false;
  UmatrixDone = false;
  lambdaDone = false;
}
//--------------------------------------------------------------
double SaintRun::BasicAngleD(const double& aa)
// Reduce angle in degrees to range -180 - +180
{
  double a = aa;
  while (a <=  -180.) {a+=360.;}
  while (a >  +180.) {a-=360.;}
  return a;
}
//--------------------------------------------------------------
bool SaintRun::SameAngle(const double& a1, const double& a2)
// true is angles a1 & a2 (degrees) are the same within tolerance (modulo 360)
//  a2 is assumed to be already "Basic" ie in range 0-360
{
  const double TOL = 0.01;
  // True if angles are the same modulo 360
  return (std::abs(BasicAngleD(a1) -a2) < TOL);  
}
//--------------------------------------------------------------
bool SaintRun::SameRun(const int& Ibatch, const int& Iaxis,
		       const double& Chi, const double& Other,
		       const int& Icryst, const double& Swing)

{
  if (Ibatch != ibatch) return false;
  if (Icryst != icryst) return false;
  if (Iaxis != iaxis) return false;
  if (!SameAngle(Chi, chi)) return false;
  if (!SameAngle(Other, other)) return false;
  if (!SameAngle(Swing, swing)) return false;
  return true;
}
//--------------------------------------------------------------
  void SaintRun::StoreReflection(const IVect3& hkl, const int& batch,
			       const DVect3& s0r, const DVect3& s2r,
			       const double& rot, const int& istl)
{
  hklv.push_back(hkl);
  s0rv.push_back(s0r);
  s2rv.push_back(s2r);
  istlv.push_back(istl);
  if (iaxis == 2) {
    // Omega scan
    omegar = clipper::Util::d2rad(rot);
  } else {
    // Phi scan
    phir = clipper::Util::d2rad(rot);
  }
  // Goniostat rotation matrix
  Rv.push_back(Rmatrix(omegar,chir,phir));
  // do we have this batch already?
  int ib = lookup.lookup(batch);
  if (ib < 0) {
    // New batch, add to list
    nbatches++;
    if (nbatches > 0.8*table_size) {
      Message::message(Message_fatal
		       ("SCAIO: increase TABLE_SIZE"));
    }
    lookup.add(batch, nbatches-1);
    rotRange.push_back(scala::Range());
    ASSERT (int(rotRange.size()) == nbatches);
    ib = nbatches-1;
  }
  rotRange.at(ib).update(rot);
  //  std::cout << "Batch, rot " << batch << " " << rot << "\n";
}
//--------------------------------------------------------------
DMat33 SaintRun::Rmatrix(const double& omegar, const double& chir, const double& phir)
{
  // Total goniostat rotation
  scala::Rotation R(scala::Euler_explicit
		      (e1, omegar,
		       e2, chir,
		       e3, phir));
  return R.matrix();
}
//--------------------------------------------------------------
int SaintRun::MaxBatch() const
// maximum batch number including offset
{
  int maxbatch = -1;
  for (int i=0;i<nbatches;++i) {
    maxbatch = Max(maxbatch, lookup.number(i));
  }
  return maxbatch;
}
//--------------------------------------------------------------
std::string SaintRun::ScanAxis() const
{
  return (iaxis==2) ? "OMEGA" : "PHI";
}
//--------------------------------------------------------------
double SaintRun::RotAxistoBeamAngle() const
// return angle (degrees) between rotation axis (omega or phi) and beam
{
  if (iaxis == 3) {
    // Phi scan
    // Angle between rotation axis e3 Phi and s0
    DVect3 e3r = Rmatrix(omegar,chir,0.0) * e3; // [Omega][Chi] e3
    return clipper::Util::rad2d(acos(e3r * s0));  // unit vectors
  } else {
    // Omega scan
    return clipper::Util::rad2d(acos(e1 * s0));  // unit vectors
  }
}
//--------------------------------------------------------------
scala::Scell SaintRun::CalcCell()
// calculate cell (returned) from dstar information, linear LSQ
//
// For each reflection, we have an observational equation
//   e = hT [M] h - q
// where q = |d*|^2  = 4 (sin theta/lambda)^2 (from the file)
//       h = (h k l)
//       [M] is the reciprocal metrix tensor = [B]T [B]
//
// We can solve for [M] as a 6-vector m,
//    elements of m indexed as m(l) = m(ij) for l = 0,5 or j=0,2; i=0,j
//    ie l = j(j+1)/2 + i
// Then  m = [H]^-1 r
// where [H]kl = Sum(n=1,N)  h(ik)n h(jk)n h(il)n h(jl)n
//          ik,jk are i,j corresponding to k
//          il,jl are i,j corresponding to l
//       r(ij) = Sum(n=1,N)   h(ik)n h(jk)n qn = r(l)
{
  if (cellDone) {return cell;}

  ASSERT ((s0rv.size() == Rv.size()) &&
	  (s2rv.size() == hklv.size()) &&
	  (s2rv.size() == istlv.size()) &&
	  (s0rv.size() == hklv.size()));
  clipper::Matrix<double> H(6,6,0.0);
  std::vector<double> r(6,0.0);

  for (size_t n=0;n<hklv.size();++n){
    IVect3 hkl = hklv[n];
    double qn = double(istlv[n])/5000.; // 2 sin theta / lambda = |d*|
    qn = qn*qn;  // |d*|^2
    for (int j=0;j<3;++j){
      for (int i=0;i<=j;++i){ // loop for vector elements & matrix column
	int l = (j*(j+1))/2 + i;
	int hij = (hkl[i]*hkl[j]);
	r[l] += qn * double(hij);
	for (int jj=0;jj<3;++jj){
	  for (int ii=0;ii<=jj;++ii){ // loop for matrix row
	    int k = (jj*(jj+1))/2 + ii;
	    H(k,l) += double(hij*hkl[ii]*hkl[jj]);
	  }
	}
      }
    }
  }
  std::vector<double> m = H.solve(r);  // solve for m
  DMat33 M; // reciprocal metric tensor
  for (int j=0;j<3;++j){
    for (int i=0;i<=j;++i){ 
      int l = (j*(j+1))/2 + i;
      if (i != j) {
	M(i,j) = 0.5*m[l];
	M(j,i) = 0.5*m[l];
      } else {
	M(i,j) = m[l];
      }
    }}

  cell = scala::Scell(M,false);
  cellDone = true;
  return cell;
}
//--------------------------------------------------------------
  DMat33 SaintRun::Umatrix(DVect3& Utrn)
  //  calculate [U] from incident beam cosines
  //  returns Utrn  translational part (should = 000)
  {
    if (UmatrixDone) {Utrn = Ut; return U;}

    // unit s0 in r frame from back-rotated [R]^-1 s0
    std::vector<clipper::Coord_orth> s0xv(s0rv.size());
    for (size_t i=0;i<s0rv.size();++i){
      s0xv[i] = clipper::Coord_orth(Rv[i].inverse() * s0);  // s0r = [R]^-1 s0
    }

    std::vector<double> reccell = cell.ReciprocalCell();
    // [C] = diag(a*,b*,c*)
    DMat33 C = MVutil::SetDiagCMat33(DVect3(reccell[0],reccell[1],reccell[2]));
    // [K] = [B] [C^-1] = [F]T
    

    // [F]^-1 = ([B][C^-1])T ^-1
    DMat33 Fi =  (cell.Bmat()*C.inverse()).transpose().inverse();
    // Orthogonalised incident beam vector s0xr = [F]^-1 s0r 
    std::vector<clipper::Coord_orth> s0xrv(s0rv.size());
    for (size_t i=0;i<s0rv.size();++i){
      s0xrv[i] = clipper::Coord_orth(Fi * s0rv[i]);  // s0x = [F]^-1 s0r
    }

    // Add in origin vector to both lists (a few times!)
    for (int i=0;i<30;i++)  {
      s0xv.push_back(clipper::Coord_orth(0.,0.,0.));
      s0xrv.push_back(clipper::Coord_orth(0.,0.,0.));
    }
    
    clipper::RTop_orth UU(s0xrv, s0xv);  // [U]
    U = UU.rot();   // rotational part
    Ut = UU.trn();  // translational part, should be zero
    Utrn = Ut;
    UmatrixDone = true;
    return U;
  }
//--------------------------------------------------------------
  double SaintRun::Wavelength()
  {
    if (lambdaDone) {return wavelength;}
    if (!cellDone) CalcCell();
    if (!UmatrixDone) Umatrix(Ut);

    DMat33 B = cell.Bmat();
    // [UB]
    DMat33 UB = U * B;
    //^
    //    std::cout << "[B]\n" << B.format() << "\n";
    //    std::cout << "[UB]\n" << UB.format() << "\n";
    //^-
    // Calculate wavelength
    // When the diffraction vector s is in the diffracting position, it makes an angle
    // of (90 - theta) with the direct beam
    scala::MeanSD mnL;
    DVect3 hkl;
    for (size_t i=0;i<hklv.size();++i) {
      hkl = DVect3(double(hklv[i][0]),double(hklv[i][1]),double(hklv[i][2]));
      DMat33 RUB = Rv[i] * UB;
      DVect3 s = RUB * hkl;
      double sintheta = s * s0 /sqrt(s*s);
      double lambda = sintheta/(double(istlv[i])/10000.);
      mnL.Add(lambda);
    }
    wavelength = mnL.Mean();
    sdwavelength = mnL.SD();
    lambdaDone = true;
    return wavelength;
  }
//--------------------------------------------------------------
  std::vector<CMtz::MTZBAT> SaintRun::MakeBatch(const float& wavelength)
  // Make batch data
{
  //  jumpax is the index (1,2,3) of the reciprocal axis closest to the
  //  rotation axis. We can get this by back-rotating the rotation axis e1
  //  into crystal space with the inverse orientation matrix [U]^-1.
  //  Then jumpax corresponds to the largest component
  DVect3 UinvE;
  if (iaxis == 2) {
    // Omega scan
    UinvE = U.inverse() * e1;
  } else {
    // Phi scan
    UinvE = U.inverse() * e3;
  }
  int jumpax=0;
  double temp = 0.0;
  for (int i=0;i<3;i++) {
    if (std::abs(UinvE[i]) > temp)  {
      jumpax = i+1;
      temp = std::abs(UinvE[i]);
    }
  }
  CMtz::MTZBAT BHeader;   // header structure
  //+  BHeader.num = BatchNumber;        /**< batch number */
  strcpy(BHeader.title, "SAINT data");   /**< batch title */	      
  strcpy(BHeader.gonlab[0], "OMEGA   ");    /**< names of the three axes */
  strcpy(BHeader.gonlab[1], "CHI     ");    // 8-characters only! 
  strcpy(BHeader.gonlab[2], "PHI     ");
  BHeader.iortyp = 0;                  /**< type of orientation block (for 
					    possible future use, now = 0) */
  /**< refinement flags for cell (constraint flags)  no constraints*/
  for (int i=0;i<6;i++) {BHeader.lbcell[i] = -1;}
  BHeader.misflg = 0;                  /**< number of phixyz used (0, 1, or 2) */
  BHeader.jumpax = jumpax;             /**< reciprocal axis closest to rotation
					  axis */
  BHeader.ncryst=icryst;               /**< crystal number */
  BHeader.lcrflg = 0;                  /**< mosaicity model: 0 = isotropic, 
					  1 = anisotropic */
  BHeader.ldtype = 2;                  /**< type of data: 2D (1), 3D (2), or 
					  Laue (3) */
  /**< datum values of goniostat axes */
    BHeader.datum[1] = chi;  // e2 = Chi
  if (iaxis == 2) {
    // Omega scan
    BHeader.jsaxs = 1;                 /**< goniostat scan axis number */
    for (int i=0;i<3;i++) {BHeader.scanax[i] = e1[i];}     /**< rotation axis in lab frame */
    BHeader.datum[0] = 0.0;    // Omega (scan)
    BHeader.datum[2] = other;  // Phi
  } else {
    BHeader.jsaxs = 3;                 /**< goniostat scan axis number */
    for (int i=0;i<3;i++) {BHeader.scanax[i] = e3[i];}     /**< rotation axis in lab frame */
    BHeader.datum[2] = 0.0;    // Phi (scan)
    BHeader.datum[0] = other;  // Omega
  }
  BHeader.nbscal = 0;                  /**< number of batch scales & Bfactors
					  (0 if unset) */
  BHeader.ngonax = 3;                  /**< number of goniostat axes */
  BHeader.lbmflg = 0;                  /**< flag for type of beam info:
					  = 0 for alambd, delamb
					  = 1 also delcor, divhd, divvd */
  BHeader.ndet = 1;                    /**< number of detectors (current maximum
					  2) */
  BHeader.nbsetid = 1;                 /**< dataset id - should be pointer? */
  for (int i=0;i<6;i++) {BHeader.cell[i] = cell[i];}   /**< cell dimensions */
  /**< orientation matrix U in Fortranic order, i.e. U(1,1), U(2,1) ... */
  int k=0;
    for (int i=0;i<3;i++) {
      for (int j=0;j<3;j++) {
	BHeader.umat[k++] = U(j,i);
      }}
    /**< missetting angles at beginning and end of oscillation */
    for (int i=0;i<2;i++) {
      for (int j=0;j<3;j++) {
	BHeader.phixyz[i][j] = 0.0;
      }}
    BHeader.crydat[0] = 0.0;              /**< mosaicity */
    for (int i=1;i<12;i++) {BHeader.crydat[i] = 0.0;}
    //+    BHeader.phistt = rotRange.min();         /**< start of phi relative to datum */
    //+    BHeader.phiend = rotRange.max();         /**< end of phi relative to datum */
    BHeader.time1 = 0.0;                /**< start time */
    BHeader.time2 = 0.0;                /**< stop time */
    BHeader.bscale = 0.0;               /**< batch scale */
    BHeader.bbfac = 0.0;                /**< batch temperature factor */
    BHeader.sdbscale = 0.0;             /**< sd bscale */
    BHeader.sdbfac = 0.0;               /**< sd bbfac */
    //+    BHeader.phirange = osc_range;       /**< phi range */
    /**< vector 1 ("Cambridge" laboratory axes defining ngonax goniostat axes */
    for (int i=0;i<3;i++) {BHeader.e1[i] = e1[i];}  /**< vector 1 */
    for (int i=0;i<3;i++) {BHeader.e2[i] = e2[i];}    /**< vector 2 */
    for (int i=0;i<3;i++) {BHeader.e3[i] = e3[i];}    /**< vector 3 */
    for (int i=0;i<3;i++) {BHeader.source[i] = s0[i];} /**< idealised source vector */
    for (int i=0;i<3;i++) {BHeader.so[i] = s0[i];}     /**< source vector (unit length)*/
    BHeader.alambd = wavelength;        /**< wavelength (A) */
    BHeader.delamb = 0.0;               /**< dispersion (deltalambda / lambda) */
    BHeader.delcor = 0.0;               /**< correlated component */
    BHeader.divhd = 0.0;                /**< horizontal beam divergence */
    BHeader.divvd = 0.0;                /**< vertical beam divergence */
    BHeader.dx[0] = 0.0;                /**< xtal to detector distance */
    BHeader.dx[1] = 0.0;                /**< xtal to detector distance */
    BHeader.theta[0] = swing;           /**< detector tilt angle */
    BHeader.theta[1] = 0.0;             /**< detector tilt angle */
    /**< min & max values of detector coords (pixels) */
    BHeader.detlm[0][0][0] = 1.;
    BHeader.detlm[0][0][1] = nxpix;
    BHeader.detlm[0][1][0] = 1.;
    BHeader.detlm[0][1][1] = nypix;
    BHeader.detlm[1][0][0] = 0.;
    BHeader.detlm[1][0][1] = 0.;
    BHeader.detlm[1][1][0] = 0.;
    BHeader.detlm[1][1][1] = 0.;
    // End of MTZBAT

    std::vector<CMtz::MTZBAT> batchHeaders;
    for (int i=0;i<nbatches;++i) {
      BHeader.num = lookup.number(i);        /**< batch number */
      ASSERT (BHeader.num > 0);
      BHeader.phistt = rotRange[i].min();         /**< start of phi relative to datum */
      BHeader.phiend = rotRange[i].max();         /**< end of phi relative to datum */
      float osc_range = rotRange[i].max() - rotRange[i].min();
      BHeader.phirange = osc_range;       /**< phi range */  
      batchHeaders.push_back(BHeader);
    }
    ASSERT (int(batchHeaders.size()) == nbatches);

    // sort on batch number
    std::sort(batchHeaders.begin(), batchHeaders.end(), CompareMtzBat());
    // Force end rot value of each batch to match start rot of next
    // ie batches should be contiguous
    if (nbatches > 1) {
      // check the first two to see if rot if ascending or descending
      bool decreasing = false;
      if (batchHeaders[1].phistt < batchHeaders[0].phistt) {
	decreasing = true;
      }
      if (decreasing) {
	// swap phi limits if rot is running backwards
	for (int i=0;i<nbatches;++i) {
	  float temp = batchHeaders[i].phiend;
	  batchHeaders[i].phiend = batchHeaders[i].phistt;
	  batchHeaders[i].phistt = temp;
	}
      }
      for (int i=0;i<nbatches-1;++i) {
	float rotLast = batchHeaders[i].phiend;
	float rotNext = batchHeaders[i+1].phistt;
	float avrot = 0.5*(rotLast+rotNext);
	batchHeaders[i].phiend = avrot;
	batchHeaders[i+1].phistt = avrot;
	batchHeaders[i].phirange = batchHeaders[i].phiend -
	  batchHeaders[i].phistt;
      }
      batchHeaders[nbatches-1].phirange =
	batchHeaders[nbatches-1].phiend - batchHeaders[nbatches-1].phistt;
    }
    return batchHeaders;
  }
//--------------------------------------------------------------
//--------------------------------------------------------------
}
