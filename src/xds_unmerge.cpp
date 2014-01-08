// xds_unmerge.cpp
//
// Read XDS data into hkl_unmerge object

#include <iostream>
#include <fstream>
#include <string>

#include "xds_unmerge.hh"
#include "controls.hh"
#include "observationflags.hh"
#include "util.hh"
#include "scala_util.hh"
#include "lattice.hh"
#include "string_util.hh"
#include "range.hh"
#define ASSERT assert


using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

namespace XDSIO
{
  //--------------------------------------------------------------
  double Vec3Length(const Vec3<double>& v)
  {
    return sqrt(Vec3<double>::dot(v, v));
  }
  //--------------------------------------------------------------
  std::string Line(const char* buf)
  // Extract line from buffer, removing any trailing Cr or Lf characters
  {
    std::string line(buf);
    size_t ll = line.size();
    while (ll >= 0) {
      if (line[ll] != '\n' && line[ll] != '\r') {break;}
      ll--;
    }
    return line.substr(0,ll-1);
  }
  //--------------------------------------------------------------
  std::vector<KeyValues>  SSplit(const std::string& line)
  // Split into keys (ending in "=") and values
  {
    std::string token;
    KeyValues keyval;
    std::vector<KeyValues> keyvalues;

    size_t tokbeg = 0, tokend = 0, tokend2 = 0;
    // state = -1 before key, 0 key found, +1 key & values found
    int state = -1;

    while (1) {
      // First first non-blank
      tokbeg = line.find_first_not_of(" ", tokend);
      if (tokbeg == std::string::npos) break;
      tokend = line.find_first_of(" ", tokbeg);
      tokend2 = line.find_first_of("=", tokbeg);
      if (tokend2 < tokend) tokend = tokend2+1;
      // Token found
      token = line.substr(tokbeg, tokend-tokbeg);
      // This is a key if the last character is "="
      if (token[token.length()-1] == '=') {
	// Store previous set if present
	if (state >=0)
	  {keyvalues.push_back(keyval);}
	// Start a new KeyValues object
	keyval = KeyValues(token);
	state = 0;
	//^
	//	  std::cout << "KeyToken:" << token << "\n";
      } else {
	// ... otherwise it is a value
	//	  if (state >= 0) {
	keyval.values.push_back(token);
	if (state >= 0)	state = +1;  // if we have had a key, then flag
	      //	      std::cout << "ValToken:" << token << "\n";
      }
      if (tokend == std::string::npos) break;
    }
    if (state < 0) {
      if (token != "") {
	//	  !!!!
	if (keyval.values.size() > 0) {
	  // No key but multiple values, make first value the key
	  keyval.key = keyval.values[0];
	  // .. and remove it from the values list
	  keyval.values.erase(keyval.values.begin());
	  keyvalues.push_back(keyval);
	} else {
	  // Only a key, no values
	  keyvalues.push_back(KeyValues(token));
	}
      }
      // Store last set if not yet done
    } else if (state > 0) {
      keyvalues.push_back(keyval);
    }
    return keyvalues;
  }
//--------------------------------------------------------------
bool GetDataLine(FILE* fo, const int& nitem, std::vector<double>& data)
// Read one line of nitem doubles into data from file fo
// Returns false if end of file or line beginning "!"
{
  const int MAXLINE=256;
  char	buf[MAXLINE+1];
  char* cptr;
  char* cptr1;

  if (fgets(buf, MAXLINE, fo) == NULL) {return false;}
  buf[MAXLINE] = '\0';
  if (buf[0] == '!') {return false;}

  cptr = buf;
  // Loop items
  for (int i=0;i<nitem;i++) {
    // Skip leading spaces
    while (*cptr == ' ') {cptr++;}
    cptr1 = cptr;
    // find terminator NULL or space
    while (*cptr != '\0' && *cptr != ' ') {cptr++;}
    // cptr points to field terminator, force it to be NULL
    *cptr = '\0';  // add field terminator
    data[i] = atof(cptr1);
    cptr++;
  }
  return true;
}
//--------------------------------------------------------------
  XDSunmergeFile::XDSunmergeFile(const std::string& xdsName,
				 const int& verbose,
				 std::string& output)
  //  xdsname            name of XDS file (or logical name)
  //
  // Open file & check first line for file type
  {
    init(xdsName, verbose, output);
  }
//--------------------------------------------------------------
  void XDSunmergeFile::init(const std::string& xdsName,
			    const int& verbose,
			    std::string& output)
  //  xdsname            name of XDS file (or logical name)
  //
  // Open file & check first line for file type
  {
    output = "";
    xdsname = xdsName;
    // Open input stream from file
    xdsin = fopen(xdsname.c_str(), "r");
    if (xdsin == NULL) {
      Message::message(Message_fatal
		       ("XDSio: cannot open file "+xdsname));
    }
    const int MAXLINE=256;
    char	buf[MAXLINE+1];
    std::string line = "";
    bool OK = false;

    if (fgets(buf, MAXLINE, xdsin) != NULL) {
      // first line
      buf[MAXLINE] = '\0';
      // Remove trailing "\n"
      line = Line(buf);
      // Header lines are distinguished by the first character "!"
      if (line[0] == '!') {
	// !FORMAT=XDS_ASCII    MERGE=FALSE    FRIEDEL'S_LAW=FALSE
	//  or
	// !OUTPUT_FILE=INTEGRATE.HKL    DATE= 6-Oct-2007
	// !Generated by INTEGRATE    (VERSION  July 5, 2007)
	if (line.find("FORMAT=XDS_ASCII") != std::string::npos) {
	  if (line.find("MERGE=FALSE") != std::string::npos) {
	    // OK unmerged & scaled
	    OK = true;
	  }
	} else if (line.find("OUTPUT_FILE=INTEGRATE.HKL")
		   != std::string::npos) {
	  // OK INTEGRATE
	  OK = true;
	}
      }
    }
    if (!OK) {
      Message::message(Message_fatal
		       ("XDSio: file "+xdsname+" has wrong header"));
    }
    rewind(xdsin);

    GeomInfo = true;   // default true for all geometrical information present in file
    nxpix = 0;
    nypix = 0;
    Umat = Mat33<double>::identity();
    mosaicrange = 0.0;

    if (verbose > 0) {
      output += FormatOutput::logTab(0,
     "\nReading XDS ascii file from file "+xdsname+"\n\n");
      output += FormatOutput::logTab(0, "Header lines:\n\n");
    }

    line = "";
    bool HeaderOK = true;
    atdata = false;
    nitem = 0;
    ClearColumnNumbers();
    int status = 0;
    FirstImage = 0;
    LastImage  = 9999999;
    osc_range = 0.0;
    phi0 = 0.0;
    image0 = 1;

    // read and parse input lines
    while (fgets(buf, MAXLINE, xdsin) != NULL) {
      buf[MAXLINE] = '\0';
      // Remove trailing "\n"
      line = Line(buf);
      // Header lines are distinguished by the first character "!"
      if (line[0] == '!') {
	// Header line
	// echo to log file
	if (verbose > 0) output += FormatOutput::logTab(0,  line);

	// Header lines should not contain "*" character
	// This will occur if things are horribly wrong, format overflow
	if (line.find("*") != std::string::npos) {
	  Message::message(Message_fatal
			   ("XDS: Header line corrupt:\n"+line));
	}

	// remove leading "!"
	line = line.substr(1, line.length()-1);
	//	    line = std::string(line, 1, line.length()-1);


	std::vector<KeyValues> keyvalues = SSplit(line);
	
	  

	for (size_t ik=0;ik<keyvalues.size();ik++) {
	  status = ProcessKeyValues(keyvalues[ik], output);
	  if (status < 0) {HeaderOK = false;}
	  if (status != 0) break;
	}
      }
      if (status != 0) break;
    }
    // finished reading header
    if (verbose > 0) output += FormatOutput::logTab(0,  "\n");

    if (!HeaderOK) 
      {Message::message(Message_fatal
			("XDSio: error in header"));}

    // Check number of items
    if (ftype == XDS_INTEGRATE) {
      //  Set number of items for INTEGRATE case
      nitem = SetItems();
    }
    if (nitem != nitem_rec) {
      // For XSCALE file, allow one extra
      if (!(ftype == XSCALE && nitem_rec == nitem+1)) {
	output += FormatOutput::logTabPrintf(0,
	    "\nERROR: number of specified items disagrees with specified number %3d %3d\n",
		    nitem, nitem_rec);
      }}

    if (ftype == XSCALE) {GeomInfo = false;}  // XSCALE output is missing geometrical information

    // Process header data to get coordinate frame transformations
    if (GeomInfo) {XDStoCambridgeFrame(output, verbose);}


    // ----  Space group
    // Is this a rhombohedral spacegroup?
    //  if so, is this a rhombohedral setting?
    char RLatTyp = 'H';
    // Might be rhombohedral, false unless angles are 90,90,120
    if (scala::RhombohedralAxes(cell.UnitCell())) {RLatTyp = 'R';}
    // get extended Hermann-Mauguin symbol
    // RLatTyp does not matter except for rhombohedral cells
    spacegroup = scala::SpaceGroup(spacegroup_number);
    xHM = spacegroup.symbol_hm();
    xHM = scala::SetRlatticetype(xHM, RLatTyp);
    atdata = true;
  }
  //--------------------------------------------------------------
  scala::FileRead XDSunmergeFile::FillHklList (std::string& output,
					      const int& verbose,
					      scala::hkl_unmerge_list& hkl_list)
  {
    // File selection flags (resolution, datasets, batches etc)
    scala::file_select file_sel;
    // Scala control classes (default settings)
    //  run controls
    //  partials controls
    scala::all_controls controls;
    scala::PxdName InputPxdName;
    scala::Scell cell;
    int fileSeries = 0;

    return AddHklList(file_sel, controls, InputPxdName, cell, fileSeries,
		      output, verbose, hkl_list);
  }
  //--------------------------------------------------------------
  scala::FileRead XDSunmergeFile::AddHklList (const scala::file_select& file_sel, 
					      const scala::all_controls& controls,
					      const scala::PxdName& InputPxdName,
					      const scala::Scell& newcell,
					      const int& fileSeries,
					      std::string& output,
					      const int& verbose,
					      scala::hkl_unmerge_list& hkl_list)
  // On entry:
  //  file_sel           flags for general selection
  //                     - dataset selection
  //                     - batch selection
  //                     - resolution limits
  //                     - detector coordinate rejection ranges
  //                     - input scale factor (MULTIPLY)
  //  controls           run controls, partial controls
  //  InputPxdName       PXD name (none in XDS file)
  //  newcell            if non-null, override cell
  //  output             output string for printing
  //  verbose            set verbosity level
  //                      = 0 silent, = +1 usual summary
  //                      >= +2 debug
  // 
  // On exit:
  //  hkl_list  has been filled, but not organised:
  //  this needs a call to "prepare" or "change_symmetry"
  {
    if (!atdata) {
      {Message::message(Message_fatal
			("XDSio: not positioned after header"));}
    }
    const double TOLERANCE = 3.0;   // angular difference limit for warning
    // Override cell if required
    if (!newcell.null()) {
      if (!cell.equalsTol(newcell, TOLERANCE)) {
	output += FormatOutput::logTab(0, 
          "**** WARNING: input CELL is significantly different from cell from XDSIN file");
	output += FormatOutput::logTabPrintf(1,"XDSIN cell: ");
	for (int i=0;i<6;i++) output += FormatOutput::logTabPrintf(0,"%6.1f",cell[i]);
	output += FormatOutput::logTab(0,"\n");
	output += FormatOutput::logTabPrintf(1,"Input cell: ");
	for (int i=0;i<6;i++) output += FormatOutput::logTabPrintf(0,"%6.1f",newcell[i]);
	output += FormatOutput::logTab(0,"\n\n");
      }
      else {
	output += FormatOutput::logTabPrintf(1,"XDSIN cell overridden, input values: ");
	for (int i=0;i<6;i++) output += FormatOutput::logTabPrintf(0,"%6.1f",newcell[i]);
      }
      cell = newcell;
    }
    //  check for compulsory columns, set DataFlags
    scala::data_flags dataflags = CheckColumns(output);

    // reserve space for some plausible number of observations, not critical
    int NreflReserve = 10000;
    // Initialise hkl list
    // Cook up some sort of title
    std::string Title("From XDS file "+scala::BaseFileName(xdsname,false)+", XDS run on "+xdsdate+
		      " from images "+image_template);
    hkl_list.init(Title, NreflReserve, scala::hkl_symmetry(xHM), controls);

    // Store list of columns present in file
    hkl_list.StoreDataFlags(dataflags);
    int Nobs = ReadObservations(xdsin, file_sel, controls, fileSeries, output, hkl_list);
    Nobs = Nobs;

    // Reconstruct orientation matrix [U] from data
    if (GeomInfo) {ReconstructOrientation(hkl_list, output);}
    
    // Make dataset information (just one dataset)
    int setid = MakeDataset(InputPxdName);
    setid = setid;
    // Make batch information, only for batches containing data
    int dataset_index = 0;  // only one dataset
    int nbatches_used = MakeBatches(hkl_list.symmetry(), dataset_index);
    nbatches_used = nbatches_used;

    // Put dataset & batch information into hkl_list
    hkl_list.StoreDatasetBatch(datasets, batches);

    if (ftype == XDS_INTEGRATE) {
      // For data from INTEGRATE, update polarisation correction for synchrtron data
      UpdatePolarisation(controls.polarisationcontrol, hkl_list, output);
    }
    return scala::FileRead(true, true, true, 0);
  }
//--------------------------------------------------------------
  int XDSunmergeFile::ProcessKeyValues(const KeyValues& keyval,
					std::string& output)
  {
    // Go through all the possible keywords & store things away
    //  Returns:
    //    0   OK
    //   +1   End of Header
    //   -1   Error

    //^    std::cout << keyval.key;
    //^    if (keyval.values.size() > 0) std::cout << " " << keyval.values[0];
    //^    std::cout << "\n";

    if (keyval.key == "END_OF_HEADER")
      {
	return +1;
      }
    // !FORMAT=XDS_ASCII    MERGE=FALSE    FRIEDEL'S_LAW=TRUE
    //   only check FORMAT
    else if (keyval.key == "FORMAT=")
      {
	if (keyval.values[0] == "XDS_ASCII") {
	  ftype = XDS_ASCII;
	} else {
	  output += FormatOutput::logTab(0,
			"ERROR in XDSIN: FORMAT must be XDS_ASCII");
	  return -1;	    
	}
      }
    else if (keyval.key == "OUTPUT_FILE=")
      {
	// !OUTPUT_FILE=INTEGRATE.HKL    DATE= 6-Oct-2007
	if (keyval.values[0] == "INTEGRATE.HKL") {
	  ftype = XDS_INTEGRATE;
	} else {
	  if (ftype == XDS_ASCII) {
	    if (keyval.values[0] != "XDS_ASCII.HKL") {
	      output += FormatOutput::logTab(0,
			"ERROR in XDSIN: OUTPUT_FILE must be XDS_ASCII.HKL");
	    }
	  } else if (ftype == XSCALE) {
	    // Anything allowed
	  } else if (ftype == XDS_INTEGRATE) {
	    output += FormatOutput::logTab(0,
			  "ERROR in XDSIN: OUTPUT_FILE must be INTEGRATE.HKL");
	    return -1;
	  }  
	}
      }
    else if (keyval.key == "GENERATED")
      {
	if (keyval.values.at(1) == "XSCALE") {
	  ftype = XSCALE;
	} else if (keyval.values.at(1) == "INTEGRATE") {
	  ftype = XDS_INTEGRATE;
	} else if (keyval.values.at(1) == "CORRECT") {
	  ftype = XDS_ASCII;
	} else {
	  Message::message(Message_fatal
			   ("XDSio: invalid field on header line GENERATED"));
	}
      }
    else if (keyval.key == "MERGE=")
      {
	if (keyval.values[0] != "FALSE")
	  {
	    output += FormatOutput::logTab(0,
			  "ERROR in XDSIN: only unmerged data accepted");
	    return -1;	    
	  }
      }
    else if (keyval.key == "DATE=")
      {xdsdate = keyval.values[0];}
    else if (keyval.key == "SPACE_GROUP_NUMBER=")
      {spacegroup_number = keyval.values[0].i();}
    else if (keyval.key == "UNIT_CELL_CONSTANTS=")
      {
	float c[6];
	for (int i=0;i<6;i++)
	  {c[i] = keyval.values[i].f64();}
	cell = scala::Scell(c);
      }
    else if (keyval.key == "REFLECTING_RANGE_E.S.D.=")
      mosaicrange = keyval.values[0].f64();  // mosaicity
    else if (keyval.key == "NAME_TEMPLATE_OF_DATA_FRAMES=")
      {image_template = keyval.values[0];}
    else if (keyval.key == "DATA_RANGE=")
      {
	FirstImage = keyval.values[0].i();
	LastImage  = keyval.values[1].i();
      }
    else if (keyval.key == "X-RAY_WAVELENGTH=")
      {wavelength = keyval.values[0].f64();}
    else if (keyval.key == "INCIDENT_BEAM_DIRECTION=")
      {
	for (int i=0;i<3;i++)
	  {beam[i] = keyval.values[i].f64();}
	beam = beam.unit();
      }	
    else if (keyval.key == "ROTATION_AXIS=")
      {
	for (int i=0;i<3;i++)
	  {rotaxis[i] = keyval.values[i].f64();}
      }	
    else if (keyval.key == "OSCILLATION_RANGE=")
      {osc_range = keyval.values[0].f64();}
    else if (keyval.key == "STARTING_ANGLE=")
      {
	phi0 = keyval.values[0].f64();
      }
    else if (keyval.key == "STARTING_FRAME=")
      {image0 = keyval.values[0].i();}
    else if (keyval.key == "DIRECTION_OF_DETECTOR_X-AXIS=")
      {
	for (int i=0;i<3;i++)
	  {dxaxis[i] = keyval.values[i].f64();}
      }	
    else if (keyval.key == "DIRECTION_OF_DETECTOR_Y-AXIS=")
      {
	for (int i=0;i<3;i++)
	  {dyaxis[i] = keyval.values[i].f64();}
      }	
    else if (keyval.key == "DETECTOR_DISTANCE=")
      {det_dist = keyval.values[0].f64();}
    else if (keyval.key == "ORGX=")
      {orgx = keyval.values[0].f64();}
    else if (keyval.key == "ORGY=")
      {orgy = keyval.values[0].f64();}
    else if (keyval.key == "NX=")
      {nxpix = keyval.values[0].i();}
    else if (keyval.key == "NY=")
      {nypix = keyval.values[0].i();}
    else if (keyval.key == "QX=")
      {qxpixsize = keyval.values[0].f64();}
    else if (keyval.key == "QY=")
      {qypixsize = keyval.values[0].f64();}
    else if (keyval.key == "NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD=")
      {nitem_rec = keyval.values[0].i();}
    else if (keyval.key == "ITEM_H=")
      {col_H = keyval.values[0].i()-1;nitem++;}
    else if (keyval.key == "ITEM_K=")
      {col_K = keyval.values[0].i()-1;nitem++;}
    else if (keyval.key == "ITEM_L=")
      {col_L = keyval.values[0].i()-1;nitem++;}
    else if (keyval.key == "ITEM_IOBS=")
      {col_IOBS = keyval.values[0].i()-1;nitem++;}
    else if (keyval.key == "ITEM_SIGMA(IOBS)=")
      {col_SIGMAI = keyval.values[0].i()-1;nitem++;}
    else if (keyval.key == "ITEM_XD=")
      {col_XD = keyval.values[0].i()-1;nitem++;}
    else if (keyval.key == "ITEM_YD=")
      {col_YD = keyval.values[0].i()-1;nitem++;}
    else if (keyval.key == "ITEM_ZD=")
      {col_ZD = keyval.values[0].i()-1;nitem++;}
    else if (keyval.key == "ITEM_RLP=")
      {col_RLP = keyval.values[0].i()-1;nitem++;}
    else if (keyval.key == "ITEM_PEAK=")
      {col_PEAK = keyval.values[0].i()-1;nitem++;}
    else if (keyval.key == "ITEM_CORR=")
      {col_CORR = keyval.values[0].i()-1;nitem++;}
    else if (keyval.key == "ITEM_PSI=")
      {col_PSI = keyval.values[0].i()-1;nitem++;}
    else if (keyval.key == "ITEM_ISET=")
      {col_ISET = keyval.values[0].i()-1;nitem++;}
    else if (keyval.key == "ITEM_DCY=")
      {col_DCY = keyval.values[0].i()-1;nitem++;}

    return 0;

  }
//--------------------------------------------------------------
  void XDSunmergeFile::XDStoCambridgeFrame(std::string& output,
					   const int& verbose)
  // Process header data to get coordinate frame transformations
  {
    //
    // In both XDS and CCP4 mtz files, the cordinate frame is fully described
    // by vectors defining the rotation axes and the incident beam, so in
    // principle any orthogonal frame may be used. Nevertheless it is convenient
    // and conventional in the CCP4 context to work in the so-called "Cambridge"
    // frame, so we will work out the conversions here
    //
    // Definitions:
    //  1) "Cambridge" frame used by Mosflm & CCP4 
    //     a) principal rotation axis e is along z
    //     b) Xray beam is ~along x (note that s0 is anti-parallel to beam ~= -x)
    //     c) y orthogonal to z & x
    //
    // Using the rotation axis & beam definitions in the XDS header, we can work
    // out a conversion matrix [Q] such that x(Cambridge) = [Q] x(XDS)
    //

    // (i)  Orthogonal base vector set in Cambridge frame is the identity matrix

    // (ii) Orthogonal base vector set in XDS frame
    //      XDS header defines direction of rotation axis & direction of beam
    Vec3<double> b3Xds = rotaxis.unit();  // xc direction is along rotation axis
    Vec3<double> b2Xds = Vec3<double>::cross(b3Xds, beam).unit(); // yc  
    Vec3<double> b1Xds = Vec3<double>::cross(b2Xds, b3Xds).unit(); // zc

    // Conversion matrix [Q] = (b1 b2 b3)^-1
    Qxds2cam = Mat33<double>(b1Xds[0], b2Xds[0], b3Xds[0],
			     b1Xds[1], b2Xds[1], b3Xds[1],
			     b1Xds[2], b2Xds[2], b3Xds[2]).inverse();

    // XDS detector frame:
    //   Header defines direction of detctor XD & YD in laboratory frame
    Vec3<double> e1d = dxaxis.unit();
    Vec3<double> e2d = dyaxis.unit();
    Vec3<double> e3d = Vec3<double>::cross(e1d, e2d).unit();
    // Ed = (e1 e2 e3)
    Mat33<double> Ed =  Mat33<double>(e1d[0], e2d[0], e3d[0],
				      e1d[1], e2d[1], e3d[1],
				      e1d[2], e2d[2], e3d[2]);
    // Conversion from pixel position in mm in XDS frame t(xds) to position
    // in Cambridge frame = [Q][Ed]
    QEd = Qxds2cam * Ed;

    if (verbose > 0)
      {
	output += FormatOutput::logTab(0,
		      "\nMatrix to transform XDS axis system to CCP4 frame:\n"+
		      Qxds2cam.format()+"\n");
	output += FormatOutput::logTab(0,
		      "\nMatrix to transform XDS detector coordinates to CCP4 frame:\n"+
		      QEd.format()+"\n");

	// Rotation axis in Cambridge frame
	Vec3<double> raxC = Qxds2cam * rotaxis;
	// Beam in Cambridge frame
	Vec3<double> beamC = Qxds2cam * beam;
	output += FormatOutput::logTabPrintf(0,
			    "\nRotation axis in CCP4 frame: (%6.3f %6.3f %6.3f)\n",
			    raxC[0],raxC[1],raxC[2]);
	output += FormatOutput::logTabPrintf(0,
			    "\nIncident beam in CCP4 frame: (%6.3f %6.3f %6.3f)\n\n",
			    beamC[0],beamC[1],beamC[2]);
      }

  }
//--------------------------------------------------------------
  int XDSunmergeFile::ReadObservations(FILE* xdsin,
				       const scala::file_select& file_sel, 
				       const scala::all_controls& controls,
				       const int& fileSeries,
				       std::string& output,
				       scala::hkl_unmerge_list& hkl_list)
  // Read all observations into list
  // Returns number of observations
  {
    if (xdsin == NULL) {
      Message::message(Message_fatal
		       ("XDSio: file not open "+xdsname));
    }



    // Rest of line
    std::vector<double> data(nitem_rec);

    // Stores for observation part
    scala::Hkl hkl;
    int isym;
    int batch_serial;
    int batch;
    Rtype I;  Rtype sigI;
    Rtype Ipr;  Rtype sigIpr;
    Rtype Xdet;  Rtype Ydet;
    Rtype phi;  Rtype time;
    Rtype fraction_calc=0.0;  Rtype width=0.0;
    Rtype LP=0.0; 
    int Npart;  int Ipart=1;
    scala::ObservationFlag ObsFlag;
    // Flag for misfit: flag as if PKratio too high, no value given (or output)
    scala::ObservationFlag ObsFlagMisfit(scala::ObservationFlag::FLAG_PKRATIO, 0.0);

    int nfract_low = 0; // PART < 100, keep but count
    int nmisfit = 0;    // misfit flagged with negative sigma, flag these
    int nreso_rej = 0;  // count observations outside resolution limits

    scala::ResoRange rrange;
    rrange.clear();
    scala::Range xdrange, ydrange;


    const double FRAC_MIN = 0.98;
    double min_fract = 10000.;
    //    Rtype  AcceptFractMin = controls.partials.accept_fract_min();

    // "BATCH" numbers will be image numbers and may range from
    //    FirstImage to LastImage
    Nbatches = LastImage - FirstImage + 1;
    // count number of observations in each batch
    NobsInBatch.assign(Nbatches, 0); //clear counts
    int maxbatchserial = 0;

    int ncount = 0;

    // Read data lines into array data until
    // either line beginning "!" or eof
    while (GetDataLine(xdsin, nitem_rec, data)) {
      ncount++;
      ObsFlag = scala::ObservationFlag();
      //^
      //^	for (int i=0;i<nitem_rec;i++) {std::cout << " " << data[i];}
      //^	std::cout << "\n";   //^

      // Reduce to asymmetric unit
      hkl = hkl_list.symmetry().put_in_asu
	(scala::Hkl(Nint(data[col_H]), Nint(data[col_K]), Nint(data[col_L])), isym);
      Dtype sSqr = hkl.invresolsq(cell);
      if (!file_sel.in_reslimits(sSqr))	{
	nreso_rej++;
	continue;
      }
      rrange.update(sSqr);
      
      // ImageNumber should be in range (FirstImage-1) to LastImage
      batch = Max(FirstImage, Min(int(data[col_ZD])+1, LastImage));
      // Do we want this batch? (batch exclusions etc)
      bool accept = file_sel.accept_batch(batch, fileSeries); 
      if (accept) {
	batch_serial = Min(Max(1, batch - FirstImage + 1), Nbatches);
	maxbatchserial = Max(maxbatchserial, batch_serial);
	// count them
	NobsInBatch[batch_serial-1]++;

	if (osc_range > 0.000001) {
	  // ZD is image number for reflection
	  //  get phi. note that for image 1, ZD is in range 0->1
	  //cv: follow the explanation on 
	  //cv: http://www.mpimf-heidelberg.mpg.de/~kabsch/xds/html_doc/xds_parameters.html#STARTING_ANGLE=
	  phi = phi0 + osc_range*(data[col_ZD] - double(image0-1));
	} else if (col_ZD >= 0) {
	  phi = data[col_ZD];
	} else {
	  phi = 0.0;
	}
	// store time = phi
	time = phi;
	
	// I & sigmaI
	I = data[col_IOBS];
	sigI = data[col_SIGMAI];
	if (sigI < 0.0) {
	  // Negative sigma means misfit in XSCALE, flag these
	  nmisfit++;
	  sigI = std::abs(sigI);
	  // Set observation flag
	  ObsFlag = ObsFlagMisfit;
	}
	
	// Apply input scale (MULTIPLY)
	I *= file_sel.InputScale();
	sigI *= file_sel.InputScale();
	//   store as IPR as well
	Ipr = I;
	sigIpr = sigI;
	// XD, YD
	Xdet = data[col_XD];
	Ydet = data[col_YD];
	xdrange.update(Xdet);
	ydrange.update(Ydet);
	
	// ITEM_PEAK  A value less than 100.0 indicates that the reflection
	// was either incompletely recorded or overlapping with neighbouring
	// reflections or included untrusted pixels in the profile region.
	Npart = 1;  // always fully recorded
	if (col_PEAK >= 0) {
	  fraction_calc = 0.01*data[col_PEAK];
	  if (fraction_calc < FRAC_MIN) {
	    nfract_low++;
	    min_fract = Min(min_fract, fraction_calc );
	  }
	}
	// LP
	if (col_RLP >= 0)
	  { LP = data[col_RLP];}
	
	// width = 0, Npart = Ipart = 1, ObsFlag = 0
	hkl_list.store_part(hkl, isym, batch, I, sigI, Ipr, sigIpr,
			  Xdet, Ydet, phi, time,
			  fraction_calc, width, LP, 
			    Npart, Ipart, ObsFlag);
      } // end accept
    } // end data read

    rrange.Set(); // resolution range is set
    int nobs = hkl_list.close_part_list(rrange, false);
    
    output += FormatOutput::logTabPrintf(0,
			"\n%8d observations accepted\n",
			nobs);
    output += FormatOutput::logTabPrintf(0,
			"         Resolution range %8.3f %8.3f\n",
			rrange.ResLow(), rrange.ResHigh());
    output += FormatOutput::logTabPrintf(0,
	"%8d accepted incomplete observations with PART < %4.2f, minimum %4.2f\n",
			nfract_low, FRAC_MIN, min_fract);
    if (file_sel.reslimits().isSet()) {
      output += FormatOutput::logTabPrintf(0,
	   "%8d observations rejected because outside resolution limits (%5.1f to %4.1fA)\n",
			   nreso_rej, file_sel.reslimits().ResLow(),
			   file_sel.reslimits().ResHigh());
      }
    output += FormatOutput::logTabPrintf(0,
	"%8d observations flagged as MISFITS in XDS, flagged as bad PKratio\n",
			nmisfit);

    bool xydreset = false;
    if (xdrange.max() > nxpix) {
      xydreset = true;
      nxpix = int(xdrange.max()) + 1;
    }
    if (ydrange.max() > nypix) {
      xydreset = true;
      nypix = int(ydrange.max()) + 1;
    }
    if (xydreset) {
      output += FormatOutput::logTabPrintf(0,"\nDetector limits reset to %8d %8d\n",
			  nxpix, nypix);
    }
    // Update batch range info if not in headers
    if (!GeomInfo) {
      Nbatches = maxbatchserial;  // serial from 1
      NobsInBatch.resize(Nbatches);
    }

    return nobs;
  }

//--------------------------------------------------------------
  void  XDSunmergeFile::ClearColumnNumbers()
  {
    // Zero all column numbers
    col_H = -1;
    col_K = -1;
    col_L = -1;
    col_IOBS = -1;
    col_SIGMAI = -1;
    col_XD = -1;
    col_YD = -1;
    col_ZD = -1;
    col_RLP = -1;
    col_PEAK = -1;
    col_CORR = -1;
    col_PSI = -1;
    col_ISET = -1;
    col_DCY = -1;
  }
//--------------------------------------------------------------
  scala::data_flags XDSunmergeFile::CheckColumns(std::string& output) const
  {
    // Check that compulsory columns are present, & set DataFlags
    // to indicate which data items are present
    bool present = true;
    if (col_H < 0) present = false;
    if (col_K < 0) present = false;
    if (col_L < 0) present = false;
    if (col_IOBS < 0) present = false;
    if (col_SIGMAI < 0) present = false;
    if (col_XD < 0) present = false;
    if (col_YD < 0) present = false;
    if (col_ZD < 0) present = false;
    if (!present) {
      output += FormatOutput::logTab(0, 
     "\n**** ERROR **** not all compulsory data items are present in XDS file ****\n");
      output += FormatOutput::logTab(0, 
	    "Compulsory items are H, K, L, IOBS, SIGMA(IOBS), XD, YD, ZD\n");
      Message::message(Message_fatal
			 ("XDSio: missing compulsory items"));
    }
    scala::data_flags dataflags;  // set MTZ compulsory columns true, others false
    // XDS compulsory so they are present
    dataflags.is_Xdet = true;
    dataflags.is_Ydet = true;
    dataflags.is_Rot = true;
    // optional columns
    if (col_PEAK >= 0) dataflags.is_fractioncalc = true;
    if (col_RLP >= 0) dataflags.is_LP = true;
    dataflags.is_ObsFlag = true;    // bad data flag for MISFIT (column FLAG)
    dataflags.is_BgPkRatio = false; // no BGPKRATIO column 
    return dataflags;
  }
//--------------------------------------------------------------
  void XDSunmergeFile::ReconstructOrientation(scala::hkl_unmerge_list& hkl_list,
					      std::string& output)
  // Reconstruct orientation matrix [U] from geometry information on the
  // observations
  // Not all observations are needed, just enough
  // Sets Umat
  {
    const int NUSED = 200;  // Number of observations to use
    const double TRNTOL = 0.0001; // Tolerance on translation length^2

    // Number of observations
    int Nobs = hkl_list.num_parts();
    // Use every Nsample'th observation to get ~ NUSED
    int Nsample = 1;
    if (Nobs > 2*NUSED) Nsample = Nobs/NUSED + 1;

    output += FormatOutput::logTabPrintf(0, 
	"\nReconstructing orientation matrix [U] from %4d observations\n",
			Nobs/Nsample);


    // Orthogonalisation matrix [B]
    Mat33<Dtype> B = cell.Bmat();

    int nused = 0; // count observations used
    scala::observation_part part;
    clipper::Coord_orth Bh;
    Vec3<double> tdXds;
    Vec3<double> s;
    s0C = -(Qxds2cam*beam);                   // s0, Cambridge frame, unit length
    Vec3<double> s0 = (1./wavelength)*s0C;   // s0, Cambridge frame, length 1/wavelength
    clipper::Coord_orth f;  // difraction vector

    // Two corresponding vector lists, one of orthogonalised hkl = Bh,
    // the other of diffraction vectors in diffracting position
    std::vector<clipper::Coord_orth> hklOrth;
    std::vector<clipper::Coord_orth> hklDiff;
    // Add in origin vector to both lists (a few times!)
    for (int i=0;i<10;i++)
      {
	hklOrth.push_back(clipper::Coord_orth(0.,0.,0.));
	hklDiff.push_back(clipper::Coord_orth(0.,0.,0.));
      }

    // Loop all observations in list
    for (int ip=0;ip<Nobs;ip++)
      {
	if (ip%Nsample == 0) // subset of observations if necessary
	  {
	    part = hkl_list.find_part(ip);
	    // Orthogonalise hkl = [B]h
	    scala::Hkl horig = hkl_list.symmetry().get_from_asu(part.hkl(), part.isym());
	    Vec3<double> v = B * Vec3<Dtype>(horig);
	    Bh = clipper::Coord_orth(v);
	    // -Phi matrix, rotation along z
	    double phi = clipper::Util::d2rad(part.phi());
	    Mat33<double> PhiInv = clipper::Rotation(clipper::Polar_ccp4(0.0,0.0,-phi)).matrix();

	    // Secondary beam in mm, XDS frame
	    tdXds[0] = qxpixsize * (part.Xdet() - orgx);
	    tdXds[1] = qypixsize * (part.Ydet() - orgy);
	    tdXds[2] = det_dist;
	    // to Cambridge frame
	    s = (QEd * tdXds).unit();
	    // diffraction vector in diffracting position = (s0 + s)
	    // s0 = negated beam direction
	    // Back rotate by -phi to phi = 0 position
	    Vec3<double> fd = (s0 + (1./wavelength)*s);
	    f = clipper::Coord_orth(PhiInv*fd);
	    hklOrth.push_back(Bh);
	    hklDiff.push_back(f);
	    nused++;
	    //^
	    //	    std::cout << "fd: " << fd.format() << "  phi " << part.phi() << "\n";
	    //	    std::cout << "|Bh| " << Vec3Length(Bh) << "  |f| " << Vec3Length(f)
	    //		      << "  ratio " << Vec3Length(Bh)/Vec3Length(f) << "\n";
	    //	    std::cout << part.hkl().format() << " " << part.isym() << " " << horig.format() << "\n";
	    //	    	    printf
	    //	    	      ("Bh %6.3f %6.3f %6.3f td %6.1f %6.1f %6.1f s %6.3f %6.3f %6.3f f %6.3f %6.3f %6.3f\n",
	    //	    		   Bh[0],Bh[1],Bh[2],
	    //	    		   tdXds[0],tdXds[1],tdXds[2],
	    //	    		   s[0],s[1],s[2],
	    //	    		   f[0],f[1],f[2]);
	    //^-
	  }
      }

    // Make rotation matrix [U] from the two lists, as a LSQ superposition
    //  see Clipper documentation
    // The algorithm employed is that of Kearsley, S.K. (1989) '
    // On the orthogonal transformation used for structural comparisons'.
    // Acta Cryst. A45, 208-210.
    // Ut includes a translation component which should be 0,0,0
    ASSERT(hklOrth.size() == hklDiff.size());
    clipper::RTop_orth Ut(hklOrth, hklDiff);
    Umat = Ut.rot();

    output += FormatOutput::logTabPrintf(1,
			"\nOrientation matrix [U]:\n%s\n      Determinant = %5.3f\n",
			Umat.format().c_str(), Umat.det());
    if ((Ut.trn()*Ut.trn()) > TRNTOL) {
      output += FormatOutput::logTabPrintf(0,
	  "WARNING: Translation component: (%6.3 f%6.3f %6.3f) (should be 0,0,0)\n\n",
			  Ut.trn()[0],Ut.trn()[1],Ut.trn()[2]);
    }
  }
//--------------------------------------------------------------
  int XDSunmergeFile::MakeDataset(const scala::PxdName& InputPxdName)
  // Enter the single dataset into list
  // Returns dataset ID number
  {
    // Default to dummy Project, Crystal, Dataset names
    // These are updated from the user input
    int setid = 1; // just one dataset, numbered 1
    // Default names
    scala::PxdName pxd("XDSproject", "XDScrystal", "XDSdataset");
    //  incorporate any non-blank elements
    pxd.update(InputPxdName);
    datasets.clear();
    datasets.push_back(scala::Xdataset(pxd, cell, wavelength, setid));
    return setid;
  }
  //--------------------------------------------------------------
  int XDSunmergeFile::MakeBatches(const scala::hkl_symmetry& Symm,
				  const int& Idataset)
  // Make list of batch data for all batches with observations
  // Returns number actually used
  {

    // Some geometrical calculations
    
    // The rotation axis (Phi) is along z
    Vec3<double> e1(0., 0., 1.0);
    Vec3<double> s0ideal(-1.0, 0., 0.);
    //  jumpax is the index (1,2,3) of the reciprocal axis closest to the
    //  rotation axis. We can get this by back-rotating the rotation axis e1
    //  into crystal space with the inverse orientation matrix [U]^-1.
    //  Then jumpax corresponds to the largest component
    Vec3<double> e1c = Umat.inverse() * e1;
    int jumpax=0;
    double temp = 0.0;
    for (int i=0;i<3;i++) {
      if (std::abs(e1c[i]) > temp) {
	jumpax = i+1;
	temp = std::abs(e1c[i]);
      }
    }

    if (!GeomInfo) { // No geometry information, set some indicators
      s0C = s0ideal;
      jumpax = 0;
      //      float nullmat[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      //      Umat = MVutil::SetCMat33(nullmat);
    }

    CMtz::MTZBAT BHeader;   // header structure
    // Store all common data, the bits that are the same for all batches
    // Items marked "//+" here are individual to each batch, but are listed
    // here for completeness
    //+     BHeader.int num;        /**< batch number */
    strcpy(BHeader.title, "XDS data");   /**< batch title */	      
    strcpy(BHeader.gonlab[0], "PHI    ");    /**< names of the three axes */
    strcpy(BHeader.gonlab[1], "       ");    // 8-characters only! 
    strcpy(BHeader.gonlab[2], "       ");
    BHeader.iortyp = 0;                  /**< type of orientation block (for 
					    possible future use, now = 0) */
    /**< refinement flags for cell (constraint flags)*/
    for (int i=0;i<6;i++) {BHeader.lbcell[i] = Symm.CellConstraint()[i];}
    BHeader.misflg = 0;                  /**< number of phixyz used (0, 1, or 2) */
    BHeader.jumpax = jumpax;             /**< reciprocal axis closest to rotation
					    axis */
    BHeader.ncryst=1;                    /**< crystal number */
    BHeader.lcrflg = 0;                  /**< mosaicity model: 0 = isotropic, 
					    1 = anisotropic */
    BHeader.ldtype = 2;                  /**< type of data: 2D (1), 3D (2), or 
					    Laue (3) */
    BHeader.jsaxs = 1;                   /**< goniostat scan axis number */
    BHeader.nbscal = 0;                  /**< number of batch scales & Bfactors
					    (0 if unset) */
    BHeader.ngonax = 1;                  /**< number of goniostat axes */
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
	BHeader.umat[k++] = Umat(j,i);
      }}
    /**< missetting angles at beginning and end of oscillation */
    for (int i=0;i<2;i++) {
      for (int j=0;j<3;j++) {
	BHeader.phixyz[i][j] = 0.0;
      }}
    BHeader.crydat[0] = mosaicrange;              /**< mosaicity */
    for (int i=1;i<12;i++) {BHeader.crydat[i] = 0.0;}
    /**< datum values of goniostat axes */
    for (int i=0;i<3;i++) {BHeader.datum[i] = 0.0;}
    //+    BHeader.phistt;         /**< start of phi relative to datum */
    //+    BHeader.phiend;         /**< end of phi relative to datum */
    for (int i=0;i<3;i++) {BHeader.scanax[i] = e1[i];}     /**< rotation axis in lab frame */
    BHeader.time1 = 0.0;                /**< start time */
    BHeader.time2 = 0.0;                /**< stop time */
    BHeader.bscale = 0.0;               /**< batch scale */
    BHeader.bbfac = 0.0;                /**< batch temperature factor */
    BHeader.sdbscale = 0.0;             /**< sd bscale */
    BHeader.sdbfac = 0.0;               /**< sd bbfac */
    BHeader.phirange = osc_range;       /**< phi range */
    /**< vector 1 ("Cambridge" laboratory axes defining ngonax goniostat axes */
    for (int i=0;i<3;i++) {BHeader.e1[i] = e1[i];}  /**< vector 1 */
    for (int i=0;i<3;i++) {BHeader.e2[i] = 0.0;}    /**< vector 2 */
    for (int i=0;i<3;i++) {BHeader.e3[i] = 0.0;}    /**< vector 3 */
    for (int i=0;i<3;i++) {BHeader.source[i] = s0ideal[i];} /**< idealised source vector */
    for (int i=0;i<3;i++) {BHeader.so[i] = s0C[i];}         /**< source vector (unit length)*/
    BHeader.alambd = wavelength;        /**< wavelength (A) */
    BHeader.delamb = 0.0;               /**< dispersion (deltalambda / lambda) */
    BHeader.delcor = 0.0;               /**< correlated component */
    BHeader.divhd = 0.0;                /**< horizontal beam divergence */
    BHeader.divvd = 0.0;                /**< vertical beam divergence */
    BHeader.dx[0] = det_dist;           /**< xtal to detector distance */
    BHeader.dx[1] = 0.0;                /**< xtal to detector distance */
    BHeader.theta[0] = 0.0;             /**< detector tilt angle */
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

    int nbatches_used = 0;
    for (int ib=0;ib<Nbatches;ib++) {
      if (NobsInBatch[ib] > 0) {
	// skip empty batches
	//	    int image = ib+1;  // numbered from 1  (batch_serial)
	int ibatch = ib + FirstImage;  // numbered from FirstImage  (batch)
	BHeader.num = ibatch;   /**< batch number */
	if (osc_range > 0.0001) {
	  /**< start of phi relative to datum */
	  BHeader.phistt = phi0 + osc_range*(ibatch - double(image0));
	  /**< end of phi relative to datum */
	  BHeader.phiend = BHeader.phistt + osc_range;
	} else {
	  /**< start of phi relative to datum */
	  BHeader.phistt = 0.0;
	  /**< end of phi relative to datum */
	  BHeader.phiend = 0.0;
	}
	// Create and store batch header
	batches.push_back(scala::Batch(BHeader, true, Idataset));
	nbatches_used++;
	//^
	//^	    MtzPrintBatchHeader(&BHeader);
	
      }
    }
    return nbatches_used;
  }
//--------------------------------------------------------------
  int XDSunmergeFile::SetItems()
  // Set up item list for INTEGRATE option
  //^H,K,L,IOBS,SIGMA,XCAL,YCAL,ZCAL,RLP,PEAK,CORR,MAXC,
  //^             XOBS,YOBS,ZOBS,ALF0,BET0,ALF1,BET1,PSI

  // Returns number of items
  {
    col_H = 0;
    col_K = 1;
    col_L = 2;
    col_IOBS=3;
    col_SIGMAI=4;
    col_XD=5;
    col_YD=6;
    col_ZD=7;
    col_RLP=8;
    col_PEAK=9;
    col_CORR=10;

    int nitems = 20;
    return nitems;
  }
//--------------------------------------------------------------
  void XDSunmergeFile::UpdatePolarisation(const scala::PolarisationControl& polarisationcontrol,
					  scala::hkl_unmerge_list& hkl_list,
					  std::string& output)
  // Update polarisation if INTEGRATE file
  // The INTEGRATE step in XDS applies a polarisation correction for an unpolarised
  // beam, so if the data are from a synchrotron we need to apply the additional correction
  // for a polarised incident beam
  //
  // P0 = 1/2 [ 1 + cos^2 2theta)  for unpolarised beam
  // P' = 1/2 Xsi cos 2rho sin^2 2theta for polarised beam
  // P = P0 + P'
  // rho is angle between the spindle projected on to an orthogonal detector plane
  //    and the projection of the diffracted beam
  {
    if (ftype != XDS_INTEGRATE) return;
    output += "\n";
    double polarisationfactor =  polarisationcontrol.Factor();
    if (polarisationcontrol.IsSet() && polarisationfactor == 0.0) {
      // explicitly set to no polarisation, do nothing
      output +=
	"No update to to polarisation correction for XDS-INTEGRATE file, assuming unpolarised beam";
      return;
    }
    bool inhouse = false;
    if (polarisationcontrol.IsSet()) {
      output += FormatOutput::logTabPrintf(0,
			     "Polarisation factor from input: %7.4f\n", polarisationfactor);
    } else {
      // try to set defaults
      if (wavelength > 0.001) {
	// Is this in-house data? check wavelength for Cu, Mo & Cr
	if (Close<double, double>(wavelength, 1.5418, 0.0019)) {
	  inhouse = true;
	  output += FormatOutput::logTabPrintf(0,
	     "Data assumed to be from unpolarised in-house Cu source, wavelength: %8.5f\n", wavelength);
	} else if (Close<double, double>(wavelength, 0.7107, 0.0002)) {
	  inhouse = true;
	  output += FormatOutput::logTabPrintf(0,
	     "Data assumed to be from unpolarised in-house Mo source, wavelength: %8.5f\n", wavelength);
	} else if (Close<double, double>(wavelength, 2.29, 0.01)) {
	  inhouse = true;
	  output += FormatOutput::logTabPrintf(0,
	     "Data assumed to be from unpolarised in-house Cr source, wavelength: %8.5f\n", wavelength);
	}
      } else {
	output == "No wavelength available, assume synchrotron source";
      }
      if (inhouse) {
	output +=
  "No update to to polarisation correction for XDS-INTEGRATE file, assuming unpolarised beam";
	return;
      }
      // default for synchotron
      polarisationfactor = polarisationcontrol.Default();
      output += FormatOutput::logTabPrintf(0,
        "Default polarisation factor used: %7.4f\n", polarisationfactor);
    }
    bool neg = (polarisationfactor < 0.0);
    double pfxds = 0.5*std::abs(polarisationfactor)+0.5;
    if (neg) pfxds = -pfxds;
    output += FormatOutput::logTabPrintf(1,
			 "equivalent to XDS  FRACTION_OF_POLARIZATION= %7.4f\n",
					 pfxds);
    // Apply updates to all parts in list
    scala::Range corrections = hkl_list.UpdatePolarisationCorrections(false, polarisationfactor);
    output += FormatOutput::logTabPrintf(0,
			 "Range of corrections to polarisation: %8.4f to %8.4f\n",
					 corrections.min(), corrections.max());
  }
//--------------------------------------------------------------

}
