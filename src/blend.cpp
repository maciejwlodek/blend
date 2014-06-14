///************************************************************************************************************/
/************************************************************************************************************/
/********* blend.cpp                                                                                *********/
/*********                                                                                          *********/
/********* Copyright (C) 2014 Diamond Light Source & Imperial College London                        *********/
/*********                                                                                          *********/
/********* Authors: James Foadi & Gwyndaf Evans                                                     *********/
/*********                                                                                          *********/
/********* This code is distributed under the BSD license, a copy of which is                       *********/
/********* included in the root directory of this package.                                          *********/
/************************************************************************************************************/
/************************************************************************************************************/
// CHANGES IN VERSION 0.5.3 - 14/06/2014
// - Added condition on images and resolution to skip radiation damage procedure: the procedure is not carried
//   out if there are less than 6 images and / or if the highest resolution is less than 5 angstroms. Radiation
//   damage procedure is also skipped if keyword RADFRAC is less than 0.000001.
// CHANGES IN VERSION 0.5.2 - 20/05/2014
// - Small correction / addition to html documentation.
// - Fixed a bug in the calculation and plotting of LCV values when some of the original data sets are discarded.
// - New syntax for making combination mode easier. Now the input command line can accept the following items:
//    a) d or d1-d2 or d1,d2             : means add data set d, or all data sets from d1 to d2, or data sets d1
//                                         and d2
//    b) [n] or [n1-n2] or [n1,n2]       : means add all data sets composing cluster n, or all data sets
//                                         composing clusters from number n1 to number n2, or all data sets composing
//                                         clusters n1 and n2;
//    c) [[d]] or [[d1-d2]] or [[d1,d2]] : means discard data set d, or all data sets from d1 to d2, or data sets d1
//                                         and d2.
// - Added keyword "RUN" for AIMLESS.
//
// CHANGES IN VERSION 0.5.1 - 13/03/2014
// - Computing LCV for all nodes in the dendrogram. Top 5 nodes display LCV values in the dendrogram.
// - Added in log file proper CCP4 banner stuff and tags in log file to have tables and plots properly
//   handled by the new CCP4 log viewer.
// - Temporarily mended a problem in R function nparWilson, when matrix wpar has Na's. The present solution
//   is to interpolate these values using a 10-degree polynomial. This is not satisfactory and I will later
//   research more on missing data handling. 
// - Now the "merged_*.mtz" files are called, more appropriately, "unscaled_*.mtz"
// - Added 2 more keywords for AIMLESS: ANOMALOUS and SCALES. They have to be used similarly to those in AIMLESS.
// - Final sorting in MERGING_STATISTICS.info is on decreasing Completeness first, and increasing Rmeas second.
// - Modified / Added a few lines of code to allow code to run smoothly in both Windows and Linux / Os X systems.
// - All files including list of data (mtz or hkl) are now reported with absolute paths.
// - Added 1 keyword for AIMLESS: EXCLUDE BATCH. This has to be used similarly to the one in AIMLESS.
// - Blocked a few lines of code for the combination mode which stopped execution every time the user wanted
//   to run a same combination.
// - Added final merging statistics to log file (synthesis and combination modes) so that they may be viewed
//   with logview.
//
// CHANGES IN VERSION 0.5.0 - 12/02/2014
// - A new keyword, BLEND LAUEGROUP  [space group or laue group, POINTLESS style], has been added. This
//   allows the laue group of input XDS files to be decided by the user. If no LAUEGROUP line is used,
//   then a simple run of POINTLESS with the -c option will convert XDS files into MTZ files. If keyword
//   'LAUEGROUP AUTO' is used, then POINTLESS will decide which laue group to assign to the original XDS
//   files and convert these into MTZ files with the new space group (and eventually modified cell
//   parameters). Alternatively the user can impose his/her own choice for LAUEGROUP.
// - POINTLESS log file are saved in xds_files directory, so that one can check what went on at this stage.
// - LCV values for all merging nodes are tabulated in "CLUSTERS.txt" file.
//
// CHANGES IN VERSION 0.4.3 - 08/02/2014
// - Now keywords are given via stdin, like other CCP4 programs. The three different sections for keywords
//   are highlighted via presence of three keywords at the beginning of each line, BLEND, POINTLESS and
//   AIMLESS. For the first time, starting from this version, only 2 POINTLESS keywords (CHOOSE SPACEGROUP
//   and TOLERANCE) and 2 AIMLESS keywords (RESO and SDCORR) are allowed.
// - Added the "warn = -1" line to all R scripts. This should avoid warning issues contrasting with ccp4i.
// - Program does not crash anymore if alternative indexing leads to issues of differing space groups.
// CHANGES IN VERSION 0.4.2 - 25/01/2014
// - Fixed problem with execution of R in batch mode. This should now avoid the program to stop while is
//   being executed from within CCP4 graphical interface. Execution of R scripts is now carried out using
//   Rscript. Inside each R script q() with appropriate parameters have been included.
// - Added plotting of dendrogram and Rmeas_vs_Completeness plots in postscript formats to be viewed from
//   within the CCP4 graphical interface.
// - Now it (hopefully) catches all AIMLESS execution errors and does not stop. AIMLESS log dumped in merged_files
//   or combined_files directories even when AIMLESS fails.
// - Fixed a bug in "merge_mtzs" (modules blend2.R and blend3.R) needed to take care of alternate indexing. Now
//   also reference file goes through same procedure. POINTLESS won't stop because reference file is named
//   differently from current file in cluster or selection.
// - Eliminated cluster re-numbering in synthesis mode.
// - End-of-program messages eliminated.
// - python execution replaced with ccp4-python
// - Fixed a bug in the execution of POINTLESS to convert XDS files into MTZ files
//
// CHANGES IN VERSION 0.4.1  -  05/05/2013
// - Linear Cell Variation (LCV) parameter introduced and printed inside dendrogram picture
//
// NEW IN VERSION 0.4.0  -  13/11/2012
// - Option "-c", followed by a few integers, indicating serial numbers of datasets loaded. This option allows
//   users to create arbitrary clusters not included in the dendrogram.
// - POINTLESS log files are now stored in both merged_files and combined_files directories. If POINTLESS crashes
//   they will not be created.
//
//
// BLEND
// Multi-crystals pre-processing program. Reads several unmerged mtz files and returns one or more single merged and
// sorted files, ready to be processed by "aimless" (SCALA)
// J. Foadi & G. Evans
// Diamond Light Source and Imperial College London
// June 2012
//
// The program runs in 3 modes:
//    1) ANALYSIS - reads all mtz files and produces global statistics;
//
//       USAGE:		blend -a datasets_list.dat
//			blend -a /path/to/directory/with/integrated/data
//
//    2) SYNTHESIS - reads one or two heights chosen after looking at the dendrogram and produces merged data; 
//
//	 USAGE:		blend -s l1			if only one cutting level is chosen
// 			blend -s l1 l2			if two levels are chosen
//
//    3) COMBINATION - reads serial number of datasets to be included in merged dataset. Used to create groups of
//                     datasets not included in the dendrogram. 
//
//	 USAGE:		blend -c d1 d2 d3 d4 d5 ...	d1, d2, d3, etc are the serial number of the datasets that 
//                                                      the user would like to merge.
// 
// In analysis mode BLEND produces a dendrogram (image file of format PNG) and a summary table. BLEND only needs to be run once in analysis mode.
// It can be run as many times as wished in synthesis mode. After having run BLEND in analysis mode, the user decides on which merging nodes 
// in the dendrogram he/she is interested. This can be done in either synthesis or combination mode.
// Suppose all nodes below an height of l1 = 8.2 are deemed to be interesting. To reproduce merged files
// for all these nodes it is sufficient to type:
//
//			blend -s 8.2
//
// If, instead, all nodes between heights, say, 7.2 and 9.6 are deemed to be useful, then two values will be included in the command line:
//
//			blend -s 7.2 9.6
//
// All merged files and other material will, then , be produced and stored in directories with sequential numbers, like result_00, result_01, etc.
// Same goes for the combination mode. The only difference here is that the selection of datasets to be merged is decided by the user, and can be
// different rom what is found in the dendrogram.
//
/************************************************************************************************************/
/************************************************************************************************************/

/******************************************************/
// Standard C++ stuff
/******************************************************/
// <iostream> is included in Output.hh
// <string> is included in Output.hh
// <sstream> is included in Output.hh
#include <cstdlib>
#include <boost/algorithm/string.hpp>
using namespace boost::algorithm;


/******************************************************/
// Stuff inherited from Phil Evans
/******************************************************/
#include "Output.hh"
#include "version.hh"
#include "hkl_unmerge.hh"

/******************************************************/
// Other crystallographic libraries
/******************************************************/
#include <clipper/clipper.h>
#include <ccp4/ccp4_program.h>

/******************************************************/
// Stuff created by James Foadi
/******************************************************/
#include "aux_blend.hh"
#ifdef _WIN32
// for windows
 int Rscp = 1;
#else
// for other systems
 int Rscp = 0;
#endif



/******************************************************/
// MAIN
/******************************************************/
int main(int argc, char* argv[])
{
 // Initial stuff (CCP4 initialization, banner, etc)
 int nargc=1;
 CCP4::ccp4fyp(nargc, argv);      // To avoid ccp4fyp to complain about stuff like "-a" we only feed program name in ccp4 iniialization
 CCP4::ccp4ProgramName(BLEND_PROGRAM_NAME.c_str());
 std::string rcsdate = "$Date: "+std::string(BLEND_PROGRAM_DATE2)+"$";
 CCP4::ccp4RCSDate(rcsdate.c_str());
 CCP4::ccp4_prog_vers(BLEND_PROGRAM_VERSION.c_str());

 // Banner with HTML style
 std::cout << "<B><FONT COLOR='#FF0000'><!--SUMMARY_BEGIN-->" << std::endl;
 std::cout << std::endl;
 CCP4::ccp4_banner();
 std::cout << std::endl;
 std::cout << "<!--SUMMARY_END--></FONT></B>" << std::endl;

 // References section

 std::cout << std::endl;
 std::cout << "$TEXT:Reference1: $$ Main reference $$" << std::endl; 
 std::cout << "   'Clustering procedures for the optimal selection of data sets" << std::endl;
 std::cout << "    from multiple crystals in macromolecular crystallography':" << std::endl;
 std::cout << "    J. Foadi, P. Aller, Y. Alguel, A. Cameron, D. Axford, R.L. Owen, W. Armour, D.G. Waterman, S. Iwata and G. Evans, (2013)" << std::endl;
 std::cout << "    Acta Crystallogr. D69, 1617-1632" << std::endl;

 std::cout << std::endl;
 std::cout << "$$" << std::endl;
 std::cout << "$SUMMARY :Reference1:  $$ Blend: $$" << std::endl;
 std::cout << ":TEXT:Reference1: $$" << std::endl;
 std::cout << std::endl;
 std::cout << "$$" << std::endl;

 try
 {
  // Get BLEND_HOME environment variable (used to find R and Python code)
  const char *val;
  std::string sval="BLEND_HOME";
  val=std::getenv(sval.c_str());
  if (val == NULL && !std::getenv("CCP4"))
  {
   int nerr=14;
   throw nerr;
  }
  if (val != NULL && !std::getenv("CCP4"))
  {
   int nerr=15;
   throw nerr;
  }
  std::string home;
  if (val != NULL)
    home = val;
  else
    home = std::string(std::getenv("CCP4")) + "/share/blend";
  std::string R_program1 = home+"/R/blend1.R";
  std::string R_program2 = home+"/R/blend2.R";
  std::string R_program3 = home+"/R/blend3.R";
  std::string Python_program1 = home+"/python/create_file_for_BLEND.py";
  std::string Python_program2 = home+"/python/merge_clusters.py";
  std::string Python_program3 = home+"/python/xds_to_mtz_list.py";

  // Check comand line input is well formatted
  int runmode=0;
  float dlevel_top,dlevel_bottom;
  std::string mode_string,filename,cl_mode,cl_height;
  //std::vector<int> arbitrary_datasets;
  std::vector<std::string> arbitrary_datasets;
  if (argc == 1)
  {
   int nerr=1;
   throw nerr;
  }
  mode_string=argv[1];
  //if (argc < 3 | argc > 4)
  if (argc < 3)
  {
   int nerr=1;
   throw nerr;
  }
  //if (mode_string != "-a" & mode_string != "-s")     // First argument after program name "blend" needs to be either "-a" or "-s"
  if (mode_string != "-a" && mode_string != "-s" && mode_string != "-c")     // First argument after program name "blend" needs to be either "-a" 
                                                                           // or "-s", or "-c"
  {
   int nerr=1;
   throw nerr;
  }
  if (mode_string == "-a")     // Analysis pass
  {
   // In order to line up BLEND with the way ccp4i works (with stdin passed keywords) this is what has been added

   // Accepted keywords
   std::vector<std::string> akeys, dvkeywdline;
   akeys.push_back("NBIN");
   dvkeywdline.push_back("NBIN      20");
   akeys.push_back("RADFRAC");
   dvkeywdline.push_back("RADFRAC   0.750");
   akeys.push_back("ISIGI");
   dvkeywdline.push_back("ISIGI     1.500");
   akeys.push_back("CPARWT");
   dvkeywdline.push_back("CPARWT    1.000");
   akeys.push_back("DATAREF");
   akeys.push_back("LAUEGROUP");
   akeys.push_back("TOLERANCE");
   akeys.push_back("CHOOSE");
   akeys.push_back("RESO");
   akeys.push_back("SDCOR");
   akeys.push_back("ANOMALOUS");
   akeys.push_back("SCALES");
   akeys.push_back("EXCLUDE");
   akeys.push_back("RUN");

   // Load in keywords from standard input
   // This is a way of reading lines from standard input and storing them into a string variable, line,
   // until a "/n" (carriage return) is entered
   int jflag;
   std::string keywdline;
   std::vector<std::string> vkeywdline;
   std::string::size_type Idx;
   std::cout << ">>>>> Input command lines <<<<<" << std::endl;
   std::cout << std::endl;
   getline(std::cin, keywdline);
   trim(keywdline);
   std::cout << keywdline << std::endl;
   if (keywdline.substr(0,3) != "END" && keywdline.substr(0,3) != "GO" && keywdline != "")
   {
    jflag = 0;
    for (unsigned int i = 0; i < akeys.size(); ++i)
    {
     Idx = keywdline.find(akeys[i]);
     //if (Idx != std::string::npos)
     if (int(Idx) == 0)
     {
      vkeywdline.push_back(keywdline.substr(Idx));
      jflag = 1;
     }
    }
    if (jflag == 0) std::cout << "'" << keywdline << "' is not a valid keyword entry." << std::endl;
   }
   //while (keywdline != "END" && keywdline != "GO" && keywdline != "" && keywdline[keywdline.size() - 1] == 0)
   std::string old_keywdline = "";
   while (keywdline != "END" && keywdline != "GO" && keywdline != old_keywdline && keywdline != "")
   {
    old_keywdline = keywdline;
    getline(std::cin, keywdline);
    trim(keywdline);
    std::cout << keywdline << std::endl;
    if (keywdline.substr(0,3) != "END" && keywdline.substr(0,3) != "GO" && keywdline != "")
    {
     jflag = 0;
     for (unsigned int i = 0; i < akeys.size(); ++i)
     {
      Idx = keywdline.find(akeys[i]);
      //if (Idx != std::string::npos)
      if (int(Idx) == 0)
      {
       vkeywdline.push_back(keywdline.substr(Idx));
       jflag = 1;
      }
     }
    if (jflag == 0) std::cout << "'" << keywdline << "' is not a valid keyword entry." << std::endl;
    }
   }
   std::cout << ">>>>>     End of input    <<<<<" << std::endl;
   std::cout << std::endl;

   // Output to old-format BLEND_KEYWORDS.dat
   std::ofstream keywd_ostream("BLEND_KEYWORDS.dat",std::ios::out);
   keywd_ostream << "BLEND KEYWORDS" << std::endl;
   for (int i = 0; i < 4; ++i)
   {
    int k = 0;
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos)
     if (int(Idx) == 0)
     {
      k = 1;
      keywd_ostream << vkeywdline[j] << std::endl;
     }
    }
    if (k == 0) keywd_ostream << dvkeywdline[i] << std::endl;
   }
   for (int i = 4; i < 6; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream << "POINTLESS KEYWORDS" << std::endl;
   for (int i = 6; i < 8; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream << "AIMLESS KEYWORDS" << std::endl;
   for (int i = 8; i < 14; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream.close();

   // Recover LAUEGROUP if present (std::string lauegroup defined earlier)
   std::string lauegroup;
   lauegroup = "";
   for (unsigned i = 0; i < vkeywdline.size(); ++i)
   {
    Idx = vkeywdline[i].find(akeys[5]);
    //if (Idx != std::string::npos) lauegroup = vkeywdline[i];
    if (int(Idx) == 0) lauegroup = vkeywdline[i];
   }

   // Carry on checking correct command-line input
   if (argc != 3) {int nerr=1; throw nerr;}
   runmode=1;
   std::string tmpstring=argv[2];
   int icheck=isdir(tmpstring.c_str());
   if (icheck == 0)                 // Input mtz or xds files are listed in a file
   {
    int Python_status;
    std::ostringstream Python_command_line;
    std::cout << std::endl;
    std::cout << "Checking if there are xds files to be converted in input list " << tmpstring << " ..." << std::endl;
    Python_command_line << "ccp4-python " << Python_program3 << " " << tmpstring << " '" << lauegroup << "'";
    Python_status=std::system((Python_command_line.str()).c_str());
    if (Python_status != 0)
    {
     int nerr=16;
     throw nerr;
    }
    filename="mtz_names.dat";
    std::cout << "Done!" << std::endl;
   }
   if (icheck == 1)                 // Input mtz or xds files are included in a single directory
   {
    int Python_status;
    std::ostringstream Python_command_line;
    std::cout << std::endl;
    std::cout << "Checking if there are xds files to be converted in directory " << tmpstring << " ..." << std::endl;
    Python_command_line << "ccp4-python " << Python_program1 << " " << tmpstring << " '" << lauegroup << "'";
    Python_status=std::system((Python_command_line.str()).c_str());
    if (Python_status != 0)
    {
     int nerr=16;
     throw nerr;
    }
    filename="mtz_names.dat";
    std::cout << "Done!" << std::endl;
   }
   if (icheck == -1)
   {
    int nerr=10;
    throw nerr;
   }
  }
  if (mode_string == "-s")     // Synthesis pass
  {
   // In order to line up BLEND with the way ccp4i works (with stdin passed keywords) this is what has been added

   // Accepted keywords
   std::vector<std::string> akeys, dvkeywdline;
   akeys.push_back("NBIN");
   dvkeywdline.push_back("NBIN      20");
   akeys.push_back("RADFRAC");
   dvkeywdline.push_back("RADFRAC   0.750");
   akeys.push_back("ISIGI");
   dvkeywdline.push_back("ISIGI     1.500");
   akeys.push_back("CPARWT");
   dvkeywdline.push_back("CPARWT    1.000");
   akeys.push_back("DATAREF");
   akeys.push_back("LAUEGROUP");
   akeys.push_back("TOLERANCE");
   akeys.push_back("CHOOSE");
   akeys.push_back("RESO");
   akeys.push_back("SDCOR");
   akeys.push_back("ANOMALOUS");
   akeys.push_back("SCALES");
   akeys.push_back("EXCLUDE");
   akeys.push_back("RUN");

   // Load in keywords from standard input
   // This is a way of reading lines from standard input and storing them into a string variable, line,
   // until a "/n" (carriage return) is entered
   int jflag;
   std::string keywdline;
   std::vector<std::string> vkeywdline;
   std::string::size_type Idx;
   std::cout << ">>>>> Input command lines <<<<<" << std::endl;
   std::cout << std::endl;
   getline(std::cin, keywdline);
   trim(keywdline);
   std::cout << keywdline << std::endl;
   if (keywdline.substr(0,3) != "END" && keywdline.substr(0,3) != "GO" && keywdline != "")
   {
    jflag = 0;
    for (unsigned int i = 0; i < akeys.size(); ++i)
    {
     Idx = keywdline.find(akeys[i]);
     //if (Idx != std::string::npos)
     if (int(Idx) == 0)
     {
      vkeywdline.push_back(keywdline.substr(Idx));
      jflag = 1;
     }
    }
    if (jflag == 0) std::cout << "'" << keywdline << "' is not a valid keyword entry." << std::endl;
   }
   //while (keywdline != "END" && keywdline != "GO" && keywdline != "" && keywdline[keywdline.size() - 1] == 0)
   std::string old_keywdline = "";
   while (keywdline != "END" && keywdline != "GO" && keywdline != old_keywdline && keywdline != "")
   {
    old_keywdline = keywdline;
    getline(std::cin, keywdline);
    trim(keywdline);
    std::cout << keywdline << std::endl;
    if (keywdline.substr(0,3) != "END" && keywdline.substr(0,3) != "GO" && keywdline != "")
    {
     jflag = 0;
     for (unsigned int i = 0; i < akeys.size(); ++i)
     {
      Idx = keywdline.find(akeys[i]);
      //if (Idx != std::string::npos)
      if (int(Idx) == 0)
      {
       vkeywdline.push_back(keywdline.substr(Idx));
       jflag = 1;
      }
     }
    if (jflag == 0) std::cout << "'" << keywdline << "' is not a valid keyword entry." << std::endl;
    }
   }
   std::cout << ">>>>>     End of input    <<<<<" << std::endl;
   std::cout << std::endl;

   // Output to old-format BLEND_KEYWORDS.dat
   std::ofstream keywd_ostream("BLEND_KEYWORDS.dat",std::ios::out);
   keywd_ostream << "BLEND KEYWORDS" << std::endl;
   for (int i = 0; i < 4; ++i)
   {
    int k = 0;
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos)
     if (int(Idx) == 0)
     {
      k = 1;
      keywd_ostream << vkeywdline[j] << std::endl;
     }
    }
    if (k == 0) keywd_ostream << dvkeywdline[i] << std::endl;
   }
   for (int i = 4; i < 6; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream << "POINTLESS KEYWORDS" << std::endl;
   for (int i = 6; i < 8; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream << "AIMLESS KEYWORDS" << std::endl;
   for (int i = 8; i < 14; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream.close();

   // Carry on checking correct command-line input
   runmode=2;
   if (argc == 3)
   {
    if (std::atof(argv[2]) == 0) {int nerr=11; throw nerr;}
    dlevel_top=std::atof(argv[2]);
    dlevel_bottom=-100;
   }
   if (argc == 4)
   {
    if (std::atof(argv[3]) == 0) {int nerr=11; throw nerr;}
    dlevel_top=std::atof(argv[2]);
    dlevel_bottom=std::atof(argv[3]);
   }
   if (argc > 4)
   {
    int nerr=11;
    throw nerr;
   }
   if (dlevel_bottom > dlevel_top) {int nerr=12; throw nerr;}
   //std::cout << "DLEVEL_TOP = " << dlevel_top << ". DLEVEL_BOTTOM = " << dlevel_bottom << std::endl;
  }
  if (mode_string == "-c")     // Combination pass
  {
   // In order to line up BLEND with the way ccp4i works (with stdin passed keywords) this is what has been added

   // Accepted keywords
   std::vector<std::string> akeys, dvkeywdline;
   akeys.push_back("NBIN");
   dvkeywdline.push_back("NBIN      20");
   akeys.push_back("RADFRAC");
   dvkeywdline.push_back("RADFRAC   0.750");
   akeys.push_back("ISIGI");
   dvkeywdline.push_back("ISIGI     1.500");
   akeys.push_back("CPARWT");
   dvkeywdline.push_back("CPARWT    1.000");
   akeys.push_back("DATAREF");
   akeys.push_back("LAUEGROUP");
   akeys.push_back("TOLERANCE");
   akeys.push_back("CHOOSE");
   akeys.push_back("RESO");
   akeys.push_back("SDCOR");
   akeys.push_back("ANOMALOUS");
   akeys.push_back("SCALES");
   akeys.push_back("EXCLUDE");
   akeys.push_back("RUN");

   // Load in keywords from standard input
   // This is a way of reading lines from standard input and storing them into a string variable, line,
   // until a "/n" (carriage return) is entered
   int jflag;
   std::string keywdline;
   std::vector<std::string> vkeywdline;
   std::string::size_type Idx;
   std::cout << ">>>>> Input command lines <<<<<" << std::endl;
   std::cout << std::endl;
   getline(std::cin, keywdline);
   trim(keywdline);
   std::cout << keywdline << std::endl;
   if (keywdline.substr(0,3) != "END" && keywdline.substr(0,3) != "GO" && keywdline != "")
   {
    jflag = 0;
    for (unsigned int i = 0; i < akeys.size(); ++i)
    {
     Idx = keywdline.find(akeys[i]);
     //if (Idx != std::string::npos)
     if (int(Idx) == 0)
     {
      vkeywdline.push_back(keywdline.substr(Idx));
      jflag = 1;
     }
    }
    if (jflag == 0) std::cout << "'" << keywdline << "' is not a valid keyword entry." << std::endl;
   }
   //while (keywdline != "END" && keywdline != "GO" && keywdline != "" && keywdline[keywdline.size() - 1] == 0)
   std::string old_keywdline = "";
   while (keywdline != "END" && keywdline != "GO" && keywdline != old_keywdline && keywdline != "")
   {
    old_keywdline = keywdline;
    getline(std::cin, keywdline);
    trim(keywdline);
    std::cout << keywdline << std::endl;
    if (keywdline.substr(0,3) != "END" && keywdline.substr(0,3) != "GO" && keywdline != "")
    {
     jflag = 0;
     for (unsigned int i = 0; i < akeys.size(); ++i)
     {
      Idx = keywdline.find(akeys[i]);
      //if (Idx != std::string::npos)
      if (int(Idx) == 0)
      {
       vkeywdline.push_back(keywdline.substr(Idx));
       jflag = 1;
      }
     }
    if (jflag == 0) std::cout << "'" << keywdline << "' is not a valid keyword entry." << std::endl;
    }
   }
   std::cout << ">>>>>     End of input    <<<<<" << std::endl;
   std::cout << std::endl;

   // Output to old-format BLEND_KEYWORDS.dat
   std::ofstream keywd_ostream("BLEND_KEYWORDS.dat",std::ios::out);
   keywd_ostream << "BLEND KEYWORDS" << std::endl;
   for (int i = 0; i < 4; ++i)
   {
    int k = 0;
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos)
     if (int(Idx) == 0)
     {
      k = 1;
      keywd_ostream << vkeywdline[j] << std::endl;
     }
    }
    if (k == 0) keywd_ostream << dvkeywdline[i] << std::endl;
   }
   for (int i = 4; i < 6; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream << "POINTLESS KEYWORDS" << std::endl;
   for (int i = 6; i < 8; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream << "AIMLESS KEYWORDS" << std::endl;
   for (int i = 8; i < 14; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream.close();

   // Carry on checking correct command-line input
   runmode=3;

   // Turn arguments into integer numbers (in vector "arbitrary_datasets") to be later used in R code
   //for (int i=2;i < argc;i++) arbitrary_datasets.push_back(atoi(argv[i])); 
   // Turn arguments into a vector, "arbitrary_datasets", of characters to be later used in R code
   for (int i=2;i < argc;i++) arbitrary_datasets.push_back(argv[i]); 
  }

  // I like well-formatted output
  std::cout.setf(std::ios::fixed);

  // Amend command line arguments for CCP4 initialization
  //int ncorrect=0;
  //if (runmode == 1) ncorrect=3;
  //if (runmode == 2)
  //{
  // if (dlevel_bottom < 0) ncorrect=3;
  // if (dlevel_bottom >= 0) ncorrect=4;
  //}
  //if (runmode == 3) ncorrect=arbitrary_datasets.size()+2;
  //argc-=ncorrect;
  //for (int i=0;i < argc;i++) argv[i]=argv[i+1]; 

  // Analysis mode
  if (runmode == 1)
  {
   // Name of ascii file containing list of mtz files
   std::string tmpstringa;
   std::vector<std::string> mtzin_name;

   // Start CCP4 before anything else.
   //CCP4::ccp4fyp(argc,argv);
   //CCP4::ccp4_banner();
   std::cout << std::endl;
   std::cout << "<B><FONT COLOR='#FF0000'><!--SUMMARY_BEGIN-->" << std::endl;
   std::cout << std::endl;
   std::cout << "You are now running BLEND in analysis mode." << std::endl; 
   std::cout << std::endl;
   std::cout << "<!--SUMMARY_END--></FONT></B>" << std::endl;
   std::cout << std::endl;
   
   // Load crystals in unmerged data structures
   std::vector<scala::hkl_unmerge_list> hkl_list=load_crystals(filename,runmode);   // This is the correct expression for a copy constructor. Defining hkl_list first
                                                                                    // and then using hkl_list=load_crystals(filename) doesn't work. Ultimately this is
                                                                                    // due to the way the hkl_unmerged_list class has been implemented by Phil.
   // Label crystals according to dataset validity.
   // crystal_flag = 0:   crystal is accepted for further processing
   // crystal_flag = 1:   crystal is rejected because data file contains no reflections
   // crystal_flag = 2:   crystal is rejected because data file is made up of multiple runs
   // crystal_flag = 3:   crystal is rejected because data file is made up of multiple datasets
   
   // Initially all crystals are assumed to have valid datasets
   std::vector<int> crystal_flag(hkl_list.size(),0);

   // Label crystals
   crystal_flag=label_crystals(hkl_list,crystal_flag);
   //for (int i=0;i < crystal_flag.size();i++) std::cout << "crystal flag " << i+1 << ": " << crystal_flag[i] << std::endl;
   
   // Maps to store types of crystals belonging to different bravais lattices. These maps contain the "sg" bit as initially we were
   // handling space groups, rather than bravais lattices. Bravais lattice numbers are as follows:
   //
   // Triclinic   :     aP       1
   // Monoclinic  :     mP       2
   //                   mS       3
   // Orthorombic :     oP       4
   //                   oS       5
   //                   oF       6
   //                   oI       7
   // Tetragonal  :     tP       8
   //                   tI       9
   // Hexagonal   :     hP      10
   // Rhombohedral:     hH      11   (it is also possible hR, but we haven't implemented this way. To be sorted later)
   // Cubic       :     cP      12
   //                   cF      13
   //                   cI      14
   //
   // sg_to_crystal:       bravais lattice number  ->  crystal serial number
   // crystal_to_sg:    crystal serial number  ->  bravais lattice number
   std::cout << "Partitioning datasets into groups with same bravais lattice ........." << std::endl;
   std::map<std::string,int> bl_symbol_to_number;
   std::map<std::string,int>::iterator pos_bl;
   bl_symbol_to_number.insert(std::make_pair<std::string,int>("aP",1));
   bl_symbol_to_number.insert(std::make_pair<std::string,int>("mP",2));
   bl_symbol_to_number.insert(std::make_pair<std::string,int>("mS",3));
   bl_symbol_to_number.insert(std::make_pair<std::string,int>("oP",4));
   bl_symbol_to_number.insert(std::make_pair<std::string,int>("oS",5));
   bl_symbol_to_number.insert(std::make_pair<std::string,int>("oF",6));
   bl_symbol_to_number.insert(std::make_pair<std::string,int>("oI",7));
   bl_symbol_to_number.insert(std::make_pair<std::string,int>("tP",8));
   bl_symbol_to_number.insert(std::make_pair<std::string,int>("tI",9));
   bl_symbol_to_number.insert(std::make_pair<std::string,int>("hP",10));
   //bl_symbol_to_number.insert(std::make_pair<std::string,int>("hR",11));
   bl_symbol_to_number.insert(std::make_pair<std::string,int>("hH",11));
   bl_symbol_to_number.insert(std::make_pair<std::string,int>("cP",12));
   bl_symbol_to_number.insert(std::make_pair<std::string,int>("cF",13));
   bl_symbol_to_number.insert(std::make_pair<std::string,int>("cI",14));
   std::map<int,int> crystal_to_sg;
   std::multimap<int,int> sg_to_crystal;
   std::multimap<int,int>::iterator pos_cs;
   scala::hkl_symmetry symmetry;
   scala::CrystalType crystal_type;
   clipper::Spacegroup cspacegroup;
   clipper::Spgr_descr cspgr_descr;
   //size_t lnum;
   int lnum;
   std::string blatsymb;
   for (unsigned int i=0;i < hkl_list.size();i++)
   {
    if (crystal_flag[i] == 0)
    {
     symmetry=hkl_list[i].symmetry();
     crystal_type=CrystalType(symmetry);
     blatsymb=crystal_type.format();
     lnum=blatsymb.find("A");
     if (lnum != -1)
     {
      blatsymb.erase(lnum,1);
      blatsymb.insert(lnum,"S");
     }
     lnum=blatsymb.find("B");
     if (lnum != -1)
     {
      blatsymb.erase(lnum,1);
      blatsymb.insert(lnum,"S");
     }
     lnum=blatsymb.find("C");
     if (lnum != -1)
     {
      blatsymb.erase(lnum,1);
      blatsymb.insert(lnum,"S");
     }
     pos_bl=bl_symbol_to_number.find(blatsymb);
     //std::cout << "Lattice type: " << blatsymb << " Number is: " << pos_bl->second << std::endl;
     cspgr_descr=clipper::Spgr_descr(symmetry.symbol_Hall());
     cspacegroup=clipper::Spacegroup(cspgr_descr);
     //sg_to_crystal.insert(std::make_pair<int,int>(cspacegroup.spacegroup_number(),i));
     //crystal_to_sg.insert(std::make_pair<int,int>(i,cspacegroup.spacegroup_number()));
     sg_to_crystal.insert(std::make_pair<int,int>(pos_bl->second,i));
     crystal_to_sg.insert(std::make_pair<int,int>(i,pos_bl->second));
    }
    else
    {
     sg_to_crystal.insert(std::make_pair<int,int>(0,i));
     crystal_to_sg.insert(std::make_pair<int,int>(i,0));
    }
   }

   // Divide crystals in groups having same spacegroup
   std::vector<int> spacegroup_class;
   for (pos_cs=sg_to_crystal.begin();pos_cs != sg_to_crystal.end();++pos_cs)
   {
    if (pos_cs == sg_to_crystal.begin())
    {
     spacegroup_class.push_back(pos_cs->first);
    } 
    else
    {
     int flag=0;
     for (unsigned int i=0;i < spacegroup_class.size();i++)
     {
      if (pos_cs->first == spacegroup_class[i]) flag=1;
     }
     if (flag == 0) spacegroup_class.push_back(pos_cs->first);
    }
   }
   if (spacegroup_class.size() == 1)
   {
    std::cout << "These datasets are partitioned in " << spacegroup_class.size() << " group." << std::endl;
   }
   else
   {
    std::cout << "These datasets are partitioned in " << spacegroup_class.size() << " groups." << std::endl;
   } 
   for (unsigned int i=0;i < spacegroup_class.size();i++)
   {
    if (spacegroup_class[i] == 0) 
    {
     std::cout << "The following crystals have invalid datasets: ";
     for (pos_cs=sg_to_crystal.begin();pos_cs != sg_to_crystal.end();++pos_cs)
     {
      if (pos_cs->first == spacegroup_class[i]) std::cout << pos_cs->second+1 << "  ";
     }
    }
    else
    {
     std::cout << "The following crystals are part of bravais lattice n. " << spacegroup_class[i] << ": ";
     for (pos_cs=sg_to_crystal.begin();pos_cs != sg_to_crystal.end();++pos_cs)
     {
      if (pos_cs->first == spacegroup_class[i]) std::cout << pos_cs->second+1 << "  ";
     }
    }
    std::cout << std::endl; 
   }

   // Output summary table
   std::cout << "Building summary table listing all crystals ........." << std::endl;
   output_summary_table(hkl_list,sg_to_crystal,spacegroup_class,crystal_flag);

   // Output ascii files for R (only if there are at least 2 crystals)
   if (hkl_list.size() > 1)
   {
    std::cout << "Building R files ........." << std::endl;
    statistics_with_R(hkl_list,sg_to_crystal,spacegroup_class,crystal_flag,R_program1,Rscp);
   }

   // Delete BLEND_KEYWORDS.dat before termination
   if (remove("BLEND_KEYWORDS.dat") != 0) std::cout << "Warning! Unable to remove file 'BLEND_KEYWORDS.dat'." << std::endl;

   // End of BLEND - Analysis mode
   std::cout << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << "## BLEND - Normal Termination                                   ##" << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << std::endl;
  }

  // Synthesis mode
  if (runmode == 2)
  {
   // Start CCP4 before anything else.
   //CCP4::ccp4fyp(argc,argv);
   //CCP4::ccp4_banner();
   std::cout << std::endl;
   std::cout << "<B><FONT COLOR='#FF0000'><!--SUMMARY_BEGIN-->" << std::endl;
   std::cout << std::endl;
   std::cout << "You are now running BLEND in synthesis mode." << std::endl; 
   std::cout << std::endl;
   std::cout << "<!--SUMMARY_END--></FONT></B>" << std::endl;
   std::cout << std::endl;

   // Run R code
   std::cout << std::endl;
   std::cout << "Running R code to read statistical analysis from a previous run of BLEND and produce information on clusters..." << std::endl; 
   std::cout << std::endl;
   int R_status;
   std::ostringstream R_command_line;
   //R_command_line << "R --vanilla --slave --quiet < " << R_program2 << " --args " << dlevel_top << " " << dlevel_bottom;
   if (Rscp == 0)
   {
    R_command_line << "Rscript " << R_program2 << " " << dlevel_top << " " << dlevel_bottom;
   }
   else
   {
    R_command_line << "Rscript.exe " << R_program2 << " " << dlevel_top << " " << dlevel_bottom;
   }
   R_status=std::system((R_command_line.str()).c_str());
   if (R_status != 0)
   {
    int nerr=13;
    throw nerr;
   }

   // Delete BLEND_KEYWORDS.dat before termination
   if (remove("BLEND_KEYWORDS.dat") != 0) std::cout << "Warning! Unable to remove file 'BLEND_KEYWORDS.dat'." << std::endl;

   // End of BLEND - Synthesis mode
   std::cout << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << "## BLEND - Normal Termination                                   ##" << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << std::endl;
  }

  // Combination mode
  if (runmode == 3)
  {
   // Start CCP4 before anything else.
   //CCP4::ccp4fyp(argc,argv);
   //CCP4::ccp4_banner();
   std::cout << std::endl;
   std::cout << "<B><FONT COLOR='#FF0000'><!--SUMMARY_BEGIN-->" << std::endl;
   std::cout << std::endl;
   std::cout << "You are now running BLEND in combination mode." << std::endl; 
   std::cout << std::endl;
   std::cout << "<!--SUMMARY_END--></FONT></B>" << std::endl;
   std::cout << std::endl;

   // Run R code
   std::cout << std::endl;
   std::cout << "Running R code to read statistical analysis from a previous run of BLEND and produce information on clusters..." << std::endl; 
   std::cout << std::endl;
   int R_status;
   std::ostringstream R_command_line,tmp;
   for (unsigned int i=0;i < arbitrary_datasets.size();i++) tmp << arbitrary_datasets[i] << " ";
   //R_command_line << "R --vanilla --slave --quiet < " << R_program3 << " --args " << tmp.str();
   if (Rscp == 0)
   {
    R_command_line << "Rscript " << R_program3 << " " << tmp.str();
   }
   else
   {
    R_command_line << "Rscript.exe " << R_program3 << " " << tmp.str();
   }
   //std::cout << R_command_line.str() << std::endl;
   R_status=std::system((R_command_line.str()).c_str());
   if (R_status != 0)
   {
    int nerr=13;
    throw nerr;
   }

   // Delete BLEND_KEYWORDS.dat before termination
   if (remove("BLEND_KEYWORDS.dat") != 0) std::cout << "Warning! Unable to remove file 'BLEND_KEYWORDS.dat'." << std::endl;

   // End of BLEND - Synthesis mode
   std::cout << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << "## BLEND - Normal Termination                                   ##" << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << std::endl;
  }
 }
 catch (int nerr)
 {
  if (nerr == 1)
  {
   std::cerr << "\n BLEND ERROR!\n"
             << "Wrong command line format: too few, or too many arguments. Correct format is:\n"
             << "                                  \n"
             << "   blend -a name_of_file.dat                                  (analysis mode)\n"
             << "                 or               \n"
             << "   blend -a /path/to/directory                                (analysis mode)\n"
             << "                 or               \n"
             << "   blend -s l1 (numeric height in dendrogram)                (synthesis mode)\n"
             << "                 or               \n"
             << "   blend -s l1 l2 (numeric heights in dendrogram)            (synthesis mode)\n"
             << "                 or               \n"
             << "   blend -c d1 d2 d3 ... (serial number of datasets)       (combination mode)\n" << std::endl;
   return EXIT_FAILURE;
  }
  if (nerr == 7)
  {
   std::cerr << "\n BLEND ERROR!\n"
             << "File ... does not exist" << std::endl;
   return EXIT_FAILURE;
  }
  if (nerr == 8)
  {
   std::cerr << "\n BLEND ERROR!\n"
             << "No valid MTZ files are listed in input file, or included in input directory." << std::endl;
   return EXIT_FAILURE;
  }
  if (nerr == 9)
  {
   std::cerr << "\n USER ERROR!\n"
             << "The specific MTZ file cannot be read." << std::endl;
  }
  if (nerr == 10)
  {
   std::cerr << "\n USER ERROR!\n"
             << "After -a you need to input a valid file name or directory name." << std::endl;
  }
  if (nerr == 11)
  {
   std::cerr << "\n USER ERROR!\n"
             << "After -s you need to input one or two valid numbers indicating levels in the dendrogram." << std::endl;
  }
  if (nerr == 12)
  {
   std::cerr << "\n USER ERROR!\n"
             << "Bottom level selected in dendrogram needs to be smaller than top level." << std::endl;
  }
  if (nerr == 13)
  {
   std::cerr << "\n EXECUTION ERROR!\n"
             << "An error occurred in the execution of R code associated with BLEND. Please, report this to J. Foadi" << std::endl;
  }
  if (nerr == 14)
  {
   std::cerr << "\n USER ERROR!\n"
             << "Environment variables BLEND_HOME and CCP4 have not been set up." << std::endl;
  }
  if (nerr == 15)
  {
   std::cerr << "\n USER ERROR!\n"
             << "Environment variable CCP4 has not been set up." << std::endl;
  }
  if (nerr == 16)
  {
   std::cerr << "\n EXECUTION ERROR!\n"
             << "An error occurred in the execution of Python code associated with BLEND. Please, report this to J. Foadi" << std::endl;
  }
 }

 return 0;
}
