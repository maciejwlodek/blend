/************************************************************************************************************/
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
// CHANGES IN VERSION 0.6.23  -  09/01/2018
// - Set crystal_flag=3 only when the MTZ file contains multiple *non-empty*
//   datasets. This works around a xia2 issue in which some MTZ files have an
//   extra empty dataset that can be ignored.
// CHANGES IN VERSION 0.6.22  -  16/09/2016
// - Corrected a bug in function "pruning_plan", affecting runs in combination mode.
//   (modules blend.cpp,version.hh,blend3.R)
// CHANGES IN VERSION 0.6.21  -  18/04/2016
// - Keyword DATAREF has been chenged into DREF.
// - DREF is now only used if the user wants to assign same space group as the one of
//   the reference file. Alternative indexing, when a reference file is not present, is
//   done with respect to the sweep-longest file within a cluster.
// - Fixed a difference in writing AIMLESS keyword RESO LOW HIGH that made synthesis output
//   slightly different from combination output. Now, when RESO keyword is not assigned by
//   user, only RESO HIGH is written out.
// - Eliminated the spacing between characters N:1,2,3 for plotting tables. Now BLEND log
//   can also be viewed with loggraph.
//   (modules blend.cpp,aux_blend.cpp,version.hh,blend1.R,blend2.R,blend3.R)
// CHANGES IN VERSION 0.6.20  -  04/04/2016
// - Fixed annoying halting of the program in endless loop when asking for dendrogram-type
//   graphics (graphics mode). Now the deepest possible level is always selected even when
//   the user chooses an higher level. Output in "graphics" directory has also level index
//   limited to deepest level, no matter what user decides.
// - Increased formatting space for I and sigI in "refs_*_*.dat" files.
//   (modules blend4.R,aux_blend.cpp)
// CHANGES IN VERSION 0.6.19  -  29/02/2016
// - Sorted out problem with recognition of the I2 setting of space group C2.
// - Now LAUEGROUP with specific symmetry symbol calls CHOOSE LAUEGROUP for POINTLESS.
// - Fixed a bug that stopped program from runnning in both cF and cP variants because one
//   of AIMLESS log tables changed from 13 to 11 columns if data have too few images for
//   proper smoothing.
// - Selection of the interested line in AIMLESS log with "Rmerge" was failing in some cases.
//   The selection string is now more robust and does not seem to be failing anymore.
// - Added the ability to add two columns with furthest datasets in CLUSTERS.txt, when esecution
//   is with the -aDO variant. Also, now the furthest datasets are in increasing order.
// - Fixed a few typos and included list of discarded files who are not unmerged mtz.
//   (modules blend.cpp,aux_blend.cpp,version.hh,xds_to_mtz_list.py,create_file_for_BLEND.py,
//            blend0.R,blend1.R,blend3.R)
// CHANGES IN VERSION 0.6.18  -  30/11/2015
// - Eliminated "pruning" word from log sentence "pruning cycles".
//   (module blend3.R)
// CHANGES IN VERSION 0.6.17  -  28/11/2015
// - Fixed a bug responsible for wrong output of the level of certain clusters.
//   (module blend4.R)
// CHANGES IN VERSION 0.6.16  -  20/11/2015
// - Added missing code for displaying correctly aLCV annotated dendrograms in PS format.
//   (module blend4.R)
// CHANGES IN VERSION 0.6.15  -  17/11/2015
// - The graphics mode now includes a "DO" option (in addition to the "D") option. This is meant to produce
//   an annotated dendrogram with only the aLCV values even when the synthesis mode has been executed.
//   Furthermore, the dendrogram with the aLCV annotations is called "aLCV_annotated*", while the one
//   including all statistics is called "stats_annotated*".
//   (module blend4.R)
// CHANGES IN VERSION 0.6.14  -  12/08/2015
// - Fixed a bug connected to mode "-cF" that wrote wrong cominations in GROUPS.info
//   (module blend3.R)
// CHANGES IN VERSION 0.6.13  -  30/07/2015
// - Fixed inconsistencies in keywords length. Now all keywords can also be typed as 4-initial-letters.
// - Added variant "-cP" to combination mode. Input groups of data follow the usual syntax. Based on the first
//   run of AIMLESS on these group, individual runs are progressively cut to obtain better statistics in an
//   automated way.
// - Added variant "-cF" to combination mode. Whole data sets are filtered out one after the other until
//   default or user-selected completeness is reached, or until a maximum number of cycles is reached. In the
//   end the result with lowest Rpim is chosen for output.
//   (modules blend.cpp blend1.R blend2.R and blend3.R)
// CHANGES IN VERSION 0.6.12  -  16/07/2015
// - When only 1 dataset is to be scaled now use EXCLUDE BATCH ...", rather than "EXCLUDE FILE 1 BATCH ... "
//   in POINTLESS. This cause POINTLESS to ignore the cuts.
//   (modules blend2.R and blend3.R)
// CHANGES IN VERSION 0.6.11  -  15/07/2015
// - Added aLCV values to annotated dendrogram, if no "BLEND.RMergingStatistics" is found in working directory.
//   In this case the nodes in the annotated dendrogram have square, rather than circular shape.
//   (modules affected blend4.R)
// CHANGES IN VERSION 0.6.10  -  14/07/2015
// - Now aLCV has priority over LCV as the first one is directly connected to resolution (both in A).
//   (modules blend0.R and blend1.R)
// CHANGES IN VERSION 0.6.9  -  28/05/2015
// - Fixed a bug concerning interpolation of NA's in dynamic wilson plots for secondary descriptors
// - Added print message while loading datasets to find out which one is not including data
//   (modules aux_blend.cpp, blend1.R).
// CHANGES IN VERSION 0.6.8  -  20/05/2015
// - Fixed a "nasty" (and wrong piece of code!) bug in routine okRadDam that would make data sets with
//   less than 6 images still undergo radiation damage routine (module blend1.R).
// CHANGES IN VERSION 0.6.7  -  15/04/2015
// - Fixed a bug responsible for the deletion of any file whose name includes "reference". This is not
//   necessary (modules blend2.R and blend3.R).
// CHANGES IN VERSION 0.6.6  -  01/04/2015
// - Added two columns in file "CLUSTERS.txt" which tell two datasets having current values of 
//   LCV and aLCV (module blend1.R).
// - Fixed bug for mode "-g D". Now data sets at the bottom are printed correctly (module blend4.R).
// CHANGES IN VERSION 0.6.5  -  24/03/2015
// - Included a few print lines in log file to allow users to understand why a file has been rejected
//   (modules blend.cpp and aux_blend.cpp).
// - Fixed a bug in the radiation damage routine. Also now the estimate is possible even in those cases
//   where intensities are negative (simply raise temporarily plot to positive values) (module blend1.R).
// - Fixed a bug in the resolution interpolation procedure (module blend1.R).
// CHANGES IN VERSION 0.6.4  -  16/03/2015
// - Added error catch code to module "blend2.R" and "blend3.R", acting when BLEND.RData has not
//   been found.
// CHANGES IN VERSION 0.6.3  -  12/03/2015
// - Fixed bug in "blend4.R" to compute annotated dendrogram when only part of dendrogram's nodes
//   have been merged and scaled.
// - Added PS graphs and amended PS dendrogram (now it is vertical, rather than horizontal).
// - Added a new, quick "dendrogram-only" mode (-aDO) to simply read cell parameters and return
//   dendrogram and a few files.
// CHANGES IN VERSION 0.6.2  -  09/03/2015
// - Graphics (D) mode can now also be used before running synthesis. It will zoom in dendrogram
//   without returning statistics, completeness and resolution.
// CHANGES IN VERSION 0.6.1  -  03/03/2015
// - Added new mode "graphical", with the purpose of producing annotated dendrograms (-g D).
//   New PNG annotated dendrograms are stored in directory "graphics". All options of graphics
//   mode are included in R code "blend4.R". Module "blend2.R" has also been slightly modified.
// CHANGES IN VERSION 0.5.14 - 13/01/2015
// - Fixed a silly code-syntax error accidentally inserted in previous changes (module blend3.R).
// CHANGES IN VERSION 0.5.13 - 12/01/2015
// - Changed the "No result..." message in both modules blend2.R and blend3.R so that now stopping can
//   be due to either POINTLESS or AIMLESS and one can find out by exploring content of directory "merged_files"
//   or "combined_files".
// CHANGES IN VERSION 0.5.12 - 23/12/2014
// - Added check to ascertain whether dataref file exists. If not the program exits with error code (files
//   blend2.R and blend3.R).
// - Bug fixes: DATAREF keyword line was all made uppercase, which resulted in reference file path to be
//   changed. Now only keyword DATAREF (and not all line) is turned to uppercase.
// - Reference to J. Foadi for support has been eliminated. Now users should report or contact CCP4-related
//   people or bullettin board.
// CHANGES IN VERSION 0.5.11 - 14/12/2014
// - Fixed a bug in keywords management. Now accepting both uppercase and lowercase.
// - Fixed a bug for the calculation of aLCV (module blend1.R)
// CHANGES IN VERSION 0.5.10 - 12/11/2014
// - Now only one POINTLESS run and one AIMLESS run are needed for each cluster or group. All input
//   and output files are retained (no more need of "new_" and "final_" files). Input, output and
//   keywords for both POINTLESS and AIMLESS runs can be seen inside each POINTLESS or AIMLESS log.
//   Thus, there is no more need of an aimless_keywords.dat file. Easier to reproduce POINTLESS or
//   AIMLESS runs now
// - Fixed a bug responsible for re-writing RESOLUTION keyword twice in aimless_keywords.dat
// - Fixed order with which data set are fed to POINTLESS and AIMLESS inside "blend2.R". Now
//   synthesis output gives same results of combination output corresponding to same cluster
// - Replaced REBATCH with POINTLESS for images removal
// - Data sets in clusters or groups are now always listed in increasing order 
// - Detecting when POINTLESS fails for cell discrepancy (TOLERANCE) and enabling the program to output
//   the log file
// - New possibilities in synthesis mode: -sLCV and -saLCV which allows to select based on LCV
//   or aLCV values, rather than based on cluster height
// CHANGES IN VERSION 0.5.9 - 09/10/2014
// - Small change in code for both blend2.R and blend3.R to reflect new resolution criteria
//   in AIMLESS for CC1/2 and Mn(I/sd)
// CHANGES IN VERSION 0.5.8 - 08/10/2014
// - Further bug fixes for aLCV calculation
// CHANGES IN VERSION 0.5.7 - 07/10/2014
// - Fixed bug for aLCV calculation
// CHANGES IN VERSION 0.5.6 - 03/10/2014
// - Added the new aLCV, i.e. the linear cell variation in absolute (angstroms) value
// - Chenged part of the documentation and tutorials to reflect aLCV addition
// CHANGES IN VERSION 0.5.5 - 07/08/2014
// - Fixed a bug in both blend2.R and blend3.R in relation to altidx; added "TOLERANCE 1000" so that now
//   POINTLESS does not fail just for alternative indexing.
// CHANGES IN VERSION 0.5.4 - 04/07/2014
// - Fixed a bug in blend3.R. Now if [[n]] is not contained in the group of data sets to combine, it just
//   ignores it.
// CHANGES IN VERSION 0.5.3 - 14/06/2014
// - Added condition on images and resolution to skip radiation damage procedure: the procedure is not carried
//   out if there are less than 6 images and / or if the highest resolution is less than 5 angstroms. Radiation
//   damage procedure is also skipped if keyword RADFRAC is less than 0.000001.
// - Fixed a bug which crashed the program when POINTLESS failed to assign SG. It relates to CC12 and Mn2.
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
  std::string R_program0 = home+"/R/blend0.R";
  std::string R_program0NC = home+"/R/blend0NC.R";
  std::string R_program0SF = home+"/R/blend0SF.R";
  std::string R_program0S6 = home+"/R/blend0S6.R";
  std::string R_programK = home+"/R/blendK.R";
  std::string R_programKNC = home+"/R/blendKNC.R";
  std::string R_programKS6 = home+"/R/blendKS6.R";
  std::string R_programKSF = home+"/R/blendKSF.R";
  std::string R_program1 = home+"/R/blend1.R";
  std::string R_program1NC = home+"/R/blend1NC.R";
  std::string R_program1SF = home+"/R/blend1SF.R";
  std::string R_program1S6 = home+"/R/blend1S6.R";
  std::string R_program2 = home+"/R/blend2.R";
  std::string R_program3 = home+"/R/blend3.R";
  std::string R_program4 = home+"/R/blend4.R";
  std::string Python_program1 = home+"/python/create_file_for_BLEND.py";
  std::string Python_program2 = home+"/python/merge_clusters.py";
  std::string Python_program3 = home+"/python/xds_to_mtz_list.py";

  // Check comand line input is well formatted
  int runmode=0;
  int addrunmode2 = 0; // mode_string "-sLCV" gives addrunmode2 = 1; mode_string "-saLCV" gives addrunmode2 = 2
  int addrunmode3 = 0; // mode_string "-c" gives addrunmode3 = 0; mode_string "-cP" gives addrunmode3 = 1
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
  if (argc < 3)
  {
   int nerr=1;
   throw nerr;
  }
  if (mode_string != "-a" && mode_string != "-aDO" &&                              // First argument after program name "blend" needs to be either "-a" or "-aDO"
      mode_string != "-aNC" && mode_string != "-aDONC" &&                          // possibly with NC appended
      mode_string != "-aSF" && mode_string != "-aDOSF" &&                          // possibly with SF appended
      mode_string != "-aS6" && mode_string != "-aDOS6" &&                          // possibly with S6 appended
      mode_string != "-k" && mode_string != "-kNC" &&                              // or "-k" for kmeans mode
      mode_string != "-kS6" && mode_string != "-kSF" &&                            // possibly with NC, SF, or S6 appended
      mode_string != "-s" && mode_string != "-sLCV" && mode_string != "-saLCV" &&  // or "-s" or "-sLCV" or "-saLCV"
      mode_string != "-sNC" && mode_string != "-sLCVNC" && mode_string != "-saLCVNC" &&  // or "-sNC" or "-sLCVNC" or "-saLCVNC"
      mode_string != "-sSF" && mode_string != "-sLCVSF" && mode_string != "-saLCVSF" &&  // or "-sSF" or "-sLCVSF" or "-saLCVSF"
      mode_string != "-sS6" && mode_string != "-sLCVS6" && mode_string != "-saLCVS6" &&  // or "-sS6" or "-sLCVS6" or "-saLCVS6"
      mode_string != "-c" && mode_string != "-cP" && mode_string != "-cF" &&       // or "-c" or "-cP" or "-cF"
      mode_string != "-cNC" && mode_string != "-cPNC" && mode_string != "-cFNC" && // or "-cNC" or "-cPNC" or "-cFNC"
      mode_string != "-cSF" && mode_string != "-cPSF" && mode_string != "-cFSF" && // or "-cSF" or "-cPSF" or "-cFSF"
      mode_string != "-cS6" && mode_string != "-cPS6" && mode_string != "-cFS6" && // or "-cS6" or "-cPS6" or "-cFS6"
      mode_string != "-gNC" &&                                                     // or "-gNC
      mode_string != "-gSF" &&                                                     // or "-gSF
      mode_string != "-gS6" &&                                                     // or "-gS6
      mode_string != "-g")                                                         // or "-g"
  {
   int nerr=1;
   throw nerr;
  }
  if (mode_string == "-a" 
      || mode_string == "-aDO" 
      || mode_string == "-aNC" 
      || mode_string == "-aDONC" 
      || mode_string == "-aSF" 
      || mode_string == "-aDOSF"
      || mode_string == "-aS6" 
      || mode_string == "-aDOS6"
      || mode_string == "-k"
      || mode_string == "-kNC"
      || mode_string == "-kS6"
      || mode_string == "-kSF")     // Analysis pass
  {
   // In order to line up BLEND with the way ccp4i works (with stdin passed keywords) this is what has been added

   // Accepted keywords
   std::vector<std::string> akeys, dvkeywdline;
   akeys.push_back("NBIN");
   dvkeywdline.push_back("NBIN      20");
   akeys.push_back("RADF");
   dvkeywdline.push_back("RADFRAC   0.750");
   akeys.push_back("ISIG");
   dvkeywdline.push_back("ISIGI     1.500");
   akeys.push_back("CPAR");
   dvkeywdline.push_back("CPARWT    1.000");
   if(mode_string == "-k" || mode_string == "-kNC" || mode_string == "-kS6" || mode_string == "-kSF") {
        akeys.push_back("K");
        dvkeywdline.push_back("K         3");
   }
   //akeys.push_back("DATA");
   akeys.push_back("DREF");
   akeys.push_back("LAUE");
   akeys.push_back("MAXC");
   akeys.push_back("COMP");
   akeys.push_back("CUTF");
   akeys.push_back("TOLE");
   akeys.push_back("CHOO");
   akeys.push_back("RESO");
   akeys.push_back("SDCO");
   akeys.push_back("ANOM");
   akeys.push_back("SCAL");
   akeys.push_back("EXCL");
   akeys.push_back("RUN");

   // Load in keywords from standard input
   // This is a way of reading lines from standard input and storing them into a string variable, line,
   // until a "/n" (carriage return) is entered
   int jflag;
   std::string keywdline, keywdline2;
   std::vector<std::string> split_vec;
   std::vector<std::string> vkeywdline;
   std::string::size_type Idx;
   std::cout << ">>>>> Input command lines <<<<<" << std::endl;
   std::cout << std::endl;
   getline(std::cin, keywdline);
   trim(keywdline);
   split(split_vec,keywdline,is_any_of(" "));
   keywdline2 = split_vec[0];
   to_upper(keywdline2);
   //if (keywdline2.substr(0,4) == "DATA")
   if (keywdline2.substr(0,4) == "DREF")
   {
    keywdline = keywdline2 + " ";
    for (unsigned int i = 1; i < split_vec.size(); ++i)
    {
     keywdline = keywdline + split_vec[i];
    }
   }
   else
   {
    to_upper(keywdline);
   }
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
    split(split_vec,keywdline,is_any_of(" "));
    keywdline2 = split_vec[0];
    to_upper(keywdline2);
    //if (keywdline2.substr(0,4) == "DATA")
    if (keywdline2.substr(0,4) == "DREF")
    {
     keywdline = keywdline2 + " ";
     for (unsigned int i = 1; i < split_vec.size(); ++i)
     {
      keywdline = keywdline + split_vec[i];
     }
    }
    else
    {
     to_upper(keywdline);
    }
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
   for (int i = 4; i < 9; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream << "POINTLESS KEYWORDS" << std::endl;
   for (int i = 9; i < 11; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream << "AIMLESS KEYWORDS" << std::endl;
   for (int i = 11; i < 17; ++i)
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
  if (mode_string == "-s" 
      || mode_string == "-sLCV" 
      || mode_string == "-saLCV"
      || mode_string == "-sNC" 
      || mode_string == "-sLCVNC" 
      || mode_string == "-saLCVNC"
      || mode_string == "-sSF" 
      || mode_string == "-sLCVSF" 
      || mode_string == "-saLCVSF"
      || mode_string == "-sS6" 
      || mode_string == "-sLCVS6" 
      || mode_string == "-saLCVS6"
  )     // Synthesis pass
  {
   // In order to line up BLEND with the way ccp4i works (with stdin passed keywords) this is what has been added

   // Accepted keywords
   std::vector<std::string> akeys, dvkeywdline;
   akeys.push_back("NBIN");
   dvkeywdline.push_back("NBIN      20");
   akeys.push_back("RADF");
   dvkeywdline.push_back("RADFRAC   0.750");
   akeys.push_back("ISIG");
   dvkeywdline.push_back("ISIGI     1.500");
   akeys.push_back("CPAR");
   dvkeywdline.push_back("CPARWT    1.000");
   //akeys.push_back("DATA");
   akeys.push_back("DREF");
   akeys.push_back("LAUE");
   akeys.push_back("MAXC");
   akeys.push_back("COMP");
   akeys.push_back("CUTF");
   akeys.push_back("TOLE");
   akeys.push_back("CHOO");
   akeys.push_back("RESO");
   akeys.push_back("SDCO");
   akeys.push_back("ANOM");
   akeys.push_back("SCAL");
   akeys.push_back("EXCL");
   akeys.push_back("RUN");

   // Load in keywords from standard input
   // This is a way of reading lines from standard input and storing them into a string variable, line,
   // until a "/n" (carriage return) is entered
   int jflag;
   std::string keywdline,keywdline2;
   std::vector<std::string> split_vec;
   std::vector<std::string> vkeywdline;
   std::string::size_type Idx;
   std::cout << ">>>>> Input command lines <<<<<" << std::endl;
   std::cout << std::endl;
   getline(std::cin, keywdline);
   trim(keywdline);
   split(split_vec,keywdline,is_any_of(" "));
   keywdline2 = split_vec[0];
   to_upper(keywdline2);
   //if (keywdline2.substr(0,4) == "DATA")
   if (keywdline2.substr(0,4) == "DREF")
   {
    keywdline = keywdline2 + " ";
    for (unsigned int i = 1; i < split_vec.size(); ++i)
    {
     keywdline = keywdline + split_vec[i];
    }
   }
   else
   {
    to_upper(keywdline);
   }
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
    split(split_vec,keywdline,is_any_of(" "));
    keywdline2 = split_vec[0];
    to_upper(keywdline2);
    //if (keywdline2.substr(0,4) == "DATA")
    if (keywdline2.substr(0,4) == "DREF")
    {
     keywdline = keywdline2 + " ";
     for (unsigned int i = 1; i < split_vec.size(); ++i)
     {
      keywdline = keywdline + split_vec[i];
     }
    }
    else
    {
     to_upper(keywdline);
    }
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
   for (int i = 4; i < 9; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream << "POINTLESS KEYWORDS" << std::endl;
   for (int i = 9; i < 11; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream << "AIMLESS KEYWORDS" << std::endl;
   for (int i = 11; i < 17; ++i)
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
   if (mode_string == "-sLCV")  addrunmode2 = 1;
   if (mode_string == "-saLCV") addrunmode2 = 2;
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
  if (mode_string == "-c" || mode_string == "-cP" || mode_string == "-cF")     // Combination pass
  {
   // In order to line up BLEND with the way ccp4i works (with stdin passed keywords) this is what has been added

   // Accepted keywords
   std::vector<std::string> akeys, dvkeywdline;
   akeys.push_back("NBIN");
   dvkeywdline.push_back("NBIN      20");
   akeys.push_back("RADF");
   dvkeywdline.push_back("RADFRAC   0.750");
   akeys.push_back("ISIG");
   dvkeywdline.push_back("ISIGI     1.500");
   akeys.push_back("CPAR");
   dvkeywdline.push_back("CPARWT    1.000");
   //akeys.push_back("DATA");
   akeys.push_back("DREF");
   akeys.push_back("LAUE");
   akeys.push_back("MAXC");
   akeys.push_back("COMP");
   akeys.push_back("CUTF");
   akeys.push_back("TOLE");
   akeys.push_back("CHOO");
   akeys.push_back("RESO");
   akeys.push_back("SDCO");
   akeys.push_back("ANOM");
   akeys.push_back("SCAL");
   akeys.push_back("EXCL");
   akeys.push_back("RUN");

   // Load in keywords from standard input
   // This is a way of reading lines from standard input and storing them into a string variable, line,
   // until a "/n" (carriage return) is entered
   int jflag;
   std::string keywdline,keywdline2;
   std::vector<std::string> split_vec;
   std::vector<std::string> vkeywdline;
   std::string::size_type Idx;
   std::cout << ">>>>> Input command lines <<<<<" << std::endl;
   std::cout << std::endl;
   getline(std::cin, keywdline);
   trim(keywdline);
   split(split_vec,keywdline,is_any_of(" "));
   keywdline2 = split_vec[0];
   to_upper(keywdline2);
   //if (keywdline2.substr(0,4) == "DATA")
   if (keywdline2.substr(0,4) == "DREF")
   {
    keywdline = keywdline2 + " ";
    for (unsigned int i = 1; i < split_vec.size(); ++i)
    {
     keywdline = keywdline + split_vec[i];
    }
   }
   else
   {
    to_upper(keywdline);
   }
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
    split(split_vec,keywdline,is_any_of(" "));
    keywdline2 = split_vec[0];
    to_upper(keywdline2);
    //if (keywdline2.substr(0,4) == "DATA")
    if (keywdline2.substr(0,4) == "DREF")
    {
     keywdline = keywdline2 + " ";
     for (unsigned int i = 1; i < split_vec.size(); ++i)
     {
      keywdline = keywdline + split_vec[i];
     }
    }
    else
    {
     to_upper(keywdline);
    }
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
   for (int i = 4; i < 9; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream << "POINTLESS KEYWORDS" << std::endl;
   for (int i = 9; i < 11; ++i)
   {
    for (unsigned int j = 0; j < vkeywdline.size(); ++j)
    {
     Idx = vkeywdline[j].find(akeys[i]);
     //if (Idx != std::string::npos) keywd_ostream << vkeywdline[j] << std::endl;
     if (int(Idx) == 0) keywd_ostream << vkeywdline[j] << std::endl;
    }
   }
   keywd_ostream << "AIMLESS KEYWORDS" << std::endl;
   for (int i = 11; i < 17; ++i)
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
   if (mode_string == "-cP")  addrunmode3 = 1;
   if (mode_string == "-cF")  addrunmode3 = 2;

   // Turn arguments into integer numbers (in vector "arbitrary_datasets") to be later used in R code
   //for (int i=2;i < argc;i++) arbitrary_datasets.push_back(atoi(argv[i])); 
   // Turn arguments into a vector, "arbitrary_datasets", of characters to be later used in R code
   for (int i=2;i < argc;i++) arbitrary_datasets.push_back(argv[i]); 
  }
  if (mode_string == "-g"
      || mode_string == "-gNC"
      || mode_string == "-gSF"
      || mode_string == "-gS6")     // Graphics mode
  {
   runmode = 4;

   // Turn arguments into a vector, "arbitrary_datasets", of characters to be later used in R code
   for (int i=2;i < argc;i++) arbitrary_datasets.push_back(argv[i]); 
  }

  // I like well-formatted output
  std::cout.setf(std::ios::fixed);

  // Analysis mode or Dendrogram-Only mode
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
   if (mode_string == "-a")
   {
    std::cout << "You are now running BLEND in analysis mode." << std::endl; 
   }
   if (mode_string == "-aDO")
   {
    std::cout << "You are now running BLEND in dendrogram-only mode." << std::endl; 
   }
   if (mode_string == "-aNC")
   {
    std::cout << "You are now running BLEND in analysis mode using NCDist." << std::endl; 
   }
   if (mode_string == "-aDONC")
   {
    std::cout << "You are now running BLEND in dendrogram-only mode using NCDist." << std::endl; 
   }
   if (mode_string == "-aSF")
   {
     std::cout << "You are now running BLEND in analysis mode using SFDist." << std::endl;
   }
   if (mode_string == "-aDOSF")
   {
     std::cout << "You are now running BLEND in dendrogram-only mode using SFDist." << std::endl;
   }
   if (mode_string == "-aS6")
   {
     std::cout << "You are now running BLEND in analysis mode using S6Dist." << std::endl;
   }
   if (mode_string == "-aDOS6")
   {
     std::cout << "You are now running BLEND in dendrogram-only mode using S6Dist." << std::endl;
   }
   if (mode_string == "-k")
   {
     std::cout << "You are now running BLEND in k-means mode." << std::endl;
   }
   if (mode_string == "-kNC")
   {
     std::cout << "You are now running BLEND in k-means mode using NCDist." << std::endl;
   }
   if (mode_string == "-kS6")
   {
     std::cout << "You are now running BLEND in k-means mode using S6Dist." << std::endl;
   }
   if (mode_string == "-kSF")
   {
     std::cout << "You are now running BLEND in k-means mode using SFDist." << std::endl;
   }

   std::cout << std::endl;
   std::cout << "<!--SUMMARY_END--></FONT></B>" << std::endl;
   std::cout << std::endl;
   
   // Rejection information (if any) for users
   std::cout << "TYPES OF FLAGS ASSIGNED BY BLEND TO DATA SETS:" << std::endl;
   std::cout << std::endl;
   std::cout << "   crystal_flag = 0:   crystal is accepted for further processing" << std::endl;
   std::cout << "   crystal_flag = 1:   crystal is rejected because data file contains no reflections" << std::endl;
   std::cout << "   crystal_flag = 2:   crystal is rejected because data file is made up of multiple runs" << std::endl;
   std::cout << "   crystal_flag = 3:   crystal is rejected because data file is made up of multiple non-empty datasets" << std::endl;
   std::cout << std::endl;
   // Load crystals in unmerged data structures
   std::vector<scala::hkl_unmerge_list> hkl_list=load_crystals(filename,runmode);   // This is the correct expression for a copy constructor. Defining hkl_list first
                                                                                    // and then using hkl_list=load_crystals(filename) doesn't work. Ultimately this is
                                                                                    // due to the way the hkl_unmerged_list class has been implemented by Phil.
   // Label crystals according to dataset validity.
   // crystal_flag = 0:   crystal is accepted for further processing
   // crystal_flag = 1:   crystal is rejected because data file contains no reflections
   // crystal_flag = 2:   crystal is rejected because data file is made up of multiple runs
   // crystal_flag = 3:   crystal is rejected because data file is made up of multiple non-empty datasets
   
   // Initially all crystals are assumed to have valid datasets
   std::vector<int> crystal_flag(hkl_list.size(),0);

   // Label crystals
   crystal_flag=label_crystals(hkl_list,crystal_flag);
   std::cout << "FLAGS have been assigned as follows:" << std::endl;
   std::cout << std::endl;
   for (int i=0;i < crystal_flag.size();i++) std::cout << "   crystal flag for data set " << std::setw(4) << i+1 << ": " << crystal_flag[i] << std::endl;
   std::cout << std::endl;
   
   // Maps to store types of crystals belonging to different bravais lattices. These maps contain the "sg" bit as initially we were
   // handling space groups, rather than bravais lattices. Bravais lattice numbers are as follows:
   //
   // Triclinic   :     aP         1
   // Monoclinic  :     mP         2
   //                   mS         3
   // Orthorombic :     oP         4
   //                   oS         5
   //                   oF         6
   //                   oI         7
   // Tetragonal  :     tP         8
   //                   tI         9
   // Hexagonal   :     hP        10
   // Rhombohedral:     hH        11   (it is also possible hR, but we haven't implemented this way. To be sorted later)
   // Cubic       :     cP        12
   //                   cF        13
   //                   cI        14
   //
   // sg_to_crystal:       bravais lattice number  ->  crystal serial number
   // crystal_to_sg:    crystal serial number  ->  bravais lattice number
   // For NCDist we need the actual centering symbol: P, A, B, C, ..., rather than replacing A , B or C with S

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
   std::map<int,std::string> crystal_to_centering;
   std::multimap<int,int> sg_to_crystal;
   std::multimap<int,int>::iterator pos_cs;
   scala::hkl_symmetry symmetry;
   scala::CrystalType crystal_type;
   clipper::Spacegroup cspacegroup;
   clipper::Spgr_descr cspgr_descr;
   //size_t lnum;
   int lnum;
   std::string blatsymb;
   std::string centsymb;
   for (unsigned int i=0;i < hkl_list.size();i++)
   {
    if (crystal_flag[i] == 0)
    {
     symmetry=hkl_list[i].symmetry();
     crystal_type=CrystalType(symmetry);
     blatsymb=crystal_type.format();
     centsymb = blatsymb.substr(1,1);
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
     lnum=blatsymb.find("I");
     if (lnum != -1)
     {
      // Replace with S only for I 1 2 1
      if (symmetry.symbol_xHM() == "I 1 2 1")
      {
       blatsymb.erase(lnum,1);
       blatsymb.insert(lnum,"S");
      }
     }
     pos_bl=bl_symbol_to_number.find(blatsymb);
     //std::cout << "Lattice type: " << blatsymb << " Number is: " << pos_bl->second << std::endl;
     cspgr_descr=clipper::Spgr_descr(symmetry.symbol_Hall());
     cspacegroup=clipper::Spacegroup(cspgr_descr);
     //sg_to_crystal.insert(std::make_pair<int,int>(cspacegroup.spacegroup_number(),i));
     //crystal_to_sg.insert(std::make_pair<int,int>(i,cspacegroup.spacegroup_number()));
     sg_to_crystal.insert(std::make_pair<int,int>(pos_bl->second,i));
     crystal_to_sg.insert(std::make_pair<int,int>(i,pos_bl->second));
     crystal_to_centering[i] = centsymb;
    }
    else
    {
     sg_to_crystal.insert(std::make_pair<int,int>(0,i));
     crystal_to_sg.insert(std::make_pair<int,int>(i,0));
     crystal_to_centering[i] = std::string("?");
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
   output_summary_table(hkl_list,sg_to_crystal,crystal_to_centering,spacegroup_class,crystal_flag,mode_string);

   // Output ascii files for R (only if there are at least 2 crystals)
   if (hkl_list.size() > 1)
   {
    std::cout << "Building R files ........." << std::endl;
    if (mode_string == "-a") statistics_with_R(hkl_list,sg_to_crystal,crystal_to_centering,
                                                                        spacegroup_class,crystal_flag,R_program1,Rscp,mode_string);
    if (mode_string == "-aDO") statistics_with_R2(hkl_list,sg_to_crystal,crystal_to_centering,
                                                                        spacegroup_class,crystal_flag,R_program0,Rscp,mode_string);
    if (mode_string == "-aNC") statistics_with_R(hkl_list,sg_to_crystal,crystal_to_centering,
                                                                        spacegroup_class,crystal_flag,R_program1NC,Rscp,mode_string);
    if (mode_string == "-aDONC") statistics_with_R2(hkl_list,sg_to_crystal,crystal_to_centering,
                                                                        spacegroup_class,crystal_flag,R_program0NC,Rscp,mode_string);
    if (mode_string == "-aSF") statistics_with_R(hkl_list,sg_to_crystal,crystal_to_centering,
                                                    spacegroup_class,crystal_flag,R_program1SF,Rscp,mode_string);
    if (mode_string == "-aDOSF") statistics_with_R2(hkl_list,sg_to_crystal,crystal_to_centering,
                                                       spacegroup_class,crystal_flag,R_program0SF,Rscp,mode_string);
    if (mode_string == "-aS6") statistics_with_R(hkl_list,sg_to_crystal,crystal_to_centering,
                                                    spacegroup_class,crystal_flag,R_program1S6,Rscp,mode_string);
    if (mode_string == "-aDOS6") statistics_with_R2(hkl_list,sg_to_crystal,crystal_to_centering,
                                                       spacegroup_class,crystal_flag,R_program0S6,Rscp,mode_string);
    if (mode_string == "-k") statistics_with_R2(hkl_list, sg_to_crystal, crystal_to_centering,
                                                        spacegroup_class, crystal_flag, R_programK, Rscp, mode_string);
    if (mode_string == "-kNC") statistics_with_R2(hkl_list, sg_to_crystal, crystal_to_centering,
                                                        spacegroup_class, crystal_flag, R_programKNC, Rscp, mode_string);
    if (mode_string == "-kS6") statistics_with_R2(hkl_list, sg_to_crystal, crystal_to_centering,
                                                        spacegroup_class, crystal_flag, R_programKS6, Rscp, mode_string);
    if (mode_string == "-kSF") statistics_with_R2(hkl_list, sg_to_crystal, crystal_to_centering,
                                                        spacegroup_class, crystal_flag, R_programKSF, Rscp, mode_string);
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
    R_command_line << "Rscript " << R_program2 << " " << addrunmode2 << " " << dlevel_top << " " << dlevel_bottom;
   }
   else
   {
    R_command_line << "Rscript.exe " << R_program2 << " " << addrunmode2 << " " << dlevel_top << " " << dlevel_bottom;
   }
   //std::cout << R_command_line.str() << std::endl;
   R_status=std::system((R_command_line.str()).c_str());
   if (R_status != 0)
   {
    // Delete BLEND_KEYWORDS.dat before termination
    if (remove("BLEND_KEYWORDS.dat") != 0) std::cout << "Warning! Unable to remove file 'BLEND_KEYWORDS.dat'." << std::endl;
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
    //R_command_line << "Rscript " << R_program3 << " " << tmp.str();
    R_command_line << "Rscript " << R_program3 << " " << addrunmode3 << " " << tmp.str();
   }
   else
   {
    //R_command_line << "Rscript.exe " << R_program3 << " " << tmp.str();
    R_command_line << "Rscript.exe " << R_program3 << " " << addrunmode3 << " " << tmp.str();
   }
   //std::cout << R_command_line.str() << std::endl;
   R_status=std::system((R_command_line.str()).c_str());
   if (R_status != 0)
   {
    // Delete BLEND_KEYWORDS.dat before termination
    if (remove("BLEND_KEYWORDS.dat") != 0) std::cout << "Warning! Unable to remove file 'BLEND_KEYWORDS.dat'." << std::endl;
    int nerr=13;
    throw nerr;
   }

   // Delete BLEND_KEYWORDS.dat before termination
   if (remove("BLEND_KEYWORDS.dat") != 0) std::cout << "Warning! Unable to remove file 'BLEND_KEYWORDS.dat'." << std::endl;

   // End of BLEND - Combination mode
   std::cout << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << "## BLEND - Normal Termination                                   ##" << std::endl;
   std::cout << "##################################################################" << std::endl;
   std::cout << std::endl;
  }

  // Graphics mode
  if (runmode == 4)
  {
   // If less than 2 characters stop
   if (arbitrary_datasets.size() < 1)
   {
    int nerr=17;
    throw nerr;
   }
   std::cout << std::endl;
   std::cout << "<B><FONT COLOR='#FF0000'><!--SUMMARY_BEGIN-->" << std::endl;
   std::cout << std::endl;
   std::cout << "You are now running BLEND in graphics mode." << std::endl; 
   std::cout << std::endl;
   std::cout << "<!--SUMMARY_END--></FONT></B>" << std::endl;
   std::cout << std::endl;

   // Run R code
   std::cout << std::endl;
   std::cout << "Running R code to read information from previous runs of BLEND and produce some graphics..." << std::endl; 
   std::cout << std::endl;
   int R_status;
   std::ostringstream R_command_line,tmp;
   for (unsigned int i=0;i < arbitrary_datasets.size();i++) tmp << arbitrary_datasets[i] << " ";
   if (Rscp == 0)
   {
    R_command_line << "Rscript " << R_program4 << " " << tmp.str();
   }
   else
   {
    R_command_line << "Rscript.exe " << R_program4 << " " << tmp.str();
   }
   R_status=std::system((R_command_line.str()).c_str());
   if (R_status != 0)
   {
    int nerr=13;
    throw nerr;
   }

   // End of BLEND - Graphics mode
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
             << "In each of the following xx may be omitted, NC, SF or S6\n"
             << "Wrong command line format. Correct format is:\n"
             << "                                  \n"
             << "   blend -aDOxx name_of_file.dat                                               (dendrogram-only mode)\n"
             << "                 or               \n"
             << "   blend -aDOxx /path/to/directory                                             (dendrogram-only mode)\n"
             << "                 or               \n"
             << "   blend -axx name_of_file.dat                                                        (analysis mode)\n"
             << "                 or               \n"
             << "   blend -axx /path/to/directory                                                      (analysis mode)\n"
             << "                 or               \n"
             << "   blend -kxx name_of_file.dat                                                         (k-means mode)\n"
             << "                 or               \n"
             << "   blend -kxx /path/to/directory                                                       (k-means mode)\n"
             << "                 or               \n"
             << "   blend -sxx l1 (numeric height in dendrogram)                                      (synthesis mode)\n"
             << "   blend -sLCVxx l1 (LCV value) \n"
             << "   blend -saLCVxx l1 (aLCV value) \n"
             << "                 or               \n"
             << "   blend -sxx l1 l2 (numeric heights in dendrogram)                                  (synthesis mode)\n"
             << "   blend -sLCVxx l1 l2 (LCV values) \n"
             << "   blend -saLCVxx l1 (aLCV values) \n"
             << "                 or               \n"
             << "   blend -cxx  d1 d2 d3 ... (serial number of datasets)                            (combination mode)\n" 
             << "   blend -cPxx d1 d2 d3 ... (serial number of datasets)                                              \n" 
             << "   blend -cFxx d1 d2 d3 ... (serial number of datasets)                                              \n" 
             << "                 or               \n"
             << "   blend -gxx DO clN (cluster number) lN (level)           (graphics mode: aLCV annotated dendrogram)\n"
             << "   blend -gxx D clN (cluster number) lN (level)  (graphics mode: merging stats. annotated dendrogram)\n"
             //<< "   blend -gxx RM d1 d2 d3 ... (serial number of datasets   )             (graphics mode: Rmerge means\n"
             << std::endl;
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
             << "No valid unmerged MTZ files are listed in input file, or included in input directory." << std::endl;
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
             << "An error occurred in the execution of R code associated with BLEND." << std::endl;
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
             << "An error occurred in the execution of Python code associated with BLEND." << std::endl;
  }
  if (nerr == 17)
  {
   std::cerr << "\n EXECUTION ERROR!\n"
             << "Minimum number of command-line parameters in graphics mode is 1." << std::endl;
  }
 }

 return 0;
}
