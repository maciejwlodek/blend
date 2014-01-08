/************************************************************************************************************/
/************************************************************************************************************/
// !!!!!! REMEMBER: IF USING INTERACTIVELY WITH R UNCOMMENT LINE IN FUNCTION (statistics_with_R) in "aux_blend.cpp"
//
// BLEND
// Multi-crystals pre-processing program. Reads several unmerged mtz files and returns one or more single merged and
// sorted files, ready to be processed by "aimless" (SCALA)
// J. Foadi & G. Evans
// Diamond Light Source and Imperial College London
// September 2009
//
// The program runs in 2 modes:
//    1) reads all mtz files and produces global statistics;
//    2) reads height at which to cut dendrogram (tree from cluster analysis). Then produces as many mtz
//       files as clusters created by the cutting; 
/************************************************************************************************************/
/************************************************************************************************************/

/******************************************************/
// Standard C++ stuff
/******************************************************/
// <iostream> is included in Output.hh
// <string> is included in Output.hh
// <sstream> is included in Output.hh
#include <cstdlib>


/******************************************************/
// Stuff inherited from Phil Evans
/******************************************************/
#include "Output.hh"
#include "hkl_unmerge.hh"

/******************************************************/
// Other crystallographic libraries
/******************************************************/
#include <clipper/clipper.h>

/******************************************************/
// Stuff created by James Foadi
/******************************************************/
#include "aux_blend.hh"
std::string R_program1="/home/james/Dropbox/source_code/BLEND/ubuntu/src_R/blend1.R";       // Global variable
std::string R_program2="/home/james/Dropbox/source_code/BLEND/ubuntu/src_R/blend2.R";       // Global variable
std::string Python_program1="/home/james/Dropbox/source_code/BLEND/ubuntu/src_python/create_file_for_BLEND.py";       // Global variable
std::string Python_program2="/home/james/Dropbox/source_code/BLEND/ubuntu/src_python/merge_clusters.py";       // Global variable
std::string Python_program3="/home/james/Dropbox/source_code/BLEND/ubuntu/src_python/xds_to_mtz_list.py";       // Global variable



/******************************************************/
// MAIN
/******************************************************/
int main(int argc, char* argv[])
{
 try
 {
  // Check comand line input is well formatted
  int runmode=0;
  std::string filename,cl_height,cl_mode;
  if (argc > 2)
  {
   std::string tmpstring=argv[1];
   if (argc == 3 & tmpstring == "-d")
   {
    int Python_status;
    std::ostringstream Python_command_line;
    tmpstring=argv[2];    // Using this temporarily. Then I'll have to set value back to "-d", for code that follows
    std::cout << "Checking if there are xds files to be converted in directory " << tmpstring << " ..." << std::endl;
    Python_command_line << "python " << Python_program1 << " " << tmpstring;
    Python_status=std::system((Python_command_line.str()).c_str());
    filename="mtz_names.dat";
    runmode=1;
    tmpstring="-d";
   }
   if (argc == 3 & tmpstring != "-d")
   { 
    if (tmpstring != "single" & tmpstring != "ward" & tmpstring != "WARD" & tmpstring != "SINGLE")
    {
     int nerr=1;
     throw nerr;
    }
    if (tmpstring == "single" | tmpstring == "SINGLE") runmode=2;
    if (tmpstring == "ward" | tmpstring == "WARD") runmode=3;
    cl_height=argv[2];
   }
  }
  if (argc == 2)
  {
   if (std::atof(argv[1]) == 0)
   {
    // List of unmerged data could have mixed xds, mtz types. Convert xds into mtz
    std::string tmpstring=argv[1];
    std::cout << "Checking if there are xds files to be converted in input list " << tmpstring << " ..." << std::endl;
    int Python_status;
    std::ostringstream Python_command_line;
    Python_command_line << "python " << Python_program3 << " " << tmpstring;
    Python_status=std::system((Python_command_line.str()).c_str());
    filename="mtz_names.dat";
    runmode=1;
   }
   if (std::atof(argv[1]) != 0)
   {
    int nerr=1;
    throw nerr;
   }
  }
  
  if (argc < 2 | argc > 3)
  {
   int nerr=1;
   throw nerr;
  }

  // I like well-formatted output
  std::cout.setf(std::ios::fixed);

  // Amend command line arguments for CCP4 initialization
  int ncorrect=0;
  if (runmode == 1) ncorrect=2;
  if (runmode == 2) ncorrect=3;
  argc-=ncorrect;
  for (int i=0;i < argc;i++) argv[i]=argv[i+1]; 

  // Mode 1.
  if (runmode == 1)
  {
   // Name of ascii file containing list of mtz files
   //std::string tmpstringa,filename;
   std::string tmpstringa;
   std::vector<std::string> mtzin_name;
   //filename=argv[1];

   // Start CCP4 before anything else.
   CCP4::ccp4fyp(argc,argv);
   CCP4::ccp4_banner();
   std::cout << "You are running BLEND in Mode " << runmode << std::endl; 
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
   size_t lnum;
   std::string blatsymb;
   for (int i=0;i < hkl_list.size();i++)
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
     for (int i=0;i < spacegroup_class.size();i++)
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
   for (int i=0;i < spacegroup_class.size();i++)
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
    statistics_with_R(hkl_list,sg_to_crystal,spacegroup_class,crystal_flag,R_program1);
   }
  }

  // Mode 2 (runmode=2 | runmode=3)
  if (runmode == 2 | runmode == 3)
  {
   std::cout << std::endl;
   std::cout << "You are running BLEND in Mode 2" << std::endl; 
   std::cout << std::endl;

   // Associate single clustering to 1 and ward clustering to 2
   if (runmode == 2) cl_mode="1";
   if (runmode == 3) cl_mode="2";

   // Run R code
   std::cout << std::endl;
   std::cout << "Running R code to read statistical analysis from a previous run of BLEND and produce information on clusters..." << std::endl; 
   std::cout << std::endl;
   int R_status;
   std::ostringstream R_command_line;
   R_command_line << "R --vanilla --slave --quiet < " << R_program2 << " --args " << cl_height << " " << cl_mode;
   R_status=std::system((R_command_line.str()).c_str());

   // Run python code for merging
   int Python_status;
   std::ostringstream Python_command_line,argsline;
   if (runmode == 2) argsline << " ./single_" << cl_height << "/CLUSTERS_SINGLE.dat"; 
   if (runmode == 3) argsline << " ./ward_" << cl_height << "/CLUSTERS_WARD.dat"; 
   //if (runmode == 2) Python_command_line << "python " << Python_program2 << " CLUSTERS_SINGLE.dat";
   //if (runmode == 3) Python_command_line << "python " << Python_program2 << " CLUSTERS_WARD.dat";
   Python_command_line << "python " << Python_program2 << argsline.str();
   Python_status=std::system((Python_command_line.str()).c_str());

   // Exit message for user
   std::cout << std::endl;
   if (runmode == 2) std::cout << "All clusters are shown in the file CLUSTERS_SINGLE.dat" << std::endl;
   if (runmode == 3) std::cout << "All clusters are shown in the file CLUSTERS_WARD.dat" << std::endl;
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
             << "   blend mtz_names.dat                                                    (mode 1)\n"
             << "                 or               \n"
             << "   blend mtz_names.dat [height] CLUSTERS_SINGLE                           (mode 2)\n"
             << "                 or               \n"
             << "   blend mtz_names.dat [height] CLUSTERS_WARD                             (mode 2)\n"
             << "    (where [height] is the selected distance in one of the two cluster dendrograms)" << std::endl;
   return EXIT_FAILURE;
  }
  if (nerr == 2)
  {
   std::cerr << "\n BLEND ERROR!\n"
             << "Wrong command line format: first command line argument should be the name\n"
             << "of a file listing the name of mtz files. Correct format is:\n"
             << "                                  \n"
             << "   blend mtz_names.dat                                                    (mode 1)\n"
             << "                 or               \n"
             << "   blend mtz_names.dat c1 c2 c3 ... cn                                    (mode 2)\n"
             << "                 or               \n"
             << "   blend mtz_names.dat HKLOUT name_of_output.mtz c1 c2 c3 ... cm          (mode 3)" << std::endl;
   return EXIT_FAILURE;
  }
  if (nerr == 3)
  {
   std::cerr << "\n BLEND ERROR!\n"
             << "You appear to be running BLEND in mode 3, but your input\n"
             << "should consist of at least 4 arguments. Correct format is:\n" 
             << "                                  \n"
             << "   blend mtz_names.dat                                                    (mode 1)\n"
             << "                 or               \n"
             << "   blend mtz_names.dat c1 c2 c3 ... cn                                    (mode 2)\n"
             << "                 or               \n"
             << "   blend mtz_names.dat HKLOUT name_of_output.mtz c1 c2 c3 ... cm          (mode 3)" << std::endl;
   return EXIT_FAILURE;
  }
  if (nerr == 4)
  {
   std::cerr << "\n BLEND ERROR!\n"
             << "You appear to be running BLEND in mode 3, but your 3rd argument is not\n"
             << "the name of an mtz file. Correct format is:\n"
             << "                                  \n"
             << "   blend mtz_names.dat                                                    (mode 1)\n"
             << "                 or               \n"
             << "   blend mtz_names.dat c1 c2 c3 ... cn                                    (mode 2)\n"
             << "                 or               \n"
             << "   blend mtz_names.dat HKLOUT name_of_output.mtz c1 c2 c3 ... cm          (mode 3)" << std::endl;
   return EXIT_FAILURE;
  }
  if (nerr == 5)
  {
   std::cerr << "\n BLEND ERROR!\n"
             << "You appear to be running BLEND in mode 3, but one of your command arguments from the\n"
             << "fourth onward is not a number. Correct format is:\n"
             << "                                  \n"
             << "   blend mtz_names.dat                                                    (mode 1)\n"
             << "                 or               \n"
             << "   blend mtz_names.dat c1 c2 c3 ... cn                                    (mode 2)\n"
             << "                 or               \n"
             << "   blend mtz_names.dat HKLOUT name_of_output.mtz c1 c2 c3 ... cm          (mode 3)" << std::endl;
   return EXIT_FAILURE;
  }
  if (nerr == 6)
  {
   std::cerr << "\n BLEND ERROR!\n"
             << "You appear to be running BLEND in mode 2, but one of your command arguments from the\n"
             << "second onward is not a number. Correct format is:\n"
             << "                                  \n"
             << "   blend mtz_names.dat                                                    (mode 1)\n"
             << "                 or               \n"
             << "   blend mtz_names.dat c1 c2 c3 ... cn                                    (mode 2)\n"
             << "                 or               \n"
             << "   blend mtz_names.dat HKLOUT name_of_output.mtz c1 c2 c3 ... cm          (mode 3)" << std::endl;
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
             << "No valid MTZ files are listed in ... file" << std::endl;
   return EXIT_FAILURE;
  }
  if (nerr == 9)
  {
   std::cerr << "\n USER ERROR!\n"
             << "The specific MTZ file cannot be read." << std::endl;
  }
 }

 return 0;
}
