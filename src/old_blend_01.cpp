/************************************************************************************************************/
/************************************************************************************************************/
// !!!!!! REMEMBER: IF USING INTERACTIVELY WITH R UNCOMMENT LINE IN FUNCTION (statistics_with_R) in "aux_blend.cpp"
// BLEND
// Multi-crystals pre-processing program. Reads several unmerged mtz files and returns a single merged and
// sorted file, ready to be processed by "ainless" (SCALA)
// J. Foadi & G. Evans
// Diamond Light Source and Imperial College London
// September 2009
//
// The program runs in 3 modes:
//    1) reads all mtz files and produces global statistics;
//    2) reads user choice (based on previous global statistics) and carries on with overlap matrix
//    3) reads user choice (based on previous global statistics and overlap matrix) and carries on with data weighting and merging 
/************************************************************************************************************/
/************************************************************************************************************/

/******************************************************/
// Standard C++ stuff
/******************************************************/
// <iostream> is included in Output.hh
// <string> is included in Output.hh
// <sstream> is included in Output.hh
#include <cstdlib>
#include <iomanip>


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
std::string R_program="/Users/james/Dropbox/source_code/BLEND/src_R/blend.R";       // Global variable



/******************************************************/
// MAIN
/******************************************************/
int main(int argc, char* argv[])
{
 try
 {
  std::string mtzout_name;

  // Check comand line input is well formatted
  std::vector<int> crystals_choice(read_command_line(argc,argv));   // crystals_choice stores crystals serial number starting from 0
  int runmode=crystals_choice[crystals_choice.size()-1];
  crystals_choice.pop_back();

  // I like well-formatted output
  std::cout.setf(std::ios::fixed);

  // Amend command line arguments for CCP4 initialization
  int ncorrect=0;
  if (runmode == 1) ncorrect=2;
  if (runmode == 2)
  {
   ncorrect=2+crystals_choice.size();
   mtzout_name="HKLOUT";
  } 
  if (runmode == 3)
  {
   ncorrect=4+crystals_choice.size();
   mtzout_name=argv[3];
  }
  argc-=ncorrect;
  for (int i=0;i < argc;i++) argv[i]=argv[i+1]; 

  // Mode 1.
  if (runmode == 1)
  {
   // Name of ascii file containing list of mtz files
   std::string tmpstringa,filename;
   std::vector<std::string> mtzin_name;
   filename=argv[1];

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
   for (int i=0;i < crystal_flag.size();i++) std::cout << "crystal flag " << i << ": " << crystal_flag[i] << std::endl;
   
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
    statistics_with_R(hkl_list,sg_to_crystal,spacegroup_class,crystal_flag,R_program);
   }

  }

  // Mode 2.
  if (runmode == 2)
  {
   // Name of ascii file containing list of mtz files
   std::string filename;
   std::vector<std::string> mtzin_name;
   filename=argv[1];

   // Change input file for runmode 2 and 3
   //std::ostringstream ofilename;
   //ofilename << "NEW_" << filename;
   //std::string new_filename=ofilename.str();
   std::string new_filename=filename;

   // Start CCP4 before anything else.
   CCP4::ccp4fyp(argc,argv);
   CCP4::ccp4_banner();
   std::cout << "You are running BLEND in Mode " << runmode << std::endl; 
   
   // Load crystals in unmerged data structures
   std::vector<scala::hkl_unmerge_list> hkl_list=load_crystals(new_filename,runmode);   // This is the correct expression for a copy constructor. Defining hkl_list first
                                                                                        // and then using hkl_list=load_crystals(filename) doesn't work. Ultimately this is
                                                                                        // due to the way the hkl_unmerged_list class has been implemented by Phil.

   // Label crystal_flag = 4 if it hasn't been selected by user
   std::vector<int> crystal_flag(hkl_list.size(),4);
   for (int i=0;i < crystals_choice.size();i++)
   {
    crystal_flag[crystals_choice[i]]=0;
   }

   // Label crystals
   crystal_flag=label_crystals(hkl_list,crystal_flag);
   
   // Overlaps matrix among selected crystals only
   std::cout << "Computing overlaps matrix ........." << std::endl;
   const int msize=crystals_choice.size();
   std::vector< std::vector<float> > overlaps_matrix=build_overlaps_matrix(hkl_list,crystals_choice);

   // Output overlaps matrix
   std::ostringstream olinea;
   std::string linea;
   std::ofstream overlaps_file("OVERLAPS.txt",std::ios::out);
   overlaps_file.setf(std::ios::fixed);   // Formatted output for overlaps matrix file
   overlaps_file << std::endl;
   overlaps_file << "  Overlaps matrix for crystals ";
   for (int i=0;i < crystals_choice.size();i++) overlaps_file << "  " << crystals_choice[i]+1;
   overlaps_file << std::endl;
   overlaps_file << "         " << std::endl;
   overlaps_file << "         " << std::endl;
   overlaps_file << "  ";
   for (int i=0;i < crystals_choice.size();i++) overlaps_file << std::setw(6) << crystals_choice[i]+1 << "  ";
   overlaps_file << std::endl; 
   overlaps_file << "  -";
   for (int i=0;i < crystals_choice.size();i++) overlaps_file << "--------";
   overlaps_file << std::endl;
   for (int irow=0;irow < msize;irow++)
   {
    overlaps_file << "  |";
    for (int icol=0;icol < msize;icol++)
    {
     overlaps_file << " " << std::setw(5) << std::setprecision(3) << overlaps_matrix[irow][icol] << " |";
    }
    overlaps_file << std::setw(3) << crystals_choice[irow]+1 << std::endl;
    overlaps_file << "  -";
    for (int i=0;i < crystals_choice.size();i++) overlaps_file << "--------";
    overlaps_file << std::endl;
   }
   overlaps_file.close();
   std::cout << "Matrix of overlapping reciprocal-space zones between all crystals is contained in file OVERLAPS.txt." << std::endl;

  }

  // Mode 3.
  if (runmode == 3)
  {
   // Name of ascii file containing list of mtz files
   std::string tmpstringa,filename;
   std::vector<std::string> mtzin_name;
   filename=argv[1];

   // Change input file for runmode 2 and 3
   //std::ostringstream ofilename;
   //ofilename << "NEW_" << filename;
   //std::string new_filename=ofilename.str();
   std::string new_filename=filename;

   // Start CCP4 before anything else.
   CCP4::ccp4fyp(argc,argv);
   CCP4::ccp4_banner();
   std::cout << "You are running BLEND in Mode " << runmode << std::endl; 
   
   // Load crystals in unmerged data structures
   std::vector<scala::hkl_unmerge_list> hkl_list=load_crystals(new_filename,runmode);   // This is the correct expression for a copy constructor. Defining hkl_list first
                                                                            // and then using hkl_list=load_crystals(filename) doesn't work. Ultimately this is
                                                                            // due to the way the hkl_unmerged_list class has been implemented by Phil.

   // Label crystal_flag = 4 if it hasn't been selected by user
   std::vector<int> crystal_flag(hkl_list.size(),4);
   for (int i=0;i < crystals_choice.size();i++)
   {
    crystal_flag[crystals_choice[i]]=0;
   }

   // Label crystals
   crystal_flag=label_crystals(hkl_list,crystal_flag);
   
   // Do something (e.g. weighting) to reflections for all crystals

   // Merge all datasets into one file
  }

 }
 catch (int nerr)
 {
  if (nerr == 1)
  {
   std::cerr << "\n BLEND ERROR!\n"
             << "Wrong command line format: too few arguments. Correct format is:\n"
             << "                                  \n"
             << "   blend mtz_names.dat                                                    (mode 1)\n"
             << "                 or               \n"
             << "   blend mtz_names.dat c1 c2 c3 ... cn                                    (mode 2)\n"
             << "                 or               \n"
             << "   blend mtz_names.dat HKLOUT name_of_output.mtz c1 c2 c3 ... cm          (mode 3)" << std::endl;
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
