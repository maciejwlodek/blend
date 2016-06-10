/************************************************************************************************************/
/************************************************************************************************************/
/********* aux_blend.cpp                                                                            *********/
/*********                                                                                          *********/
/********* Copyright (C) 2014 Diamond Light Source & Imperial College London                        *********/
/*********                                                                                          *********/
/********* Authors: James Foadi & Gwyndaf Evans                                                     *********/
/*********                                                                                          *********/
/********* This code is distributed under the BSD license, a copy of which is                       *********/
/********* included in the root directory of this package.                                          *********/
/************************************************************************************************************/
/************************************************************************************************************/
//

#include "aux_blend.hh"

std::vector<scala::hkl_unmerge_list> load_crystals(std::string filename,int runmode)
{
 int nwrong_mtz=0,nlines=0;
 std::ifstream file(filename.c_str(),std::ios::in);
 std::vector<std::string> mtz_file,no_unmerged_files;
 std::string tmpstring;
 if (!file.is_open())
 {
  int nerr=7;
  throw nerr;
 }
 else
 {
  std::cout << "Reading from file " << filename << " ........." << std::endl;
  while (!file.eof())
  {
   file >> tmpstring;
   if (tmpstring.size() == 0)    // Stop if list-of-mtz-files file is empty
   {
    int nerr=8;
    throw nerr;
   }
   nlines++;
   if (tmpstring.substr(tmpstring.length()-4,4) == ".mtz")   // check all filenames end up by ".mtz"
   {
    // Check this is really an unmerged mtz file
    MtzIO::MtzUnmrgFile tmpmtzin;
    bool fileread=tmpmtzin.open_read(tmpstring);
    tmpmtzin.close_read();
    if (fileread)
    {
     mtz_file.push_back(tmpstring);
     if (file.eof()) mtz_file.pop_back();
    }
    else
    {
     no_unmerged_files.push_back(tmpstring);
     nwrong_mtz++;
     if (file.eof()) nwrong_mtz--;
    } 
   }
   else
   {
    no_unmerged_files.push_back(tmpstring);
    nwrong_mtz++;
    if (file.eof()) nwrong_mtz--;
   }
   if (file.eof()) nlines--;
  }
 }
 file.close();
 if (nwrong_mtz > 0)
 {
  std::cout << std::endl;
  std::cout << "********* ********* *********" << std::endl;
  std::cout << "Your input list includes " << nwrong_mtz  << " missing or not properly-formatted unmerged mtz files" << std::endl;
  std::cout << "They are: " << std::endl;
  for (unsigned int i=0; i < no_unmerged_files.size(); i++)
  {
   std::cout << "          " << no_unmerged_files[i] << std::endl;
  }
  std::cout << "********* ********* *********" << std::endl;
  std::cout << std::endl;
 }
  
 // Stop execution if no MTZ files are found
 if (mtz_file.size() == 0)
 {
  // Error coded as 8
  int nerr=8;
  throw nerr;
 } 

 // Vector of unmerged data sets
 const int ncrystals=mtz_file.size();
 std::vector<scala::hkl_unmerge_list> hkl_list(ncrystals); 

 // Output object, needed in MtzIO::MtzUnmrgFile method AddHklList. Not really needed for anything else here.
 //phaser_io::Output output;
 std::string output;

 // File selection flags (resolution, datasets, batches etc)
 scala::file_select file_sel;

 // Intensity or profile selection ( =0 means we use intensities)
 scala::col_controls column_selection;
 column_selection.SetIcolFlag(0,0.0);

 // Make list of required columns in unmerged HKLIN file
 MtzIO::column_labels column_list=MtzIO::setup_columns();

 // Scala control classes (default settings)
 scala::all_controls controls;

 // Loop over input mtz files.
 std::cout << std::endl;
 std::cout << "Loading unmerged mtz files ........." << std::endl;
 for (unsigned int i=0;i < mtz_file.size();i++)
 {
  // MTZ valid or not
  bool fileread;

  // Open mtz file
  MtzIO::MtzUnmrgFile mtzin;

  // Make sure hkl_unmerge_list object is empty (to avoid compatibility-check errors)
  hkl_list[i].clear();

  // Fill hkl_list[i]
  const int verbose=+1;
  std::cout << "   Trying with dataset " << mtz_file[i] << " ......... ";
  fileread=(mtzin.AddHklList(i,mtz_file[i],file_sel,column_selection,column_list,controls,PxdName(),Scell(),double(2.0),output,verbose,hkl_list[i])).Read();
  std::cout << "Done." << std::endl;
  if (!fileread) nwrong_mtz++;
  //std::cout << "Partial controls: " << controls.partials.accept_fract_min() << "        " << controls.partials.accept_fract_max() << std::endl;
 }
 std::cout << "Loading complete." << std::endl;

 // File with working mtz files is written in "New_[name of old file]"
 if (runmode == 1)
 {
  //std::ostringstream ofilename;
  //ofilename << "NEW_" << filename;
  //std::string new_filename=ofilename.str();
  std::string new_filename="NEW_list_of_files.dat";
  std::ofstream new_file(new_filename.c_str(),std::ios::out);
  for (unsigned int i=0;i < mtz_file.size();i++)
  {
   new_file << mtz_file[i] << std::endl;
  }
  new_file.close();
  //std::cout << "List of valid files to be further analysed is included in file " << new_filename << "." << std::endl;
  std::cout << std::endl;
 }

 return hkl_list;
}

// Label crystals according to dataset validity.
// crystal_flag = 0:   crystal is accepted for further processing
// crystal_flag = 1:   crystal is rejected because data file contains no reflections
// crystal_flag = 2:   crystal is rejected because data file is made up of multiple runs
// crystal_flag = 3:   crystal is rejected because data file is made up of multiple datasets
// crystal_flag = 4;   crystal is rejected by user (Mode 2 and Mode 3)
std::vector<int> label_crystals(const std::vector<scala::hkl_unmerge_list>& hkl_list,std::vector<int> crystal_flag_in)
{
 std::vector<int> crystal_flag_out(crystal_flag_in);   // Way to copy all elements of a vector container

 // Loop over all datasets
 for (unsigned int i=0;i < hkl_list.size();i++)
 {
  // Does dataset contain data?
  if (hkl_list[i].IsEmpty() && crystal_flag_in[i] != 4) crystal_flag_out[i]=1;

  // Does dataset consists of multiple datasets or runs? (At the moment we do not accept this. For us 1 crystal = 1 dataset = 1 run)
  //std::cout << "Dataset " << i+1 << " number of runs " << hkl_list[i].num_runs() << std::endl;
  //std::cout << "Dataset " << i+1 << " number of datasets " << hkl_list[i].num_datasets() << std::endl;
  if (hkl_list[i].num_runs() != 1 && crystal_flag_in[i] != 4) crystal_flag_out[i]=2;
  if (hkl_list[i].num_datasets() != 1 && crystal_flag_in[i] != 4) crystal_flag_out[i]=3;
 } 

 return crystal_flag_out;
}

// Function to output summary table for all crystals (produced in Mode 1)
void output_summary_table(std::vector<scala::hkl_unmerge_list>& hkl_list,std::multimap<int,int> sg_to_crystal,
                          std::vector<int> spacegroup_class,std::vector<int> crystal_flag)
{
 // Build map container for bravais lattice number to symbol correspondence
 std::map<int,std::string> bl_number_to_symbol;
 std::map<int,std::string>::iterator pos_bl;
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(1,"aP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(2,"mP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(3,"mS"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(4,"oP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(5,"oS"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(6,"oF"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(7,"oI"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(8,"tP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(9,"tI"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(10,"hP"));
 //bl_number_to_symbol.insert(std::make_pair<int,std::string>(11,"hR"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(11,"hH"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(12,"cP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(13,"cF"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(14,"cI"));

 // Objects needed to output information from each crystal dataset
 scala::Xdataset dataset;
 scala::Scell unitcell;
 //scala::hkl_symmetry symmetry;
 scala::ResoRange resorange;
 clipper::Cell_descr ccell_descr;
 clipper::Cell ccell;
 

 // Objects needed for getting information from individual reflections
 scala::reflection reflection;
 scala::observation observation;
 scala::observation_part obspart;

 // Open file
 std::ofstream summary_file("BLEND_SUMMARY.txt",std::ios::out);

 // Header
 summary_file.setf(std::ios::fixed);   // Permanent until reset
 summary_file << "\n" << "********* ********* ********* SUMMARY INFORMATION ON ALL CRYSTALS ********* ********* *********" << "\n\n" << std::endl;

 // Iterators
 std::multimap<int,int>::iterator pos_cs;
 std::pair<std::multimap<int,int>::iterator,std::multimap<int,int>::iterator> pos_pair;

 // Loop over different bravais lattices
 for (unsigned int i=0;i < spacegroup_class.size();i++)
 {
  if (spacegroup_class[i] != 0)
  {
  //symmetry=scala::hkl_symmetry(spacegroup_class[i]);

  // Header for each space group
  pos_bl=bl_number_to_symbol.find(spacegroup_class[i]);
  summary_file << "         BRAVAIS LATTICE NUMBER " << spacegroup_class[i] << ": (" << pos_bl->second << ")\n" << std::endl;
  summary_file << "==============|=============================================================|=============|=++==========|================|=====================|=============|"
               << std::endl;
  summary_file << "              |                              CELL                           |             |             |   RESOLUTION   | CRYSTAL TO DETECTOR |             |"
               << std::endl;
  summary_file << "  CRYSTAL N.  |                                                             | CELL VOLUME | MOSAICITY   |     RANGE      |                     | WAVELENGTH  |"
               << std::endl;
  summary_file << "              |       a         b         c      alpha    beta     gamma    |             |             | Low       High |      DISTANCE       |             |"
               << std::endl;
  summary_file << "==============|=============================================================|=============|=============|================|=====================|=============|"
               << std::endl;

  // Suff for loggraph, etc
  std::cout << std::endl;
  std::cout << "$TABLE: Summary of type I parameters :" << std::endl;
  //std::cout << "$GRAPHS:        Cell Parameters   : N : 1, 2, 3, 4, 5, 6, 7 :" << std::endl;
  std::cout << "$GRAPHS:        Cell Parameters (a,     b,    c)     :N:1,2,3,4:" << std::endl;
  std::cout << "       :        Cell Parameters (alpha, beta, gamma) :N:1,5,6,7:" << std::endl;
  std::cout << "       :        Cell Volume                          :N:1,8:" << std::endl; 
  std::cout << "       :        Mosaicity                            :N:1,9:" << std::endl; 
  std::cout << "$$" << std::endl;
  std::cout << "Crystal  a  b  c  alpha  beta  gamma  Cell_Volume  Mosaicity  Resolution(low)  Resolution(high)  Wavelength  $$" << std::endl;    
  std::cout << "$$" << std::endl;      // Part of this table carries on further down in code


  // Pair of iterators pointing at beginning and end of values (crystals) with same key (spacegroup number)
  pos_pair=sg_to_crystal.equal_range(spacegroup_class[i]);
  for (pos_cs=pos_pair.first;pos_cs != pos_pair.second;++pos_cs)
  {
   //std::cout << "pos_cs->first " << pos_cs->first << " pos_cs->second " << pos_cs->second << "  crystal_flag " << crystal_flag[pos_cs->second] << std::endl;
   if (crystal_flag[pos_cs->second] == 0)
   {
    resorange=hkl_list[pos_cs->second].ResRange();
    dataset=hkl_list[pos_cs->second].xdataset(0);   // At present we only consider crystals with 1 dataset
    unitcell=dataset.cell();
    //ccell_descr=clipper::Cell_descr(unitcell.VCell()[0],unitcell.VCell()[1],unitcell.VCell()[2],unitcell.VCell()[3],unitcell.VCell()[4],unitcell.VCell()[5]);
    ccell_descr=clipper::Cell_descr(unitcell.UnitCell()[0],unitcell.UnitCell()[1],unitcell.UnitCell()[2],unitcell.UnitCell()[3],unitcell.UnitCell()[4],unitcell.UnitCell()[5]);
    ccell=clipper::Cell(ccell_descr);
    scala::Batch batch=hkl_list[pos_cs->second].batch(0);  // I'm using any batch (for instance number 0) just to recover crystal-to-detector distance
    CMtz::MTZBAT batchmtz=batch.batchdata();
    float ctoddist=batchmtz.dx[0];
    float wlength=dataset.wavelength();

    // Organise data (Phil's code. Needed to count unique reflections, etc)
    int nrefs=hkl_list[pos_cs->second].prepare();
    int  nobs=hkl_list[pos_cs->second].sum_partials();

    // Completeness in resolution shells (low, middle, high)
    //int nbins=resorange.Nbins();   // If SetNbins has not been used, Nbins() returns an optimal number of bins
    //std::cout << "Number of bins: " << resorange.Nbins() << std::endl;
    //int nbins=3;   // 3 bins for now: low, middle, high
    //resorange.SetNbins(nbins);

    // Find completeness in full reciprocal-space sphere

    // 1) Theoretical number of reflections in full-sphere shells
    //std::vector<int> Nrefres(nbins,0);
    //std::vector<int> Nrefacen(nbins,0);
    //scala::NumberComplete(resorange,symmetry,unitcell,Nrefres,Nrefacen);

    // 2) Actual number of reflections in full-sphere shells (cycle through all reflections)
    //std::vector<int> NumRefSphere(nbins,0);
    //int NumSymm=symmetry.Nsym();
    //while (hkl_list[pos_cs->second].next_reflection(this_refl) >= 0)
    //{
     //Rtype invresolsq=this_refl.invresolsq();
     //bool Centric=symmetry.is_centric(this_refl.hkl());
     //float epsiln=symmetry.get_multiplicity(this_refl.hkl());
     //int multcy=Nint(float(NumSymm)/epsiln);
     //if (!Centric) multcy=multcy*2;
     //int mres=resorange.bin(invresolsq);   // Specific bin this reflection belongs to
     //NumRefSphere[mres]+=multcy;   // Total number of reflections in resolution shell
    //}

    // 3) Completeness in resolution shells and global
    //int genint1=0;
    //std::vector<float> bin_completeness(nbins,0);
    //for (int j=0;j < nbins;j++) genint1+=Nrefres[j];
    //int genint2=0;
    //for (int j=0;j < nbins;j++)
    //{
    // genint2+=NumRefSphere[j];
    // bin_completeness[j]=100*float(NumRefSphere[j])/float(Nrefres[j]);
    //}
    //float completeness=100*float(genint2)/float(genint1);

    // Output to summary file
    summary_file << "  " << std::setw(10) << pos_cs->second+1 << "  |" << std::setw(10) << std::setprecision(3) << unitcell.UnitCell()[0]
                                                              << std::setw(10) << std::setprecision(3) << unitcell.UnitCell()[1]
                                                              << std::setw(10) << std::setprecision(3) << unitcell.UnitCell()[2]
                                                              << std::setw(9) << std::setprecision(2) << unitcell.UnitCell()[3]
                                                              << std::setw(9) << std::setprecision(2) << unitcell.UnitCell()[4]
                                                              << std::setw(9) << std::setprecision(2) << unitcell.UnitCell()[5] << "    |"
                                                              << " " << std::setw(11) << std::setprecision(2) << ccell.volume() << " |"
                                                              << std::setw(11) << std::setprecision(5) << dataset.Mosaicity() << "  |"
                                                              << std::setw(7) << std::setprecision(3) << resorange.ResLow() << " "
                                                              << std::setw(7) << std::setprecision(3) << resorange.ResHigh() << " |       "
                                                              << std::setw(7) << std::setprecision(2) << ctoddist << "       |"
                                                              << std::setw(10) << std::setprecision(5) << wlength << "   |"
                                                              //<< "    " << std::setw(5) << std::setprecision(1) << completeness << "     |"
                                                              //<< " " << std::setw(5) << std::setprecision(1) << bin_completeness[0] << "  "
                                                              //<< " " << std::setw(5) << std::setprecision(1) << bin_completeness[1] << "  "
                                                              //<< " " << std::setw(5) << std::setprecision(1) << bin_completeness[2] << "  |"
                                                              << std::endl;

    // For loggraph, etc
    std::cout << "  " << std::setw(10) << pos_cs->second+1 << std::setw(10) << std::setprecision(3) << unitcell.UnitCell()[0]
                                                           << std::setw(10) << std::setprecision(3) << unitcell.UnitCell()[1]
                                                           << std::setw(10) << std::setprecision(3) << unitcell.UnitCell()[2]
                                                           << std::setw(9) << std::setprecision(2) << unitcell.UnitCell()[3]
                                                           << std::setw(9) << std::setprecision(2) << unitcell.UnitCell()[4]
                                                           << std::setw(9) << std::setprecision(2) << unitcell.UnitCell()[5]
                                                           << " " << std::setw(11) << std::setprecision(2) << ccell.volume()
                                                           << std::setw(11) << std::setprecision(5) << dataset.Mosaicity()
                                                           << std::setw(9) << std::setprecision(3) << resorange.ResLow()
                                                           << std::setw(9) << std::setprecision(3) << resorange.ResHigh()
                                                           << std::setw(10) << std::setprecision(5) << wlength
                                                           << std::endl;
   }
  }

  // For loggraph, etc
  std::cout << "$$" << std::endl;
  std::cout << std::endl;

  summary_file << "==============|=============================================================|=============|=============|================|=====================|=============|"
               << std::endl;
  summary_file << "\n\n" << std::endl;
  }
 }
 summary_file.close();
 std::cout << "Summary information is included in file BLEND_SUMMARY.txt" << std::endl;
}

// Function to output ascii files to be read by R (produced in Mode 1)
void statistics_with_R(std::vector<scala::hkl_unmerge_list>& hkl_list,std::multimap<int,int> sg_to_crystal,
                          std::vector<int> spacegroup_class,std::vector<int> crystal_flag,std::string R_program,int Rscp)
{
 // Build map container for bravais lattice number to symbol correspondence
 std::map<int,std::string> bl_number_to_symbol;
 std::map<int,std::string>::iterator pos_bl;
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(1,"aP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(2,"mP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(3,"mS"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(4,"oP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(5,"oS"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(6,"oF"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(7,"oI"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(8,"tP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(9,"tI"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(10,"hP"));
 //bl_number_to_symbol.insert(std::make_pair<int,std::string>(11,"hR"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(11,"hH"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(12,"cP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(13,"cF"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(14,"cI"));

 // Objects needed to output information from each crystal dataset
 scala::Xdataset dataset;
 scala::Scell unitcell;
 scala::ResoRange resorange;
 scala::reflection reflection;
 scala::observation observation;
 clipper::Cell_descr ccell_descr;
 clipper::Cell ccell;

 // File "forR_raddam.dat" containing names of all refs_ files
 std::ofstream forR_file("forR_raddam.dat",std::ios::out);

 // File "forR_macropar.dat" containing cell parameters and mosaicity for all crystals
 std::ofstream forR_file2("forR_macropar.dat",std::ios::out);
 forR_file2.setf(std::ios::fixed);   // Permanent until reset

 // Ascii files with tables to be read by R. Names change from bravais lattice to bravais lattice
 std::string filename;
 std::ostringstream ofilename;

 // Iterators
 std::multimap<int,int>::iterator pos_cs;
 std::pair<std::multimap<int,int>::iterator,std::multimap<int,int>::iterator> pos_pair;

 // Loop over different space groups
 for (unsigned int i=0;i < spacegroup_class.size();i++)
 {
  if (spacegroup_class[i] != 0)
  {
   // Pair of iterators pointing at beginning and end of values (crystals) with same key (spacegroup number)
   pos_pair=sg_to_crystal.equal_range(spacegroup_class[i]);
   for (pos_cs=pos_pair.first;pos_cs != pos_pair.second;++pos_cs)
   {
    if (crystal_flag[pos_cs->second] == 0)
    {
     resorange=hkl_list[pos_cs->second].ResRange();
     dataset=hkl_list[pos_cs->second].xdataset(0);   // At present we only consider crystals with 1 dataset
     unitcell=dataset.cell();
     ccell_descr=clipper::Cell_descr(unitcell.UnitCell()[0],unitcell.UnitCell()[1],unitcell.UnitCell()[2],unitcell.UnitCell()[3],unitcell.UnitCell()[4],unitcell.UnitCell()[5]);
     ccell=clipper::Cell(ccell_descr);
     scala::Batch batch=hkl_list[pos_cs->second].batch(0);  // I'm using any batch (for instance number 0) just to recover crystal-to-detector distance
     CMtz::MTZBAT batchmtz=batch.batchdata();
     float ctoddist=batchmtz.dx[0];
     float wlength=dataset.wavelength();
     //std::cout << "WAVELENGTH " << wlength << std::endl;

     // Write to file "forR_macropar.dat"
     forR_file2 << "  " << std::setw(10) << pos_cs->second+1 << std::setw(10) << std::setprecision(3) << unitcell.UnitCell()[0]
                                                             << std::setw(10) << std::setprecision(3) << unitcell.UnitCell()[1]
                                                             << std::setw(10) << std::setprecision(3) << unitcell.UnitCell()[2]
                                                             << std::setw(9) << std::setprecision(2) << unitcell.UnitCell()[3]
                                                             << std::setw(9) << std::setprecision(2) << unitcell.UnitCell()[4]
                                                             << std::setw(9) << std::setprecision(2) << unitcell.UnitCell()[5]
                                                             << std::setw(11) << std::setprecision(5) << dataset.Mosaicity()
                                                             << std::setw(12) << std::setprecision(2) << ctoddist
                                                             << std::setw(10) << std::setprecision(5) << wlength
                                                             << std::endl;

     // Organise data (Phil's code. Needed to count unique reflections, etc)
     int nrefs=hkl_list[pos_cs->second].prepare();
     int  nobs=hkl_list[pos_cs->second].sum_partials();
     //std::cout << "NREFS:  " << nrefs << std::endl;
     //std::cout << "NOBS:  " << nobs << std::endl;

     // Create and output files for R, to determine radiation damage
     ofilename << "refs_" << spacegroup_class[i] << "_" << (pos_cs->second+1) << ".dat";
     filename=ofilename.str();
     forR_file << filename << std::endl;
     std::ofstream infofile(filename.c_str(),std::ios::out);   // Remember, first argument is a string literal, not a string!
     infofile.setf(std::ios::fixed);   // Permanent until reset
     for (int j=0;j < nrefs;j++)
     {
      reflection=hkl_list[pos_cs->second].get_reflection(j);
      int nobsnow=reflection.num_observations();
      for (int jj=0;jj < nobsnow;jj++)
      {
       observation=reflection.get_observation(jj);
       // Extract batch object  for each observation
       int nBatch=observation.Batch();
       int nBatch_serial=hkl_list[pos_cs->second].batch_serial(nBatch);
       scala::Batch batch=hkl_list[pos_cs->second].batch(nBatch_serial);
       std::string rispo="N";

       // Find average coordinates corresponding to different parts on detector   
       if(observation.IsFull()) rispo="S";
       infofile << std::setw(6) << reflection.hkl().h() << std::setw(6) << reflection.hkl().k() << std::setw(6) << reflection.hkl().l() 
                                                        << std::setw(10) << nBatch
                                                        << "   " << rispo << "   "
                                                        << std::setw(15) << std::setprecision(5) << observation.TotalFraction()
                                                        << std::setw(25) << std::setprecision(5) << observation.I()
                                                        << std::setw(25) << std::setprecision(5) << observation.sigI()
                                                        << std::setw(20) << std::setprecision(10) << std::sqrt(reflection.invresolsq())
                                                        << std::endl;
      }
     }
     ofilename.str("");  // Clear filename string before using it (to avoid chaining subsequent strings together)
     infofile.close();
    }
   }
  }
 }
 forR_file2.close();
 forR_file.close();

 // Run R script to perform statistical analysis.
 std::cout << "Performing statistical analysis (this might take a while!)........." << std::endl;
 int R_status;
 std::ostringstream R_command_line;
 if (Rscp == 0)
 {  
  R_command_line << "Rscript " << R_program;
 }
 else
 {
  R_command_line << "Rscript.exe " << R_program;
 }
 R_status=std::system((R_command_line.str()).c_str());   // The ostringstream object is first turned into a string and this into a cstring
 if (R_status != 0)
 {
  int nerr=13;
  throw nerr;
 }
 
 return;    // Remove this when using standard BLEND
}

// Function to output ascii files to be read by R (produced in Mode 1, dendrogram-only version)
void statistics_with_R2(std::vector<scala::hkl_unmerge_list>& hkl_list,std::multimap<int,int> sg_to_crystal,
                          std::vector<int> spacegroup_class,std::vector<int> crystal_flag,std::string R_program,int Rscp)
{
 // Build map container for bravais lattice number to symbol correspondence
 std::map<int,std::string> bl_number_to_symbol;
 std::map<int,std::string>::iterator pos_bl;
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(1,"aP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(2,"mP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(3,"mS"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(4,"oP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(5,"oS"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(6,"oF"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(7,"oI"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(8,"tP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(9,"tI"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(10,"hP"));
 //bl_number_to_symbol.insert(std::make_pair<int,std::string>(11,"hR"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(11,"hH"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(12,"cP"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(13,"cF"));
 bl_number_to_symbol.insert(std::make_pair<int,std::string>(14,"cI"));

 // Objects needed to output information from each crystal dataset
 scala::Xdataset dataset;
 scala::Scell unitcell;
 scala::ResoRange resorange;
 scala::reflection reflection;
 scala::observation observation;
 clipper::Cell_descr ccell_descr;
 clipper::Cell ccell;

 // File "forR_macropar.dat" containing cell parameters and mosaicity for all crystals
 std::ofstream forR_file2("forR_macropar.dat",std::ios::out);
 forR_file2.setf(std::ios::fixed);   // Permanent until reset

 // Ascii files with tables to be read by R. Names change from bravais lattice to bravais lattice
 std::string filename;
 std::ostringstream ofilename;

 // Iterators
 std::multimap<int,int>::iterator pos_cs;
 std::pair<std::multimap<int,int>::iterator,std::multimap<int,int>::iterator> pos_pair;

 // Loop over different space groups
 for (unsigned int i=0;i < spacegroup_class.size();i++)
 {
  if (spacegroup_class[i] != 0)
  {
   // Pair of iterators pointing at beginning and end of values (crystals) with same key (spacegroup number)
   pos_pair=sg_to_crystal.equal_range(spacegroup_class[i]);
   for (pos_cs=pos_pair.first;pos_cs != pos_pair.second;++pos_cs)
   {
    if (crystal_flag[pos_cs->second] == 0)
    {
     resorange=hkl_list[pos_cs->second].ResRange();
     dataset=hkl_list[pos_cs->second].xdataset(0);   // At present we only consider crystals with 1 dataset
     unitcell=dataset.cell();
     ccell_descr=clipper::Cell_descr(unitcell.UnitCell()[0],unitcell.UnitCell()[1],unitcell.UnitCell()[2],unitcell.UnitCell()[3],unitcell.UnitCell()[4],unitcell.UnitCell()[5]);
     ccell=clipper::Cell(ccell_descr);
     scala::Batch batch=hkl_list[pos_cs->second].batch(0);  // I'm using any batch (for instance number 0) just to recover crystal-to-detector distance
     CMtz::MTZBAT batchmtz=batch.batchdata();
     float ctoddist=batchmtz.dx[0];
     float wlength=dataset.wavelength();
     //std::cout << "WAVELENGTH " << wlength << std::endl;

     // Write to file "forR_macropar.dat"
     forR_file2 << "  " << std::setw(10) << pos_cs->second+1 << std::setw(10) << std::setprecision(3) << unitcell.UnitCell()[0]
                                                             << std::setw(10) << std::setprecision(3) << unitcell.UnitCell()[1]
                                                             << std::setw(10) << std::setprecision(3) << unitcell.UnitCell()[2]
                                                             << std::setw(9) << std::setprecision(2) << unitcell.UnitCell()[3]
                                                             << std::setw(9) << std::setprecision(2) << unitcell.UnitCell()[4]
                                                             << std::setw(9) << std::setprecision(2) << unitcell.UnitCell()[5]
                                                             << std::setw(11) << std::setprecision(5) << dataset.Mosaicity()
                                                             << std::setw(12) << std::setprecision(2) << ctoddist
                                                             << std::setw(10) << std::setprecision(5) << wlength
                                                             << std::endl;
    }
   }
  }
 }
 forR_file2.close();

 // Run R script to perform statistical analysis.
 int R_status;
 std::ostringstream R_command_line;
 if (Rscp == 0)
 {  
  R_command_line << "Rscript " << R_program;
 }
 else
 {
  R_command_line << "Rscript.exe " << R_program;
 }
 R_status=std::system((R_command_line.str()).c_str());   // The ostringstream object is first turned into a string and this into a cstring
 if (R_status != 0)
 {
  int nerr=13;
  throw nerr;
 }
 
 return;    // Remove this when using standard BLEND
}

// Overlaps matrix
//std::vector< std::vector<float> > build_overlaps_matrix(std::vector<scala::hkl_unmerge_list>& hkl_list,std::vector<int> crystals_choice)
//{
// int msize=crystals_choice.size();
// std::vector<std::vector<float> > overlaps_matrix;
//
// // Initialise matrix (so that all diagonal elements will be 1)
// for (int irow=0;irow < msize;irow++)
// {
//  overlaps_matrix.push_back(std::vector<float> ());
//  for (int icol=0;icol < msize;icol++)
//  {
//   overlaps_matrix[irow].push_back(1.0);
//  }
// }
//
// // Vector containing number of reflections for chosen crystals
// int nobs;
// scala::hkl_symmetry symmetry;
// std::vector<int> nrefs;
// for (int i=0;i < msize;i++)
// {
//  // Symmetry info extracted only once (all crystal belongs to same symmetry) 
//  if (i == 0) symmetry=hkl_list[crystals_choice[i]].symmetry();
//
//  // Number of unique reflections obtained while preparing data
//  nrefs.push_back(hkl_list[crystals_choice[i]].prepare());   // prepare() put all reflections in CCP4 asymmetric unit
//  nobs=hkl_list[crystals_choice[i]].sum_partials();   // Needed as part of data preparation
// }
// 
// // Work out off-diagonal elements of the overlaps matrix; upper triangle
// //int isym;
// scala::reflection this_refl1,this_refl2;
// scala::Hkl hkl1;
// for (int irow=0;irow < msize-1;irow++)
// {
//  for (int icol=irow+1;icol < msize;icol++)
//  {
//   // Insert here bit to compute reflections overlap. Use hkl_list[i] and hkl_list[j]
//   int noverlaps=0;
//   for (int iref=0;iref < nrefs[irow];iref++)
//   {
//    this_refl1=hkl_list[crystals_choice[irow]].get_reflection(iref);
//    hkl1=this_refl1.hkl();
//    int k=hkl_list[crystals_choice[icol]].get_reflection(this_refl2,hkl1);
//    if (k != -1) noverlaps++;
//   }
//   overlaps_matrix[irow][icol]=float(noverlaps)/float(std::max(nrefs[irow],nrefs[icol]));
//  }
// }
//
// // Copy lower triangle off-diagonal elements of the overlaps matrix
// for (int i=1;i < msize;i++)
// {
//  for (int j=0;j < i;j++)
//  {
//   overlaps_matrix[i][j]=overlaps_matrix[j][i];
//  }
// }
//
// return overlaps_matrix;
//}

// ISDIR: Returns 1 if the string is a directory,
//                0 if it's a file,
//               -1 otherwise
int isdir(const char *in_string)
{
DIR *dir; // <dirent.h>

// Access function 0 determines if an object exits, whether file or directory
if (access(in_string,0) != 0) // <io.h>
   return -1; // Function failed. No clean up, just return the "fail" flag.

// The object, whatever it is, exists

// Try to open the object as a directory stream
dir = opendir(in_string);
if (dir != NULL) 
   {
   // Object successfully opened as a directory stream
   // Close the open directory stream and return the "is a directory" flag
   closedir(dir);
   return 1;
   }

// All that's left is for it to be is a file. Return the "is a file" flag
return 0;
}
