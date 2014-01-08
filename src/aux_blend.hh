// Module containing functions to be used in conjunction with blend.cpp

#ifndef AUX_BLEND
#define AUX_BLEND

#include "Output.hh"
#include "hkl_unmerge.hh"
#include "hkl_controls.hh"
#include "hkl_symmetry.hh"
#include "mtz_unmerge_io.hh"
#include "controls.hh"
#include "openinputfile.hh"
#include "numbercomplete.hh"
#include "util.hh"               // This also includes cmath
#include <clipper/clipper.h>
#include <iomanip>
#include <dirent.h>

// Load crystals into a vector container of hkl_unmerge_list objects
std::vector<scala::hkl_unmerge_list> load_crystals(std::string,int);
 
// Check which crystal in a list is suitable for further processing
std::vector<int> label_crystals(const std::vector<scala::hkl_unmerge_list>&,std::vector<int>);

// Create and output summary table
void output_summary_table(std::vector<scala::hkl_unmerge_list>&,std::multimap<int,int>,std::vector<int>,std::vector<int>);

// Create ascii data files for R
void statistics_with_R(std::vector<scala::hkl_unmerge_list>&,std::multimap<int,int>,std::vector<int>,std::vector<int>,std::string);

// Create matrix of reciprocal space overlaps between couples of crystal datasets
std::vector< std::vector<float> > build_overlaps_matrix(std::vector<scala::hkl_unmerge_list>&,std::vector<int>);

// Check if in_string is directory, file or other
int isdir(const char*);

#endif
