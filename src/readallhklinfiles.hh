// readallhklinfiles.hh

#ifndef READALLHKLINFILES_HEADER
#define READALLHKLINFILES_HEADER

#include "io_files.hh"
#include "hkl_controls.hh"
#include "controls.hh"
#include "mtz_unmerge_io.hh"
#include "globalcontrols.hh"
#include "hkl_datatypes.hh"
#include "hkl_merged_list.hh"
#include "hkl_unmerge.hh"
#include "Output.hh"


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
		    hkl_unmerge_list& hkl_list);
  // Add this file into hkl_list, if compatible with whatever is already there
  // Returns true if file was read
  //
  // fileSeries is index number for file or file-series (from 1) for batch
  //  exclusion in file_sel (MTZ files only)
  // non-null cell overrides file cell, or adds it to SCALEPACK files
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
			const RefListType& HklinIsMerged,
			phaser_io::Output& output,
			hkl_unmerge_list& hkl_list);
  //--------------------------------------------------------------
  // Extract maximum resolution from HKLIN files
  // Only works for MTZ input: in other cases returns a default value
  double GetResolutionFromHKLIN(const IO_files& AllFiles);
}

#endif
