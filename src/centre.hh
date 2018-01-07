// centre.hh

#ifndef CENTRE_HEADER
#define CENTRE_HEADER

#include "hkl_unmerge.hh"
#include "normalise.hh"

namespace scala
{
  // On entry:
  //  hkl_list            data list
  //  LGname              Laue group name
  //  reindex             reindex operator
  //  MinIsigratio        for resolution cut-off
  //  GridMax             maximum indices for grid search
  //  hklout_filename     filename for optional output file
  //  output              output object
  void FindCentre(hkl_unmerge_list& hkl_list,
		  const std::string& LGname,
		  const ReindexOp& reindex,
		  const double& MinIsigRatio,
		  const clipper::Vec3<int>& GridMax,
		  const std::string& hklout_filename,
		  phaser_io::Output& output);
  // Returns reindexing operation to fix up data
  //
  // On entry:
  //  hkl_list            data list
  //  MinIsigratio        for resolution cut-off
  //  GridMax             maximum indices for grid search
  //  output              output object
  ReindexOp Centre(hkl_unmerge_list& hkl_list,
		   const double& MinIsigRatio,
		   const clipper::Vec3<int>& GridMax,
		   phaser_io::Output& output);
}
#endif
