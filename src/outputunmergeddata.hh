// outputunmergeddata.hh

#ifndef OUTPUTUNMERGEDDATA_HEADER
#define OUTPUTUNMERGEDDATA_HEADER

#include <string>

#include "hkl_unmerge.hh"
#include "hkl_symmetry.hh"
#include "hkl_datatypes.hh"
#include "Output.hh"

namespace scala{
//--------------------------------------------------------------
void OutputUnmergedData(hkl_unmerge_list& hkl_list,
                        const std::string& hklout_filename,
                        const hkl_symmetry& SPGsymm,
                        const ReindexOp& reindex_sg,
			const bool& reindexedTestData,
                        const std::string& Title,
                        const std::string& groupflag,
			phaser_io::Output& output);
}

#endif
