cmake_minimum_required(VERSION 2.8)
project(blend C CXX)
find_package(CCP4 REQUIRED mmdb ccp4c clipper-core clipper-mmdb clipper-minimol
                           clipper-contrib clipper-ccp4 cctbx)

include_directories(${CCP4_INCLUDE_DIRS} ${CCP4C_INCLUDE_DIR}/ccp4)

#file(GLOB BLEND_SOURCES src/*.cpp)
set(BLEND_SOURCES
src/blend.cpp                   src/aux_blend.cpp
src/Output.cpp                  src/hkl_unmerge.cpp           src/openinputfile.cpp
src/controls.cpp                src/icering.cpp               src/range.cpp
src/hash.cpp                    src/lattice.cpp               src/rotation.cpp
src/hkl_controls.cpp            src/matvec_utils.cpp          src/runthings.cpp
src/hkl_datatypes.cpp           src/mtz_unmerge_io.cpp        src/scala_util.cpp
src/hkl_symmetry.cpp            src/observationflags.cpp      src/mtz_utils.cpp
src/jiffy.cpp                   src/writeunmergedmtz.cpp      src/cellgroup.cpp
src/normalise.cpp               src/pointgroup.cpp            src/linearlsq.cpp
src/zone.cpp                    src/average.cpp               src/globalcontrols.cpp
src/latsym.cpp                  src/getsubgroups.cpp          src/probfunctions.cpp
src/score_datatypes.cpp         src/spline.cpp                src/numbercomplete.cpp
src/CCP4base.cpp                src/cctbx_utils.cpp           src/centre.cpp
src/checkcompatiblesymmetry.cpp src/columnlabels.cpp          src/copymergedmtz.cpp
src/datatypes.cpp               src/eprob.cpp                 src/Errors.cpp
src/expandfilename.cpp          src/fileread.cpp              src/hkl_merged_list.cpp
src/InputAll.cpp                src/intensitystatistics.cpp   src/interpretcommandline.cpp
src/keywords.cpp                src/listdir.cpp               src/mtz_merge_io.cpp
src/pgscore.cpp                 src/Preprocessor.cpp          src/scsignificance.cpp
src/setscores.cpp               src/spacegroupreindex.cpp     src/string_util.cpp
src/sysabstest.cpp              src/outputunmergeddata.cpp    src/printthings.cpp
src/io_files.cpp                src/filetype.cpp              src/readallhklinfiles.cpp
src/comparestring.cpp           src/sysabszones.cpp           src/sysabsstats.cpp
src/makehkl_listscompatible.cpp src/xds_unmerge.cpp           src/sca_unmerge.cpp
src/testlauegroup.cpp           src/checklatticecentering.cpp src/checknullbatches.cpp
src/normaliseunmerge.cpp        src/testindex.cpp             src/fix_merged_scores.cpp
src/getCCsigfac.cpp             src/pointless.cpp             src/choosesolution.cpp
src/writesummedlist.cpp         src/completeness.cpp          src/timer.cpp
src/twintest.cpp                src/tablegraph.cpp)

add_executable(blend ${BLEND_SOURCES})
target_link_libraries(blend ${CCP4_LIBRARIES})
install(TARGETS blend DESTINATION bin)
install(DIRECTORY python R DESTINATION share/blend/)
install(FILES license.txt DESTINATION share/blend/)

