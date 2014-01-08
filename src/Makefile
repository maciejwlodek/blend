# Compilers
CC = cc
CXX= g++

# Compiler flags
#OPT = -O2 -Wall
OPT = -O2 -Wall -Wl,-rpath=${CCP4}/lib
CFLAGS= ${OPT}
CXXFLAGS=  ${OPT}

# Link flags
LFLAGS= ${OPT} -static-libgcc -static-libstdc++

# CCP4 locations
CCP4loc=${CCP4}/include/ccp4

# CCTBX locations
CCTBX=${CCP4}/include/cctbx

# Clipper location
CLPR=${CCP4}/include

#<<<<< Includes
# CCP4
ICCP4loc=-I${CCP4loc}

# Clipper
ICLPR=-I${CLPR}

# cctbx
ICCTBX=-I${CCTBX}

INCLUDE=-I. ${ICCTBX} ${ICLPR} ${ICCP4loc}

#<<<<< Libraries

# CCP4
LCCP4=-L${CLIB} -lrfftw -lfftw -lccp4f -lccp4c -lmmdb -lm

# Clipper
LCLPR=-L${CLIB} -lclipper-ccp4 -lclipper-contrib -lclipper-minimol -lclipper-mmdb -lclipper-core 

# cctbx
LCCTBX= -L${CLIB} -lcctbx 

LDLIBS=${LCCTBX} ${LCLPR} ${LCCP4}


CPPFLAGS= ${INCLUDE}
CCFLAGS = ${INCLUDE} 

EXE = blend

OBJECTS = ${EXE}.o \
Output.o                    hkl_unmerge.o               openinputfile.o \
controls.o                  icering.o                   range.o \
hash.o                      lattice.o                   rotation.o \
hkl_controls.o              matvec_utils.o              runthings.o \
hkl_datatypes.o             mtz_unmerge_io.o            scala_util.o \
hkl_symmetry.o              observationflags.o          mtz_utils.o \
jiffy.o                     writeunmergedmtz.o          cellgroup.o \
normalise.o                 pointgroup.o                linearlsq.o \
zone.o                      average.o                   globalcontrols.o \
latsym.o                    getsubgroups.o              probfunctions.o \
score_datatypes.o           spline.o                    numbercomplete.o \
CCP4base.o                  cctbx_utils.o               centre.o \
checkcompatiblesymmetry.o   columnlabels.o              copymergedmtz.o \
datatypes.o                 eprob.o                     Errors.o \
expandfilename.o            fileread.o                  hkl_merged_list.o \
InputAll.o                  intensitystatistics.o       interpretcommandline.o \
keywords.o                  listdir.o                   mtz_merge_io.o \
pgscore.o                   Preprocessor.o              scsignificance.o \
setscores.o                 spacegroupreindex.o         string_util.o \
sysabstest.o                outputunmergeddata.o        printthings.o \
io_files.o                  filetype.o                  readallhklinfiles.o \
comparestring.o             sysabszones.o               sysabsstats.o \
makehkl_listscompatible.o   xds_unmerge.o               sca_unmerge.o \
testlauegroup.o             checklatticecentering.o     checknullbatches.o \
normaliseunmerge.o          testindex.o                 fix_merged_scores.o \
getCCsigfac.o               pointless.o                 choosesolution.o \
writesummedlist.o           completeness.o              timer.o \
twintest.o                  tablegraph.o                aux_blend.o


${EXE}: ${OBJECTS}
	${CXX} ${LFLAGS} -o ../${EXE} ${OBJECTS} ${LDLIBS}


#.SUFFIXES: .cpp .o
