// writeunmergedmtz.cpp
//


// Clipper
#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"
using clipper::Message;
using clipper::Message_fatal;

#define ASSERT assert
#include <assert.h>


// CCP4
#include "csymlib.h"    // CCP4 symmetry stuff
#include "ccp4_general.h"
#include "cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)

#include "hkl_datatypes.hh"
#include "hkl_controls.hh"
#include "scala_util.hh"
#include "mtz_utils.hh"
#include "hkl_symmetry.hh"
#include "writeunmergedmtz.hh"
#include "openinputfile.hh"
#include "observationflags.hh"
#include "cellgroup.hh"
#include "version.hh"

using namespace CMtz;

namespace MtzIO
{
  //--------------------------------------------------------------
  void OptAddCol(const bool& coln,
		 MTZCOL* col[], int& ic, MTZ* mtzout, MTZSET* baseset,
		 const char* label, const char* type)
    // Conditional column addition, only if coln > 0
  {
    if (coln)
	col[ic++] = MtzAddColumn(mtzout, baseset, label, type);
  }
  //--------------------------------------------------------------
  bool WriteUnmergedMTZ(scala::hkl_unmerge_list& hkl_list,
			const scala::hkl_symmetry& NewSymm,
			const std::string& filename_out,
			const std::string& title)
    // Write unmerged MTZ file from hkl_unmerge_list object
    // returns true if OK
  {
    // Organise list if necessary
    int Nrefl = hkl_list.prepare();
    Nrefl = Nrefl;

    // Initialise MTZ data structure and output file
    MTZ* mtzout = MtzMalloc(0,0);
    mtzout->refs_in_memory = 0;   // not in memory
    mtzout->fileout = MtzOpenForWrite(filename_out.c_str());
    if (mtzout->fileout == NULL)
      {Message::message(Message_fatal("Can't open file "+filename_out));}


    // Write title
    ccp4_lwtitl(mtzout, title.c_str(), 0);

    // Symmetry
    // Check for rhombohedral lattice (H or R)
    char HorR = 'H'; // default H setting
    if (NewSymm.lattice_type() == 'H' || NewSymm.lattice_type() == 'R') {
      if (RhombohedralAxes(hkl_list.Cell().UnitCell())) { // true if not H
	HorR = 'R';
      }}
    mtzout->mtzsymm = spg_to_mtz(NewSymm.GetSpaceGroup(), HorR);
    // for cell constraints
    CCtbxSym::CellGroup CG(NewSymm.GetSpaceGroup());

    // Crystals & datasets
    int Nxd = hkl_list.num_datasets();
    std::vector<std::string> xtl;
    float ucell[6];
    std::string lastxname = "";
    MTZXTAL* xtal;
    MTZSET* set;

    // Put in HKL_base dataset explicitly
    scala::Scell HKLcell = CG.constrain(hkl_list.cell());  // constrained cell
    for (int i=0;i<6;i++)
      {ucell[i] = HKLcell[i];}
    xtal = MtzAddXtal(mtzout, "HKL_base", "HKL_base",
		      ucell);
    MTZSET* baseset = MtzAddDataset(mtzout, xtal, "HKL_base", 0.0);

    for (int jxd=0; jxd<Nxd;jxd++)
      // loop xdatasets
      // NB assume that datasets belonging to a crystal are contiguous
      {
	std::string xname = hkl_list.xdataset(jxd).pxdname().xname();
	if (xname != lastxname)
	  // New crystal
	  {
	    lastxname = xname;
	    HKLcell = CG.constrain(hkl_list.xdataset(jxd).cell());
	    for (int i=0;i<6;i++) {
	      ucell[i] = HKLcell[i];
	    }
	    xtal = MtzAddXtal(mtzout, xname.c_str(),
			      hkl_list.xdataset(jxd).pxdname().pname().c_str(),
			      ucell);
	  }
	// Dataset
	set = MtzAddDataset(mtzout, xtal,
			    hkl_list.xdataset(jxd).pxdname().dname().c_str(),
			    hkl_list.xdataset(jxd).wavelength());
      }

    // Columns, all in base dataset
    data_flags  col_sel = hkl_list.DataFlags();

    // How many columns?
    MTZCOL* col[MAXNCOLUMNS];  // MAXNCOLUMNS defined in openinputfile.hh

    int ic=0;
    col[ic++] = MtzAddColumn(mtzout, baseset, "H", "H");
    col[ic++] = MtzAddColumn(mtzout, baseset, "K", "H");
    col[ic++] = MtzAddColumn(mtzout, baseset, "L", "H");
    col[ic++] = MtzAddColumn(mtzout, baseset, "M/ISYM", "Y");
    col[ic++] = MtzAddColumn(mtzout, baseset, "BATCH", "B");
    col[ic++] = MtzAddColumn(mtzout, baseset, "I", "J");
    col[ic++] = MtzAddColumn(mtzout, baseset, "SIGI", "Q");
    // Optional columns, add in only if present (ie read from input file)
    // also increment ic
    OptAddCol(col_sel.is_Ipr, col, ic, mtzout, baseset, "IPR", "J");
    OptAddCol(col_sel.is_sigIpr, col, ic, mtzout, baseset, "SIGIPR", "Q");
    OptAddCol(col_sel.is_fractioncalc, col, ic, mtzout, baseset, "FRACTIONCALC", "R");
    OptAddCol(col_sel.is_Xdet, col, ic, mtzout, baseset, "XDET", "R");
    OptAddCol(col_sel.is_Ydet, col, ic, mtzout, baseset, "YDET", "R");
    OptAddCol(col_sel.is_Rot, col, ic, mtzout, baseset, "ROT", "R");
    OptAddCol(col_sel.is_Width, col, ic, mtzout, baseset, "WIDTH", "R");
    OptAddCol(col_sel.is_LP, col, ic, mtzout, baseset, "LP", "R");
    OptAddCol(col_sel.is_Mpart, col, ic, mtzout, baseset, "MPART", "I");
    OptAddCol(col_sel.is_ObsFlag, col, ic, mtzout, baseset, "FLAG", "I");
    OptAddCol(col_sel.is_BgPkRatio, col, ic, mtzout, baseset, "BGPKRATIOS", "R");
    OptAddCol(col_sel.is_time, col, ic, mtzout, baseset, "TIME", "R");
    OptAddCol(col_sel.is_scale, col, ic, mtzout, baseset, "SCALE", "R");
    OptAddCol(col_sel.is_sigscale, col, ic, mtzout, baseset, "SIGSCALE", "R");
    int NumCol = ic;
    
    // List is sorted on the first 5 columns
    MtzSetSortOrder(mtzout, col);

    
    // History    
    char history[MTZRECORDLENGTH];
    char date[11];
    char time[9];
    CCP4::ccp4_utils_date(date);
    CCP4::ccp4_utils_time(time);
    std::string text = "Pointless, version "+std::string(PROGRAM_VERSION)+
      ", run on "+std::string(date)+" at "+
      std::string(time);
    strcpy(history, text.c_str()); 
    int Nhist = MtzAddHistory(mtzout, &history, 1);

    // Count observation parts in each batch
    std::vector<int> nobsbatch(hkl_list.num_batches(),0);
    // Write all data
    float data[MAXNCOLUMNS];
    for (int i=0;i<hkl_list.num_parts();i++)  {
      scala::observation_part part = hkl_list.find_part(i);
      scala::Hkl hkl = part.hkl();
      data[0] = hkl.h();
      data[1] = hkl.k();
      data[2] = hkl.l();
      // Packed M/ISYM, M = 1 for partial
      int isym = part.isym();
      int M_Isym = isym;
      if (part.Npart() != 1) M_Isym = 256 + isym;
      data[3] = M_Isym;
      data[4] = part.batch();
      data[5] = part.Ic();
      data[6] = part.sigIc();
      // count observation parts by batch serial
      nobsbatch.at(hkl_list.batch_serial(Nint(data[4])))++;

     // Optional columns
      ic = 7;
      if (col_sel.is_Ipr) data[ic++] = part.Ipr();
      if (col_sel.is_sigIpr) data[ic++] = part.sigIpr();
      if (col_sel.is_fractioncalc) data[ic++] = part.fraction_calc();
      if (col_sel.is_Xdet) data[ic++] = part.Xdet();
      if (col_sel.is_Ydet) data[ic++] = part.Ydet();
      if (col_sel.is_Rot) data[ic++] = part.phi();
      if (col_sel.is_Width) data[ic++] = part.width();
      if (col_sel.is_LP) data[ic++] = part.LP();
      if (col_sel.is_Mpart) {
	int Mpart = 0;
	if (part.Npart() > 1) Mpart = 100*part.Npart() + part.Ipart();
	data[ic++] = Mpart;
      }
      if (col_sel.is_ObsFlag) data[ic++] = part.ObsFlag().Flags();
      if (col_sel.is_BgPkRatio) data[ic++] = part.ObsFlag().BgPk();
      if (col_sel.is_time) data[ic++] = part.time();
      // Optional scale: dummy for now
      if (col_sel.is_scale) data[ic++] = 1.0;
      if (col_sel.is_sigscale) data[ic++] = 0.0;
      ASSERT (ic == NumCol);
      
      ccp4_lwrefl(mtzout, data, col, NumCol, i+1);
    }

    // Batches: construct linked list in mtzout object
    MTZBAT* batch;
    MTZBAT* previous_batch = 0;
    int nbat = 0;

    for (int jbat=0;jbat<hkl_list.num_batches();jbat++)  {
      // Only output accepted batches
      if (hkl_list.batch(jbat).Accepted() && nobsbatch[jbat] > 0) {
	batch = MtzMallocBatch(); // make space for batch data
	
	if (nbat == 0) {
	  mtzout->batch = batch; // pointer to first batch
	} else {
	  // Link previous batch to this one
	  previous_batch->next = batch;
	}
	nbat++;
	previous_batch = batch;

	*batch = hkl_list.batch(jbat).batchdata(); // copy data
	// reset time limits if no time data
	if (!col_sel.is_time) {
	  batch->time1 = 0.0;
	  batch->time2 = 0.0;
	}
	// Fix up cell, constrain to symmetry
	scala::Scell Bcell = scala::Scell(batch->cell);
	// Is this a valid cell?
	bool valid = true;
	for (int i=0;i<6;i++) {
	  if (Bcell[i] <= 0.001) {valid = false;}
	}
	if (valid) {
	  scala::Scell Bcell = CG.constrain(scala::Scell(batch->cell));
	  for (int i=0;i<6;i++) {batch->cell[i] = Bcell[i];}
	}
	// Update NBsetid if required, ie index in dataset list
	// look it up in new mtzout structure
	std::string path = "/"+
	  hkl_list.xdataset(hkl_list.batch(jbat).index()).pxdname().xname()+"/"+
	  hkl_list.xdataset(hkl_list.batch(jbat).index()).pxdname().dname();
	batch->nbsetid = MtzSetLookup(mtzout, path.c_str())->setid;  // setid
	batch->next = NULL;  // for last one
      }
    }

    ccp4_lhprt(mtzout,1);
    if (!MtzPut(mtzout, " ")) {
      Message::message(Message_fatal("Can't write file "+filename_out));
    }
    MtzFree(mtzout);

    return true;
  }
}
