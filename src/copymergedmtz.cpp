// copymergedmtz.cpp
//


// Clipper
#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"
using clipper::Message;
using clipper::Message_fatal;


// CCP4
#include "csymlib.h"    // CCP4 symmetry stuff
#include "ccp4_general.h"
#include "cmtzlib.h"    // CCP4 MTZlib headers (namespace CMtz)

#include "copymergedmtz.hh"
#include "scala_util.hh"
#include "mtz_utils.hh"
#include "hkl_symmetry.hh"
#include "cellgroup.hh"
#include "pointgroup.hh"
#include "printthings.hh"

using namespace CMtz;
using phaser_io::LOGFILE;

typedef std::pair<CMtz::MTZCOL*,CMtz::MTZCOL*> MtzColPair;

namespace MtzIO
{
  bool CopyMergedMTZ(const std::string& filename_in,
		     const std::string& filename_out,
		     const std::string& SpaceGroupName,
		     const scala::ReindexOp& reindex,
		     const bool& reduce,
		     const bool& verbose,
		     phaser_io::Output& output)
  // Copy all MTZ file changing hkl according to reindex
  // SpaceGroupName new SpaceGroupName
  //            if blank or "HKLIN" don't change
  // returns true if OK
  {
    // Sanity checks
    if (filename_in == "" || filename_out == "" ||
	filename_in == filename_out)
      Message::message(Message_fatal
		       ("CopyMergedMTZ invalid filenames"));

    // Don't read input file into memory
    CMtz::MTZ* mtzin = MtzGet(filename_in.c_str(), 0);

    if (verbose) {
      output.logTab(0,LOGFILE,
		    "\n===============================================================");
      output.logTab(0,LOGFILE,
		    "\nCopying merged MTZ file from "+filename_in+
		    "\n                          to "+filename_out);
      output.logTabPrintf(1,LOGFILE,"\n");
      output.logTabPrintf(1,LOGFILE,
			  "Reindexing operator         %s\n\n",
			  reindex.as_hkl().c_str());
      output.logTabPrintf(1,LOGFILE,
			  "Real space transformation   %s\n\n",
			  CCtbxSym::Reindex_as_xyz(reindex).c_str());
      
      if (!reduce) {
	output.logTab(0,LOGFILE,
		      "\nReflections are not reduced to asymmetric unit\n");
      }
    }

    // New symmetry
    std::string spacegroupname = SpaceGroupName; // local copy
    bool sgnamegiven = true;
    if (spacegroupname == "" || spacegroupname == "HKLIN") {
      spacegroupname = std::string(mtzin->mtzsymm.spcgrpname);
      sgnamegiven = false; // space group name not specified explicitly
      //^      std::cout << "MTZ sgname: " << spacegroupname << "\n";
    }

    bool reindexed = !reindex.IsIdentity();  // true if reindex op is not h,k,l
    bool translation = reindex.IsTranslation(); // true if translation component

    scala::hkl_symmetry NewSymm= scala::hkl_symmetry(spacegroupname);

    if (!sgnamegiven) { // space group name not given, so derive from MTZ name + reindex
      NewSymm.ChangeBasis(reindex);
      if (verbose) {
	output.logTabPrintf(1,LOGFILE,
			    "Reindexed space group : %s\n\n",
			    NewSymm.symbol_xHM().c_str());
      }
    } else {
      if (verbose) {
	output.logTabPrintf(1,LOGFILE,
			    "Space group : %s\n\n",
			    NewSymm.symbol_xHM().c_str());
      }
    }
    
    char history[MTZRECORDLENGTH];
    char date[11];
    char time[9];
    CCP4::ccp4_utils_date(date);
    CCP4::ccp4_utils_time(time);
    std::string text = "POINTLESS, "+std::string(date)+" "+std::string(time);
    if (reindexed) {
      text += " reindexed "+reindex.as_hkl();
    }
    text += " Space group "+spacegroupname;
    strcpy(history, text.c_str()); 

    int Nhist = MtzAddHistory(mtzin, &history, 1);

    // vectors of pointers to interesting columns
    std::vector<CMtz::MTZCOL*> DanoCol;
    //  ... and pairs
    std::vector<MtzColPair> PairCol;
    //  ... phase columns
    std::vector<CMtz::MTZCOL*> PhaseCol;
    //  ... ABCD columns (1st of set)
    std::vector<CMtz::MTZCOL*> ABCDCol;

    SpaceGroup newspacegroup =  NewSymm.GetSpaceGroup();

    // for cell constraints
    CCtbxSym::CellGroup CG(newspacegroup);

    // Process all crystals, changing cells by reindex
    // Loop crystals
    scala::Scell newcell;
    if (verbose) {
      output.logTab(0,LOGFILE, "Unit Cells: ");
    }
    for (int x=0; x < CMtz::MtzNxtal(mtzin); x++) {
      CMtz::MTZXTAL* xtl = CMtz::MtzIxtal(mtzin,x);
      newcell = CG.constrain(scala::Scell(xtl->cell).change_basis(reindex));
      for (int i=0;i<6;i++)
	xtl->cell[i] = newcell[i];
      if (verbose) {
	output.logTabPrintf(1,LOGFILE,"Crystal: %s\n",xtl->xname);
	output.logTab(1,LOGFILE, "Cell: "+scala::FormatCell(newcell));
      }
      // Loop datasets
      for (int set=0;set<xtl->nset;set++) {
	CMtz::MTZSET* pset = xtl->set[set];
	// Loop columns
	int Ncol = pset->ncol;
	std::vector<bool> ColFree(Ncol,true);
	for (int col=0;col<Ncol;col++) {  // loop columns
	  if (ColFree[col]) { // column not yet assigned
	    CMtz::MTZCOL* pcol = pset->col[col];
	    char ctype = pcol->type[0]; // column type

	    // Check for any phase columns (type P)
	    // or phase coefficient columns (type A)
	    // These cannot be reindexed by this routine if there is a translation
	    if (translation && (ctype == 'P' || ctype == 'A')) {
	      PrintError
  ("**** Can't apply translation reindex to Phase or Phase Coefficient column ****",
		 output);
	      output.logTabPrintf(0,LOGFILE,"   Column label %s Type %c\n\n",
				  pcol->label, pcol->type[0]);
	      
	      return false;
	    }

	    if (reduce) {
	      // mark phase & ABCD columns for later processing
	      if  (ctype == 'P') { // phase column
		PhaseCol.push_back(pcol);
	      } else if (ctype == 'A') {
		// Check for 4 columns of type A
		bool found = true;
		if (col <= Ncol-4) {
		  for (int col2=col+1;col2<col+4;col2++) {
		    if (pset->col[col2]->type[0] == ctype) { // OK type A
		      ColFree[col2] = false; // mark as used
		    } else {
		      found = false; // wrong type
		    }
		  }
		}
 		if (!found) {
		  Message::message(Message_fatal
				   ("CopyMergedMTZ: incomplet ABCD set"+
				    std::string(pcol->label)+" type"+ctype));
		}
		ABCDCol.push_back(pcol); // OK
	      } else if (ctype == 'D') { // Check for anomalous difference columns
		DanoCol.push_back(pcol);
	      }
	      // Check for columns of type F+, F- (G, sigma L),
	      // or I+, I- (K, sigma M)
	      else if (ctype == 'G' || ctype == 'L' || ctype == 'K' || ctype == 'M') {
		// Search forward for a matching column
		bool found = false;
		if (col < Ncol-1) {
		  for (int col2=col+1;col2<Ncol;col2++)	{
		    if (pset->col[col2]->type[0] == ctype) {
		      // xtl/set/col & col2 are a matched pair to be swapped
		      PairCol.push_back(MtzColPair(pcol, pset->col[col2]));
		      ColFree[col2] = false;
		      found = true;
		      break;
		    }
		  }
		}
		if (!found) {
		  Message::message(Message_fatal
				   ("CopyMergedMTZ: unmatched anomalous column"+
				    std::string(pcol->label)+" type"+ctype));		  
		}
	      } // end F/I+- check
	    } // end reduce
	  }  // end column not yet assigned
	}  // end column loop
      }  // end dataset loop
    }  // end crystal loop

    if (verbose) {
      if (DanoCol.size() > 0) {
	output.logTabPrintf(0,LOGFILE,
			    "\nColumns containing anomalous differences will be negated\n");
	output.logTabPrintf(0,LOGFILE,
			    " if hand of indices hkl is changed in reducing to asymmetric unit\n");
	output.logTabPrintf(0,LOGFILE,
			    "\n   Anomalous column labels: ");
	for (size_t i=0;i<DanoCol.size();i++) {
	  output.logTabPrintf(0,LOGFILE," %s",DanoCol[i]->label);
	}
	output.logTab(0,LOGFILE,"\n");
      }
      if (PairCol.size() > 0) {
	output.logTabPrintf(0,LOGFILE,
			    "\nAnomalous pairs will be swapped if hand");
	output.logTabPrintf(0,LOGFILE,
			    " of indices hkl is changed in reducing to asymmetric unit\n\n");
	output.logTabPrintf(0,LOGFILE,
			    " Anomalous Column              Label1               Label2  Type\n\n");
	for (size_t i=0;i<PairCol.size();i++) {
	  char ctype = PairCol[i].first->type[0];
	  std::string content;
	  if (ctype == 'G') {content = "      F+, F- ";}
	  else if (ctype == 'L') {content = "sigma(F+, F-)";}
	  else if (ctype == 'K')  {content = "      I+, I- ";}
	  else if (ctype == 'M') {content = "sigma(I+, I-)";}
	  output.logTabPrintf(0,LOGFILE,
			      "%16s %20s %20s %c\n",
			      content.c_str(),
			      PairCol[i].first->label,
			      PairCol[i].second->label, ctype);
	}
	output.logTab(0,LOGFILE,"\n");
      } // end paircol
      if (PhaseCol.size() > 0) {
	output.logTabPrintf(0,LOGFILE,
	  "\nPhases may be changed if indices hkl are changed in reducing to asymmetric unit\n");
	output.logTabPrintf(0,LOGFILE,
			    "\n   Phase column labels: ");
	for (size_t i=0;i<PhaseCol.size();i++) {
	  output.logTabPrintf(0,LOGFILE," %s",PhaseCol[i]->label);
	}
	output.logTab(0,LOGFILE,"\n");
      }
      if (ABCDCol.size() > 0) {
	output.logTabPrintf(0,LOGFILE,
	  "\nABCDs may be changed if indices hkl are changed in reducing to asymmetric unit\n");
	output.logTabPrintf(0,LOGFILE,
			    "\n   ABCD column labels: ");
	for (size_t i=0;i<ABCDCol.size();i++) {
	  output.logTabPrintf(0,LOGFILE," %s %s %s %s",
			      ABCDCol[i]->label,(ABCDCol[i]+1)->label,
			      (ABCDCol[i]+2)->label,(ABCDCol[i]+3)->label);
	}
	output.logTab(0,LOGFILE,"\n");
      }

    } // end verbose
    
    // Update file cell
    // ?? don't need to, first crystal cell is used

        // Check for rhombohedral lattice (H or R)
    char HorR = 'H'; // default H setting
    if (NewSymm.lattice_type() == 'H' || NewSymm.lattice_type() == 'R') {
      if (RhombohedralAxes(newcell.UnitCell())) { // true if not H
	HorR = 'R';
      }}

    mtzin->mtzsymm = spg_to_mtz(newspacegroup, HorR);

    int Ncolumns = CMtz::MtzNcol(mtzin);
    // Column offsets ( = number-1) for H,K,L columns
    int mtzH = MtzColLookup(mtzin, "H")->source - 1;
    int mtzK = MtzColLookup(mtzin, "K")->source - 1;
    int mtzL = MtzColLookup(mtzin, "L")->source - 1;
  
    mtzin->fileout = CMtz::MtzOpenForWrite(filename_out.c_str());

    int isym;
    int Nneg = 0;
    int Nswap = 0;
    int Nphasechange = 0;
    int Nref = 0;
    int NfractIdx = 0;
    std::vector<float> cols(Ncolumns);
    
    // Loop all reflections reindexing
    for (unsigned mtzr = 0; mtzr < unsigned(mtzin->nref); mtzr++) {
      // NB the construction &*vector.begin() gives a pointer to
      //    the data part of a vector, which must have already been
      //    set up with the correct length
      CMtz::MtzRrefl(mtzin->filein, Ncolumns, &*cols.begin());
      
      scala::Hkl oldhkl(Nint(cols[mtzH]), Nint(cols[mtzK]), Nint(cols[mtzL]));
      // Reindex
      scala::Hkl newhkl;
      // returns false if non-integral indices
      bool HklOK = oldhkl.change_basis(newhkl, reindex);
      if (!HklOK) {
	NfractIdx++;
	continue;
      }
	
      if (reduce) {
	// Reduce to asymmetric unit
	oldhkl = newhkl;
	newhkl = NewSymm.put_in_asu(oldhkl, isym);
	bool hklprint = true;
	// If hand has changed (isym even), do any necessary swaps
	if (isym%2 == 0 && !NewSymm.is_centric(newhkl)) {
	  // Anomalous difference swaps, acentric reflections only
	  //  NB Data has to be swapped as there is no reliable way just to
	  //  change the column names
	  if (DanoCol.size() > 0) {  // DANO column
	    for (size_t i=0;i<DanoCol.size();i++) {
	      // Dano, negate
	      int j = DanoCol[i]->source - 1;
	      cols[j] = -cols[j];
	      if (verbose) {
		if (Nneg < 20) {
		  if (hklprint) {
		    output.logTab(0,LOGFILE,newhkl.format());
		    hklprint = false;
		  }
		  output.logTabPrintf(0,LOGFILE,
				      "Negating Danom column %s\n",
				      DanoCol[i]->label);
		}
	      }
	    }
	    if (verbose && Nneg == 20) {
	      output.logTab(0,LOGFILE, "   ... more negated");
	    }
	    Nneg++;
	  } // end DANO column
	  if (PairCol.size() > 0) { // Paired +/- column
	    for (size_t i=0;i<PairCol.size();i++) {
	      int j = PairCol[i].first->source - 1;
	      int k = PairCol[i].second->source - 1;
	      float temp = cols[j];
	      cols[j] = cols[k];
	      cols[k] = temp;
	      if (verbose) {
		if (Nswap < 20) {
		  if (hklprint) {
		    output.logTab(0,LOGFILE,newhkl.format());
		    hklprint = false;
		  }
		  output.logTabPrintf(0,LOGFILE,
				      "Swapping columns %20s %20s\n",
				      PairCol[i].first->label,
				      PairCol[i].second->label);
		}
	      }
	    }
	    if (verbose &&Nswap == 20) {
	      output.logTab(0,LOGFILE, "   ... more swapped");
	    }
	    Nswap++;
	  } // end paired column
	} // end swap
	// --- 
	if (PhaseCol.size() > 0 || ABCDCol.size() > 0) { // Phase columns
	  // Phase shift for +hkl, degrees
	  double phaseshift = oldhkl.HKL().sym_phase_shift(newspacegroup.SymopFromIsym(isym));
	  if (verbose) {
	    if (std::abs(phaseshift) > 0.001) {
	      if (Nphasechange < 20) {
		if (hklprint) {
		  output.logTab(0,LOGFILE,newhkl.format());
		  hklprint = false;
		}
		std::string s = "Phase shift "+clipper::String(clipper::Util::rad2d(phaseshift));
		if (isym%2 == 0) s += " negated";
		output.logTab(0,LOGFILE, s);
		Nphasechange++;
	      }
	    }
	  }
	  if (PhaseCol.size() > 0) { // phase columns
	    for (size_t i=0;i<PhaseCol.size();i++) {
	      // Phase column, shift
	      int j = PhaseCol[i]->source - 1;
	      cols[j] = cols[j] + clipper::Util::rad2d(phaseshift);
	      if (isym%2 == 0) { // change of hand, negate phase
		cols[j] = -cols[j];
	      }
	    }
	  }  // end phase columns
	  
	  if (ABCDCol.size() > 0) { // ABCD columns
	    for (size_t i=0;i<ABCDCol.size();i++) {
	      // ABCD column, shift
	      int j = ABCDCol[i]->source - 1;
	      // Make ABCD object
	      clipper::datatypes::ABCD<float> HL(cols[j], cols[j+1], cols[j+2], cols[j+3]); 
	      HL.shift_phase(phaseshift);  // shift phase
	      if (isym%2 == 0) {HL.friedel();} // Friedel symmetry if required
	      cols[j]   = HL.a();
	      cols[j+1] = HL.b();
	      cols[j+2] = HL.c();
	      cols[j+3] = HL.d();
	    }
	  } // end ABCD
	} // end phase columns
      }  // end reduce
      
      cols[mtzH] = newhkl[0];
      cols[mtzK] = newhkl[1];
      cols[mtzL] = newhkl[2];
      
      CMtz::MtzWrefl(mtzin->fileout, Ncolumns, &*cols.begin());
      Nref++;
    } // end loop reflections
    
    // Reset global resolution limits
    float minres, maxres;
    CMtz::MtzResLimits(mtzin, &minres, &maxres);
    mtzin->resmax_out = maxres;
    mtzin->resmin_out = minres;

    if (! MtzPut(mtzin, filename_out.c_str()))
      Message::message(Message_fatal
		       ("CopyMergedMTZ failed to write file"+filename_out));
    
    MtzFree(mtzin);

    output.logTabPrintf(0,LOGFILE,
			 "\n%8d reflections copied to output file\n",Nref);
    
    if (NfractIdx > 0) {
      output.logTabPrintf(0,LOGFILE,
			  "\n%8d reflections with fractional indices discarded\n",
			  NfractIdx);
    }
    output.logTabPrintf(0,LOGFILE,"\n");

    return true;
  }
}
