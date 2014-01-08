// io_files.cpp
//

#include "io_files.hh"
#include "scala_util.hh"
#include "string_util.hh"
#define ASSERT assert
#include <assert.h>

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;

namespace scala
{
  HKLIN_data::HKLIN_data(const std::string& filename_in,
			 const ReflectionFileType& type_in,  // MTZ, XDS, SCA
			 const PxdName& pxdname_in,
			 const int& series_in)
    : filename(filename_in), type(type_in), pxdname(pxdname_in),
      series(series_in)
  {
    reindex_set = false;
    reindex = ReindexOp();
    cc = correl_coeff();
    unmerged = +1;  // Assume unmerged by default
}
  //--------------------------------------------------------------
  IO_files::IO_files() 
    : HKLREF_filename(""), HKLOUT_filename(""),
      XMLOUT_filename(""), XYZIN_filename(""), seriescount(0)
  {
    HKLIN_things.clear();
  }
  //--------------------------------------------------------------
  int IO_files::AddHKLINfilename(const std::string& Filename,
				 const ReflectionFileType& Type)
  {
    return AddHKLINfilename(Filename, Type, PxdName());
  }
  //--------------------------------------------------------------
  int IO_files::AddHKLINfilename(const std::string& Filename,
				 const ReflectionFileType& Type,
				 const PxdName& Pxd)
  {
    if (Filename == "") return 0;
    int nf = 0;
    seriescount++;    // count file series
    std::string fn = StringUtil::Trim(Filename);
    if (fn.size() > 0) {
      HKLIN_things.push_back(HKLIN_data(fn, Type, Pxd, seriescount));
      seriesnames.push_back(fn);
      nf++;
    }
    return nf;
  }
  //--------------------------------------------------------------
  int IO_files::AddHKLINfilename(const std::vector<std::string>& Filenames,
				 const std::vector<PxdName> Pxds,
				 const ReflectionFileType& Type,
				 const std::vector<std::string>& Seriesnames,
				 const std::vector<int>& fileSeries)
  {    
    if (Filenames.size() == 0) return 0;
    int nf = 0;
    for (size_t i=0;i<Filenames.size();i++) {
      // remove leading & trailing spaces
      std::string fn = StringUtil::Trim(Filenames[i]);
      if (fn.size() > 0) {
	HKLIN_things.push_back(HKLIN_data(fn, Type, Pxds[i], fileSeries[i]));
	seriescount = Max(seriescount, fileSeries[i]);
	nf++;
	//^
	//	  std::cout << "Store HKLIN filename " << fn << " type " << Type
	//		    << " series " << fileSeries[i] << "\n";  //^
      }
    }
    // Series names
    for (size_t i=0;i<Seriesnames.size();++i) {
      seriesnames.push_back(Seriesnames[i]);
    }
    return nf;
  }
  //--------------------------------------------------------------
  int IO_files::AddHKLINfilename(const std::vector<std::string>& Filenames,
				 const std::vector<PxdName> Pxds,
				 const ReflectionFileType& Type)
  {    
    if (Filenames.size() == 0) return 0;
    seriescount++;
    std::vector<int> fileSeries(Filenames.size());
    for (size_t i=0;i<Filenames.size();++i) {
      // if not MTZ, each files in its own series
      if (!Type.IsTypeMTZ() && i>0) seriescount++;
      fileSeries[i] = seriescount;
    }
    std::vector<std::string> seriesnames(1,Filenames[0]);
    return AddHKLINfilename(Filenames, Pxds,
			    Type, seriesnames, fileSeries);
  }
  //--------------------------------------------------------------
  int IO_files::AddHKLINfilename(const std::vector<std::string>& Filenames,
				 const ReflectionFileType& Type)
  {    
    if (Filenames.size() == 0) return 0;
    return AddHKLINfilename(Filenames,
	    std::vector<PxdName>(Filenames.size(), PxdName()),
	    Type);
  } 
  //--------------------------------------------------------------
  // Store an HKLOUT, HKLREF, XMLOUT or XYZIN filename
  // if non-blank. Strip leading & trailing spaces
  void IO_files::StoreFilename(const std::string& Streamname,
			       const std::string& Filename)
  {
    std::string fn = StringUtil::Trim(Filename);
    if (fn.size() > 0)
      {    
	if (Streamname == "HKLREF")
	  {HKLREF_filename = Filename;}
	else if (Streamname == "HKLOUT")
	  {HKLOUT_filename = Filename;}
	else if (Streamname == "XMLOUT")
	  {XMLOUT_filename = Filename;}
	else if (Streamname == "XYZIN")
	  {XYZIN_filename = Filename;}
	else
	  {
	    Message::message(Message_fatal
			     ("IO_files::StoreFilename: unrecognised stream "+Streamname));
	  }
      }
  }
  //--------------------------------------------------------------
  void IO_files::StoreHKLINpxdname(const int& index, const PxdName& Pxd)
  // Store PXD name for file
  {
    CheckIndex(index);
    HKLIN_things[index].pxdname = Pxd;
  }
  //--------------------------------------------------------------
  std::string IO_files::Filename(const std::string& Streamname) const
  {
    if (Streamname == "HKLREF")
      {return HKLREF_filename;}
    else if (Streamname == "HKLOUT")
      {return HKLOUT_filename;}
    else if (Streamname == "XMLOUT")
      {return XMLOUT_filename;}
    else if (Streamname == "XYZIN")
      {return XYZIN_filename;}
    else
      {
	Message::message(Message_fatal
		 ("IO_files::Filename: unrecognised stream "+Streamname));
      }
    return "";  // dummy
  }
  //--------------------------------------------------------------
  bool IO_files::IsFile(const std::string& Streamname) const
  {
    if (Streamname == "HKLREF")
      {return HKLREF_filename != "";}
    else if (Streamname == "HKLOUT")
      {return HKLOUT_filename != "";}
    else if (Streamname == "XMLOUT")
      {return XMLOUT_filename != "";}
    else if (Streamname == "XYZIN")
      {return XYZIN_filename != "";}
    else
      {
	Message::message(Message_fatal
		 ("IO_files::Filename: unrecognised stream "+Streamname));
      }
    return false;  // dummy
  }
  //--------------------------------------------------------------
  int IO_files::NumberOfFilesInSeries(const int& kseries) const
  // Number of files in series
  //  = 0 if series not found 
  {
    int n = 0;
    for (size_t i=0;i<HKLIN_things.size();i++) {
      if (HKLIN_things[i].series == kseries) {n++;}
    }
    return n;
  }
  //--------------------------------------------------------------
  std::string IO_files::HKLINfilename(const int& index,
				      ReflectionFileType& Type,
				      PxdName& Pxd, int& Series) const
// Return index'th HKLIN filename and its type, PxdName & series
  {
    CheckIndex(index);
    Type = HKLIN_things[index].type;
    Pxd = HKLIN_things[index].pxdname;
    Series = HKLIN_things[index].series;
    return HKLIN_things[index].filename;
  }
  //--------------------------------------------------------------
  std::string IO_files::HKLINfilename(const int& index) const
// Return index'th HKLIN filename only
  {
    CheckIndex(index);
    return HKLIN_things[index].filename;
  }
  //--------------------------------------------------------------
  std::string IO_files::HKLINfilename(const int& index,
				      ReflectionFileType& Type) const
// Return index'th HKLIN filename and its type only
  {
    CheckIndex(index);
    Type = HKLIN_things[index].type;
    return HKLIN_things[index].filename;
  }
  //--------------------------------------------------------------
  void IO_files::RemoveHKLINfile(const int& index)
  // Remove index'th file from list
  {
    CheckIndex(index);
    HKLIN_things.erase(HKLIN_things.begin()+index);
  }
  //--------------------------------------------------------------
  void IO_files::StoreHKLINreindex(const int& fileSeries,
				   const ReindexOp& reindex,
				   const correl_coeff& CC,
				   const float& likelihood,
				   const float& confidence)
  {
    CheckSeries(fileSeries);
    for (size_t index=0;index<HKLIN_things.size();index++) {
      if (HKLIN_things.at(index).series == fileSeries) {
	HKLIN_things.at(index).reindex = reindex;
	HKLIN_things.at(index).cc = CC;
	HKLIN_things.at(index).likelihood = likelihood;
	HKLIN_things.at(index).confidence = confidence;
	HKLIN_things.at(index).reindex_set = true;
      }
    }
  }
  //--------------------------------------------------------------
  ReindexOp IO_files::HKLINreindex(const int& index) const
  {
    CheckIndex(index);
    if (!HKLIN_things.at(index).reindex_set) {
      Message::message(Message_fatal
		       ("IO_files::HKLINreindex: reindex operator not set "+
			  clipper::String(index)));
    }
    return HKLIN_things.at(index).reindex;
  }
  //--------------------------------------------------------------
  correl_coeff IO_files::HKLIN_CC(const int& index) const
  {
    CheckIndex(index);
    if (!HKLIN_things.at(index).reindex_set) {
      Message::message(Message_fatal
		       ("IO_files::HKLIN_CC: CC not set "+
			  clipper::String(index)));
    }
    return HKLIN_things.at(index).cc;
  }
  //--------------------------------------------------------------
  HKLIN_data IO_files::HKLIN_Data(const int& index) const
  {
    CheckIndex(index);
    return HKLIN_things.at(index);
  }
  //--------------------------------------------------------------
  bool IO_files::IsHKLINreindex() const
  // True if any reindex operators have been stored
  {
    bool set = false;
    for (size_t i=0;i<HKLIN_things.size();i++) {
      if (HKLIN_things.at(i).reindex_set) set = true;
    }
    return set;
  }
  //--------------------------------------------------------------
  bool IO_files::IsHKLINreindex(const int& index) const
  // True if this reindex operators has been stored
  {
    CheckIndex(index);
    return HKLIN_things.at(index).reindex_set;
  }
  //--------------------------------------------------------------
  void IO_files::CheckIndex(const int& index) const
  {
    if (index < 0 || index >= int(HKLIN_things.size())) {
      Message::message(Message_fatal
		       ("IO_files::HKLIN: file index out of range "+
			clipper::String(index)+", "+
			clipper::String(int(HKLIN_things.size())-1)));
    }
  }
  //--------------------------------------------------------------
  void IO_files::CheckSeries(const int& fileSeries) const
  {
    if (fileSeries <= 0 || fileSeries > seriescount) {
      Message::message(Message_fatal
		       ("IO_files::HKLIN: file series number out of range "+
			clipper::String(fileSeries)+", "+
			clipper::String(seriescount)));
    }
  }
  //--------------------------------------------------------------
  void IO_files::StoreHKLINsymmetry
    (const int& index, const hkl_symmetry& symm_in)
    // Store symmetry from input file
  {
    CheckIndex(index);
    HKLIN_things.at(index).symmetry = symm_in;    
  }
  //--------------------------------------------------------------
  hkl_symmetry IO_files::HKLINsymmetry(const int& index) const
    // Retrieve symmetry
  {
    CheckIndex(index);
    return HKLIN_things.at(index).symmetry;
  }    
  //--------------------------------------------------------------
  void IO_files::StoreHKLINcell
    (const int& index, const Scell& cell)
    // Store cell from input file
  {
    CheckIndex(index);
    HKLIN_things.at(index).cell = cell;    
  }
  //--------------------------------------------------------------
  Scell IO_files::HKLINcell(const int& index) const
    // Retrieve cell
  {
    CheckIndex(index);
    return HKLIN_things.at(index).cell;
  }    
  //--------------------------------------------------------------
  void IO_files::SetMergedFlag
    (const int& index, const bool& Unmerged)
  // Set [un]merged, true if file is unmerged (default)
  {
    CheckIndex(index);
    if (Unmerged) {
      HKLIN_things.at(index).unmerged = +1;
    } else {
      HKLIN_things.at(index).unmerged = -1;
    }
  }
  //--------------------------------------------------------------
  bool IO_files::MergedFlag(const int& index) const
  // Fetch unmerged flag, true if file is unmerged
  {
    CheckIndex(index);
    return (HKLIN_things.at(index).unmerged > 0);
  }
  //--------------------------------------------------------------
  bool IO_files::AreThereFileSeries() const
  // Test if any file series are defined
  // (ie with wild-cards as opposed to single files)
  {
    bool isSeries = false;
    for (int i=0;i<seriescount;i++) {
      // Note series are numbered from 1
      if (NumberOfFilesInSeries(i+1) > 1) isSeries = true;
    }
    return isSeries;
  }
  //--------------------------------------------------------------
  int IO_files::HKLIN_fileSeries(const int& index) const
  // return file series for file number
  {
    CheckIndex(index);
    return HKLIN_things.at(index).series;
  }
  //--------------------------------------------------------------
  void IO_files::CheckHklinFileType()
  // check file type for all hklin files, reassign them if necessary
  {
    for (size_t ifl=0;ifl<HKLIN_things.size();++ifl) { // loop HKLIN files
      ReflectionFileType rftype(HKLIN_things[ifl].filename);
      if (rftype.FileType() != HKLIN_things[ifl].type.FileType()) {
	// Need to change it
	//^	std::cout << "HKLIN wrong type " << ifl
	//^		  << " type " << HKLIN_things[ifl].type.format()
	//^		  << " real type " << rftype.format() << "\n"; //^
	HKLIN_things[ifl].type = rftype;
	if (rftype.Unmerged() != 0) { // reset unmerged if known
	  HKLIN_things[ifl].unmerged = rftype.Unmerged();
	}
      }
    }
  }
  //--------------------------------------------------------------
} // namespace scala

