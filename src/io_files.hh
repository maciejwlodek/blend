// io_files.hh
//
// Object to contain all io filenames etc
// principally HKLIN or XDSIN (possibly multiple), HKLOUT, HKLREF XMLOUT
//

#ifndef IO_FILES_HEADER
#define IO_FILES_HEADER

#include <vector>
#include <string>
#include "hkl_datatypes.hh"
#include "score_datatypes.hh"
#include "hkl_symmetry.hh"
#include "filetype.hh"

namespace scala
{
  ///  enum ReflectionFileType {MTZ, XDS, SCA};
  enum RefListType {NONE, MERGED, UNMERGED}; // NONE means unknown

  class HKLIN_data {
  public:
    HKLIN_data(const std::string& filename_in,
	       const ReflectionFileType& type_in,  // MTZ, XDS, SCA
	       const PxdName& pxdname_in,
	       const int& series_in);

    std::string filename;
    ReflectionFileType type;  // MTZ, XDS..., SCA... etc
    hkl_symmetry symmetry;
    PxdName pxdname;
    bool reindex_set;    // true is reindex & CC 
    ReindexOp reindex;
    Scell cell;
    correl_coeff cc;
    float likelihood;
    float confidence;
    int series;          // internal series count
			 // a "series" is either a file or a set of files
			 // specified with wild-cards
    int unmerged;    //  +1  unmerged, -1  merged,  0 unknown
  };

  class IO_files
  {
  public:
    IO_files();

    // Store one or more HKLIN filenames (this includes XDSIN) if non-blank.
    //   Strip leading & trailing spaces
    //   Type is MTZ or XDS
    // Return number of filenames added
    int AddHKLINfilename(const std::string& Filename,
			 const ReflectionFileType& Type);

    int AddHKLINfilename(const std::string& Filename,
			 const ReflectionFileType& Type,
			 const PxdName& Pxd);

    int AddHKLINfilename(const std::vector<std::string>& Filenames,
			 const ReflectionFileType& Type);

    int AddHKLINfilename(const std::vector<std::string>& Filenames,
			 const std::vector<PxdName> Pxds,
			 const ReflectionFileType& Type);

    int AddHKLINfilename(const std::vector<std::string>& Filenames,
			 const std::vector<PxdName> Pxds,
			 const ReflectionFileType& Type,
			 const std::vector<std::string>& Seriesnames,
			 const std::vector<int>& fileSeries);

    // Store an HKLOUT, HKLREF, XYZIN or XMLOUT filename
    //   if non-blank. Strip leading & trailing spaces
    void StoreFilename(const std::string& Streamname, const std::string& Filename);

    // Store PXD name for file
    void StoreHKLINpxdname(const int& index, const PxdName& Pxd);

    // Retrieve or test filenames
    std::string Filename(const std::string& Streamname) const;
    bool IsFile(const std::string& Streamname) const; // returns true if non-blank

    int NumberOfHKLINnames() const {return HKLIN_things.size();}

    // Number of HKLIN file series
    int NumberOfFileSeries() const {return seriescount;}
    // Number of HKLIN files in series
    int NumberOfFilesInSeries(const int& kseries) const;
    // Test if any HKLIN file series are defined
    // (ie with wild-cards as opposed to single files)
    bool AreThereFileSeries() const;
    // Series names
    std::vector<std::string> SeriesNames() const {return seriesnames;}


    //  true if index is valid filename index
    bool IsHKLINfile(const int& index) const
    {return (index >= 0 && index < int(HKLIN_things.size()));}
  
    // Return index'th HKLIN filename and its type, PxdName & series if valid
    std::string HKLINfilename(const int& index,
			      ReflectionFileType& Type,
			      PxdName& Pxd, int& Series) const;
    // Return index'th HKLIN filename only
    std::string HKLINfilename(const int& index) const;
    // Return index'th HKLIN filename and its type
    std::string HKLINfilename(const int& index,
			      ReflectionFileType& Type) const;
    // Remove index'th file from list
    void RemoveHKLINfile(const int& index);

    // Store results for all files in file series
    void StoreHKLINreindex(const int& fileSeries,
			   const ReindexOp& reindex,
			   const correl_coeff& CC,
			   const float& likelihood,
			   const float& confidence);
    // True if any reindex operators have been stored
    bool IsHKLINreindex() const;
    bool IsHKLINreindex(const int& index) const;

    ReindexOp HKLINreindex(const int& index) const;
    correl_coeff HKLIN_CC(const int& index) const;
    // return data for file index
    HKLIN_data HKLIN_Data(const int& index) const;
    // return file series for file number
    int HKLIN_fileSeries(const int& index) const;

    // Store symmetry from input file
    void StoreHKLINsymmetry(const int& index, const hkl_symmetry& symm_in);
    // Retrieve symmetry
    hkl_symmetry HKLINsymmetry(const int& index) const;

    // Store cell
    void StoreHKLINcell(const int& index, const Scell& cell);
    // Retrieve cell
    Scell HKLINcell(const int& index) const;

    // Set [un]merged, true if file is unmerged (default if not set here)
    void SetMergedFlag(const int& index, const bool& Unmerged=false);
    // Fetch unmerged flag, true if file is unmerged
    bool MergedFlag(const int& index) const;
    
    // check file type for all hklin files, reassign them if necessary
    void CheckHklinFileType();

  private:
    void CheckIndex(const int& index) const;
    void CheckSeries(const int& fileSeries) const;

    std::vector<HKLIN_data> HKLIN_things;
    std::string HKLREF_filename;
    std::string HKLOUT_filename;
    std::string XMLOUT_filename;
    std::string XYZIN_filename;
    int seriescount;   // count from 1 for first file or series
    std::vector<std::string> seriesnames;
  };

}
#endif
