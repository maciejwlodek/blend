// xds_unmerge.hh
//
// Read XDS data into hkl_unmerge object
//

#ifndef XDS_UNMERGE_HEADER
#define XDS_UNMERGE_HEADER

#include <string>
#include "hkl_controls.hh"
#include "controls.hh"
#include "hkl_unmerge.hh"
#include "hkl_datatypes.hh"
#include "range.hh"

// Clipper
#include <clipper/clipper.h>
using clipper::Vec3;
using clipper::Mat33;


namespace XDSIO
{
  //--------------------------------------------------------------
  class KeyValues
  {
    // A key & a list of values
  public:
    KeyValues() {}
    KeyValues(const clipper::String& Key) 
    {
      key = Key;
      values.clear();
    }
    
    clipper::String key;
    std::vector<clipper::String> values;
  };
  //--------------------------------------------------------------
  enum XDStype {NONE, XDS_INTEGRATE, XDS_ASCII, XSCALE};
  //--------------------------------------------------------------
  class XDSunmergeFile
  {
  public:
    XDSunmergeFile(): ftype(NONE) {}  // Dummy constructor
    //! constructor from file name
    XDSunmergeFile(const std::string& xdsName,
		   const int& verbose,
		   std::string& output);
    //! init from file name
    void init(const std::string& xdsName,
	      const int& verbose,
	      std::string& output);

    //! Fill an hkl_unmerge_list object from this file with default options
    scala::FileRead FillHklList (std::string& output,
		   const int& verbose,
		   scala::hkl_unmerge_list& hkl_list);

    //! Fill an hkl_unmerge_list object from this file
    scala::FileRead AddHklList (const scala::file_select& file_sel, 
		   const scala::all_controls& controls,
		   const scala::PxdName& InputPxdName,
		   const scala::Scell& newcell,
		   const int& fileSeries,
 		   std::string& output,
		   const int& verbose,
		   scala::hkl_unmerge_list& hkl_list);
    //!< On entry:
    //!<  xdsname            name of XDS file (or logical name)
    //!<  file_sel           flags for general selection
    //!<                     - dataset selection
    //!<                     - batch selection
    //!<                     - resolution limits
    //!<                     - detector coordinate rejection ranges
    //!<             only resolution limits used at present
    //!<  controls           run controls
    //!<  InputPxdName       PXD name (none in XDS file)
    //!<  newcell            if non-null, override cell
    //!<  output             output string for printing
    //!<  verbose            set verbosity level
    //!<                      = 0 silent, = +1 usual summary
    //!<                      >= +2 debug
    //!< 
    //!< On exit:
    //!<  hkl_list  has been filled, but not organised:
    //!<  this needs a call to "prepare" or "change_symmetry"

    int Ncols() const {return nitem;}
    //! return SpaceGroup object if known, or blank
    scala::SpaceGroup Spacegroup() const {return spacegroup;}
    //! return cell
    scala::Scell Cell() const {return cell;}

  private:
    // Private member functions
    void ClearColumnNumbers();
    int ProcessKeyValues(const KeyValues& keyval, std::string& output);
    void XDStoCambridgeFrame(std::string& output,
			     const int& verbose);
    int ReadObservations(FILE* xdsin,
			 const scala::file_select& file_sel,
			 const scala::all_controls& controls,
			 const int& fileSeries,
			 std::string& output,
			 scala::hkl_unmerge_list& hkl_list);
    scala::data_flags CheckColumns(std::string& output) const;
    void ReconstructOrientation(scala::hkl_unmerge_list& hkl_list,
				std::string& output);
    int MakeDataset(const scala::PxdName& InputPxdName);
    int MakeBatches(const scala::hkl_symmetry& Symm,
		    const int& Idataset);
    int SetItems();
    void UpdatePolarisation(const scala::PolarisationControl& polarisationcontrol,
			    scala::hkl_unmerge_list& hkl_list,
			    std::string& output);
    

    // Private data
    std::vector<scala::Xdataset> datasets;
    std::vector<scala::Batch> batches;

    std::string xdsname;
    FILE* xdsin;
    bool atdata; // true if file is positioned after header
    std::string xHM; // space group name
    bool GeomInfo; // true if all geometrical information present in file
    scala::SpaceGroup spacegroup;

    // Things from XDS header
    // FORMAT or OUTPUT_FILE
    XDStype ftype;
    // SPACE_GROUP_NUMBER
    int spacegroup_number;
    // UNIT_CELL_CONSTANTS
    scala::Scell cell;
    // DATE
    std::string xdsdate;
    // NAME_TEMPLATE_OF_DATA_FRAMES
    std::string image_template;
    // DATA_RANGE
    int FirstImage, LastImage;
    double wavelength;
    // INCIDENT_BEAM_DIRECTION
    clipper::Vec3<double> beam;
    // ROTATION_AXIS
    clipper::Vec3<double> rotaxis;
    // OSCILLATION_RANGE
    double osc_range;
    // STARTING_ANGLE
    double phi0;
    // STARTING_FRAME  (may be different from FirstImage)
    int image0;
    // DIRECTION_OF_DETECTOR_X-AXIS
    clipper::Vec3<double> dxaxis;
    // DIRECTION_OF_DETECTOR_X-AXIS
    clipper::Vec3<double> dyaxis;
    // DETECTOR_DISTANCE
    double det_dist;
    // ORGX, ORGY
    double orgx, orgy;
    // NX, NY   number of pixels along detector x & y
    int nxpix, nypix;
    // QX, QY    pixel size along detector x & y
    double qxpixsize, qypixsize;
    // REFLECTING_RANGE_E.S.D.  (== mosaicity)
    double mosaicrange;
    // NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD
    int nitem_rec;
    int nitem;

    int col_H, col_K, col_L, col_IOBS, col_SIGMAI, col_XD, col_YD, col_ZD;
    int col_RLP, col_PEAK, col_CORR, col_PSI, col_ISET, col_DCY;

    // Conversion matrix from XDS frame to Cambridge frame
    Mat33<double> Qxds2cam;
    // Conversion from pixel position in mm in XDS frame p(xds) to
    // position in Cambridge frame = [Q][Ed]
    //
    //  p(xds)  =  ( qx*(xd-origx) )
    //             ( qy*(yd-origy) )
    //             ( F )
    //    where F is crystal to detector distance
    //
    Mat33<double> QEd;

    // Orientation matrix [U]
    Mat33<double> Umat;
    Vec3<double> s0C;     // s0, Cambridge frame, unit length

    // Batch information
    int Nbatches;
    std::vector<int> NobsInBatch;
  };
}
#endif
