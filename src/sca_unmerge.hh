// sca_unmerge.hh
//
// Read scalepack data into hkl_unmerge object
// Also reads scalepack merged files & ShelX intensity files
// For ShelX files, there is no simple way to tell if they are merged or unmerged
//

#ifndef SCA_UNMERGE_HEADER
#define SCA_UNMERGE_HEADER

#include <fstream>

// Clipper
#include "rotation.hh"
#include <clipper/clipper.h>
using clipper::Vec3;
using clipper::Mat33;

#include <string>
#include "hkl_controls.hh"
#include "controls.hh"
#include "hkl_unmerge.hh"
#include "hkl_datatypes.hh"
#include "hash.hh"


namespace SCAIO
{
  enum ScaFileType {SCA_NONE, SCA_MERGED, SCA_UNMERGED, SCA_SHELX,
		    SCA_SAINT}; // SCA_NONE means unknown

  class KeyValues;

  class SCAunmergeFile
  //! Class to read files in Scalepack, Shelx or Saint formats into unmerged object
  {
  public:
    SCAunmergeFile(){}  // Dummy constructor

    //! Open file & check file type
    /*!  \param scaname            name of file (or logical name) */
    SCAunmergeFile(const std::string& scaname,
		   std::string& output);

    ~SCAunmergeFile() {if (scain.is_open()) scain.close();} // close if open

    //! Return file type
    ScaFileType FileType() const {return sca_type;}

    //! Fill an hkl_unmerge_list object from this file with default options
    scala::FileRead FillHklList(std::string& output,
				scala::hkl_unmerge_list& hkl_list);
    //! Fill an hkl_unmerge_list object from this file with cell & default options
    scala::FileRead FillHklList(std::string& output,
				const scala::Scell& input_cell,
				scala::hkl_unmerge_list& hkl_list);

    //! Fill an hkl_unmerge_list object from this file with selection options
    // Recognised formats:
    //   SCALEPACK unmerged
    //   SCALEPACK merged file
    //   ShelX HKL file (h k l I sigI) format 3I4,2F8.2
    //   Saint format
    scala::FileRead AddHklList(const scala::file_select& file_sel, 
			const scala::all_controls& controls,
			const scala::PxdName& InputPxdName,
			const scala::Scell& input_cell,
			const double& input_wavelength,
			std::string& output,
			scala::hkl_unmerge_list& hkl_list);
    //!< On entry:
    //!<  file_sel           flags for general selection
    //!<                     - dataset selection
    //!<                     - batch selection
    //!<                     - resolution limits
    //!<                     - detector coordinate rejection ranges
    //!<             only resolution limits used at present
    //!<  controls           run controls
    //!<  InputPxdName       PXD name (none in SCA file)
    //!<  input_cell         cell from input
    //!<  input_wavelength   wavelength from input = -1.0 is not set
    //!<  output             output object for printing
    //!< 
    //!< On exit:
    //!<  hkl_list  has been filled, but not organised:
    //!<  this needs a call to "prepare" or "change_symmetry"


  private:
    // Private member functions
    int ReadObservations(std::ifstream& scain,
			 const scala::file_select& file_sel,
			 const scala::all_controls& controls,
			 std::string& output,
			 scala::hkl_unmerge_list& hkl_list);
    int MakeDataset(const scala::PxdName& InputPxdName);
    int MakeBatches(const scala::hkl_symmetry& Symm,
		    const int& Idataset);

    int ReadSaint(std::ifstream& scain,
		  const scala::file_select& file_sel,
		  const scala::all_controls& controls,
		  std::string& output,
		  scala::hkl_unmerge_list& hkl_list);
    void ZeroRecVector(std::vector<clipper::Coord_orth>& vc);
    scala::data_flags SetSaintDataFlags();
    bool CheckCells(const scala::Scell& input_cell, const scala::Scell& file_cell,
		    const std::string& label, std::string& output);


    // Private data
    std::vector<scala::Xdataset> datasets;
    std::vector<scala::Batch> batches;

    // File type
    ScaFileType sca_type; // SCA_NONE, SCA_MERGED, SCA_UNMERGED, SCA_SHELX

    std::string scaname;    
    std::ifstream scain;
    clipper::String line;  // current line
    bool lineFull;        // true if line has something in it (eg from header)
    std::string spaceGroupName;
    scala::Scell cell;
    double wavelength;

    scala::Scell inputcell;  // from control input, if given
    double inputwavelength;

    // Batch information
    int Nbatches;
    hash_table lookup;
  };
  //  SCAunmergeFile
//--------------------------------------------------------------
class SaintRun {
  //  A "run" or "sweep" from a Saint file, ie information about a set of reflections
  //  with the same scan axis and non-scan goniostat axes
  //
  // The SAINT file contains information which can be used to extract the unit cell,
  // & crystal orientation information. The only external information needed in the
  // wavelength which cannot be inferred from the file.
  //
  // Each reflection record contains the following information on cell & orientation:
  //   h,k,l        indices
  //   s0r(3)       direction cosines of incident beam s0, projected on reciprocal
  //                lattice axes. This in the reciprocal axis frame "r", not
  //                necessarily orthogonal
  //   s2r(3)       direction cosines of secondary beam in reciprocal axis frame
  //   rot          rotation angle (omega or phi) degrees
  //   iaxis        scan axis 2=omega, 3=phi
  //   istl         sin theta/lambda * 10000
  //   chi          chi angle degrees
  //   other        phi or omega, non-scan angle degrees
  //
  //  s = [R] [U] [F^-1] [F] [K] [C] h   = [R] [U] [B] h 
  //  l      z   x      r         frames
  //
  //  orthogonalisation matrix [B] = [K][C] where [C] = diag(a*,b*,c*)
  //  [F] = [K]T 
  //  ie [F] takes a unit orthogonal vector v -> f (cosines on a*, b*, c*)
  //     f = (v.a*/|a*|  v.b*/|b*|  v.c*/|c*|)
  //       =  [F] v
  //
  //     [F] = (a*/|a*|) = [K]T = [B][C^-1]
  //           (b*/|b*|)
  //           (c*/|c*|)
  //  
  //  Coordinate frames:
  //    r   reciprocal axis frame
  //    x   crystal frame
  //    z   zero angle frame (all goniostat angles zero)
  //    l   laboratory frame
  //
  // Determination of parameters
  //  1) the cell is determined from sin theta/ lambda values
  //     This may be overridden by input values
  //  2) the orientation matrix [U] is determined from the incident beam cosines
  //  3) the wavelength is determined from {U][B]: this may be overridden from input
  //  
public:
  SaintRun(){}

  SaintRun(const int& Ibatch, const int& Iaxis, const double& Chi, const double& Other,
	    const int& Icryst, const double& Swing);

  void init(const int& Ibatch, const int& Iaxis, const double& Chi, const double& Other,
	    const int& Icryst, const double& Swing);

  bool SameRun(const int& Ibatch, const int& Iaxis,
	       const double& chi, const double& other,
	       const int& icryst, const double& swing);

  void StoreReflection(const IVect3& hkl, const int& batch,
		       const DVect3& s0r, const DVect3& s2r,
		       const double& rot, const int& istl);

  // calculate cell (returned) from dstar information, linear LSQ
  scala::Scell CalcCell();

  void SetBatchOffset(const int& batchOffset) {batchoffset = batchOffset;}
  int BatchOffset() const {return batchoffset;}
  int MaxBatch() const;  // maximum batch number including offset
  int Nbatches() const { return nbatches;}

  std::string ScanAxis() const;
  double Chi() const {return chi;}
  double Other() const {return other;}
  double Swing() const {return swing;}
  // return angle (degrees) between rotation axis (omega or phi) and beam
  double RotAxistoBeamAngle() const;

  //  calculate [U] from incident beam cosines
  //  returns Utrn  translational part (should = 000)
  DMat33 Umatrix(DVect3& Utrn);

  double Wavelength();
  double SdWavelength() const {return sdwavelength;}
  scala::Scell Cell() const {return cell;}

  void SetWavelength(const double& wvl)
  {wavelength = wvl; sdwavelength=0.0;lambdaDone=true;}
  void SetCell(const scala::Scell& Cell) {cell = Cell;cellDone=true;}

  static DMat33 Rmatrix(const double& omegar, const double& chir, const double& phir);
  static DVect3 s0v() {return s0;}

  std::vector<CMtz::MTZBAT> MakeBatch(const float& wavelength);

private:
  double BasicAngleD(const double& aa);
  bool SameAngle(const double& a1, const double& a2);

  struct CompareMtzBat {
    // This construct seems to be a way of getting pointer-to-function
    // into the argument for sort. Copied from the Web.
    bool operator()(const CMtz::MTZBAT& a, const CMtz::MTZBAT& b);
  };


  int ibatch;  // batch number from file
  int iaxis;   //  scan axis 2=omega, 3=phi
  double chi;  // chi angle, degrees
  double other; // non-scan axis (phi or omega), degrees
  int icryst;
  double swing;
  std::vector<scala::Range> rotRange;

  double omegar, chir, phir; // radians
  //  Lists of all reflections from file
  std::vector<DVect3> s0rv;   // unit s0 in reciprocal axis (r) frame from file
  std::vector<DVect3> s2rv;   // unit s2 in reciprocal axis frame from file
  std::vector<IVect3> hklv;    // hkl
  std::vector<int> istlv;  // d* = 2sin theta/lambda
  std::vector<DMat33> Rv;     // Goniostat rotation matrix [Omega][Chi][Phi]


  // results
  bool cellDone;     // true if cell has been calculated
  bool UmatrixDone;  // true if [U] has been calculated
  bool lambdaDone;   // true if wavelength has been calculated
  DMat33 U;          // orientation matrix [U]
  DVect3 Ut;         //  translational part, should be 0
  scala::Scell cell;
  double wavelength;
  double sdwavelength;


  // Batch information
  int nbatches;
  hash_table lookup;
  int table_size;
  int batchoffset;

  //  Saint goniostat things: axis directions
  // Note: "Cambridge" axis system
  //   z along principle rotation axis (e1 = Omege in ths case)
  //   x in plane containing z & s0
  //   y along s0 x e1
  //   source vector s0 anti-parallel to beam ie approximately along -x
  static DVect3 e1;  // Omega along z
  static DVect3 e2;  // Chi along x
  static DVect3 e3;  // Phi along z
  static DVect3 s0;  // s0 along x

  static int nxpix, nypix;
};  // SaintRun
}
#endif
