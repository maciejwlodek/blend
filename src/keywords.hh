#ifndef KEYWORDS_HEADER
#define KEYWORDS_HEADER

#include "CCP4base.hh"
#include "hkl_datatypes.hh"
#include "range.hh"
#include "globalcontrols.hh"
#include "listdir.hh"
#include "expandfilename.hh"

using scala::ResoRange;
using scala::Chirality;

namespace phaser_io {

  //--------------------------------------------------------------
class CHIRAL : public InputBase, virtual public CCP4base
{
public:
  Chirality Chiral;
  
public:
  CHIRAL();
  virtual ~CHIRAL() {}
  Token_value parse(std::istringstream&);
  void setCHIRAL(Chirality);
  Chirality getCHIRAL() const;
  void analyse();
};

  //--------------------------------------------------------------
  class RESO : public InputBase, virtual public CCP4base
  {
    // Syntax:
    //  RESOlution  <high> [<low>]  either order
    //  RESOlution  HIGH <high>
    //  RESOlution  LOW  <low>
  public:
    RESO();
    virtual ~RESO() {}
    Token_value parse(std::istringstream&);

    void setRESO(const ResoRange& reso_range);
    ResoRange getRESO() const;
    void analyse(){}

    bool ValidRESO() const {return Valid;}

  private:
    ResoRange input_resolution_range;
    bool Valid;
  };

  //--------------------------------------------------------------
  class ISIG : public InputBase, virtual public CCP4base
  {

    // Syntax:
    //  ISIGLIMIT <minimum<I>/<sigmaI>>
  public:
    ISIG();
    virtual ~ISIG() {}
    Token_value parse(std::istringstream&);

    void setISIG(const float& MinIsig);
    float getISIG() const;
    void analyse(){}

  private:
    float MinIsigRatio;
  };

  //--------------------------------------------------------------
  class LABI : public InputBase, virtual public CCP4base
  {
    // Syntax: LABIN [F|I = ] <F|Ilabel> [[SIGF|I = ] <sigF|Ilabel>]
  public:
    LABI();
    virtual ~LABI() {}
    Token_value parse(std::istringstream&);

    void setLABI_I(const std::string& lab_I) {FIlabel = lab_I;}
    void setLABI_sigI(const std::string& lab_sigI) {sigFIlabel = lab_sigI;}
    std::string getLABI_I() const {return FIlabel;}
    std::string getLABI_sigI() const {return sigFIlabel;}
    void analyse(){}

  private:
    std::string FIlabel;
    std::string sigFIlabel;
  };

  //--------------------------------------------------------------
  class LABR : public InputBase, virtual public CCP4base
  {
    // Syntax: LABREF [F|I = ] <F|Ilabel> [[SIGF|I = ] <sigF|Ilabel>]
  public:
    LABR();
    virtual ~LABR() {}
    Token_value parse(std::istringstream&);

    void setLABR_I(const std::string& lab_I) {FIlabel = lab_I;}
    void setLABR_sigI(const std::string& lab_sigI) {sigFIlabel = lab_sigI;}
    std::string getLABR_I() const {return FIlabel;}
    std::string getLABR_sigI() const {return sigFIlabel;}
    void analyse(){}

  private:
    std::string FIlabel;
    std::string sigFIlabel;

  };
  //--------------------------------------------------------------
  class ORIG : public InputBase, virtual public CCP4base
  {
    // Syntax: ORIGinalLattice
  public:
    ORIG();
    virtual ~ORIG() {}
    Token_value parse(std::istringstream&);


    void setORIGINAL(const bool& OrigLat) {OriginalLat = OrigLat;}
    bool getOriginalLatFlag() const {return OriginalLat;}
    void analyse(){}

  private:
    bool OriginalLat;
  };
  //--------------------------------------------------------------
  class LAUE : public InputBase, virtual public CCP4base
  {
    // Syntax: LAUEgroup <GroupSymbol> | HKLIN | ALL
  public:
    LAUE();
    virtual ~LAUE() {}
    Token_value parse(std::istringstream&);


    void setLAUEGROUP(const std::string& LG) {Lauegroup = LG;}
    std::string getLAUEGROUP() const {return Lauegroup;}
    void analyse(){}

  private:
    std::string Lauegroup;
  };

  //--------------------------------------------------------------
  class SPAC : public InputBase, virtual public CCP4base
  {
    // Syntax: SPACegroup <GroupSymbol> | HKLIN
  public:
    SPAC();
    virtual ~SPAC() {}
    Token_value parse(std::istringstream&);


    void setSPACEGROUP(const std::string& SG) {SPGname = SG;}
    std::string getSPACEGROUP() const {return SPGname;}
    void analyse(){}

  private:
    std::string SPGname;
  };

  //--------------------------------------------------------------
  class REIN : public InputBase, virtual public CCP4base
  {
    // Syntax: REINdex [LEFTHANDED] <operator>
    //   eg reindex 2h+k,k,l
    //  LEFTHANDED allows inversion operator
  public:
    REIN();
    virtual ~REIN() {}
    Token_value parse(std::istringstream&);


    void setREINDEX(const std::string& Operator);
    scala::ReindexOp getREINDEX() const {return  Reindex;}
    bool REINDEXisSet() const {return set;}
    void analyse(){}

  private:
    scala::ReindexOp Reindex; 
    bool set;   // if it has been set
    bool left;  // true to allow lefthand operator (negative determinant)
  };

  //--------------------------------------------------------------
  class TOLE : public InputBase, virtual public CCP4base
  {
    // Syntax: TOLERANCE  <LatticeTolerance>
    //   tolerance in degrees for test lattice to fit input cell
  public:
    TOLE();
    virtual ~TOLE() {}
    Token_value parse(std::istringstream&);


    void setTOLERANCE(const double& Tol) {Tolerance = Tol;}
    double getTOLERANCE() const {return  Tolerance;}
    void analyse(){}

  private:
    double Tolerance;
  };
  //--------------------------------------------------------------
  class ACCE : public InputBase, virtual public CCP4base
  {
    // Syntax: ACCEPT <AcceptanceFraction> 
    //   
  public:
    ACCE();
    virtual ~ACCE() {}
    Token_value parse(std::istringstream&);


    void setACCEPTANCELIMIT(const double& Acc) {acceptlimit = Acc;}
    double getACCEPTANCELIMIT() const {return  acceptlimit;}
    void analyse(){}

  private:
    double acceptlimit;
  };
  //--------------------------------------------------------------
  class  CENT: public InputBase, virtual public CCP4base
  {
    // Syntax: CENTre
    //   
  public:
    CENT();
    virtual ~CENT() {}
    Token_value parse(std::istringstream&);
    clipper::Vec3<int> CentreGridMax() const {return gridmax;}
    bool Centre() const {return centre;}

    void analyse(){}

  private:
    bool centre;
    clipper::Vec3<int> gridmax;
  };
  //--------------------------------------------------------------
  class  SYSA: public InputBase, virtual public CCP4base
  {
    // Syntax: SYSAbsence <ProbabilityThreshold>
    //   
  public:
    SYSA();
    virtual ~SYSA() {}
    Token_value parse(std::istringstream&);


    void setProbabilityThreshold(const double& Pt) {ProbabilityThreshold = Pt;}
    double getProbabilityThreshold() const {return ProbabilityThreshold;}
    void analyse();

  private:
    double ProbabilityThreshold;
  };
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  class  PARTIALS: public InputBase, virtual public CCP4base
  {
    // Syntax: PARTials [TEST <lower_limit> <upper_limit>] [CORRECT <minimum_fraction>] [NOCHECK] [[NO]GAP <gap limit>]
    //   
  public:
    PARTIALS();
    virtual ~PARTIALS() {}
    Token_value parse(std::istringstream&);


    void setFracLimMin(const double& MinLim) 
    {minfraclim = MinLim;
}
    void setFracLimMax(const double& MaxLim) 
    {maxfraclim = MaxLim;}

    void setFracScaleLim(const double& MinLim)
    {minscalefrac = MinLim;}

    double getFracLimMin() const  {return minfraclim;}
    double getFracLimMax() const  {return maxfraclim;}

    double getSclMinLim() const {return minscalefrac;}
    bool getCheck() const {return check;}

    void setMaxGap(const int& MaxGap) {maxgap = MaxGap;}  // set maxgap
    int getMaxGap() const {return maxgap;} // return maxgap

    void analyse(){}

  private:
    double minfraclim, maxfraclim;
    double minscalefrac;
    bool check;
    int maxgap;
  };
  //--------------------------------------------------------------
  class  SYSABS: public InputBase, virtual public CCP4base
    // SYSTematicabsence  ON | OFF
  {
  public:
    SYSABS();
    virtual ~SYSABS() {}
    Token_value parse(std::istringstream&);

    void setSYSABSCHECK(const bool& flag) {check = flag;}

    bool getSYSABSCHECK() const {return check;}

    void analyse(){}

  private:
    bool check;
  };
  //--------------------------------------------------------------
  class HKLIN : public InputBase, virtual public CCP4base
  {
    // Syntax: HKLIN <filename>
    //  filename may contain wild-cards "*" or "?"
  public:
    HKLIN();
    virtual ~HKLIN() {}
    Token_value parse(std::istringstream&);

    // Stores file names from command line
    // returns number of files added
    int AddHKLINfilename(const std::vector<std::string>& Names);

    int HKLINnumber() const {return names.size();}
    std::vector<std::string> getHKLINnames() const;
    std::vector<scala::PxdName> getHKLINpxdNames() const {return pxdnames;}
    std::vector<int> getHKLINfileseries() const {return seriesnumber;}
    std::vector<std::string> getHKLINseriesNames() const {return seriesnames;}
    void analyse(){}
    // true if there were wild-cards in any filename
    bool HKLIN_wild() const {return wild;}

    // Strip out files in all series which are out of time order
    // ie where "later" files in order have earlier modification times
    // returns list of filenames rejected
    std::vector<std::string> RejectOutofSequenceHKLINFiles();

  private:
    bool wild;
    // these vector are the same length ie number of files -->
    std::vector<FileNameTime> names;   // actual file names & mod time
    std::vector<scala::PxdName> pxdnames;
    std::vector<int> seriesnumber;
    // <--
    std::vector<std::string> seriesnames;  // names including wild cards
    int fileseries;   // Number of file or wild-card series, from 1

    // Store names etc
    void storeNames(const scala::ExpandFileName& filenames);

    // reject & remove from list files with indices i1 to i2 (inclusive)
    // return list of filenames rejected
    std::vector<std::string> Reject(const int& i1, const int& i2);

    // return fname modified to include wildcard characters "?"
    // Match fname to template from end,
    //  replace characters equivalent to "#" in template by "?" in fname
    static std::string MakeWildTemplate(const std::string& fname,
					const std::string& ftemplate);
  };
  //--------------------------------------------------------------
  class XDSIN : public InputBase, virtual public CCP4base
  {
    // Syntax: XDSIN <filename>
  public:
    XDSIN();
    virtual ~XDSIN() {}
    Token_value parse(std::istringstream&);

    // Stores file names from command line
    // returns number of files added
    int AddXDSINfilename(const std::string& Name);

    //    void setXDSIN(const std::string& Name) {names.push_back(Name);}
    int XDSINnumber() const {return names.size();}
    std::vector<std::string> getXDSINnames() const {return names;}
    std::vector<scala::PxdName> getXDSINpxdNames() const {return pxdnames;}
    void analyse(){}

  private:
    std::vector<std::string> names;
    std::vector<scala::PxdName> pxdnames;
  };
  //--------------------------------------------------------------
  class SCAIN : public InputBase, virtual public CCP4base
  {
    // Syntax: SCAIN <filename>
  public:
    SCAIN();
    virtual ~SCAIN() {}
    Token_value parse(std::istringstream&);

    // Stores file names from command line, series 1
    // returns number of files added
    int AddSCAINfilename(const std::string& Name);

    //    void setSCAIN(const std::string& Name) {names.push_back(Name);}
    int SCAINnumber() const {return names.size();}
    std::vector<std::string> getSCAINnames() const {return names;}
    std::vector<scala::PxdName> getSCAINpxdNames() const {return pxdnames;}
    void analyse(){}

  private:
    std::vector<std::string> names;
    std::vector<scala::PxdName> pxdnames;
  };
  //--------------------------------------------------------------
  class MERGED : public InputBase, virtual public CCP4base
  {
    // Syntax: MERGED <spacegroupname>
    // Tell the program that the input file is merged, and in what space (point) group
    // Only needed for ShelX files
  public:
    MERGED();
    virtual ~MERGED() {}
    Token_value parse(std::istringstream&);

    // Returns blank string in not merged, "UNKNOWN" if unset
    std::string getHklinMergedGroup() const {return spacegroupName;}

    void analyse(){}

  private:
    std::string spacegroupName;
  };
  //--------------------------------------------------------------
  class HKLOUT : public InputBase, virtual public CCP4base
  {
    // Syntax: HKLOUT <filename>
  public:
    HKLOUT();
    virtual ~HKLOUT() {}
    Token_value parse(std::istringstream&);


    void setHKLOUT(const std::string& Name) {name = Name;}
    std::string getHKLOUT() const {return name;}
    void analyse(){}

  private:
    std::string name;
  };
  //--------------------------------------------------------------
  class HKLREF : public InputBase, virtual public CCP4base
  {
    // Syntax: HKLREF <filename>
  public:
    HKLREF();
    virtual ~HKLREF() {}
    Token_value parse(std::istringstream&);


    void setHKLREF(const std::string& Name) {name = Name;}
    std::string getHKLREF() const {return name;}
    void analyse(){}

  private:
    std::string name;
  };
  //--------------------------------------------------------------
  class XMLOUT : public InputBase, virtual public CCP4base
  {
    // Syntax: XMLOUT <filename>
  public:
    XMLOUT();
    virtual ~XMLOUT() {}
    Token_value parse(std::istringstream&);


    void setXMLOUT(const std::string& Name) {name = Name;}
    std::string getXMLOUT() const {return name;}
    void analyse(){}

  private:
    std::string name;
  };
  //--------------------------------------------------------------
  class XYZIN : public InputBase, virtual public CCP4base
  {
    // Syntax: XYZIN <filename>
  public:
    XYZIN();
    virtual ~XYZIN() {}
    Token_value parse(std::istringstream&);


    void setXYZIN(const std::string& Name) {name = Name;}
    std::string getXYZIN() const {return name;}
    void analyse(){}

  private:
    std::string name;
  };
  //--------------------------------------------------------------
  class NAME : public InputBase, virtual public CCP4base
  {
    // Syntax: NAME PROJECT <Pname> CRYSTAL <Xname> DATASET <Dname>
  public:
    NAME();
    virtual ~NAME() {}
    Token_value parse(std::istringstream&);


    static void setNAME(const scala::PxdName& PXDname);
    static scala::PxdName getNAME() {return pxdname;}
    void analyse(){}

  private:
    static scala::PxdName pxdname;
  };
  //--------------------------------------------------------------
  class COPY : public InputBase, virtual public CCP4base
  {
    // Syntax: COPY 
    //  Just copy file
  public:
    COPY();
    virtual ~COPY() {}
    Token_value parse(std::istringstream&);


    void setCOPY(const bool& Copy) {copyflag = Copy;}
    bool getCOPY() const {return copyflag;}
    void analyse(){}

  private:
    bool copyflag;
  };
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  class CELL : public InputBase, virtual public CCP4base
  {
    // Syntax: CELL   a  b  c  alpha  beta  gamma 
    //   Supply cell to override all cells in input files
  public:
    CELL();
    virtual ~CELL() {}
    Token_value parse(std::istringstream&);


    void setCELL(const scala::Scell& Cell) {cell = Cell;}
    scala::Scell  getCELL() const {return cell;}
    void analyse(){}

  private:
    scala::Scell cell;
  };
  //--------------------------------------------------------------
  class TEST : public InputBase, virtual public CCP4base
  {
    // Syntax: TESTFIRSTFILE
    //   Force Laue group testing on first file if multiple HKLIN
  public:
    TEST();
    virtual ~TEST() {}
    Token_value parse(std::istringstream&);


    static void setTESTFIRSTFILE(const bool& TestFirst) {testfirstfile = TestFirst;}
    bool TestFirstFile() const {return testfirstfile;}
    void analyse(){}

  private:
    static bool testfirstfile;
  };
  //--------------------------------------------------------------
  class NOTEST : public InputBase, virtual public CCP4base
  {
    // Syntax: NOTESTFIRSTFILE
    //   Cancel Laue group testing on first file if multiple HKLIN
  public:
    NOTEST();
    virtual ~NOTEST() {}
    Token_value parse(std::istringstream&);

    void analyse(){}
  };
  //--------------------------------------------------------------
  class ASSU : public InputBase, virtual public CCP4base
  {
    // Syntax: ASSUMESAMEINDEXING  [ON || OFF]
    //   Assume all HKLIN files have the same indexing scheme
  public:
    ASSU();
    virtual ~ASSU() {}
    Token_value parse(std::istringstream&);

    static void setASSUMESAMEINDEXING(const bool& AssumeSame);
    bool IsSetASSUMESAMEINDEXING() const {return (assumesameindexing >= 0);}
    int ASSUMESAMEINDEXING() const {return assumesameindexing;}
    void analyse(){}

  private:
    // = -1 unset, = 0 false, = +1 true
    static int assumesameindexing;
  };
  //--------------------------------------------------------------
  class NOASSU : public InputBase, virtual public CCP4base
  {
    // Syntax: DONOTASSUMESAMEINDEXING  ==  ASSUMESAMEINDEXING  OFF
    //   Do not assume all HKLIN files have the same indexing scheme
  public:
    NOASSU();
    virtual ~NOASSU() {}
    Token_value parse(std::istringstream&);

    void analyse(){}
  };
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  class EXCLUDE : public InputBase, virtual public CCP4base
  {
    // Syntax:
    //   EXCLUDE DATASET <datasetname> | <crystalname>/<datasetname>
    //     exclude dataset
    //   EXCLUDE BATCH [FILE|SERIES  <Jfile>]  <b1> <b2> <b3> ... | <b1> TO <b2>
    //     exclude batch list or range: if FILE key present, b1 etc
    //     refer to original (file) batch numbers (SERIES is the equivalent for a
    //     series of files defined with a wild-card)
    //   
  public:
    EXCLUDE();
    virtual ~EXCLUDE() {}
    Token_value parse(std::istringstream&);
    scala::BatchSelection BatchExclusions() const {return batchexclude;}

    void analyse(){}
  private:
    scala::BatchSelection batchexclude;
  };
  //--------------------------------------------------------------
  class CHOOSE : public InputBase, virtual public CCP4base
  {
    // Syntax:
    //  CHOOSE [SOLUTION] <N>
    //  CHOOSE LAUEGROUP <group>
    //  CHOOSE SPACEGROUP <group>
    //  choose solution from Laue group or space group list
  public:
   CHOOSE();
    virtual ~CHOOSE() {}
    Token_value parse(std::istringstream&);

    int getSolutionNumber() const {return solution_number;}
    std::string getChosenLaueGroup() const {return lauegroup;}
    std::string getChosenSpaceGroup() const {return spacegroup;}
    void analyse(){}

  private:
    int solution_number;
    std::string lauegroup;  
    std::string spacegroup;
  };
  //--------------------------------------------------------------
  class NEIGHBOUR : public InputBase, virtual public CCP4base
  {
    // Syntax:
    //  NEIGHBOUR <neighbourfraction>
    // Set fraction of neighbouring intensities to subtract from
    // axial intensities in systematic absence analysis
  public:
   NEIGHBOUR();
    virtual ~NEIGHBOUR() {}
    Token_value parse(std::istringstream&);

    float getNeighbourFraction() const {return neighbourfraction;}
    void analyse(){}

  private:
    float neighbourfraction;
  };

  //--------------------------------------------------------------
class SETTING : public InputBase, virtual public CCP4base
  //
  // Syntax:
  //   SETTING CELL-BASED | SYMMETRY-BASED | C2
  //
  // This affects the settings of primitive orthorhombic space groups
  // and centred monoclinic lattices. CELL-BASED follows the
  // International Tables conventions
  //
  // 1. Primitive orthorhombic lattices always have a <= b <= c
  //    Space groups 17 (P 2 2 21) & 18 (P 21 21 2) may be permuted eg
  //    to P 21 2 21. SYMMETRY-BASED means they are always in the "reference"
  //    setting irrespective of axial lengths
  //
  // 2. Centred monoclinic lattices will be I-centred rather than C-centred
  //    if that leads to a smaller beta angle, ie "Laue" group I2/m is allowed.
  //    SYMMETRY-BASED or C2 forces C2 always, never I2
  //
  // Default is CELL-BASED
{
public:
  SETTING();
  virtual ~SETTING() {}
  Token_value parse(std::istringstream&);
  int getSetting() const {return ITCsetting;}
  void analyse(){}

private:
  int ITCsetting;  // 0 CELL-BASED +1 SYMMETRY-BASED -1 C2
};
  //--------------------------------------------------------------
class OUTPUT : public InputBase, virtual public CCP4base
  //
  // Output settings
  // Syntax:
  //   OUTPUT [SUMMEDLIST <filename_summed>]
  //     SUMMEDLIST    write summed partials to formatted  file filename_summed
  //
{
public:
  OUTPUT();
  virtual ~OUTPUT() {}
  Token_value parse(std::istringstream&);
  std::string getSUMMEDLIST() const {return filename_summed;}
  void analyse(){}

private:
  std::string filename_summed;
};
//--------------------------------------------------------------
class MULTIPLY : public InputBase, virtual public CCP4base
  //
  // MULTIPLY  scale factor to apply on input
  // Syntax:
  //   MULTIPLY <input_scale>
  //
{
public:
  MULTIPLY();
  virtual ~MULTIPLY() {}
  Token_value parse(std::istringstream&);
  float getMULTIPLY() const {return scale;}
  void analyse(){}

private:
  float scale;
};
//--------------------------------------------------------------
class SCORE : public InputBase, virtual public CCP4base
  //
  // SCORE  parameters for score functions
  //   MAXMULTIPLICITY   maximum multiplicity to use in scoring
  //   
  // Syntax:
  //   SCORE MAXMULTIPLICITY <maximum_multiplicity>  [default 20]
  //
{
public:
  SCORE();
  virtual ~SCORE() {}
  Token_value parse(std::istringstream&);
  int getScoreMaxMultiplicity() const {return maxmult;}
  void analyse(){}

private:
  int maxmult;
};
//--------------------------------------------------------------
class ALLOW : public InputBase, virtual public CCP4base
  //
  // Syntax:
  //  ALLOW  OUTOFSEQUENCEFILES
  //        allow files in a series to be out of time sequence
  //   
{
public:
  ALLOW();
  virtual ~ALLOW() {}
  Token_value parse(std::istringstream&);
  bool AllowOutofSequenceFiles() const {return allowoutofsequencefiles;}
  void analyse(){}

private:
  bool allowoutofsequencefiles;
};
//--------------------------------------------------------------
class WAVELENGTH : public InputBase, virtual public CCP4base
  //
  // Syntax:
  //  WAVELENGTH  <wavelength>
  //        read wavelength
  //   
{
public:
  WAVELENGTH();
  virtual ~WAVELENGTH() {}
  Token_value parse(std::istringstream&);
  double Wavelength() const {return wavelength;}
  //  bool IsSet() const {return wavelength>0.0;}

  void analyse(){}

private:
  float wavelength;
  bool isset;
};
//--------------------------------------------------------------
class BLANK : public InputBase, virtual public CCP4base
  //
  // Read parameters for detecting blank images
  // Syntax:
  //  BLANK [RESO[FRACTION] <fraction>] [NEGATIVE <negativefraction>]
  //    <fraction>   blank image test will be done out to this fraction
  //                 of the maximum resolution, default 0.5
  //    <negativefraction> fraction of reflections with negative intensity
  //                 above which the image is considered to be blank
  //                 Default 0.3
  //   
{
public:
  BLANK();
  virtual ~BLANK() {}
  Token_value parse(std::istringstream&);
  float NullResolutionfraction() const {return nullResolutionfraction;}
  float NullNegativeReject() const {return nullNegativeReject;}

  void analyse(){}

private:
  float nullResolutionfraction;
  float nullNegativeReject;
  bool isset;
};
//--------------------------------------------------------------
class WILDFILE : public InputBase, virtual public CCP4base
  //
  // Set flag to treat HKLIN filename as if it had wild-cards
  // Syntax:
  //  WILDFILE [<template>]
  // Default template "*_###.*"
  //   
{
public:
  WILDFILE();
  virtual ~WILDFILE() {}
  Token_value parse(std::istringstream&);
  static bool WildFile() {return wildfile;}
  static std::string Template() {return ftemplate;}

  void analyse(){}

private:
  static bool wildfile;
  static std::string ftemplate;
};
//--------------------------------------------------------------
class POLARISATION : public InputBase, virtual public CCP4base
  //
  // Give polarisation factor for XDS/INTEGRATE files
  // Syntax:
  //  POLARISATION [XDS | MOSFLM] <polarisation_value>
  //   
{
public:
  POLARISATION();
  virtual ~POLARISATION() {}
  Token_value parse(std::istringstream&);

  double Polarisation() const {return polarisationfactor;}
  bool IsPolarisationSet() const {return set;}

  void analyse(){}

private:
  double polarisationfactor;
  bool set;
};
//--------------------------------------------------------------
//--------------------------------------------------------------
} // phaser_io

#endif
