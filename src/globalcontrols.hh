// globalcontrols.hh

#ifndef GLOBAL_CONTROLS
#define GLOBAL_CONTROLS

#include "Output.hh"

using phaser_io::LOGFILE;


namespace scala
{
  //=================================================================

  //--------------------------------------------------------------
  static float ToleranceDefault = 2.0;  // Default value

  class GlobalControls
  {
  public:
    // Null constructor
    GlobalControls() : chiral(CHIRAL),
		       LatTol(ToleranceDefault), ReindexSet(false),
		       assumesameindexing(false) {}
    GlobalControls(const Chirality chiral_in)
      : LatTol(ToleranceDefault), ReindexSet(false),
	assumesameindexing(false)
    {chiral = chiral_in;}

    void set_Chiral(const Chirality& ChiralFlag);
    void set_MinIsig(const double& MinIsig);

    void set_OriginalLattice(const bool& OrigLat) {OriginalLat = OrigLat;}
    void set_Lauegroup(const std::string& LG) {LaueGroup = LG;}
    void set_Spacegroup(const std::string& SG) {SpaceGroup = SG;}
    void set_Reindex(const ReindexOp& H) {Reindx = H; ReindexSet=true;}
    void set_LatticeTolerance(const double& Tol) {LatTol = Tol;}
    void set_SysAbsCheck(const bool& abscheck) {CheckSysAbs = abscheck;}
    void set_TestFirstFile(const bool& TestFirst) {testfirstfile = TestFirst;}
    void set_AssumeSameIndexing(const bool& Assumesameindexing)
      {assumesameindexing = Assumesameindexing;}
    void set_AllowI2(const int& lI2) {allowI2 = lI2;}
    void set_MaxMult(const int& Maxmult) {maxmult = Maxmult;}

    // Access
    Chirality GetChiral() const {return chiral;}
    double GetMinIsig() const {return MinIsigRatio;}

    bool OriginalLattice() const {return OriginalLat;}
    std::string Lauegroup() const {return LaueGroup;}
    std::string Spacegroup() const {return SpaceGroup;}
    ReindexOp Reindex() const {return Reindx;}
    bool IsReindexSet() const {return ReindexSet;}
    double LatticeTolerance() const {return LatTol;}
    bool SysAbsCheck() const {return CheckSysAbs;}
    bool TestFirstFile() const {return testfirstfile;}
    bool AssumeSameIndexing() const {return assumesameindexing;}
    int AllowI2() const {return allowI2;}
    int MaxMult() const {return maxmult;}

  private:
    Chirality chiral;
    double MinIsigRatio;
    bool OriginalLat;    // true to use original lattice
    double LatTol;       // tolerance (degrees) on lattice dimensions
    std::string LaueGroup;   // non-blank to choose Lauegroup
    std::string SpaceGroup;  // non-blank to choose Spacegroup
    ReindexOp Reindx;       //   & reindex operator
    bool ReindexSet;
    bool CheckSysAbs;    // normally true to check systematic absences
    bool testfirstfile;  // test first file for Laue group
    bool assumesameindexing;  // assume all HKLIN files have same indexing
    int  allowI2;        // Allow I2 setting of C2
    int maxmult;         // Maximum multiplicity for scoring
  };
//=================================================================
class ScoreAccept
{
public:
  ScoreAccept() : threshold(-1.0) {};

  ScoreAccept(const double& AcceptanceFraction,
	      const double& MaximumScore);

  bool Accept(const double& score) const;

  bool IfSet() const {return (threshold >= 0.0);}

  double Threshold() const {return threshold*scoremax;}

private:
  double threshold;
  double scoremax;
};


}
#endif
