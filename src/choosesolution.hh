// choosesolution.hh

#ifndef CHOOSESOLUTION
#define CHOOSESOLUTION

#include <string>
#include <vector>
#include "pgscore.hh"
#include "zone.hh"
#include "Output.hh"

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;

namespace scala
{
  class ChooseSolution
  {
  public:
    ChooseSolution(): solutionnumber(0), lauegroup(""), spacegroup("") {}
    ChooseSolution(const int& solution,
		   const std::string& Lauegroup,
		   const std::string& Spacegroup,
		   phaser_io::Output& output);

    // Return true if Laue group choice is specified, either as
    // solution number or as Laue group
    bool ChooseLaueGroup() const;

    // Return true if space group choice is specified
    bool ChooseSpaceGroup() const {return (spacegroup != "");}

    // Mark chosen group in Laue group solution list,
    // return solution number chosen, from 1,
    // =  0 if not found as solution number
    // = -1 if not found as Laue group
    int MarkAcceptedSubgroup(std::vector<PGscore>& SGscores);    

    // Return chosen group number in space group solution list,
    // return solution index chosen, = -1 if none
    int ChooseSpacegroup(const std::vector<scala::PossibleSpaceGroup>& AllGroups) const;
    
    int SolutionNumber() const {return solutionnumber;}
    std::string Lauegroup() const {return lauegroup;}
    std::string Spacegroup() const {return spacegroup;}
    // true if either Laue group or space group is I-centred monoclinic
    bool IsI2() const;
    // true if either Laue group or space group is rhombohedral R-lattice
    bool IsRhombohedralR() const;

  private:
    int solutionnumber;
    std::string lauegroup;  
    std::string spacegroup;
  };
}
#endif
