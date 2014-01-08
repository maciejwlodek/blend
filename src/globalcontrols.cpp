
// globalcontrols.cpp

#include "scala.hh"

namespace scala
{
  //------------------------------------------------------------
  void GlobalControls::set_Chiral(const Chirality& ChiralFlag)
  {chiral = ChiralFlag;}
  
  //------------------------------------------------------------
  void GlobalControls::set_MinIsig(const double& MinIsig)
  {MinIsigRatio =  MinIsig;}
  
  //------------------------------------------------------------
  ScoreAccept::ScoreAccept(const double& AcceptanceFraction,
			   const double& MaximumScore)
    : threshold(AcceptanceFraction), scoremax(MaximumScore)
  {}
  //------------------------------------------------------------
  bool ScoreAccept::Accept(const double& score) const
  {
    if (threshold < 0.0) return true;
    return (score >= threshold*scoremax);
  }
}
//--------------------------------------------------------------}
