//  average.cpp

#include "average.hh"

namespace scala
{
  //--------------------------------------------------------------
  IsigI average_Is(const reflection& this_refl)
  {
    Rtype I;
    Rtype sigI;
    average_I(this_refl, I, sigI);
    return IsigI(I,sigI);
  }
  //--------------------------------------------------------------
  void average_I(const reflection& this_refl, Rtype& I, Rtype& sigI)
    // Average I for all observations
    // Partials must have been added first
    // Variance weighted & no scaling or selections for now
  {
    int n = 0;
    Rtype w;

    I = 0.0;
    sigI = 0.0;

    for (int lobs = 0; lobs < this_refl.num_observations(); lobs++)
      {
	observation this_obs = this_refl.get_observation(lobs);
	if (this_obs.IsAccepted())
	  {
	    n += 1;
	    w = 1./(this_obs.sigI()*this_obs.sigI());
	    I += w * this_obs.I();
	    sigI += w;
	  }
      }
    if (n > 0) 
      {
	I = I/sigI;
	sigI = sqrt(1.0/sigI);
      }
  }
  //--------------------------------------------------------------
  IsigI average_Es(const reflection& this_refl,
		   const Normalise& NormRes)
  {
    Rtype I;
    Rtype sigI;
    average_E(this_refl, I, sigI, NormRes);
    return IsigI(I,sigI);
  }
  //--------------------------------------------------------------
  void average_E(const reflection& this_refl, Rtype& I, Rtype& sigI,
		 const Normalise& NormRes)

    // Average normalised I for all observations
    // Partials must have been added first
    // Variance weighted & no scaling or selections for now
  {
    int n = 0;
    Rtype w;

    I = 0.0;
    sigI = 0.0;
    float sSqr = this_refl.invresolsq();

    for (int lobs = 0; lobs < this_refl.num_observations(); lobs++)
      {
	observation this_obs = this_refl.get_observation(lobs);
	if (this_obs.IsAccepted())
	  {
	    int irun = this_obs.run();
	    float time = this_obs.time();
	    IsigI Is = NormRes.apply(this_obs.I_sigI(), sSqr,
				     irun, time);

	    w = 1./(Is.sigI()*Is.sigI());
	    I += w * Is.I();
	    sigI += w;
	    n++;
	  }
      }
    if (n > 0) 
      {
	I = I/sigI;
	sigI = sqrt(1.0/sigI);
      }
  }
} // scala
