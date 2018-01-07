// reindexscore.hh
//

#ifndef REINDEX_SCORE_HEADER
#define REINDEX_SCORE_HEADER

#include "hkl_datatypes.hh"
#include "setscores.hh"

namespace scala
{
  //--------------------------------------------------------------
  class ReindexScore : public  ReindexOp, public SetScores
  // A reindex operator and a set of scores, also a probability
  {
  public:
    ReindexScore() : likelihood(-1.0){};
    ReindexScore(const ReindexOp& Reindex)
      : ReindexOp(Reindex), SetScores(), likelihood(-1.0) {}

    // returns true if the score has been set
    bool IsSet() const {return (OverallCC().result().count > 0);}

    // Store likelihood value
    void SetLikelihood(const float& Likelihood)
    {likelihood = Likelihood;}
    // return likelihood (-1 if unset)
    float Likelihood() const {return likelihood;}

    friend bool operator < (const ReindexScore& a,const ReindexScore& b) {
      // for sorting by rank on likelihood or correlation coefficient score
      // bigger is better!
      if ((a.likelihood >= 0.0) && (b.likelihood >= 0.0)) {
	return (a.likelihood > b.likelihood);
      } else {
	return (a.OverallCC().result().val > b.OverallCC().result().val);
      }
    }
      
  private:
    float likelihood;
  };
}
#endif
