// checknullbatches.hh
//
// Check for and eliminate batches with essentially zero intensity
// (ie from blank images): this is surprisingly difficult
//
// Control parameters from input or default:
//   ResoFraction        fraction of maximum resolution to use (ie low reso only)
//   NullNegativeReject  proportion of negative reflections above which to reject 
// 
// Algorithm: (far from perfect)
//   For each batch accumulate Sum(I), Sum(I/sigI), N(I<0), N, for low resolution
//     reflections only (within ResoFraction of maximum resolution)
//     hence Zf(batch) = fraction of negatives
//   Split batches into 1->4 parts at least 10 batches/part
//   For each part, get <Zf>, SD(Zf) over batches in the part
//   Discard the part with the largest SD(Zf) (unless only one part), average the rest
//     to get <Zf>, SD(Zf)
//   Set Zf acceptance level = Max(NullNegativeReject, <Zf> + 4 SD(Zf) )
//     ie this may increase limit if most images are very weak
//   

//

#ifndef CHECKNULLBATCHES_HEADER
#define CHECKNULLBATCHES_HEADER

#include "hkl_unmerge.hh"
#include "Output.hh"

namespace scala {
  int CheckNullBatches(hkl_unmerge_list& hkl_list,
		       const file_select& file_sel,
		       phaser_io::Output& output);
}

#endif

