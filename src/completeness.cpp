// completeness.cpp
//
// Test completeness (usually of merged file)
//

#include "completeness.hh"
#include "numbercomplete.hh"

namespace scala {

  std::vector<double> Completeness(const hkl_unmerge_list& hkl_list, const int& Nbins)
  {
  ResoRange resrange = hkl_list.ResRange();
  resrange.SetNbins(Nbins);
  std::vector<int> NumRefSphere(Nbins,0);  // Number unique in sphere

  hkl_symmetry symm = hkl_list.symmetry();

  // number of symmetry operators including lattice centering
  float NumSymm = symm.Nsym();  

  reflection this_refl;
  hkl_list.rewind();

  while (hkl_list.next_reflection(this_refl) >= 0)  {
    Rtype invresolsq = this_refl.invresolsq();
    // Resolution bin
    int mres = resrange.bin(invresolsq);

    bool Centric = symm.is_centric(this_refl.hkl());

    // Multiplicity is number of times this reflection will occur in a
    //     complete sphere of data. This is Nsym/Epsilon, multiplied by 2
    //     for acentric reflections (since Nsym symmetry operations generate
    //     a hemisphere of acentric data, but all centric reflections)
    // NOTE this is not the same as the multiplicity used for weighting
    //     <I> in Wilson plot calculations (see Iwasaki & Ito, Acta Cryst,
      //     A33,227-229(1977))
    float epsiln = symm.epsilon(this_refl.hkl());
    int multcy;
    if (epsiln == 0.0) {
      // systematic absence
      multcy = 0;
    } else {
      multcy = Nint(float(NumSymm)/epsiln);
      if (!Centric) multcy = multcy*2;
    }
    NumRefSphere[mres] += multcy;      // Total unique in sphere
  }

  // Get number of reflections in each resolution bin in complete sphere
  std::vector<int> Nrefres(Nbins,0);  // number
  std::vector<int> Nrefacen(Nbins,0);  // number of acentrics
  // Get numbers in each resolution shell
  NumberComplete(resrange, symm, hkl_list.Cell(), Nrefres, Nrefacen);

  std::vector<double> completeness(Nbins, 0.0);
  for (int i=0;i<Nbins;++i) {
    if (Nrefres[i] > 0) {
      completeness[i] = NumRefSphere[i]/double(Nrefres[i]);
    }
  }
  return completeness;
}
} // namespace scala
