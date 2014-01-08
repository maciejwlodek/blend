// checklatticecentering.cpp
//
// Check hkl list for lattice absences: P, A, B, C, I, F, R (hexagonal setting)
//


#include "checklatticecentering.hh"


namespace scala {
  //--------------------------------------------------------------
  char CheckLatticeCentering (hkl_unmerge_list& hkl_list)
  // Returns character indicating lattice type  P, A, B, C, I, F, R (hexagonal setting)
  {
    reflection this_refl;
    hkl_list.rewind();  // point to beginning of list
    int h,k,l;
    int Ain  = 0; 
    int Aout = 0; 
    int Bin  = 0; 
    int Bout = 0; 
    int Cin  = 0; 
    int Cout = 0; 
    int Iin  = 0; 
    int Iout = 0; 
    int Fin  = 0; 
    int Fout = 0; 
    int Rin  = 0; 
    int Rout = 0; 

    while (hkl_list.next_reflection(this_refl) >= 0)  {
      Hkl hkl = this_refl.hkl();
      h = hkl.h();
      k = hkl.k();
      l = hkl.l();
      //  A    k+l = 2n
      //  B    h+l = 2n
      //  C    h+k = 2n
      //  I    h+k+l = 2n
      //  F    h,k,l = 2n or h,k,l = 2n+1
      //  R:H  -h+k+l = 3n  (hexagonal axes)
      // else P
      // Accumulate all totals
      if ((k+l)%2 == 0) {Ain++;} else {Aout++;}
      if ((h+l)%2 == 0) {Bin++;} else {Bout++;}
      if ((h+k)%2 == 0) {Cin++;} else {Cout++;}
      if ((h+k+l)%2 == 0) {Iin++;} else {Iout++;}
      if ((-h+k+l)%3 == 0) {Rin++;} else {Rout++;}
      if ((h%2 == 0 && k%2 == 0 && l%2 == 0) ||
	  (h%2 != 0 && k%2 != 0 && l%2 != 0)) {Fin++;} else {Fout++;}
    } // end loop reflections

    char LT = 'P';
    if (Fin > 0 && Fout == 0) LT = 'F';
    else if (Rin > 0 && Rout == 0) LT = 'H';
    else if (Iin > 0 && Iout == 0) LT = 'I';
    else if (Ain > 0 && Aout == 0) LT = 'A';
    else if (Bin > 0 && Bout == 0) LT = 'B';
    else if (Cin > 0 && Cout == 0) LT = 'C';
    return LT;
  }
  //--------------------------------------------------------------
  std::string MinimalLatticeGroup(const char& latType)
  // returns the minimal space group symbol for a group with given lattice type
  // viz. P 1; A,B,C,I 2; F 2 2 2; R 3
  {
    if (latType == 'R' || latType == 'H') return latType+std::string(" 3");
    if (latType == 'F') return latType+std::string(" 2 2 2");
    if (latType == 'P') return latType+std::string(" 1");
    return latType+std::string(" 2");
  }
}  // namespace scala
