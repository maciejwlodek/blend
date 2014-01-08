//  sysabszones.cpp

#include "pointgroup.hh"
#include "lattice.hh"
#include "zone.hh"
#include "sysabszones.hh"
#include "hkl_symmetry.hh"

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;


#include <vector>

namespace scala
{
  // Characteristic indices for axes and zones
  // static constants
  Hkl HKLTest::Qh00 = Hkl(1,0,0);
  Hkl HKLTest::Q0k0 = Hkl(0,1,0);
  Hkl HKLTest::Q00l = Hkl(0,0,1);
  
  Hkl HKLTest::Q0kl = Hkl(0,2,2);
  Hkl HKLTest::Qh0l = Hkl(2,0,2);
  Hkl HKLTest::Qhk0 = Hkl(2,2,0);
  
  Hkl HKLTest::Qhhl = Hkl(2,2,10);
  Hkl HKLTest::Qhkk = Hkl(10,2,2);
  Hkl HKLTest::Qhkh = Hkl(2,10,2);
  //--------------------------------------------------------------
  bool IfInAsu(const Hkl& hkl, const hkl_symmetry& Symm)
    // Returns true if hkl is in asymmetric unit defined by Symm
  {
    int new_isym;
    Hkl hkl_new = Symm.put_in_asu(hkl, new_isym);
    //^    std::cout << "IfInAsu: " << hkl.format() << " " << hkl_new.format() << "\n";
    // Is changed from input?
    return (hkl_new == hkl) ? true : false;
  }
  //--------------------------------------------------------------
  std::vector<Zone> SysAbsZones(const PointGroup& LG, const Chirality& chiral)
    // Make list of zones for testing for systematic absences
    // for Laue group LG
    // A zone may be an axis corresponding to a screw axis
    // or (if chiral != CHIRAL) a plane corresponding to a glide plane
    // Don't create zones if the equivalent systematic absences
    // would be covered by lattice absences
    //
    // For conditions see International Tables Vol A, table 3.2
    //
    // In cases where symmetry-equivalent zones or axes might be tested,
    // (eg 0kl & h0l in tetragonal), use asymmetric unit test (IfInAsu)
    // to decide which one to use
    //
  {
    std::vector<Zone> zones;
    std::string LGname = LG.RefLGname();
    hkl_symmetry Symm(LGname);

    //    char LatType = LG.GetLatType();

    // List of possible Laue group names
    //   P -1
    //   P 1 2/m 1
    //   C 1 2/m 1
    //   P m m m
    //   C m m m
    //   I m m m
    //   F m m m
    //   P 4/m
    //   P 4/m m m
    //   I 4/m
    //   I 4/m m m
    //   P -3
    //   P -3 1 m
    //   P -3 m 1
    //   H -3
    //   H -3 m
    //   P 6/m
    //   P 6/m m m
    //   P m -3
    //   P m -3 m
    //   I m -3
    //   I m -3 m
    //   F m -3
    //   F m -3 m

    if (LGname == "P -1")
      // Triclinic
      {
      }
    // Monoclinic
    else if (LGname == "P 1 2/m 1")
      {
	zones.push_back(Zone("b",2,Symm));  // 2(1) b
	if (chiral != CHIRAL) 
	  {
	    zones.push_back(Zone("b","a",Symm));  // a(b) glide
	    zones.push_back(Zone("b","c",Symm));  // c(b) glide
	    zones.push_back(Zone("b","n",Symm));  // n(b) glide
	  }
      }
    else if (LGname == "C 1 2/m 1" || LGname == "I 1 2/m 1")
      {
	if (chiral != CHIRAL) 
	  {
	    zones.push_back(Zone("b","c",Symm));  // c(b) glide
	  }
      }
    // Orthorhombic
    else if (LGname == "P m m m")
      {
	zones.push_back(Zone("a",2,Symm));  // 2(1) a
	zones.push_back(Zone("b",2,Symm));  // 2(1) b
	zones.push_back(Zone("c",2,Symm));  // 2(1) c
	if (chiral != CHIRAL) 
	  {
	    zones.push_back(Zone("a","b",Symm));  // b(a) glide
	    zones.push_back(Zone("a","c",Symm));  // c(a) glide
	    zones.push_back(Zone("a","n",Symm));  // n(a) glide
	    
	    zones.push_back(Zone("b","a",Symm));  // a(b) glide
	    zones.push_back(Zone("b","c",Symm));  // c(b) glide
	    zones.push_back(Zone("b","n",Symm));  // n(b) glide
	    
	    zones.push_back(Zone("c","a",Symm));  // a(c) glide
	    zones.push_back(Zone("c","b",Symm));  // b(c) glide
	    zones.push_back(Zone("c","n",Symm));  // n(c) glide
	  }
      }
    else if (LGname == "C m m m")
      {
	zones.push_back(Zone("c",2,Symm));  // 2(1) c
	if (chiral != CHIRAL) 
	  {
	    zones.push_back(Zone("a","c",Symm));  // c(a) glide
	    zones.push_back(Zone("b","c",Symm));  // c(b) glide
	    zones.push_back(Zone("c","a",Symm));  // a(c) glide
	  }
      }
    else if (LGname == "I m m m")
      {
	if (chiral != CHIRAL) 
	  {
	    zones.push_back(Zone("a","b",Symm));  // b(a) glide
	    zones.push_back(Zone("b","c",Symm));  // c(b) glide
	    zones.push_back(Zone("c","a",Symm));  // a(c) glide
	  }
      }
    else if (LGname == "F m m m")
      {
	if (chiral != CHIRAL) 
	  {
	    // d 100 glides:
	    zones.push_back(Zone("a","d",Symm));  // d(a) glide
	    zones.push_back(Zone("b","d",Symm));  // d(b) glide
	    zones.push_back(Zone("c","d",Symm));  // d(c) glide
	  }
      }
    // Tetragonal
    else if (LGname == "P 4/m")
      {
	zones.push_back(Zone("c",4,Symm));  // 4(1) c
	if (chiral != CHIRAL) 
	  {
	    zones.push_back(Zone("c","n",Symm));  // n(c) glide
	  }
      }
    else if (LGname == "I 4/m")
      {
	// don't test 4(2) as it's covered by the I lattice (flag -4)
	zones.push_back(Zone("c",-4,Symm));  // 4(1) c
	if (chiral != CHIRAL) 
	  {
	    zones.push_back(Zone("c","a",Symm));  // a(c) glide
	  }
      }
    else if (LGname == "P 4/m m m")
      {
	zones.push_back(Zone("c",4,Symm));  // 4(1) c
	if (IfInAsu(HKLTest::Qh00, Symm)) zones.push_back(Zone("a",2,Symm));  // 2(1) a
	if (IfInAsu(HKLTest::Q0k0, Symm)) zones.push_back(Zone("b",2,Symm));  // 2(1) b
	if (chiral != CHIRAL) 
	  {
	    zones.push_back(Zone("110","c",Symm));  // c(110) glide
	    zones.push_back(Zone("c","n",Symm));  // n(c) glide

	    // Alternatives for 0kl or h0l zones
	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","n",Symm));  // n(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","n",Symm));  // n(b) glide

	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","b",Symm));  // b(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","a",Symm));  // a(b) glide

	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","c",Symm));  // c(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","c",Symm));  // c(b) glide
	  }
      }
    else if (LGname == "I 4/m m m")
      {
	// don't test 4(2) as it's covered by the I lattice (flag -4)
	zones.push_back(Zone("c",-4,Symm));  // 4(1) c
	if (chiral != CHIRAL) 
	  {
	    zones.push_back(Zone("b","a",Symm));  // a(b) glide
	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","c",Symm));  // c(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","c",Symm));  // c(b) glide

	    zones.push_back(Zone("110","d",Symm));  // d(110) glide
	  }
      }
    // Trigonal
    else if (LGname == "P -3")
      {
	zones.push_back(Zone("c",3,Symm));  // 3(1) c      
      }
    else if (LGname == "P -3 1 m")
      {
	zones.push_back(Zone("c",3,Symm));  // 3(1) c      
	if (chiral != CHIRAL) 
	  {
	    zones.push_back(Zone("110","c",Symm));  // c(110) glide
	  }
      }
    else if (LGname == "P -3 m 1")
      {
	zones.push_back(Zone("c",3,Symm));  // 3(1) c      
	if (chiral != CHIRAL) 
	  {
	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","c",Symm));  // c(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","c",Symm));  // c(b) glide
	  }
      }
    else if (LGname == "H -3")
      {
      }
    else if (LGname == "H -3 m")
      {
	if (chiral != CHIRAL) 
	  {
	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","c",Symm));  // c(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","c",Symm));  // c(b) glide
	  }
      }
    // Hexagonal
    else if (LGname == "P 6/m")
      {
	zones.push_back(Zone("c",6,Symm));  // 6(1) c 
      }
    else if (LGname == "P 6/m m m")
      {
	zones.push_back(Zone("c",6,Symm));  // 6(1) c 
	if (chiral != CHIRAL) 
	  {
	    zones.push_back(Zone("110","c",Symm));  // c(110) glide
	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","c",Symm));  // c(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","c",Symm));  // c(b) glide
	  }
      }
    // Cubic
    else if (LGname == "P m -3")
      {
	// 2(1) along one of h00, 0k0, 00l
	if (IfInAsu(HKLTest::Qh00, Symm)) zones.push_back(Zone("a",2,Symm));  // 2(1) a
	if (IfInAsu(HKLTest::Q0k0, Symm)) zones.push_back(Zone("b",2,Symm));  // 2(1) b
	if (IfInAsu(HKLTest::Q00l, Symm)) zones.push_back(Zone("c",2,Symm));  // 2(1) c
	if (chiral != CHIRAL) 
	  {
	    // Alternatives for 0kl, h0l or hk0 zones
	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","b",Symm));  // b(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","c",Symm));  // c(b) glide
	    if (IfInAsu(HKLTest::Qhk0, Symm)) zones.push_back(Zone("c","a",Symm));  // a(c) glide

	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","c",Symm));  // c(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","a",Symm));  // a(b) glide
	    if (IfInAsu(HKLTest::Qhk0, Symm)) zones.push_back(Zone("c","b",Symm));  // b(c) glide

	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","n",Symm));  // n(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","n",Symm));  // n(b) glide
	    if (IfInAsu(HKLTest::Qhk0, Symm)) zones.push_back(Zone("c","n",Symm));  // n(c) glide
	  }
      }
    else if (LGname == "I m -3")
      {
	if (chiral != CHIRAL) 
	  {
	    // Alternatives for 0kl, h0l or hk0 zones
	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","b",Symm));  // b(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","c",Symm));  // c(b) glide
	    if (IfInAsu(HKLTest::Qhk0, Symm)) zones.push_back(Zone("c","a",Symm));  // a(c) glide
	  }
      }
    else if (LGname == "F m -3")
      {
	if (chiral != CHIRAL) 
	  {
	    // Alternatives for 0kl, h0l or hk0 zones
	    // d 100 glides:
	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","d",Symm));  // d(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","d",Symm));  // d(b) glide
	    if (IfInAsu(HKLTest::Qhk0, Symm)) zones.push_back(Zone("c","d",Symm));  // d(c) glide
	  }
      }
    else if (LGname == "P m -3 m")
      {
	// 4(1) along one of h00, 0k0, 00l
	if (IfInAsu(HKLTest::Qh00, Symm)) zones.push_back(Zone("a",4,Symm));  // 4(1) a
	if (IfInAsu(HKLTest::Q0k0, Symm)) zones.push_back(Zone("b",4,Symm));  // 4(1) b
	if (IfInAsu(HKLTest::Q00l, Symm)) zones.push_back(Zone("c",4,Symm));  // 4(1) c
	if (chiral != CHIRAL) 
	  {
	    // Alternatives for 0kl, h0l or hk0 zones
	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","n",Symm));  // n(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","n",Symm));  // n(b) glide
	    if (IfInAsu(HKLTest::Qhk0, Symm)) zones.push_back(Zone("c","n",Symm));  // n(c) glide

	    if (IfInAsu(HKLTest::Qhhl, Symm)) zones.push_back(Zone("110","n",Symm));  // n(110) glide
	    if (IfInAsu(HKLTest::Qhkk, Symm)) zones.push_back(Zone("011","n",Symm));  // n(011) glide
	    if (IfInAsu(HKLTest::Qhkh, Symm)) zones.push_back(Zone("101","n",Symm));  // n(101) glide
	  }
      }
    else if (LGname == "I m -3 m")
      {
	// don't test 4(2) as it's covered by the I lattice
	if (IfInAsu(HKLTest::Qh00, Symm)) zones.push_back(Zone("a",-4,Symm));  // 4(1) a 
	if (IfInAsu(HKLTest::Q0k0, Symm)) zones.push_back(Zone("b",-4,Symm));  // 4(1) b
	if (IfInAsu(HKLTest::Q00l, Symm)) zones.push_back(Zone("c",-4,Symm));  // 4(1) c 	  {
	if (chiral != CHIRAL) 
	  {
	    // Alternatives for 0kl, h0l or hk0 zones
	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","b",Symm));  // b(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","c",Symm));  // c(b) glide
	    if (IfInAsu(HKLTest::Qhk0, Symm)) zones.push_back(Zone("c","a",Symm));  // a(c) glide

	    if (IfInAsu(HKLTest::Qhhl, Symm)) zones.push_back(Zone("110","d",Symm));  // d(110) glide
	    if (IfInAsu(HKLTest::Qhkk, Symm)) zones.push_back(Zone("011","d",Symm));  // d(011) glide
	    if (IfInAsu(HKLTest::Qhkh, Symm)) zones.push_back(Zone("101","d",Symm));  // d(101) glide
	  }
      }
    else if (LGname == "F m -3 m")
      {
	// don't test 4(2) as it's covered by the F lattice
	if (IfInAsu(HKLTest::Qh00, Symm)) zones.push_back(Zone("a",-4,Symm));  // 4(1) a 
	if (IfInAsu(HKLTest::Q0k0, Symm)) zones.push_back(Zone("b",-4,Symm));  // 4(1) b
	if (IfInAsu(HKLTest::Q00l, Symm)) zones.push_back(Zone("c",-4,Symm));  // 4(1) c 	  {
	if (chiral != CHIRAL) 
	  {
	    // Alternatives for 0kl, h0l or hk0 zones
	    if (IfInAsu(HKLTest::Q0kl, Symm)) zones.push_back(Zone("a","d",Symm));  // d(a) glide
	    if (IfInAsu(HKLTest::Qh0l, Symm)) zones.push_back(Zone("b","d",Symm));  // d(b) glide
	    if (IfInAsu(HKLTest::Qhk0, Symm)) zones.push_back(Zone("c","d",Symm));  // d(c) glide

	    if (IfInAsu(HKLTest::Qhhl, Symm)) zones.push_back(Zone("110","c",Symm));  // c(110) glide
	    if (IfInAsu(HKLTest::Qhkk, Symm)) zones.push_back(Zone("011","a",Symm));  // a(011) glide
	    if (IfInAsu(HKLTest::Qhkh, Symm)) zones.push_back(Zone("101","b",Symm));  // b(101) glide
	  }
      }
    else
      // Invalid group
      Message::message(Message_fatal
		       ("SysAbsZones: invalid Laue group "+LGname));


    // Sort to put axes first
    // don't!    std::sort(zones.begin(),zones.end());

    return zones;
  }

}
