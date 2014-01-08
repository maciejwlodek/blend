//
// writesummedlist.cpp
//
// Write out summed partials to file
//

#include "writesummedlist.hh"

namespace scala {
  void WriteSummedList(const hkl_unmerge_list& TestList,
		       const std::string& filename)
  {
    FILE* file = fopen(filename.c_str(), "w");
    if (file == NULL)
      Message::message(Message_fatal("Can't open file "+filename));

    // Loop all reflections
    TestList.rewind();  // point to beginning of list
    reflection this_refl;
    observation this_obs;
    Hkl hkl;
    TestList.rewind();

    // Maximum hkl (original) indices in list
    while (TestList.next_reflection(this_refl) >= 0) {
	while (this_refl.next_observation(this_obs) >= 0) {
	  fprintf(file,
		  "%5d%5d%5d %8d %10.1f %10.1f %8.2f %8.1f %8.1f\n",
		  this_obs.hkl_original().h(),
		  this_obs.hkl_original().k(),
		  this_obs.hkl_original().l(),
		  this_obs.Batch(),
		  this_obs.I(), this_obs.sigI(), this_obs.phi(),
		  this_obs.XYdet().first, this_obs.XYdet().second);
	}
    }
  }
}
