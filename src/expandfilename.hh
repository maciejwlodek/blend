// expandfilename.hh

// Expand file name if it contains wild-cards "*" or "?"

#ifndef EXPANDFILENAME_HEADER
#define EXPANDFILENAME_HEADER

#include <vector>
#include <string>
#include "listdir.hh"

namespace scala 
{
  class ExpandFileName {
  public:
    // Expand file name if it contains wild-cards "*" or "?"
    // Add extension if filename is unique
    ExpandFileName(const std::string& fname,
		   const std::string& ext);

    // Get times & store with names
    ExpandFileName(const std::vector<std::string>& Names);

    // returns wild = true if fname contains wild cards
    bool Wild() const {return wild;}

    std::vector<FileNameTime> NameTimes() const {return names;} 

  private:
    std::string filename;
    std::string extension;
    bool wild;
    std::vector<FileNameTime> names; 
  };
}
#endif
