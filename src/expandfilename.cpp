// expandfilename.cpp

// Expand file name if it contains wild-cards "*" or "?"
#include <algorithm>

#include "expandfilename.hh"
#include "listdir.hh"
#include "comparestring.hh"
#include "scala_util.hh"
#include "string_util.hh"

namespace scala 
{
  //--------------------------------------------------------------
  ExpandFileName::ExpandFileName(const std::string& fname,
				 const std::string& ext)
  // Expand file name if it contains wild-cards "*" or "?"
  // Add extension if filename is unique

  // sets wild = true if fname contains wild cards
  {
    if ((fname.find("*") == std::string::npos) &&
	(fname.find("?") == std::string::npos)) {
      // Nothing to expand
      wild = false;
      filename = fname;
      scala::AddFileExtension(filename, ext);
      names.push_back(FileNameTime(filename, time_t(0)));
    } else {
      // Wild cards found
      wild = true;
      // All files in this directory
      std::string dir = scala::Directory(StringUtil::Trim(fname));
      if (dir == "") {dir = ".";}
      std::vector<FileNameTime> allnames = ListDir(dir);
      // Sort into alphabetical order, in case
      std::sort(allnames.begin(), allnames.end());
      // Test each name against the target & accumulate
      std::string filename = StringUtil::Trim(fname); // remove any leading & trailing spaces
      size_t j = filename.rfind('/');
      if (j != std::string::npos) {
	filename = filename.substr(j+1); // get base name without any directory
      }
      for (size_t i=0;i<allnames.size();i++) {
	if (CompareString(allnames[i].FileName(), filename)) {
	  names.push_back(FileNameTime
			  (dir+"/"+allnames[i].FileName(),
			   allnames[i].ModTime()));
	  //^
	  //	  std::cout << "File: " << allnames[i].FileName() << " "
	  //		    <<  allnames[i].ModTime() << "\n"; //^
	}
      }
      // "names" now contains list of matching file names &
      // their modification times
    }
  }
  //--------------------------------------------------------------
  ExpandFileName::ExpandFileName(const std::vector<std::string>& Names)
  // Get times & store with names
  {
    filename = "";
    extension = "";
    wild = (Names.size() > 1);
    names.clear();
    if (Names.size() > 0) {
      for (size_t i=0;i<Names.size();++i) {
	names.push_back(FileNameTime(Names[i], ModTime(Names[i])));
      }
      // Sort into alphabetical order, in case
      std::sort(names.begin(), names.end());
    }
  }
  //--------------------------------------------------------------

}
