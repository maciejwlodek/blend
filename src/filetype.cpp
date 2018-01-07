// \file rfiletype.cpp

#include <iostream>
#include <fstream>

#include "filetype.hh"

using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

//--------------------------------------------------------------
ReflectionFileType::ReflectionFileType(const std::string& fileName)
{
  file_type_code = UNKNOWN;
  std::ifstream filein;
  filein.open(fileName.c_str());
  if (!filein) {
    Message::message(Message_fatal
		     ("FILEIO: cannot open file "+fileName));
  }
  clipper::String line;
  char mtz[3];
  filein.read(mtz, 3);
  std::string smtz(mtz,3);
  if (smtz == "MTZ") {
    file_type_code = MTZ;
  } else {
    std::getline(filein, line);  // first line
    // remove any trailing lf or cr characters
    line = StringUtil::BuftoLine(smtz+line);
    if (line.substr(0,26) == "!OUTPUT_FILE=INTEGRATE.HKL") {
      // XDS integrate
      file_type_code = XDS_INTEGRATE;
    } else if (line.substr(0,17) == "!FORMAT=XDS_ASCII") {
      // XDS integrate
      file_type_code = XDS_ASCII;
    } else if (line.size() == 5 && CheckNumber(line.substr(0,5))) {
      // SCA Merged files have first line == "    1"  (i5)
      if (line.i() == 1) {
	std::getline(filein, line);
	if (line.size() == 5) { // 2nd line 5 characters
	  std::getline(filein, line);
	  // 3rd line should have 7 fields
	  if (line.split(" ").size() == 7) {
	    file_type_code = SCA_MERGED;  // merged file
	  }}}
    } else if (line.size() > 7  && CheckNumber(line.substr(0,5)) && isalpha(line[6])) {
      // FILE unmerged 1st 5 characters are a number, 7th character is a letter
      // unmerged file starts with number of symops & spacegroup
      int Nsymop = line.i();
      std::string spaceGroupName = line.substr(5, line.length()-1);
      std::getline(filein, line);
      if (line.split(" ").size() == 9) {
	// Skip 2 * Nsymop - 1 lines
	for (int i=0;i<2*Nsymop-1;i++) {std::getline(filein, line);}
	std::getline(filein, line);
	if (line.split(" ").size() == 12) { // reflection lines 12 fields
	  file_type_code = SCA_UNMERGED;  // unmerged file
	}}
    } else if (line.size() >= 250 && CheckNumber(line)) {
      // Long line all numbers, probably Saint format
      file_type_code = SAINT;
    } else if (line.size() >= 28  && CheckNumber(line)) {
      // Probably ShelX, 3i4,2f8.2, at least 4 fields
      //   (allowing for intensity running into l index)
      if (line.split(" ").size() >= 4) {
	file_type_code = SHELX;
      }
    }
  }
  //^  std::cout << file_type_code << " type\n";
  if (file_type_code == UNKNOWN) {
    Message::message(Message_warn
		     ("Unrecognised file type for file "+fileName));
  }
}
//--------------------------------------------------------------
//! Constructor as "MTZ", "XDS" or "SCA"
void ReflectionFileType::SetType(const std::string& filetypename)
{
  if (filetypename == "MTZ") {
    file_type_code = MTZ;
  } else if (filetypename == "XDS") {
    file_type_code =  XDS_ASCII;
  } else if (filetypename == "SCA") {
    file_type_code =  SCA_UNMERGED;
  }
}
//--------------------------------------------------------------
  //! return formatted file type
std::string ReflectionFileType::format() const
{
  if (file_type_code == NOTSET) {
    return "File type not set";
  } else if (file_type_code == UNKNOWN) {
    return "Unrecognised file type";
  } else if (file_type_code == MTZ) {
    return "MTZ";
  } else if (file_type_code == SCA_MERGED) {
    return "Scalepack merged";
  } else if (file_type_code == SCA_UNMERGED) {
    return "Scalepack unmerged";
  } else if (file_type_code == SHELX) {
    return "Shelx";
  } else if (file_type_code == SAINT) {
    return "Saint";
  } else if (file_type_code == XDS_INTEGRATE) {
    return "XDS INTEGRATE";
  } else if (file_type_code == XDS_ASCII) {
    return "XDS ASCII";
  }
  return "";
}
//--------------------------------------------------------------
//! true is file is MTZ
bool ReflectionFileType::IsTypeMTZ() const {
  return (file_type_code == MTZ);
}
//--------------------------------------------------------------
//! true is file is XDS
bool ReflectionFileType::IsTypeXDS() const {
  return ((file_type_code == XDS_ASCII)||(file_type_code == XDS_INTEGRATE));
}
//--------------------------------------------------------------
//! true is file is SCA
bool ReflectionFileType::IsTypeSCA() const {
  return ((file_type_code == SCA_MERGED)||(file_type_code == SCA_UNMERGED));
}
//--------------------------------------------------------------
//! true is file is SHELX
bool ReflectionFileType::IsTypeSHELX() const {
  return (file_type_code == SHELX);
}
//--------------------------------------------------------------
//! true is file is SAINT
bool ReflectionFileType::IsTypeSAINT() const {
  return (file_type_code == SAINT);
}
//--------------------------------------------------------------
//! true is file is SCA, Shelx or Saint
bool ReflectionFileType::IsTypeSSS() const {
  return ((file_type_code == SCA_MERGED)||(file_type_code == SCA_UNMERGED)||
	  (file_type_code == SHELX)||(file_type_code == SAINT));
}
//--------------------------------------------------------------
//! +1 unmerged, -1 merged, 0 unknown
int ReflectionFileType::Unmerged() const
{
  if (file_type_code == MTZ) return 0; // MTZ unknown
  if (IsTypeXDS()) return -1;  // assume XDS unmerged
  if (file_type_code == SHELX) return 0; // Shelx unknown
  if (file_type_code == SCA_MERGED) return +1; // SCA
  if (file_type_code == SCA_UNMERGED) return -1; // SCA
  if (file_type_code == SAINT) return -1; // Saint unmerged
  return 0;
}
//--------------------------------------------------------------
bool ReflectionFileType::CheckNumber(const std::string& s) const
{
  bool number = true;
  for (size_t i=0;i<s.size();++i) {
    if (!(((s[i]>='0') && (s[i]<='9')) || s[i] == '-' ||
	  s[i] == '+' || s[i] == '.' || s[i] == ' '))
      {number = false; break;}
  }
  return number;
}

