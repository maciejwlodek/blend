// \file filetype.hh

#ifndef REFLECTIONFILETYPE_HEADER
#define REFLECTIONFILETYPE_HEADER

#include <string>

#include <clipper/clipper.h>

#include "string_util.hh"

//! Class to automatically determine the file type of a reflection file
class ReflectionFileType {
public:
  //! Reflection file type
  enum FileTypeCode {NOTSET, ABSENT, UNKNOWN,
		     MTZ, SCA_MERGED, SCA_UNMERGED, SHELX,
		     SAINT, XDS_INTEGRATE, XDS_ASCII};

  ReflectionFileType() : file_type_code(NOTSET) {}
  //! Constructor: try to determine file type, close file after testing
  ReflectionFileType(const std::string& fileName);
  //! Set type code
  void SetTypeCode(const FileTypeCode& code) {file_type_code = code;}

 //! Set type as "MTZ", "XDS" or "SCA"
  void SetType(const std::string& filetypename);

  //! return file type code
  FileTypeCode FileType() const {return file_type_code;}
  //! true is file is MTZ
  bool IsTypeMTZ() const;
  //! true is file is XDS
  bool IsTypeXDS() const;
  //! true is file is SCA
  bool IsTypeSCA() const;
  //! true is file is SHELX
  bool IsTypeSHELX() const;
  //! true is file is SAINT
  bool IsTypeSAINT() const;

  //! true is file is SCA, SHELX or SAINT
  bool IsTypeSSS() const;

  //! +1 unmerged, -1 merged, 0 unknown
  int Unmerged() const;
  
  //! return formatted file type
  std::string format() const;
  
private:
  FileTypeCode file_type_code;

  // returns true if all characters in s are numeric
  bool CheckNumber(const std::string& s) const;
};

#endif
