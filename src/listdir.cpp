// listdir.cpp
//
// Return contents of directory

#include "listdir.hh"

#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>

#ifdef _WIN32
/* Compiling for Windows */
  #include <windows.h>
#else
/* Compiling for UNIX / POSIX */
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <dirent.h>
#endif

// Compare filenames for sorting
bool operator < (const FileNameTime& a,const FileNameTime& b)
{return a.filename<b.filename;}
bool operator != (const FileNameTime& a,const FileNameTime& b)
{return (a.filename!=b.filename)||(a.modtime!=b.modtime);}

std::vector<FileNameTime> ListDir(const std::string& DirectoryName)
// Returns list of all files in directory
{
  std::vector<FileNameTime>  filenames;

#ifdef _WIN32
  /* Compiling for Windows */
  
  BOOL            fFinished;
  HANDLE          hList;
  TCHAR           szDir[MAX_PATH+1];
  WIN32_FIND_DATA FileData;
  
  // Get the proper directory path
  sprintf(szDir, "%s\\*", DirectoryName.c_str());
  
  // Get the first file
  hList = FindFirstFile(szDir, &FileData);
  if (hList == INVALID_HANDLE_VALUE)
    { 
      printf("No files found\n\n");
    }
  else
    {
      struct stat attrib;  // file attributes
      time_t mtime;
      // Traverse through the directory structure
      fFinished = FALSE;
      while (!fFinished) {
	// Check the object is a directory or not
	if (!(FileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
	  // file attributes
	  stat(FileNameTime(std::string(FileData.cFileName), &attrib));
	  mtime = attrib.st_mtime;
	  time_t mtime(0);
	  filenames.push_back
	    (FileNameTime(std::string(FileData.cFileName), mtime));
	}
	if (!FindNextFile(hList, &FileData)) {
	  if (GetLastError() == ERROR_NO_MORE_FILES)
	    fFinished = TRUE;
	}
      }
    }
  FindClose(hList);
  
#else
  /* Compiling for UNIX / POSIX */
  
  DIR *dir = opendir(DirectoryName.c_str());
  if(dir) {
    struct dirent *ent;
    struct stat attrib;  // file attributes
    time_t mtime;
    while((ent = readdir(dir)) != NULL)	{
      // Get last modified time for this file
      std::string fname = ent->d_name;
      std::string fullname = DirectoryName+"/"+fname;
      //      stat(ent->d_name, &attrib);  // file attributes
      int istat = stat(fullname.c_str(), &attrib);  // file attributes
      mtime = !istat ? attrib.st_mtime : 0;
      filenames.push_back(FileNameTime(fname, mtime));
      /*      //^
      time_t mt =  attrib.st_mtime;
      time_t at =  attrib.st_atime;
      time_t ct =  attrib.st_ctime;

      std::string smt(ctime(&mt));
      std::string sat(ctime(&at));
      std::string sct(ctime(&ct));

      std::cout << "listdir: " << fullname << " : "
		<< mt << " " << at
		<< " " << ct << "\n"
		<< smt << " " << sat << " " << sct << "\n";
      */
    }
  } else {
    fprintf(stderr, "Error opening directory %s\n",DirectoryName.c_str());
  }
#endif
  return filenames;
}
#ifdef _WIN32
//--------------------------------------------------------------
time_t ModTime(const std::string& fname)
{return 0;}
#else
//--------------------------------------------------------------
time_t ModTime(const std::string& fname)
// Get modification time for file "fname"
{
  struct stat attrib;  // file attributes
  //  time_t mtime;
  int istat = stat(fname.c_str(), &attrib);  // file attributes
  return !istat ? attrib.st_mtime : 0;
}
#endif

