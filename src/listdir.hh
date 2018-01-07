// listdir.hh

#ifndef LISTDIR_HEADER
#define LISTDIR_HEADER

#include <vector>
#include <string>

#ifdef _WIN32
/* Compiling for Windows */
  #include <windows.h>
#else
/* Compiling for UNIX / POSIX */
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <dirent.h>
#endif

class FileNameTime{
  // A file & its last modification date
public:
  FileNameTime(){}
  FileNameTime(const std::string& FileName, const time_t& ModTime)
    : filename(FileName), modtime(ModTime) {}

  std::string FileName() const {return filename;};
  time_t ModTime() const {return modtime;}
  // Compare filenames for sorting
  friend bool operator < (const FileNameTime& a,const FileNameTime& b);
  friend bool operator != (const FileNameTime& a,const FileNameTime& b);

private:
  std::string filename;
  time_t modtime;  // last modification epoch
};


// Returns list of all files in directory
std::vector<FileNameTime> ListDir(const std::string& DirectoryName);

// Get modification time for file "fname"
time_t ModTime(const std::string& fname);

#endif
