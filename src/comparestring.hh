
#ifndef COMPARESTRING_HEADER
#define COMPARESTRING_HEADER

#include <string>

#undef CompareString // it can be defined by Windows headers

  // Returns true if strings match, allowing for
  // wild-cards "*" or "?" in string b only
bool CompareString(const std::string& a, const std::string& b);

#endif
