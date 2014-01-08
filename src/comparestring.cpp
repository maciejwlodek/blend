// comparestring.cpp
//
// Compare two strings, including 
// wild-cards "*" or "?" in string b only


#include <string>
#include <iostream>
#include <vector>

#include "comparestring.hh"
#include "scala_util.hh"
#include "string_util.hh"

int TestString(const std::string& a, const int& ja, const std::string& bs) 
{
  // TestString(a,ja,bs) compares substring bs (which may contain "?") to
  // substring a[ja]..., returns index in a of match, -1 if no match

  // If either string is empty, no match (if both empty they match!)
  int la = a.size();  // length of string a
  int las = la - ja;  // length of substring a
  int lb = bs.size();     // length of b
  if (lb == 0 && las == 0) return  0;
  if (lb == 0 || las == 0) return -1;
  // substring a cannot be shorter than b
  if (las < lb) return -1;
  
  // Loop possible start points in a for substring bs
  for (int ia=ja;ia<=la-lb;ia++) {
    bool success = true;
    // Loop characters is bs to test
    for (int ib=0;ib<lb;ib++) {
      if (bs[ib] != '?') {   // '?'  matches anything
	if (bs[ib] != a[ia+ib]) {
	  success = false;    // fail
	  break;
	}
      }
    }
    if (success) {return ia;}  // string found
  }
  return -1; // fail  
}

bool CompareString(const std::string& a, const std::string& b)
{
  // Returns true if strings match, allowing for
  // wild-cards "*" or "?" in string b only

  // Split at "*" characters
  std::vector<std::string> bs = StringUtil::split(b, "*");

  if (bs.size() == 0) {return true;}  // blank or "*" matches anything
  int ia = 0;
  for (size_t i=0;i<bs.size();i++) {
    int ia1 = TestString(a, ia, bs[i]);
    if (ia1 < 0) {return false;} // not found, fail
    // string found starting at index ia1
    // This should either be == ia unless the previous
    //  test character is "*"
    if (ia1 > ia) {
      int ib = 0;
      if (i > 0) ib = bs[i-1].size();
      if (b[ib] != '*') return false;
    }
    ia = ia1 + bs[i].size();
    if (ia >= int(a.size())) return true;  
  }
  int j = bs.size() - 1;
  if (b[0] == '*' || bs.size() > 1) {
    // Check end of string
    ia = a.size() - bs[j].size();
    if (TestString(a, ia, bs[j]) == ia) {return true;}
  }
  if (b[b.size()-1] == '*') return true;
  return false;
}
