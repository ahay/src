#ifndef __RVL_IWAVE_PARSER_DECL
#define __RVL_IWAVE_PARSER_DECL

#include "usempi.h"
#include "except.hh"
#include "std_cpp_includes.hh"

namespace RVL {

  /** Class to hold key=value pairs - could do this with map, but
      this is so simple. Essentially a struct with constructors for
      easy build without a string of assignment statements
  */

  class STRING_PAIR {
  public:
    string key;
    string val;
    STRING_PAIR(string _key = "", string _val = ""): key(_key), val(_val) {}
    STRING_PAIR(STRING_PAIR const & pr): key(pr.key), val(pr.val) {}
    ~STRING_PAIR() {}
    STRING_PAIR operator=(STRING_PAIR const & rhs) {
      key=rhs.key;
      val=rhs.val;
      return *this;
    }
    
    ostream & write(ostream & str) const { 
      str<<"STRING_PAIR: key="<<this->key<<" val="<<this->val<<"\n";
      return str;
    }
  };

}



#endif
