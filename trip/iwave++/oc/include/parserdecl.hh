#ifndef __RVL_IWAVE_PARSER_DECL
#define __RVL_IWAVE_PARSER_DECL

#include "usempi.h"
#include "parser.h"
#include "except.hh"
#include "std_cpp_includes.hh"

namespace RVL {

  /** template wrapper around Igor's parser functions. NOTE: returns
      bool - TRUE if successful, FALSE otherwise, more or less opposite
      to behaviour of underlying C functions. General case always false -
      only specializations may return true.

      Eventually to be used to construct bijection between RVL::Table and
      IWAVE::PARAMARRAY.
  */

  template<typename T>
  bool parse(PARARRAY & par, std::string name, T & val);

  template<>
  bool parse<string>(PARARRAY & par, std::string name, string & val);

  template<>
  bool parse<char>(PARARRAY & par, std::string name, char & val);

  template<>
  bool parse<short>(PARARRAY & par, std::string name, short & val);

  template<>
  bool parse<int>(PARARRAY & par, std::string name, int & val);

  template<>
  bool parse<long>(PARARRAY & par, std::string name, long & val);

  template<>
  bool parse<unsigned short>(PARARRAY & par, std::string name, unsigned short & val);

  template<>
  bool parse<unsigned int>(PARARRAY & par, std::string name, unsigned int & val);

  template<>
  bool parse<unsigned long>(PARARRAY & par, std::string name, unsigned long & val);

  template<>
  bool parse<float>(PARARRAY & par, std::string name, float & val);

  template<>
  bool parse<double>(PARARRAY & par, std::string name, double & val);

  template<typename T>
  void parse_except(PARARRAY & par, std::string name, T & val);

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
