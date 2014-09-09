#ifndef __RVL_IWAVE_PARSER
#define __RVL_IWAVE_PARSER

#include "usempi.h"
#include "parser.h"
#include "parserdecl.hh"

namespace RVL {

  template<typename T>
  bool parse(PARARRAY & par, std::string name, T & val) { return false; }

  template<>
  bool parse<string>(PARARRAY & par, std::string name, string & val) {
    char * p;
    int res=ps_flcstring(par,name.c_str(),&p);
    if (!res) {
      val=p;
      userfree_(p);
      return true;
    }
    cerr<<"parser: failed on key="<<name<<endl;
    ps_printall(par,stderr);
    return false;
  }

  template<>
  bool parse<char>(PARARRAY & par, std::string name, char & val) {
    return !ps_flchar(par,name.c_str(),&val);
  }

  template<>
  bool parse<short>(PARARRAY & par, std::string name, short & val) {
    return !ps_flshort(par,name.c_str(),&val);
  }

  template<>
  bool parse<int>(PARARRAY & par, std::string name, int & val) {
    return !ps_flint(par,name.c_str(),&val);
  }

  template<>
  bool parse<long>(PARARRAY & par, std::string name, long & val) {
    return !ps_fllong(par,name.c_str(),&val);
  }

  template<>
  bool parse<unsigned short>(PARARRAY & par, std::string name, unsigned short & val) {
    return !ps_flushort(par,name.c_str(),&val);
  }

  template<>
  bool parse<unsigned int>(PARARRAY & par, std::string name, unsigned int & val) {
    return !ps_fluint(par,name.c_str(),&val);
  }

  template<>
  bool parse<unsigned long>(PARARRAY & par, std::string name, unsigned long & val) {
    return !ps_flulong(par,name.c_str(),&val);
  }

  template<>
  bool parse<float>(PARARRAY & par, std::string name, float & val) {
    return !ps_flfloat(par,name.c_str(),&val);
  }

  template<>
  bool parse<double>(PARARRAY & par, std::string name, double & val) {
    return !ps_fldouble(par,name.c_str(),&val);
  }

  /** throws exception if returns false */
  template<typename T>
  void parse_except(PARARRAY & par, std::string name, T & val) {
    if (!parse<T>(par,name,val)) {
      RVLException e;
      e<<"Error: parse_except\n";
      e<<"  failed to parse key = "<<name<<"\n";
      throw e;
    }
  }
}

#endif
