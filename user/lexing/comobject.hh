#ifndef _COMOBJECT_HPP_
#define _COMOBJECT_HPP_

#include "commoninc.hh"

using std::string;
using std::map;

class ComObject
{
protected:
  string _prefix;
public:
  ComObject(const string& prefix): _prefix(prefix) {;}
  virtual ~ComObject() {;}
  //-------------------------
  const string& prefix() { return _prefix; }
};

#endif
