#include "cdt.hh"

namespace TSOpt {


  Time & StdDiscreteTime::operator=(Time const & t) {
    if (this != &t) {
      try {
	StdDiscreteTime const & myt = dynamic_cast<StdDiscreteTime const &>(t);
	*this=myt;
	return *this;
      }
      catch (bad_cast) {
	RVLException e;
	e<<"Error: StdDiscreteTime::operator=\n";
	e<<"input not StdDiscreteTime\n";
	throw e;
      }
    }
    return *this;
  }
   
  StdDiscreteTime & StdDiscreteTime::operator=(StdDiscreteTime const & t) {
    if (this != &t) {
      if (!_it) _it=new int;
      *_it=*(t._it);
    }
    return *this;
  }

  StdDiscreteTime & StdDiscreteTime::operator=(int it) {
    if (!_it) _it=new int;
    *_it = it;
    return *this;
  }
  
  bool StdDiscreteTime::operator<(Time const & t) const {
    try {
      StdDiscreteTime const & myt = dynamic_cast<StdDiscreteTime const &>(t);
      if (myt._it) {
	return (*_it<*(myt._it));
      }
      else {
	RVLException e;
	e<<"Error: StdDiscreteTime::operator<\n";
	e<<"input StdDiscreteTime not initialized\n";
	throw e;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: StdDiscreteTime::operator<\n";
      e<<"input not StdDiscreteTime\n";
      throw e;
    }
  }
  
  bool StdDiscreteTime::operator>(Time const & t) const {
    try {
      StdDiscreteTime const & myt = dynamic_cast<StdDiscreteTime const &>(t);
      if (myt._it) {
	return (*_it>*(myt._it));
      }
      else {
	RVLException e;
	e<<"Error: StdDiscreteTime::operator>\n";
	e<<"input StdDiscreteTime not initialized\n";
	throw e;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: StdDiscreteTime::operator>\n";
      e<<"input not StdDiscreteTime\n";
      throw e;
    }
  }
  
  StdDiscreteTime & StdDiscreteTime::operator++() {
    if (_it) {
      (*(_it))++;
    }
    else {
      RVLException e;
      e<<"Error: StdDiscreteTime::operator++\n";
      e<<"input StdDiscreteTime not initialized\n";
      throw e;
    }    
    return *this;
  }

  StdDiscreteTime & StdDiscreteTime::operator--() {
    if (_it) {
      (*(_it))--;
    }
    else {
      RVLException e;
      e<<"Error: StdDiscreteTime::operator--\n";
      e<<"input StdDiscreteTime not initialized\n";
      throw e;
    }    
    return *this;
  }

  ostream & StdDiscreteTime::write(ostream & str) const {
    if (_it) {
      str<<"Discrete Time: index = "<<*_it<<"\n";
    }
    else {
      str<<"Discrete Time: uninitialized\n";
    }
    return str;
  }
  
}
