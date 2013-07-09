/*************************************************************************

Copyright Rice University, 2004, 2005, 2006, 2007
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/
#ifndef __RVL_REALTIME
#define __RVL_REALTIME

#include "write.hh"
#include <limits>

namespace TSOpt {

  using RVL::Writeable;

  /** 
      RealTime is a Time class for adaptive time stepping, internally holding
      a double. Comparison operators are defined to take into account
      machine epsilon tolerances
  */
  
  class RealTime: public Time {
  private:
    double * _myTime;

  public:
    RealTime(): _myTime(NULL)  {}
    
    RealTime(RealTime const & t): _myTime(NULL)  {
      if (t._myTime) { 
	_myTime=new double; 
	*_myTime=*(t._myTime);
      }
      
    }
      
    ~RealTime() { 
      if (_myTime) { delete _myTime; _myTime = NULL;}
    }
    
    /** assignment - post-construction initialization */
    RealTime & operator=(Time const & t) {
      if(this != &t) {
	try {
	  RealTime const & myt = dynamic_cast<RealTime const &>(t);
	  *this=myt;
	  return *this;
	}
	catch (bad_cast) {
	  RVLException e;
	  e<<"Error: RealTime::operator=\n";
	  e<<"input not RealTime\n";
	  throw e;
	}
      }
      return *this;
      
    }
    

    /** assignment operator for a Realtime object */
    RealTime & operator=(RealTime const & t) {

      if (this != &t) {
	if (!_myTime) _myTime=new double;
	*_myTime=*(t._myTime);
      }

      return *this;
   
    }
    
  
    /** assignment operator for a double */
    RealTime & operator=(double it) {
      if (!_myTime) _myTime=new double;
      *_myTime = it;
      return *this;
    }
    


    /** ordering, lesser */
    bool operator<(Time const & t) const {
        try {
	RealTime const & myt = dynamic_cast<RealTime const &>(t);
	double eps = 1000*numeric_limits<double>::epsilon();
	if (myt._myTime) {
	  return ( ( *(myt._myTime) - *_myTime ) > eps); 
	  //return (*_myTime<*(myt._myTime));
	  //return ( ((*(myt._myTime) + eps - *_myTime) > 0) ||  ( (*(myt._myTime) - eps) - *_myTime > 0)); 
	}
	else {
	  RVLException e;
	  e<<"Error: RealTime::operator>\n";
	  e<<"input RealTime not initialized\n";
	  throw e;
	}
      }
      catch (bad_cast) {
	RVLException e;
	e<<"Error: RealTime::operator>\n";
	e<<"input not RealTime\n";
	throw e;
      }
    }
    
    /** ordering, greater */
    bool operator>(Time const & t) const {
      try {
	RealTime const & myt = dynamic_cast<RealTime const &>(t);
	double eps = 1000*numeric_limits<double>::epsilon();
	if (myt._myTime) {
	  return ( (*_myTime - *(myt._myTime)) > eps);
	  //	  return (*_myTime>*(myt._myTime));
	  //return ( (*_myTime + eps - *(myt._myTime)) > 0) ||  ( (*_myTime - eps) -  *(myt._myTime) > 0 ); 
	}
	else {
	  RVLException e;
	  e<<"Error: RealTime::operator>\n";
	  e<<"input RealTime not initialized\n";
	  throw e;
	}
      }
      catch (bad_cast) {
	RVLException e;
	e<<"Error: RealTime::operator>\n";
	e<<"input not RealTime\n";
	throw e;
      }
    }
    
    /** equality - implemented */
    bool operator==(Time const & t) const {
      return (!(((*this)<t) || ((*this)>t)));
    }
    
    /** inequality - implemented */
    bool operator!=(Time const & t) const {
      return (!(*this==t));
    }

    ostream & write(ostream & str) const {
      if (_myTime) {
	str<<"RealTime Time: " << getreal() << "\n";
      }
      else {
	str<<"RealTime: uninitialized\n";
      }
      return str;
    }


    /** returns a double representing the time */
    double getreal() const { return *(_myTime);}
		  
    RealTime * clone() const {
      return new RealTime(*this);	
    }
		  
  };

}

#endif
