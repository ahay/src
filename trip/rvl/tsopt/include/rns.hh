/*************************************************************************

Copyright Rice University, 2004, 2005, 2006, 2007, 2008
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
#ifndef __RVL_RNSTATE
#define __RVL_RNSTATE

#include "terminator.hh"
#include "sim.hh"
#include "cdt.hh"


namespace TSOpt {
 
  /** Basic array type struct - stand-in for many non-object data
      containers with float data and some notion of discrete time,
      that might be encountered in simulator definitions.*/

  typedef struct {
    /** time index */
    mutable int it;     
    /** state dim */
    size_t nu;     
    /** control dim */
    size_t nc;     
    /** state samples */
    float * u;  
    /** control samples */
    float * c;  

  } rn;

  /** default constructor */
  int create_rn(rn & s);
  /** initialization from obvious data */
  int init_rn(rn & s, int _it, size_t _nu, size_t _nc);
  /** copy constructor */
  int copy_rn(rn & s, rn const & t);
  /** destructor */
  void delete_rn(rn & s);
  /** print contents */ 
  void fprint_rn(rn const & s, ostream & fp);
    
  /** Simple class wrapper around rn, implicitly wrapping time index
      in StdDiscreteTime object */
 
  class RnState {
  private:
    /** rn data */
    rn s;
    /** discrete time **/
    mutable StdDiscreteTime dt;

  public:

    /** default constructor */
    RnState();

    /** copy constructor */
    RnState(RnState const & t);

    /** destructor */
    ~RnState();
    
    /** initialization - separate from construction. Note
	synchronization of StdDiscreteTime with integer count.
    */
    void initialize(size_t nu=1, size_t nc=1, int it=0);

    /** implemented by cast to StdDiscreteTime */
    void setTime(Time const & t);
    /** synchronization of StdDiscreteTime and rn data members */
    Time const & getTime() const;

    /** copy method */
    void copy(RnState const & x);

    /** access to rn data member */
    //    rn & getrn() { return s; }
    rn const & getrn() const { return s; }

    /** element-wise access to data, all const methods. 
	All "State" classes should implement these functions
    */
    float getStateElement(size_t idx) const { return s.u[idx];}
    void setStateElement(float datum, size_t idx) { s.u[idx] = datum; }

    float getControlElement(size_t idx) const { return s.c[idx]; }
    void setControlElement(float datum, size_t idx) { s.c[idx] = datum; }

    int getTimeIndex() const {return s.it;}
    size_t getStateDim() const {return s.nu;}
    size_t getControlDim() const {return s.nc;}

    /** reveal your inner rn */
    ostream & write(ostream & str) const {
      str<<"RnState: = "<<this<<"\n";
      str<<"RnState: inner rn = "<<&s<<endl;
      str<<"RnState: u = "<<s.u<<endl;
      str<<"RnState: c = "<<s.c<<endl;
      fprint_rn(s,str);
      return str;
    }
  };


 
}

#endif
