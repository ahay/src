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
#ifndef __RVL_REAL_RNSTATE
#define __RVL_REAL_RNSTATE

#include "terminator.hh"
#include "sim.hh"
#include "rt.hh"

namespace TSOpt {
 
  typedef struct {
    /** current time */
    double time;

    /** state dim */
    size_t nu;     
    /** control dim */
    size_t nc;     
    /** state samples */
    double * u;  
    /** control samples */
    double * c;  

  } realrn;

  /** default constructor */
  inline int create_realrn(realrn & s) {
    s.time = 0.0;
    // s.it=0;
    s.nu=0;
    s.nc=0;
    s.u=NULL;
    s.c=NULL;
    return 0;

  }

  /** initialization from obvious data */
  inline int init_realrn(realrn & s, double t, size_t _nu, size_t _nc) {
    if (_nu<1 || _nc<1 || s.u || s.c) { 
      cout << "_nu = " << _nu << endl;
      cout << "_nc = " << _nc << endl;
      cout << "s.u = " << s.u[0] << endl;
      cout << "s.c =  " << s.c[0] << endl;
      return 1;
    }
    s.time = t; 
    //s.it=_it; 
    s.nu=_nu; 
    s.nc=_nc; 
    s.u = new (nothrow) double[s.nu];
    s.c = new (nothrow) double[s.nc];
    return 0;

  }


  /** copy constructor */
  inline int copy_realrn(realrn & s, realrn const & t) {

   if (t.nu<1||t.nc<1) {
      cout << "t.nu = " << t.nu << endl;
      cout << "t.nc = " << t.nc << endl;
      return 1;
    }
    
    //printf("t.it = %d, t.nu = %d,  t.nc = %d", t.it, t.nu, t.nc); 
    
   s.time=t.time; 
   //s.it=t.it; 
   s.nu=t.nu; 
   s.nc=t.nc; 
  
   //if (s.u) {delete [] s.u; s.u = NULL;}
   //if (s.c) {delete [] s.c; s.c = NULL;}

   s.u = new double[s.nu];
   s.c = new double[s.nc];

   for (size_t i=0;i<s.nu;i++) s.u[i]=t.u[i];
   for (size_t i=0;i<s.nc;i++) s.c[i]=t.c[i];

   return 0;
  }

  /** destructor */
  inline void delete_realrn(realrn & s) {
    if (s.u) {delete [] s.u;  s.u = NULL; }
    if (s.c) {delete [] s.c;  s.c = NULL; }
    
  }


  /** print contents */ 
  inline void fprint_realrn(realrn const & s, ostream & fp) {
    fp << "realrn struct: nu="<<s.nu<<" nc="<<s.nc; // 
    // <<" it="<<s.it
    fp << " time = " << s.time << "\n";
  
    if (s.u) {
      for (size_t i=0;i<s.nu;i++) fp<<"i="<<i<<" u="<<s.u[i]<<"\n";
    }
    else {
      fp<<"data uninitialized\n";
    }
    
    if (s.c) {
      for (size_t i=0;i<s.nc;i++) fp<<"i="<<i<<" c="<<s.c[i]<<"\n";
    }
    else {
      fp<<"control uninitialized\n";
    }
    
  }
    
  /** Simple class wrapper around realrn, with a RealTime object. 
      Used for adaptive time stepping
  */
 
  class RealRnState {
  private:
    /** rn data */
    realrn s;
    /** discrete time **/
    mutable RealTime dt;


  public:

    /** default constructor */
    RealRnState() {
      if (create_realrn(s)) {
	RVLException e;
	e<<"Error: RealRnState constructor from create_rn\n";
	throw e;
      }
      dt = 0.0;
      
    }

    /** copy constructor */
    RealRnState(RealRnState const & t) {

      if (copy_realrn(s,t.s)) {
	RVLException e;
	e<<"Error: RealRnState::copy from copy_rn\n";
	throw e;
      } 
      dt=t.dt; 
    }

    /** destructor */
    ~RealRnState() { delete_realrn(s); }
    
    /** initialization - separate from construction. Note
	synchronization of StdDiscreteTime with integer count.
    */
    void initialize(size_t nu=1, size_t nc=1, double t=0.0) {

      if (init_realrn(s,t,nu,nc)) {
	RVLException e;
	e<<"Error: RealRnState::initialize from init_rn\n";
	throw e;
      }
      dt = t;

    }

    /** implemented by cast to RealTime */
    void setTime(Time const & t) {
      try {
	dt=t;
	// need to set time comp of s here	
	// s.it=dt.getint();
	s.time = dt.getreal();
      }
      catch (RVLException & e) {
	e<<"\ncalled from RealRnState::setTime\n";
	throw e;
      }

      
    }

    /** synchronization of RealTime and rn data members */
    Time const & getTime() const {
      dt = s.time; 
      return dt;
      
    }

    /** copy method */
    void copy(RealRnState const & x) {
  
      //cout << "in RnState::copy" << endl;
      //fprint_rn(s, cout);
      //fprint_rn(x.s,cout);
      
      try {
	if (copy_realrn(s,x.s)) {
	  RVLException e;
	  e<<"Error: RealRnState::copy from copy_rn\n";
	  throw e;
	} 
	dt=x.dt;
      }
      catch (RVLException & e) {
	e<<"\ncalled from RealRnState::copy\n";
	throw e;
      }
    }

    /** access to rn data member */
    realrn & getrealrn() { return s; }
    realrn const & getrealrn() const { return s; }

    /** element-wise access to data, all const methods. 
	All "State" classes should implement these functions
    */
    double getStateElement(int idx) const { return s.u[idx];}
    void setStateElement(double datum, int idx) { s.u[idx] = datum; }

    double getControlElement(int idx) const { return s.c[idx]; }
    void setControlElement(double datum, int idx) { s.c[idx] = datum; }

    //int getTimeIndex() const {return s.it;}
    int getStateDim() const {return s.nu;}
    int getControlDim() const {return s.nc;}

    /** reveal your inner rn */
    ostream & write(ostream & str) const {
      str<<"RealRnState: inner rn = \n";
      fprint_realrn(s,str);
      return str;
    }
  };




}

#endif
