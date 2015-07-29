/*************************************************************************

Copyright Rice University, 2004.
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

#ifndef __RVLALG_SCALAR_TERMINATOR
#define __RVLALG_SCALAR_TERMINATOR

#include "alg.hh"
#include "functional.hh"
#include "except.hh"

/**  scalarterm.H
     All Terminators which watch the behavior of a scalar or integer
     to determine when to stop.

*/

namespace RVLAlg {

  /** This terminator contains an internal integer which it increments when called.
      After incrementing, the integer is compared to a maximum value.
      Stops when \f$count \geq maxcount\f$.
  */
class CountTerminator: public Terminator {
public:

  /** Initializes count to zero.
      When query is called, increments count by 1 and then tests to see if maxcount 
      has been reached.

      returns true when count == maxcount
      
              false when count < maxcount
  */
  CountTerminator( int maxcount ) : mc(maxcount), count(0), i(1) {}

  /** Initializes count to init.
      When query is called, increments count by inc and then tests to see if maxcount 
      has been reached.

      returns true when count >= maxcount

              false when count < maxcount

	      throws an error if inc is moving count away from maxcount
  */
  CountTerminator( int init, int maxcount, int inc ) 
    : mc(maxcount), count(init), i(inc) {
    if ( inc <= 0 ) {
      RVL::RVLException e; e << "Error in CountTerminator constructor: parameters cause infinite loop";
      throw e;
    }
  }
  
  ~CountTerminator() {}

  virtual bool query() {
    count += i;
    return (count >= mc );
  }

  int getCount() { return count; }

protected:
  int mc;
  int count;
  int i;

private:

  CountTerminator();
  CountTerminator(CountTerminator const &);
};

/** Terminator which takes a scalar and a maximum value.

    STOP: \f$x > maxval\f$
*/

template<class Scalar>
class MaxTerminator: public Terminator {
public:
  MaxTerminator(const Scalar & tx, Scalar maxval) : x(tx), mv(maxval) {}
  
  virtual bool query() {
    return (x > mv);
  }

protected:
  const Scalar & x;
  Scalar mv;
};

/** Terminator which takes a scalar and a minimum value.

    STOP: \f$x < minval\f$
*/

template<class Scalar>
class MinTerminator: public Terminator {
public:
  MinTerminator(const Scalar & tx, Scalar minval) : x(tx), mv(minval) {}
  
  virtual bool query() {
    cout<<"   MinTerminator: test scalar = "<<x<<" threshold = "<<mv<<endl;
    return (x < mv);
  }

protected:
  const Scalar & x;
  Scalar mv;
};

/** Terminator which takes a scalar and a minimum value.

    STOP: \f$x < minval\f$
*/
  
  template<class Scalar>
  class MinTerminatorFE: public Terminator {
  public:
    MinTerminatorFE(RVL::FunctionalEvaluation<Scalar> & _fx, Scalar minval) 
      : fx(_fx), mv(minval) {}
    
    virtual bool query() {
      Scalar x = fx.getValue();
      cout<<"   MinTerminator: test scalar = "<<x<<" threshold = "<<mv<<endl;
      return (x < mv);
    }
    
  protected:
    Scalar mv;
  private:
    RVL::FunctionalEvaluation<Scalar> & fx;
  };


}

#endif
