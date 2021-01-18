// alg.H
// created by ADP

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

#ifndef __RVL_ALG
#define __RVL_ALG

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <exception>
#include <stdexcept>
#include <typeinfo>
#include <string>
#include "except.hh"

namespace RVLAlg {

  using namespace std;

  /**
     Algorithm is the abstract base class for the algorithm package.  An
     Algorithm is something you can run().  We intentionally do not
     restrict the behavior any further.  This allows implementors to
     write a variety of constructors and access methods.  We recommend
     that Algorithms own very little internal data; This should mostly
     be intermediate calculations that are stored to avoid unnecessary
     calculations.  Constructors will often take references to external
     objects, and the algorithm object will then store this reference so
     that it can operate on these objects when told to run().

  */

  class Algorithm {
  public:

    Algorithm() {}
#if __cplusplus >= 201103L  // We are using C++11 or a later version
    virtual ~Algorithm() noexcept(false) {}
#else
    virtual ~Algorithm() {}
#endif

    /** 
	This is the only required member function.  When called
	the object should execute its expected behavior.
    */
    virtual void run() = 0;

  };

  /** A vacuous algorithm, for use as a placeholder
      where needed. */

  class NoAlg: public Algorithm {
  public:
    void run() {}
  };

  /**
     ListAlg behaves like a linked list of algorithms.  The
     object contains an algorithm and a ListAlg.

     ListAlg -> Algorithm 
     | Algorithm ListAlg

     When run(), the ListAlg runs its Algorithm, and then runs the trailing
     ListAlg.

     Note that since ListAlg is an Algorithm, you can in fact build 
     arbitrary trees of Algorithms by passing ListAlgs as both inputs
     to the constructor.  Such a tree would then run by running the left
     branch followed by the right branch.

     Return true if all algs in list are successful.
  */
  class ListAlg: public Algorithm {
  public:
    ListAlg(Algorithm & first): islist(false), one(first), two(*this) {}
    ListAlg(Algorithm & first, Algorithm & next)
      : islist(true), one(first), two(next) {}

    virtual void run() {
      one.run();
      if( islist ) two.run();
    }


  protected:
    bool islist;
    Algorithm & one;
    Algorithm & two;

  };

#define STOP_LOOP true
#define CONTINUE_LOOP false

  /** This is the abstract base class for a termination criterion.  
      In general, it is an object which we can query for a boolean value:

      true = stop looping
      false = criteria not satisfied, so continue iterating

      Note that the query function has no parameters.  Thus, it is imperative
      that appropriate constructors be written and references kept for specific
      termination criteria.  See the examples below.

      Also, generally speaking these are lightweight objects (i.e. don't
      own any large data members) and thus can be created, destroyed, passed
      as necessary.  The author intended for each object to be used in a
      single loop or nest of loops.

      It is intended that these objects be allowed to have side effects.
      For example, one of the desired behaviors we would like to eventually 
      include is to allow a user to interact with their problem as an
      algorithm progresses by modifying the objective function and constraints. 
      These I/O interactions and modifications would occur during a query to 
      a special terminator that had references to problem data as data members.

      A final note: the only difference between a Terminator and an
      Algorithm is that the Terminator has something to say when it's
      done.

  */
  class Terminator {
  public:
#if __cplusplus >= 201103L // We are using C++11 or a later version   
    virtual ~Terminator() noexcept(false) {}
#else
    virtual ~Terminator() {}
#endif
    virtual bool query() = 0;
  };


  /**
     LoopAlg takes an inside Algorithm and a Terminator.
     It runs the inside algorithm as long as the Terminator returns
     false, and the inside algorithm continues to run successfully.
   
     Thus, the inside alg will run zero or more times.  This mimics
     the behavior of a C while loop.
  */
  class LoopAlg: public Algorithm {
  public:
    LoopAlg(Algorithm & alg, Terminator & stop) : inside(alg), term(stop) { }

    virtual void run() {
      try {
	while( (!term.query()) ) {
	  inside.run();
	}
      }
      catch (RVL::RVLException & e) {
	e<<"\ncalled from LoopAlg::run\n";
	throw e;
      }
    }

  protected:
    Algorithm & inside;
    Terminator & term;
  };

  /**
     Behaves very much like LoopAlg, except that the inside alg
     runs AT LEAST once.  This mimics the behavior of a C do-while loop.
  */

  class DoLoopAlg: public LoopAlg {
  public:
    DoLoopAlg(Algorithm & alg, Terminator & stop) : LoopAlg(alg,stop) {}

    virtual void run() {
      inside.run();
      LoopAlg::run();
    }
  };

  /**
     CondListAlg is like a ListAlg, except that it only runs the 
     second alg conditionally, when the Terminator supplied tests false.

     This is useful for constructing backup methods if the original
     method fails.
  */
  class CondListAlg: public ListAlg {

  public:
    CondListAlg(Algorithm & first, Algorithm & next, Terminator & _stop)
      : ListAlg(first,next), stop(_stop) {}

    virtual void run() {
      one.run();
      if( !stop.query() ) {
	two.run();
      }
    }

  protected:

    Terminator & stop;
  };

  /**
     A StateAlg is an algorithm with the addition of an explicit state
     variable.  The StateAlg must know how to set its state to a specified
     value and be able to return a reference to its state when requested.  
     We do not require it to own a copy of the state at all times, as this
     may be expensive and the state can be assembled from other data when 
     necessary.  However, we expect that most StateAlgs will own a private
     data member which is the state.
  */
  template<class T>
  class StateAlg: public Algorithm {
  public:
    /*
      virtual void setState(const T & x) = 0;
    */
    virtual T & getState() = 0;
    virtual const T & getState() const = 0;
  };

  /** Uses a terminator to select a branch in an algorithm. 
   */
  class BranchAlg : public Algorithm {
    
  public:
    BranchAlg(Terminator & iftest, 
	      Algorithm & thenclause, 
	      Algorithm & elseclause )
      : thencl(thenclause), elsecl(elseclause), test(iftest) {}
    
    virtual void run() {
      if( test.query() )
	thencl.run();
      else
	elsecl.run();
    }
      
  protected:
    Algorithm & thencl;
    Algorithm & elsecl;
    Terminator & test;
  };


}

#endif
