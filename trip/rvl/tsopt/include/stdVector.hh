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
#ifndef __RVL_STDVECTOR
#define __RVL_STDVECTOR


#include "stackBase.hh"

namespace TSOpt {

  using RVL::Writeable;
  using RVL::RVLException;
  using RVLAlg::Algorithm;
  using RVLAlg::ListAlg;
  using RVLAlg::LoopAlg;
  using RVLAlg::StateAlg;
  using RVLAlg::Terminator;


  /** Wrapper for the std::vector class, to conform to the methods
      in stackbase */
  template<typename State, 
	   typename Alloc>
  class StdVector : public adaptStackBase<State> {

  private:
    std::vector<State*> _ldclist;
    
  public:
    StdVector() : _ldclist() {}

    StdVector(StdVector<State, Alloc> const & s): _ldclist(s._ldclist) 
    {}

    ~StdVector() {
      //      cerr<<"begin stdV destr\n";
     for (int i=0; i<_ldclist.size(); ++i) {
	delete _ldclist.at(i);
      }

     _ldclist.clear();
     //     cerr<<"  end stdV destr\n";

    }

    /** Place state element at the back of the stack */
    void push_back(State const & t ) {
      State * tmp;
      tmp = new State(t);

      _ldclist.push_back(tmp);
    }

    /** Pop the state element at the top of the stack */
    void pop_back() {
      delete _ldclist.back();
      _ldclist.pop_back();
    }

    /** Returns stack size */
    int size() { return _ldclist.size(); }


    /** Allows access to a specific state at index idx,
        places element idx into state r (assumes r has a
        copy method!)  */
    
    State & at(int idx) { return *(_ldclist.at(idx));  }
    State const & at(int idx) const { return *(_ldclist.at(idx)); }
    
    void getAt(State & r, int idx){  
      try {
	r.copy(*(_ldclist.at(idx)));
      }
      catch (RVLException & e) {
	e<<"\ncalled from stdVector::getAt\n";
	throw e;
      }
    }

    /** Returns reference at the head of the stack */
    State & front() { return *(_ldclist.front()); }   
    State const & front() const { return *(_ldclist.front()); }
    
    /** Returns reference at the tail of the stack */
    State & back() { return *(_ldclist.back()); } 
    State const & back() const { return *(_ldclist.back()); }
    
    
    // made up functions 
    void placeAt(int idx, State const & t) { 
      State *tmp2;
      tmp2 = new State(t);
      _ldclist.at(idx) = tmp2;
    }


    // 0-offset for time indices
    int getOrigin() {return 0;}

  };
  

  


}


#endif
