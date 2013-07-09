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
#ifndef __RVL_STACKBASE
#define __RVL_STACKBASE


namespace TSOpt {

  using RVL::Writeable;
  using RVL::RVLException;
  using RVLAlg::Algorithm;
  using RVLAlg::ListAlg;
  using RVLAlg::LoopAlg;
  using RVLAlg::StateAlg;
  using RVLAlg::Terminator;

  /** Base class for all stack types for use with RASim and CPSim 
      (or any other sim subclass that requires saving of states)
  */
  template <typename State>
  class stackBase {

  public:
    stackBase() {}
    stackBase(stackBase const & rhs) {}
    virtual ~stackBase() {}
 
    /** Place state element at the back of the stack */
    virtual void push_back(State const & t ) = 0; 

    
    /** Pop the state element at the top of the stack */
    virtual void pop_back() = 0;

    /** Returns stack size */
    virtual int size() = 0; 
    
    /** Allows access to a specific state at index idx */
    virtual State & at(int idx) = 0;

    /** Allows access to a specific state at index idx, 
	copies element idx to state r */
    virtual void getAt(State & r, int idx) = 0; 
    
    /** Returns reference at the head of the stack */
    virtual State & front() = 0;
    
    /** Returns reference at the tail of the stack */
    virtual State & back() = 0;
    
    /** Alters the offset of time indices used by checkpointing
	if necessary */
    virtual int getOrigin() = 0;

    /** Remove all elements in the stack */
    virtual void clear() {
      while (size() != 0) {
	pop_back();
      }
    }

  };
  


  /** Derived stackBase class, extending functionality for adaptive checkpointing */
  template <typename State>
  class adaptStackBase : public stackBase<State> {
  public:
    virtual void placeAt(int idx, State const & s) = 0;
    

  };



}


#endif
