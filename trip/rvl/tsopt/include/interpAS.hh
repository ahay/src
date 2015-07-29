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
#ifndef __RVL_INTERP_AS
#define __RVL_INTERP_AS

#include <list>
#include <deque>
#include <iterator>
#include "t.hh"
#include "rns.hh"

namespace TSOpt {

  using RVL::Writeable;
  using RVL::RVLException;

  /** "Moving Buffer" class, for use with meACPsim (adaptive checkpointing). 
      Uses a deque to implement the "moving buffer". Templated on a State
      type and an Interpolator type
  */
  template <typename State, typename Interpolator>
  class interpAS {
  private:
    int _numNodes;
    //std::list< State* > _interpBuffer;
    std::deque<State> _interpBuffer;

  public: 
    /** default constructor, only need a int to determine the size
	of the moving buffer. Typically, this should be (n+1), where
	n is the order of the time-stepping scheme being used
    */
    interpAS(int numNodes) : _numNodes(numNodes), _interpBuffer()   {
      if (numNodes < 2) {
	RVLException e;
	e << "Error in interpAS constructor. numNodes must be at least 2 \n";
	throw e;
      }

    }

    /** copy constructor */
    interpAS(interpAS const & s) : _numNodes(s._numNodes),
				   _interpBuffer(s._interpBuffer) {}

    /** destructor */
    ~interpAS() { 
      for (int i=0; i<_interpBuffer.size(); ++i) {
	//delete _interpBuffer.back();
	_interpBuffer.pop_back();
      }
      
  
    }

    /** Checks to see if the time is in the limits of the 
        interpolation buffer (meaning, we can interpolate
        if this is true) 
    */
    bool inInterpN(const Time & t) {

      if (! isBufferFull() ) {
	return false;
      }

      // Assumes that _interpBuffer is always ordered from least 
      //   to greatest (front to back), in terms of time
      if ( notBelowInterpN(t)  && notAboveInterpN(t) ) {
	return true;
      }

      return false;

    }


    bool notAboveInterpN(const Time & t) {
      #ifdef DEBUG
      cout << "notAboveInterpN" << endl;
      cout << "t = " << endl;
      t.write(cout);
      cout << "backTime = " << endl;
      _interpBuffer.back().getTime().write(cout);
      #endif
      if ((t < _interpBuffer.back().getTime()) || (t == _interpBuffer.back().getTime())) {       
	return true;
      }

      return false;
    }

    bool notBelowInterpN(const Time & t) {
      #ifdef DEBUG
      cout << "notBelowInterpN" << endl;
      cout << "t = " << endl;
      t.write(cout);
      cout << "frontTime = " << endl;
      _interpBuffer.front().getTime().write(cout);
      #endif
      if ( (t > _interpBuffer.front().getTime()) || (t == _interpBuffer.front().getTime())  ) {
	return true;
      }

      return false;

    }

    bool isBufferFull() { return (_interpBuffer.size() == _numNodes); }
    int getNumNodes() const {return _numNodes;}

    /** Necessary if the interpolation buffer allotted is greater than
        the total number of fwd steps taken
    */
    void overrideNumNodes(int numNodes) {
      cout << "Warning: overriding number of nodes in interpolator" << endl;
      cout << "from " << _numNodes << " to " << numNodes << endl;
      _numNodes = numNodes;
    }
    
    /** Update for moving buffer, while performing a forward evolution  */
    void updateInterpBuffer_F(State & s) {
      
      // fill buffer first
      if (_interpBuffer.size() < _numNodes) {
	_interpBuffer.push_back(s);
      }
      
      // else update buffer by deleting the latest state
      else {
	_interpBuffer.pop_front();
	_interpBuffer.push_back(s);

      }

      cout << "!!!! interpBuffer.size() = " << _interpBuffer.size() << endl;
      //typename list<State*>::iterator iter;
      typename deque<State>::iterator iter;
      for (iter=_interpBuffer.begin(); iter!=_interpBuffer.end() ; iter++) {
	(*iter).write(cout);
      }

    }

    /** Update for a moving buffer, while performing a backward evolution */
   void updateInterpBuffer_A(State & s) {
      
      // fill buffer first
      if (_interpBuffer.size() < _numNodes) {
	_interpBuffer.push_front(s);
      }
      
      // else update buffer by deleting the latest state
      else {
	_interpBuffer.pop_back();
	_interpBuffer.push_front(s);

      }

      cout << "!!!! interpBuffer.size() = " << _interpBuffer.size() << endl;
      //typename list<State*>::iterator iter;
      typename deque<State>::iterator iter;
      for (iter=_interpBuffer.begin(); iter!=_interpBuffer.end() ; iter++) {
	(*iter).write(cout);
      }

    }


    /** Update for moving buffer, asking to put the state in the middle
	of the buffer, effectively "refining" the interpolation span. 
	Should be called in the backward traversal, once Revolve
	calls "takeshot" */
    void updateInterpBuffer_M(State & s, const Time & t) {
      typename deque<State>::iterator iter;
      bool insert = false;

      for (iter=_interpBuffer.begin(); iter<(_interpBuffer.end()-1) ; iter++) {
	
	if( ( (*iter).getTime() < s.getTime() )  && 
	    (s.getTime() < (*(iter+1)).getTime()) ) {
	  cout << "found candidate for interpbuf_M!" << endl;

	  (*iter).getTime().write(cout);
	  s.getTime().write(cout);
	  (*(iter+1)).getTime().write(cout);


	  _interpBuffer.insert(iter+1,s);
	  //_interpBuffer.pop_front();
	  insert = true;
	}
	
	cout << "current list of times" << endl;
	showBufferTimes();
	  
	cout << "target time = " << endl;
	t.write(cout);
	
	   
	int blah;
	cin >> blah;

	if (insert) {
	  int counter = 1;
	  for (iter=_interpBuffer.begin(); iter<(_interpBuffer.end()-1) ; iter++) {
	   
 
	  if( ((*iter).getTime() < t)   &&  (t < (*(iter+1)).getTime()) ) {
	      
	    cout << "counter = " << counter << endl;
	    if( counter < static_cast<int>(_interpBuffer.size()/2) ) {
	      _interpBuffer.pop_back();
	    }
	    else  {
	      _interpBuffer.pop_front();
	    }
	    
	    cout << "current list of times, now" << endl;
	    showBufferTimes();
	 
	  }
	  
	  counter++;
	    
	  }	
	}
	


      }

    }


    /** Interpolates data buffers at requested time and returns 
        the corresponding state. Makes use of the interface in the 
	template argument Interpolator
    */
    void interpBuffer(const Time & t, State & s) {

      if (!isBufferFull()) {
 	RVLException e;
	e << "interpAS Error: buffer size != numNodes\n";
	e << "buffer size = " << _interpBuffer.size();
	e << "\nnumNodes = " << _numNodes;
	throw e;
      }

      if (!inInterpN(t)) {
	RVLException e;
	e << "interpAS Error: requested time not in buffer neighborhood\n";
	e << "Cannot perform interpolation\n";
	e << "Requested Time: \n";
	t.write(e);
	throw e;
      }


      Interpolator::getInterpState(s, t, _interpBuffer);
  
      s.setTime(t);
      
      cout << "'Interpolated' State" << endl;
      s.write(cout);

      
    }



    /** A debug feature: shows times in the buffer */
    void showBufferTimes() {
      cout << "List of Times in Interpolation Buffer: " << endl;
      typename deque<State>::iterator iter;
      for (iter=_interpBuffer.begin(); iter!=_interpBuffer.end() ; iter++) {
	(*iter).getTime().write(cout);
      }

    }

    
  };

}
#endif
