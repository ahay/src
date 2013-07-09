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
#ifndef __RVL_ACP_SIM
#define __RVL_ACP_SIM

#include "sim.hh"
#include "revolve.hh"
#include <vector>

//#define DEBUG

namespace TSOpt {
	
  using RVL::Writeable;
  using RVL::RVLException;
  using RVLAlg::Algorithm;
  using RVLAlg::ListAlg;
  using RVLAlg::LoopAlg;
  using RVLAlg::StateAlg;
  using RVLAlg::Terminator;

  template<typename State,
	   typename myVector>
  class ACPSim: public Sim<State> {
	     
  private:
        
    myVector _ldclist;
    std::vector<int> _checklist;
    
    // keeps track if initial checkpoint schedule was already calculated
    bool _gotsched;   
    
    // Revolve object
    Revolve *_r;      

    // Revolve-related ints
    int _steps; 
    int _snaps; 
    int _info; 
    int _capo;
    int _check; 
    int _fine;


    // number of time steps
    int _nSteps;

    int _currStepCount;
    

  public:
	     
    ACPSim(TimeStep<State> & step,
	   TimeTerm & term,
	   Algorithm & initStep,
	   unsigned int nsnaps)
      : Sim<State>(step,term, initStep), _ldclist(),  
	_checklist(),  
	_r(NULL), _steps(0),
	_snaps(nsnaps), _info(3), _capo(0), _check(0), _fine(0),
	_currStepCount(0) { }
    
    ACPSim(ACPSim<State, myVector> const & s)
      : Sim<State>(s), _ldclist(s._ldclist),  
	_checklist(s._checklist), 
	_r(s._r), 
	_steps(s._steps), _snaps(s._snaps), _info(s._info), _capo(s._capo), 
	_check(s._check), _fine(s._fine), 
	_currStepCount(0) {}
    
    ~ACPSim() {

      if (_r) delete _r;
    
    }
    
    void run() {

      //      cout << "*******Calling Run!**********" << endl;

      // push initial values to state vectors
      if (_ldclist.size() == 0) {
	/*
	State * tmp;
	tmp = new State(this->getState());
	_ldclist.push_back(tmp);
	*/
	
	_ldclist.push_back(this->step.getState());
	_checklist.push_back(0);
	
      }
        
      //      cout << "STEPS = " << _steps << endl;

      // Forward Simulation
      while (this->term.getTargetTime() > this->step.getTime() ) {

	_capo = 0;
	_steps = 0;
	int oldcapo;
	enum ACTION::action whatodo;

	// Note: snaps set upon construction

	_fine =  _steps + _capo;
	_check = -1;  

	// Declare Revolve Object
	_r = new Revolve(_snaps, _info);
	
	
	do {
	  oldcapo = _capo;
	  
	  whatodo = _r->revolve(&_check, &_capo, &_fine, _snaps, &_info);

	    #ifdef DEBUG
	    cout << "whatodo = " << whatodo << " , "; 
	    #endif

	    if ((whatodo == ACTION::takeshot) && (_info > 1)){
	      #ifdef DEBUG
	      cout << " takeshot at " << setw(6) << _capo << " in CP " << _check << endl;
	      cout << "_checklist.size() = " << _checklist.size() << endl;	      
	      cout << "_ldclist.size() = "<< _ldclist.size() << endl;
              #endif
	      //_checklist.push_back(_capo);

	      if (_ldclist.size() < _snaps) {
		_ldclist.push_back(this->step.getState());
		_checklist.push_back(_currStepCount);
	      }
	      
	      else if (_ldclist.size() >= _snaps) {		
		_ldclist.placeAt(_check, this->step.getState());
		_checklist.at(_check) = _currStepCount;
	      }
	      
#ifdef DEBUG    
	      cout << "LDC List" << endl;
	      for (int i=0; i<_ldclist.size(); ++i) {
		_ldclist.at(i).write(cout);
	      }
	      
	      cout << "CheckList " << endl;
	      for (int i=0; i<_ldclist.size(); ++i) {
		cout << _checklist.at(i) << endl;
	      }
#endif

	    }
	       
#ifdef DEBUG	    
	    if (whatodo == ACTION::firsturn) { cout << "firsturn" << endl; }
	    if (whatodo == ACTION::terminate) {cout <<"terminate" << endl;}

	    if (whatodo == ACTION::restore) { 
	      cout << " restore at " << setw(7) << _capo << endl;  
	    }
		
	    if ( whatodo == ACTION::youturn ) {cout << "youturn" << endl; }
#endif  

	    if ( whatodo == ACTION::advance ) {
#ifdef DEBUG	     
	      cout << " advance to " << setw(7) << _capo << endl;   
#endif
	      // cout << "(this->getTime()).getint() = " << (this->getTime()).getint()  << endl;
	      //	      while((this->step.getTime()).getint() < _capo) 
	      while (_currStepCount < _capo)   {  
#ifdef DEBUG
		cout << "(F) running fwd sim" << endl;
#endif
		this->step.run(); 
#ifdef DEBUG
		this->step.getState().write(cout);
#endif
		++_steps;
		++_currStepCount;
		
		// cout << "thisTime = " << (this->step.getTime()).getint() << endl;
#ifdef DEBUG
		cout << "_currStepCount = " << _currStepCount << endl;
		cout << "_capo = " << _capo << endl;
#endif
		// <ME?> This is serious stuff -- targettime needs to have the *final*
		//  simulation time, not an intermediate 'final' time
		// -> TO DO: consider passing in tf to acpsim constructor. 
		if ( (this->step.getTime() > this->term.getTargetTime()) || (this->step.getTime() == this->term.getTargetTime())) {
		  _capo = _steps;
		  _fine = _steps+1;
		  _r->turn(_fine);
		  
		}
		
	
	      }
	    
	      /*
	      cout << "LDC List" << endl;
	      for (int i=0; i<_ldclist.size(); ++i) {
		  _ldclist.at(i).write(cout);
	      }
	      */

	    }
	    
	    if (whatodo == ACTION::error) {
	      RVLException e;
	      e << "Error: GriewankStateHistory::setTime()\n";
	      e << " irregular termination of revolve \n";
	      throw e;
	    }
	    


	} while (whatodo != ACTION::firsturn);
	
	

      }
      


      // Backward Simulation
      while ( (this->term.getTargetTime() < this->step.getTime()) ) {

	int oldcapo;
	enum ACTION::action whatodo;
#ifdef DEBUG
	cout << "----------------------" << endl;
	cout << "TARGET TIME!" << endl;
	this->term.getTargetTime().write(cout);
	cout << "GETTIME!" << endl;
	this->getTime().write(cout);
	cout << "----------------------" << endl;
#endif
	// Error check here, to see if fwd sim has completed?

    
	if (_steps > 2) {
	  
	  
	  do {
	  
	    oldcapo = _capo;
	    whatodo = _r->revolve(&_check, &_capo, &_fine, _snaps, &_info); 
	    //	    cout << "whatodo = " << whatodo << " , " << endl;

	    if ((whatodo == ACTION::takeshot) && (_info > 1)) {
              #ifdef DEBUG
	      cout << "*** takeshot at " << setw(6) << _capo << " in CP " << _check << endl;
              #endif

	      _ldclist.placeAt(_check, this->step.getState());
#ifdef DEBUG	     
	      cout << "_checklist.size() " << _checklist.size() << endl;
#endif
	      _checklist.at(_check) = _currStepCount;
	      
	    } // Matches if (whatodo == ACTION::takeshot)
	  
	    

	    if ((whatodo == ACTION::advance) && (_info > 2))  {  
              #ifdef DEBUG
	      cout << "*** advance to " << setw(7) << _capo << endl;
              #endif

	      // cout << endl << "begin comparison" << endl;
	      // cout << "_ldclist.at(_check).getTime.getint = " <<  (_ldclist.at(_check).getTime()).getint() << endl;
	      // cout << "_checklist.at(_check) = " << _checklist.at(_check) << endl;
	      // cout << "end comparison" << endl << endl;
	    
	            
	      int limit;
	      //limit = (_ldclist.at(_check).getTime()).getint() + static_cast<int>(_capo-oldcapo);
	      limit = _currStepCount + static_cast<int>(_capo-oldcapo);
#ifdef DEBUG
	      cout << "int limit = " << limit << endl;
#endif  
	      //while((this->step.getTime()).getint() < limit) 
	      while (_currStepCount < limit) 
		{  
#ifdef DEBUG
		  cout << "(B) running fwd sim" << endl;
#endif
		  this->step.run(); 
#ifdef DEBUG
		  this->step.getState().write(cout);
#endif
		  ++_currStepCount;
		}

	  
	     
	    }// Matches if (whatodo == ACTION::advance)
	

	    if ((whatodo == ACTION::restore) && (_info > 2)) {
              #ifdef DEBUG
	      cout << "*** restore at " << setw(7) << _capo << endl; 
	      cout << "LDClist size" << _ldclist.size() << endl;
              #endif

	      
	      State & r = this->step.getState();
	      r.copy(_ldclist.at(_check));

	      _currStepCount = _checklist.at(_check);
	      
	      //this->setState(_ldclist.at(_check));
	      
	    } // Matches if (whatodo == ACTION::restore)

	    
	    
	    if (whatodo == ACTION::error) {
	    
	      RVLException e;
	      e << "Error: Irregular execution of Revolve\n";
	      throw e;

	    } // matches if (whatodo == ACTION::error);

	    /*
	    cout << "LDC List" << endl;
	    for (int i=0; i<_ldclist.size(); ++i) {
	      _ldclist.at(i).write(cout);
	    }
	    */
		
	  } while( (whatodo != ACTION::terminate) && (whatodo != ACTION::youturn) );

	}

	
	// If number of time steps is less than 2
	else {
	
	  _capo=0; _check=0; _info=0; _fine=0; _snaps=1; 
	  oldcapo=_capo; 
	  whatodo = static_cast<ACTION::action>(0);
	  
	  State & r = this->step.getState();
	  r.copy(_ldclist.at(_check));
	  
	  //this->setState(_ldclist.at(_check));
	  
	}

	
	
	
      }
      
    
    }



    
  };
  


}







#endif
