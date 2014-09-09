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
#ifndef __RVL_CP_SIM
#define __RVL_CP_SIM

#include "sim.hh"
#include "revolve.hh"
#include <vector>

//#define DEBUG
//#define CP_INFO
//#define ITER_INFO

namespace TSOpt {
	
  using RVL::Writeable;
  using RVL::RVLException;
  using RVLAlg::Algorithm;
  using RVLAlg::ListAlg;
  using RVLAlg::LoopAlg;
  using RVLAlg::StateAlg;
  using RVLAlg::Terminator;

  /**  
      CPSim -- Checkpointing Sim Class, which is a combination of StdSim (save none)
      and RASim (save all) strategies for handling the reference simulation states. 
      This class uses Griewank's optimal checkpointing schedule, allowing users
      to save only a subset of the reference simulation states for adjoint computations. 
      This sim class requires two unique inputs: the size of the buffer for checkpointing,
      and also a stack object, whose base class is defined in stackbase.hh. Also, 
      the class is templated on a state type and a discrete time type. 

      Note that the source file defines several debugging / info flags
      which produce varying levels of descriptive output. By default,
      all of these are turned off.
  */
  
  template<typename State,
	   typename myDiscTime>
  class CPSim: public Sim<State> {
    
  private:

    stackBase<State> & _ldclist;
    std::vector<myDiscTime> _checklist;

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

  public:

    /** standard constructor */
    CPSim(TimeStep<State> & step,
	  TimeTerm & term,
	  Algorithm & initStep,
	  unsigned int nsteps,
	  unsigned int nsnaps,
	  stackBase<State> & ldclist)
      : Sim<State>(step,term, initStep), _ldclist(ldclist),  
	_checklist(),  
	_gotsched(false), _r(NULL),
	_steps(nsteps), _snaps(nsnaps), _info(3), _capo(0), _check(0), _fine(0) { 	
      if (nsnaps <= 0 || nsnaps > 500) {
	RVLException e;
	e << "Error: nsnaps either negative, or exceeds maximum value (500) \n";
	e << "Possible explanation: nsnaps not initialized \n";
	throw e;
      }
    }
    
    /** copy constructor */
    CPSim(CPSim<State, myDiscTime> const & s)
      : Sim<State>(s), _ldclist(s._ldclist),  
	_checklist(s._checklist), 
	_gotsched(s._gotsched),  
	_r(s._r),
	_steps(s._steps), _snaps(s._snaps), _info(s._info), _capo(s._capo), 
	_check(s._check), _fine(s._fine)  {
      
      if (_snaps <= 0 || _snaps > 500) {
	RVLException e;
	e << "Error: nsnaps either negative, or exceeds maximum value (500) \n"; 
	e << "Possible explanation: nsnaps not initialized \n";
	throw e;
      }
    }
    
    /** use arg struct - for generic policy */
    CPSim(StdSimData<State> const & d)
      : Sim<State>(d.getStep(), d.getTerm(), d.getInit()), 
	_ldclist(*(d.getStack())), _checklist(), _gotsched(false), _r(NULL), 
	_steps(d.getNsteps()), _snaps(d.getNsnaps()), _info(3), 
	_capo(0), _check(0), _fine(0) {}

    /** destructor */
    ~CPSim() {
      if (_r) delete _r;
    }
    
    /** run() is split into two phases: the forward phase and the backward phase. 
	During the forward phase, the algorithm gets the checkpoint list from 
	Griewank's algorithm, and saves the appropriate reference states to the
	buffer. During the backward phase, Griewank's algorithm dictates what 
	is to be done with the saved states (restore, use as a starting point
	for evolution, etc.) */
    void run() {

      try {

#ifdef ITER_INFO
	  cerr << "--------------- call CPSIM::run ---------------" << endl;
	  cerr << "TARGET TIME = ";
	  this->getTargetTime().write(cerr);
	  cerr << "CURRENT TIME = ";
	  this->getTime().write(cerr);
	  if (this->getTargetTime() > this->getTime()) 
	    cerr<<"--- case advance\n";
	  else if (this->getTargetTime() < this->getTime()) 
	    cerr<<"--- case restore\n";
	  else 
	    cerr<<"--- case you're there - no-op\n";
	  cerr << "-----------------------------------------------" << endl;
#endif

	// push initial values to state vectors
	if (_ldclist.size() == 0) {	  
#ifdef CP_INFO
	  cerr<<"*** checkpoint loop: record initial data\n";
#endif
	  _ldclist.push_back(this->step.getState());
	}
	
	// Run Revolve to get the initial checkpoint list 
	if (!_gotsched) { this->getSched(); } 
	
	// Forward Simulation
	while ( (this->term.getTargetTime() > this->step.getTime() )  )	 {

#ifdef CP_INFO
	  cerr<<"*** checkpoint loop\n";
#endif
	  this->step.run();
	  
	  try {
	    const myDiscTime & myDT  = dynamic_cast<const myDiscTime &>(this->step.getTime());
#ifdef DEBUG
	    cerr<<"CPSim: print myDT\n";
	    myDT.write(cerr);
#endif
	    if (timeToCheck(myDT) ) {
#ifdef CP_INFO
	      cerr<<"*** checkpoint loop: push back state, top = "<<_ldclist.size()<<"\n*** time = ";
	      myDT.write(cerr);
#endif
	      _ldclist.push_back(this->step.getState());
	    }
#ifdef ITER_INFO
	  cerr << "--------------- CPSIM::run checkpoint loop bottom ----" << endl;
	  cerr << "TARGET TIME = ";
	  this->getTargetTime().write(cerr);
	  cerr << "CURRENT TIME = ";
	  this->getTime().write(cerr);
	  if (this->getTargetTime() > this->getTime()) 
	    cerr<<"--- run again\n";
	  else 
	    cerr<<"--- finished\n";
	  cerr << "----------------------------------------------------" << endl;
#endif

#ifdef CP_INFO
	    cerr<<"*** checkpoint loop: run again, top = "<<_ldclist.size()<<"\n";
#endif
	  }
	  catch (bad_cast) {
	    cerr << "Dynamic cast from step.getTime() to myDiscTime failed!" << endl;
	    cerr << "Probable error: state time type is not equal to myDiscTime" << endl;
	  }
	}

	// Backward Simulation
	while ( (this->term.getTargetTime() < this->step.getTime()) ) {

#ifdef CP_INFO
	  cerr<<"*** revolve loop\n";
#endif
	  // Error check here, to see if fwd sim has completed?
	  
	  int oldcapo; 
	  enum ACTION::action whatodo;
	  
	  if (_steps > 2) {
	    do {
	      oldcapo = _capo;
	      whatodo = _r->revolve(&_check, &_capo, &_fine, _snaps, &_info); 
	      if ((whatodo == ACTION::takeshot) && (_info > 1)) {
#ifdef CP_INFO
		cerr << "*** takeshot at " << setw(6) << _capo << " in CP " << _check << endl;
#endif
		while (_check < _ldclist.size()){
		  _ldclist.pop_back();
		}
		
		_ldclist.push_back(this->step.getState());	
#ifdef CP_INFO
		cerr<<"*** pushed back buffer "<<_ldclist.size()<<endl;
#endif
	      } // Matches if (whatodo == ACTION::takeshot)
	      
	      if ((whatodo == ACTION::advance) && (_info > 2))  {  
#ifdef CP_INFO
		cerr << "*** advance to " << setw(7) << _capo << endl;
#endif
		const myDiscTime & myDT  = dynamic_cast<const myDiscTime &>(this->step.getTime());
		myDiscTime targetDT(myDT); 
		for (int i=oldcapo;i<_capo;++i) ++targetDT;
       		while(myDT < targetDT) {		 
		  this->step.run(); 
		}
	      } // Matches if (whatodo == ACTION::advance)
	      
	      if ((whatodo == ACTION::restore) && (_info > 2)) {
#ifdef CP_INFO
		cerr << "*** restore at " << setw(7) << _capo << endl; 
		cerr << "*** check = "<<_check<<"\n";
#endif
		State & r = this->getState();
		_ldclist.getAt(r,_check);
	      } // Matches if (whatodo == ACTION::restore)
	      
	      if (whatodo == ACTION::error) {
		RVLException e;
		e << "Error: Irregular execution of Revolve\n";
		throw e;
	      } // matches if (whatodo == ACTION::error);
	      
	      if (whatodo == ACTION::terminate) {
		#ifdef CP_INFO
		cerr << endl << "called TERMINATE" << endl;
		#endif
	      }

	      if (whatodo == ACTION::youturn) {
		#ifdef CP_INFO
		cerr << endl << "called YOUTURN" << endl;
		#endif
	      }
	    } while( (whatodo != ACTION::terminate) && (whatodo != ACTION::youturn) );
	  }
	
	  // If number of time steps is less than 2
	  else {
	    _capo=0; _check=0; _info=0; _fine=0; _snaps=1; 
	    oldcapo=_capo; 
	    whatodo = static_cast<ACTION::action>(0);
	    State & r = this->getState();
	    _ldclist.getAt(r,_check);
	  }
	}

#ifdef DEBUG
	cerr<<"*** goodbye\n";
#endif
      }
      
      catch (RVLException & e) {
	e<<"\ncalled from CPSim::run.\n";
	throw e;
      }
   }

    /** Resets the sim, so that the CP algorithm may be run all over again */
    void resetSim() {
      _gotsched = false;
      _ldclist.clear();
    }
  
    /** Checks to see if the current discrete time is in Griewank's checkpointing
	list. Used as a criteria to store simulaton state
    */
    bool timeToCheck(const myDiscTime & myDT) {
	for (int k=0; k< (int)_checklist.size(); ++k){
	if (_checklist[k] == myDT )  {
#ifdef CP_INFO
	  cerr<<"timeToCheck: found time match k = "<<k<<endl;
	  cerr<<"----> time = ";
	  myDT.write(cerr);
#endif
	  return true;
	}
      }
      return false;
    }
    
    /** Gets checkpointing schedule from Griewank's algorithm, Revolve, and saves them
	in a vector of ints
    */
    void getSched() {

      enum ACTION::action whatodo;
      _capo = 0;
      _info = 3;
      
      // Note: snaps and steps are set upon construction
      
      _fine = _steps + _capo;
      _check = -1;  
      
      // Declare Revolve Object
      _r = new Revolve(_steps, _snaps, _info);
      
      // Clear checkpointed times
      _checklist.clear();
      
      do {
	whatodo = _r->revolve(&_check, &_capo, &_fine, _snaps, &_info);
	
#ifdef CP_INFO
	cerr << "whatodo = " << whatodo << " , "; 
#endif
	
	// only care if I'm saving (time) checkpoints
	
	if ((whatodo == ACTION::takeshot) && (_info > 1)){
#ifdef CP_INFO
	  cerr << " takeshot at " << setw(6) << _capo << " in CP " << _check << endl;
#endif
	  int offset = _ldclist.getOrigin();
	  
	  myDiscTime capoDT;
      	  capoDT = _capo + offset;
	  _checklist.push_back(capoDT);
#ifdef CP_INFO
	  cerr<<"push back: capo = "<<_capo<<" offset = "<<offset<<" capoDT = ";
	  capoDT.write(cerr);
#endif
	}

	if (whatodo == ACTION::error) {
	  RVLException e;
	  e << "Error: GriewankStateHistory::setTime()\n";
	  e << " irregular termination of revolve \n";
	  throw e;
	}
	
#ifdef CP_INFO
	if (whatodo == ACTION::firsturn)  { cerr << "firsturn" << endl; }
	if (whatodo == ACTION::terminate) { cerr <<"terminate" << endl;}
	if (whatodo == ACTION::restore)   { cerr << " restore at " << setw(7) << _capo << endl;  }
	if (whatodo == ACTION::youturn )  { cerr << "youturn" << endl; }
	if (whatodo == ACTION::advance )  { cerr << " advance to " << setw(7) << _capo << endl;   } 
#endif
      }  while (whatodo != ACTION::firsturn) ;
      _gotsched = true;
    }
  };
}

#endif
