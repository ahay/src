
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
#ifndef __RVL_MEACP_SIM
#define __RVL_MEACP_SIM

#include "sim.hh"
#include "revolve.hh"
#include "interpAS.hh"
#include "stackBase.hh"
#include <vector>

 #define DEBUG

namespace TSOpt {
	
  using RVL::Writeable;
  using RVL::RVLException;
  using RVLAlg::Algorithm;
  using RVLAlg::ListAlg;
  using RVLAlg::LoopAlg;
  using RVLAlg::StateAlg;
  using RVLAlg::Terminator;



  template<typename State,
	   typename Interpolator>
  class meACPSim: public Sim<State> {
	     
  private:
        
    adaptStackBase<State> &  _ldclist;
    std::vector<double> _dtlist;
    
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
    int _currStepCount;

    // moving interpolating window
    interpAS<State, Interpolator> _terpy;

    // forward/backward mode flag
    bool _forwardMode;

    // lock for revolve
    bool _revolveLock;

    // for revolve to remembe actions
    enum ACTION::action _whatodo;

    // allows user to enable/disable forward evolution in
    //  the backward traversal. (Forward traversal in 
    //  the backward mode is not necessary due to the 
    //  interpolation buffer. However, it could 
    //  be used to refine the interpolation buffer)
    bool _refineInterpBuf;

    bool _testSwitch;

  public:
	     
    meACPSim(TimeStep<State> & step,
	     TimeTerm & term,
	     Algorithm & initStep,
	     unsigned int nsnaps, 	  
	     int bufSize, 
	     adaptStackBase<State> & ldclist,
	     bool refineInterpBuf = false)
      : Sim<State>(step,term, initStep),   
	_dtlist(), 
	_r(NULL), _steps(0),
	_snaps(nsnaps), _info(3), _capo(0), _check(0), _fine(0), 
	_currStepCount(0),
	_terpy(bufSize), _ldclist(ldclist), _forwardMode(true), 
	_revolveLock(false), _refineInterpBuf(refineInterpBuf), _testSwitch(false)
    { 
    }
    
    meACPSim(meACPSim<State, Interpolator> const & s)
      : Sim<State>(s), _ldclist(s._ldclist),  
	_dtlist(s._dtlist), 
	_r(s._r), 
	_steps(s._steps), _snaps(s._snaps), _info(s._info), _capo(s._capo), 
	_check(s._check), _fine(s._fine), _currStepCount(s._currStepCount), _terpy(s._terpy),
	_ldclist(s._ldclist), _forwardMode(s._forwardMode), _revolveLock(s._revolveLock),
	_refineInterpBuf(s._refineInterpBuf), _testSwitch(s._testSwitch)
      
    {}
    
    ~meACPSim() {
      if (_r) delete _r;
    }
    

    void run() {

      cout << endl << "CurrTime: " << endl;
      this->getTime().write(cout);

      cout << endl << "TT:" << endl;
      this->term.getTargetTime().write(cout);
      

      // Forward Mode
      if (_forwardMode) {

	cout << "*******FwdSim!**********" << endl;

	_capo = 0;	
	_steps = 0;
	
	int oldcapo;
	// enum ACTION::action whatodo;


	// Note: snaps set upon construction

	_fine =  _steps + _capo;
	_check = -1;  

	// Declare Revolve Object
	_r = new Revolve(_snaps, _info);
	
	
	
	do {
	  oldcapo = _capo;
	  
	  _whatodo = _r->revolve(&_check, &_capo, &_fine, _snaps, &_info);

	    cout << "ch ca fi sn in = " << _check << " " << _capo << " " << _fine << " " << _snaps << " " << _info << endl;

	    #ifdef DEBUG
	    cout << "whatodo = " << _whatodo << " , "; 
	    #endif

	    if ((_whatodo == ACTION::takeshot) && (_info > 1)){
	      #ifdef DEBUG
	      cout << " takeshot at " << setw(6) << _capo << " in CP " << _check << endl;
	      #endif
	      

	      if (_ldclist.size() < _snaps) {
		_ldclist.push_back(this->step.getState());
		_dtlist.push_back(dynamic_cast<PeekTimeStep<State>  &>(this->step).getDt());
	      }
	      else {
		_ldclist.placeAt(_check, this->step.getState());
		_dtlist.at(_check) = 	dynamic_cast<PeekTimeStep<State>  &>(this->step).getDt();
	      }
	      
	      /*
	      cout << "LDC List" << endl;
	      for (int i=0; i<_ldclist.size(); ++i) {
		_ldclist.at(i).write(cout);
	      }
	      */

	      cout << "timeList " << endl;
	      for (int i=0; i<_ldclist.size(); ++i) {
		cout << _ldclist.at(i).getTime().getreal() << endl;
	      }

	    }
	       
	    if (_whatodo == ACTION::firsturn) { cout << "firsturn" << endl; }
	    if (_whatodo == ACTION::terminate) {cout << "terminate" << endl;}
	    
	    if (_whatodo == ACTION::restore) { cout << " restore at " << setw(7) << _capo << endl;  }
	    if (_whatodo == ACTION::youturn ) {cout << "youturn" << endl;   }


	    if (_whatodo == ACTION::advance ) {
	     
	      cout << " advance to " << setw(7) << _capo << endl;   
	
	
	      dynamic_cast<PeekTimeStep<State>  &>(this->step).peekTT(this->term.getTargetTime());
	      this->step.run(); 
	      this->step.getState().write(cout);

	      // update buffer here
	      _terpy.updateInterpBuffer_F(this->step.getState());

	      ++_steps;
	      ++_currStepCount;
	      
	      cout << "steps = " << _steps << " and currStepCount = " << _currStepCount << endl;
	      
	      if ( (this->step.getTime() == this->term.getTargetTime()) ) {
		// push last state, at t = T

		/*
		  if (_ldclist.size() < _snaps) {
		  _ldclist.push_back(this->step.getState());
		  _dtlist.push_back(dynamic_cast<PeekTimeStep<State>  &>(this->step).getDt());
		  }
		  else {
		  _ldclist.placeAt(_check, this->step.getState());
		  _dtlist.at(_check) = 	dynamic_cast<PeekTimeStep<State>  &>(this->step).getDt();
		  }
		*/
		
		/*
		  _capo = _steps-1;
		  _fine = _steps;
		  _r->turn(_fine);
		*/
		
		_fine = _capo+1;
		_r->turn(_fine);

		_forwardMode = false;
	      }


	    }
	    
	    if (_whatodo == ACTION::error) {
	      RVLException e;
	      e << "Error: GriewankStateHistory::setTime()\n";
	      e << " irregular termination of revolve \n";
	      throw e;
	    }
	    

	    if (_whatodo == ACTION::firsturn) {
	      cout << "curr step count =" << _currStepCount << endl;
	    }
	    


	} while (_whatodo != ACTION::firsturn);
	
	

      }
      
       

      // Backward Simulation
      else {

       cout << "~~~~~~BwdSim!~~~~~~" << endl;
     
       cout << "curr step count =" << _currStepCount << endl;


	int oldcapo;
	//	enum ACTION::action whatodo;
	cout << "----------------------" << endl;
	cout << "TARGET TIME!" << endl;
	this->term.getTargetTime().write(cout);
	cout << "GETTIME!" << endl;
	this->step.getTime().write(cout);
	cout << "----------------------" << endl;

      
	if (_fine > 0) {
	  
	  
	  do {

	    cout << "ch ca fi sn in = " << _check << " " << _capo << " " << _fine << " " << _snaps << " " << _info << endl;


	    oldcapo = _capo;


	    if (!_revolveLock) {
	      _whatodo = _r->revolve(&_check, &_capo, &_fine, _snaps, &_info); 
	    }
	    
	    cout << "whatodo = " << _whatodo << " , " << endl;

	    #ifdef DEBUG
	    cout << "LDC List Times:" << endl;
	    for (int i=0; i<_ldclist.size(); ++i) {
	      _ldclist.at(i).getTime().write(cout);
	    }
	    #endif



	    if ((_whatodo == ACTION::takeshot) && (_info > 1)) {
              #ifdef DEBUG
	      cout << "*** takeshot at " << setw(6) << _capo << " in CP " << _check << endl;
              #endif

	       _ldclist.placeAt(_check,this->step.getState());
	       _dtlist.at(_check) = dynamic_cast<PeekTimeStep<State>  &>(this->step).getDt();
	       
	       _testSwitch = true;
    
	    } // Matches if (whatodo == ACTION::takeshot)
	  
	    

	    if ((_whatodo == ACTION::advance) && (_info > 2))  {  
              #ifdef DEBUG
	      cout << "*** advance to " << setw(7) << _capo << endl;
              #endif
	      
	      if (_refineInterpBuf) {
	      for (int i=oldcapo; i<_capo; ++i){
		// prevents wasteful forward time-stepping
		if (this->step.getTime() < this->term.getTargetTime()) {
	
		cout << "In Backward Mode [Advance]. Running fwd sim" << endl;
		this->step.run(); 
		this->step.getState().write(cout);
		++_currStepCount;
		
		}
	     	
		_testSwitch = true;
	      }
	      

	      }
		

	    } // Matches if (whatodo == ACTION::advance)
	

	    if ((_whatodo == ACTION::restore) && (_info > 2)) {
              #ifdef DEBUG
	      cout << endl << "*** restore at " << setw(7) << _capo << endl; 
              #endif

	      
	      cout << "LDClist size = " << _ldclist.size() << endl;
	      cout << "dtlist size = " << _dtlist.size() << endl;

	      State & r = this->step.getState();

	      r.copy(_ldclist.at(_check));
	      this->step.setTime(_ldclist.at(_check).getTime());	      
	      dynamic_cast<PeekTimeStep<State>  &>(this->step).setDt(_dtlist.at(_check));
	    
	      cout << endl << "restored state : " << endl;
	      r.write(cout);

	      cout << endl << "in restore, TT = " << endl;
	      this->term.getTargetTime().write(cout);


	      _testSwitch = true;

	    } // Matches if (whatodo == ACTION::restore)

	    
	    
	    if (_whatodo == ACTION::error) {
	    
	      RVLException e;
	      e << "Error: Irregular execution of Revolve\n";
	      throw e;

	    } // matches if (whatodo == ACTION::error);

	    
	    if (_whatodo == ACTION::youturn) {
	      cout << endl << " *** youturn *** " << endl;

	      _terpy.showBufferTimes();

	      cout << "TargetTime = " << endl;
	      this->term.getTargetTime().write(cout);

	      if (_testSwitch) {
		cout << "current adjsim state" << endl;
		this->step.getState().write(cout);
	      }

	      // lock calling of revolve while targettime is in span of buffer
	      if ( (_terpy.inInterpN(this->term.getTargetTime()))  ) {
		cout << "In interp buffer, locking revolve call" << endl;
		_revolveLock = true;

		State & r = this->step.getState();
		_terpy.interpBuffer(this->term.getTargetTime(),r);
		_testSwitch = false;
	
	      }
	      else {
	      	cout << "unlocking revolve call" << endl;
	      	_revolveLock = false;

		if (_testSwitch) {
		  State & r = this->step.getState();
		  _terpy.updateInterpBuffer_A(r);
		  _testSwitch = false;		  
		  _revolveLock = true; 
		}
		
	      }		
		
	    }
	  


	  } 
	  while( ( !(this->term.getTargetTime() == this->step.getTime())  || 
		   (_whatodo != ACTION::youturn)) && (_whatodo != ACTION::terminate) 	); 
	 
	}

	
	// If fine = 0, Revolve is terminating 
	else {
      	  _capo=0; _check=0; _info=0; _fine=0; _snaps=1; 
	  oldcapo=_capo; 
	  _whatodo = static_cast<ACTION::action>(0);
	  
	  cout << "fine = 0, teminating backward mode" << endl;

	  // 0th state is always the first one (t=t0)
	  _terpy.updateInterpBuffer_A(_ldclist.at(0));
	  State & r = this->step.getState();

	  _terpy.showBufferTimes();


	  _terpy.interpBuffer(this->term.getTargetTime(),r);
	  
	  resetSim();

	}

	
     }


      //  Reset params when at the end of adj evolution
      //  <ME?> Take "end" as an input to constructor later
      RealTime end;
      end = 0.0;
      if (this->step.getTime() == end) {	
	cout << "Reached end of Backward Traversal" << endl;
	_terpy.showBufferTimes();
	resetSim();
      }
      
    }

    
    void resetSim() {
      cout << "Resetting meACPSim" << endl;
      _forwardMode = true;
      _revolveLock = false;
      _testSwitch = false;
      _dtlist.clear();
      _ldclist.clear();
      _currStepCount = 0;
    }
    
  };
  


}







#endif
