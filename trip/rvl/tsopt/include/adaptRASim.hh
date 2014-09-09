#ifndef ADAPT_RASIM_HH
#define ADAPT_RASIM_HH

#include "alg.hh"
#include "tterm.hh"

#include <stack>
#include "dtterm.hh"
#include "cdt.hh"

#include "rt.hh"

namespace TSOpt {

  /** Remember-All Simulator for adaptive time stepping. Capable of 
      interpolating state data. As usual, myStack handles how states 
      are saved */
  template<typename State, 
	   typename StateInterp >
  class adaptRASim: public Sim<State> {

  private:

    adaptRASim();

    stackBase<State> & stateHist;

    bool filledSH;
 

  public: 

    adaptRASim(PeekTimeStep<State> & step,
	       TimeTerm & term,
	       Algorithm & initStep,
	       stackBase<State> & stateHist_)
      : Sim<State>(step,term, initStep), filledSH(false),  stateHist(stateHist_)
    { }
    
    
    adaptRASim(adaptRASim<State, StateInterp> const & s)
      : Sim<State>(s), filledSH(s.filledSH), stateHist(s.stateHist)
    {}
    
    virtual ~adaptRASim() {
      while (stateHist.size() > 0) {
	stateHist.pop_back();
      }
    }

    
    void fillSH() {
      // run initialization step
      this->initstep.run();

      // push initial state to stack
      stateHist.push_back(this->step.getState());


      cout << "step.time = " << endl;
      this->step.getTime().write(cout);
      cout << "Term.ttime = " << endl;
      this->term.getTargetTime().write(cout);

      
      while ( this->step.getTime() < this->term.getTargetTime() ) {

	// pass target time to ensure that no over-stepping takes place
	dynamic_cast<PeekTimeStep<State>  &>(this->step).peekTT(this->term.getTargetTime());
	this->step.run();
	stateHist.push_back(this->step.getState());
	cout << "stateHist.size() = " << stateHist.size() << endl;
      }
      
      filledSH = true;

     }

    /** to reset the adaptRASim, the stack has to be cleared first! 
	Function fillSH() takes care of the initstep.run() part */
    void resetSim() {
      cout << "Calling aRA resetSim()" << endl;
      stateHist.clear();
      filledSH = false;
      
    }

    
    void run() {
  
      if (!filledSH) {  fillSH(); }

      while (( this->step.getTime() > this->term.getTargetTime() )) {
	
	try {
	  this->step.setTime(this->term.getTargetTime());		  
	  State & r = this->step.getState();

	  // interpolate forward states in stateHist	  
	  StateInterp::getInterpState(r, this->term.getTargetTime(), stateHist);



	  if (dynamic_cast<RealTime const &>(r.getTime()).getreal() == 0.0) {
	    resetSim();
	  }
	
	 
	}
	catch (RVLException & e) {
	  e<<"\ncalled from adaptRASim::run\n";
	  throw e;
	}
      }

	  State & r = this->step.getState();
          double rt = dynamic_cast<RealTime const &>(r.getTime()).getreal();
	  cout << "aRA rt = " << rt << endl;      

    }

  };



 
  
}


#endif
