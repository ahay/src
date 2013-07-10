
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
#ifndef __RVL_JET
#define __RVL_JET

#include "sim.hh"

namespace TSOpt {

  using RVL::RVLException;

  /** Base class for Simulator 1-Jet. Three Sims provide forward,
      linearized, and adjoint linearized simulation. Putting them
      together in a single struct is an assertion that they are
      appropriately related. */

  template<typename State>
  class Sim1Jet {
    
  public:

    Sim1Jet() {}
    Sim1Jet(Sim1Jet<State> const &) {}
    virtual ~Sim1Jet() {}

    /** returns ref to simulator */
    virtual Sim<State> & getFwd() const = 0;
    /** returns ref to linearized simulator */
    virtual Sim<State> & getLin() const = 0;
    /** returns ref to adjoint linearized simulator */
    virtual Sim<State> & getAdj() const = 0;
  };

  /** Standard realization of Sim1Jet. Note that the reference
      simulator and the TimeStep and TimeTerm defining the lin and adj
      sim s are initalized outside and referenced.

      Construction of linearized and adjoint linearized simulators:
      <ol> <li> reference simulator and linearized step combined in
      TargetCurrTimeStep; its run method advance reference simulator
      to time of linearized step.</li> <li> linearized step and
      TargetCurrTimeStep combined in TimeStepList, with latter as
      postop.</li> <li> TimeStepList and TimeTerm combined in
      StdSim.</li> <li> adjoint simulator built similarly.</li> </ol>
  */
  template<typename State, typename RefState>
  class StdSim1Jet: public Sim1Jet<State> {

  private:

    Sim<State> & fwd;
    Algorithm &fwdinitstep;

    TimeStep<State> & linstep;    
    StateAlg<State> & lininitstep;
    TimeTerm & linterm;
    mutable ListAlg *lininit;
    mutable TargetCurrTimeStep<State, RefState> * linpost;
    mutable TimeStepList<State> * linlist;
    mutable StdSim<State> * lin;
   
    mutable ListAlg *adjinit;
    TimeStep<State> & adjstep;
    StateAlg<State> & adjinitstep;
    TimeTerm & adjterm;
    mutable TargetNextTimeStep<State, RefState> * adjpre;
    mutable TimeStepList<State> * adjlist;
    mutable StdSim<State> * adj;

    StdSim1Jet();
    StdSim1Jet(StdSim1Jet<State, RefState> const &);

  public:

    /** main constructor: takes forward simulator and raw ingredients
	for linearized, adjoint simulators. */
    StdSim1Jet(Sim<State> & _fwd,

	       TimeStep<State> & _linstep,
	       StateAlg<State> &_lininitstep,
	       TimeTerm & _linterm,	

	       TimeStep<State> & _adjstep,
	       StateAlg<State> &_adjinitstep,
	       TimeTerm & _adjterm
	
	       ) 
      : fwd(_fwd),
	fwdinitstep(_fwd.getInitStep()),
	linstep(_linstep),
	lininitstep(_lininitstep),
	linterm(_linterm),
	adjstep(_adjstep),
	adjinitstep(_adjinitstep),
	adjterm(_adjterm),

	lin(NULL),linpost(NULL),linlist(NULL), lininit(NULL),
	adj(NULL),adjpre(NULL), adjlist(NULL), adjinit(NULL)
    {}

    ~StdSim1Jet() {

      //      cerr<<" begin StdSim1Jet destr\n";
      if (lin) delete lin;
      if (linlist) delete linlist;
      if (linpost) delete linpost;
      if (lininit) delete lininit;

      if (adj) delete adj;
      if (adjpre) delete adjpre;
      if (adjlist) delete adjlist;
      if (adjinit) delete adjinit;
      //      cerr<<"   end StdSim1Jet destr\n";

    }

    Sim<State> & getFwd() const { return fwd; }
    
    Sim<State> & getLin() const {
      try  {
	if (!lininit) lininit = new ListAlg(fwdinitstep, lininitstep);
	if (!linpost) linpost = new TargetCurrTimeStep<State, RefState>(fwd, linstep);
	if (!linlist) linlist = new TimeStepList<State>(linstep, *linpost, false);
	if (!lin)         lin = new StdSim<State>(*linlist, linterm, *lininit);
	if (lin)  return *lin; 
	RVLException e; 
	e<<"Error: Sim::getLin - linearized Sim member not constructed\n";
	throw e;
      }
      catch (RVLException & e) {
	e<<"\ncalled from StdSim1Jet::getLin\n";
	throw e;
      }
    }
    
    Sim<State> & getAdj() const {
      try {
	if (!adjinit) adjinit = new ListAlg(fwdinitstep, adjinitstep);
	if (!adjpre)  adjpre  = new TargetNextTimeStep<State, RefState>(fwd, adjstep);
	if (!adjlist) adjlist = new TimeStepList<State> (adjstep, *adjpre, true);
	if (!adj)         adj = new StdSim<State>(*adjlist, adjterm, *adjinit);

	//	cerr<<"********** in StdSim1Jet::getAdj"<<endl;
	//	adjstep.write(cerr);

	if (adj) return *adj;

	RVLException e; 
	e<<"Error: Sim::getAdj - adj Sim member not constructed\n";
	throw e;
      }

      catch(RVLException &e) {
	e<<"\ncalled from StdSim1Jet::getAdj\n";
	throw e;
      }

    }
    
  };


  
  /** Jet class with offset adj step -- calls getCurrTimeStep
      on the adjpre step, instead of getNextTimeStep */
  template<typename State, typename RefState>
  class StdSim1Jet_AC: public Sim1Jet<State> {

  private:

    Sim<State> & fwd;
    Algorithm &fwdinitstep;

    TimeStep<State> & linstep;    
    StateAlg<State> & lininitstep;
    TimeTerm & linterm;
    mutable ListAlg *lininit;
    mutable TargetCurrTimeStep<State, RefState> * linpost;
    mutable TimeStepList<State> * linlist;
    mutable StdSim<State> * lin;
   
    mutable ListAlg *adjinit;
    TimeStep<State> & adjstep;
    StateAlg<State> & adjinitstep;
    TimeTerm & adjterm;
    mutable TargetCurrTimeStep<State, RefState> * adjpre;
    mutable TimeStepList<State> * adjlist;
    mutable StdSim<State> * adj;

    StdSim1Jet_AC();
    StdSim1Jet_AC(StdSim1Jet_AC<State, RefState> const &);

  public:

    /** main constructor: takes forward simulator and raw ingredients
	for linearized, adjoint simulators. */
    StdSim1Jet_AC(Sim<State> & _fwd,

		  TimeStep<State> & _linstep,
		  StateAlg<State> &_lininitstep,
		  TimeTerm & _linterm,	

		  TimeStep<State> & _adjstep,
		  StateAlg<State> &_adjinitstep,
		  TimeTerm & _adjterm
	
		  ) 
      : fwd(_fwd),
	fwdinitstep(_fwd.getInitStep()),
	linstep(_linstep),
	lininitstep(_lininitstep),
	linterm(_linterm),
	adjstep(_adjstep),
	adjinitstep(_adjinitstep),
	adjterm(_adjterm),

	lin(NULL),linpost(NULL),linlist(NULL), lininit(NULL),
	adj(NULL),adjpre(NULL), adjlist(NULL), adjinit(NULL)
    {}

    ~StdSim1Jet_AC() {

      if (lin) delete lin;
      if (linlist) delete linlist;
      if (linpost) delete linpost;
      if (lininit) delete lininit;

      if (adj) delete adj;
      if (adjpre) delete adjpre;
      if (adjlist) delete adjlist;
      if (adjinit) delete adjinit;
    }

    Sim<State> & getFwd() const { return fwd; }
    
    Sim<State> & getLin() const 
    {
      try  
	{
	  if (!lininit) lininit = new ListAlg(fwdinitstep, lininitstep);
	  if (!linpost) linpost = new TargetCurrTimeStep<State, RefState>(fwd, linstep);
	  if (!linlist) linlist = new TimeStepList<State>(linstep, *linpost, false);
	  if (!lin)	lin =  new StdSim<State>(*linlist, linterm, *lininit);
	  if (lin)  return *lin; 
	  RVLException e; 
	  e<<"Error: Sim::getLin - linearized Sim member not constructed\n";
	  throw e;
	}
      catch (RVLException & e) 
	{
	  e<<"\ncalled from StdSim1Jet_AC::getLin\n";
	  throw e;
	}
    }
    
    Sim<State> & getAdj() const 
    {
      try
	{
	  if (!adjinit) adjinit =  new ListAlg(fwdinitstep, adjinitstep);
	  if (!adjpre)  adjpre	= new TargetCurrTimeStep<State, RefState>(fwd, adjstep);
	  if (!adjlist) adjlist = new TimeStepList<State> (adjstep, *adjpre, true);
	  if (!adj)	adj =  new StdSim<State>(*adjlist, adjterm, *adjinit);

	  if (adj) return *adj;

	  RVLException e; 
	  e<<"Error: Sim::getAdj - adj Sim member not constructed\n";
	  throw e;
	}

      catch(RVLException &e) 
	{
	  e<<"\ncalled from StdSim1Jet_AC::getAdj\n";
	  throw e;
	}

    }

    
  };








  /** Jet for adaptive time stepping classes */
  template<typename State, typename RefState>
  class RealSim1Jet: public Sim1Jet<State> {

  private:

    Sim<State> & fwd;
    Algorithm &fwdinitstep;

    TimeStep<State> & linstep;    
    StateAlg<State> & lininitstep;
    TimeTerm & linterm;
    mutable ListAlg *lininit;
    mutable TargetCurrTimeStep<State, RefState> * linpost;
    mutable TimeStepList<State> * linlist;
    mutable StdSim<State> * lin;
   
    mutable ListAlg *adjinit;
    TimeStep<State> & adjstep;
    StateAlg<State> & adjinitstep;
    TimeTerm & adjterm;
    mutable TargetNextTimeStep<State, RefState> * adjpre;
    mutable TimeStepList<State> * adjlist;
    mutable StdSim<State> * adj;

    RealSim1Jet();
  
  public:

    /** copy constructor for jet; needed for external use of this class 
	with the alg::umin package */
    RealSim1Jet(RealSim1Jet<State, RefState> const & _j): fwd(_j.fwd), fwdinitstep(_j.fwdinitstep),
							  linstep(_j.linstep), lininitstep(_j.lininitstep),
							  linterm(_j.linterm), lininit(_j.lininit),
							  linpost(_j.linpost), linlist(_j.linlist),
							  lin(_j.lin), adjinit(_j.adjinit), adjstep(_j.adjstep),
							  adjinitstep(_j.adjinitstep), adjterm(_j.adjterm),
							  adjpre(_j.adjpre), adjlist(_j.adjlist), adj(_j.adj) {}

    

    /** main constructor: takes forward simulator and raw ingredients
	for linearized, adjoint simulators. */
    RealSim1Jet(Sim<State> & _fwd,

		TimeStep<State> & _linstep,
		StateAlg<State> &_lininitstep,
		TimeTerm & _linterm,	

		TimeStep<State> & _adjstep,
		StateAlg<State> &_adjinitstep,
		TimeTerm & _adjterm
	
		) 
      : fwd(_fwd),
	fwdinitstep(_fwd.getInitStep()),
	linstep(_linstep),
	lininitstep(_lininitstep),
	linterm(_linterm),
	adjstep(_adjstep),
	adjinitstep(_adjinitstep),
	adjterm(_adjterm),

	lin(NULL),linpost(NULL),linlist(NULL), lininit(NULL),
	adj(NULL),adjpre(NULL), adjlist(NULL), adjinit(NULL)
    {}

    ~RealSim1Jet() {

      if (lin) delete lin;
      if (linlist) delete linlist;
      if (linpost) delete linpost;
      if (lininit) delete lininit;

      if (adj) delete adj;
      if (adjpre) delete adjpre;
      if (adjlist) delete adjlist;
      if (adjinit) delete adjinit;
    }

    Sim<State> & getFwd() const { return fwd; }
    
    Sim<State> & getLin() const 
    {

      try  
	{
	  if (!lininit) lininit = new ListAlg(fwdinitstep, lininitstep);
	  if (!linpost) linpost = new TargetCurrTimeStep<State,RefState>(fwd, linstep);
	  if (!linlist) linlist = new TimeStepList<State>(linstep, *linpost, false);
	  if (!lin)	lin =  new StdSim<State>(*linlist, linterm, *lininit);
	  if (lin)  return *lin; 
	  RVLException e; 
	  e<<"Error: Sim::getLin - linearized Sim member not constructed\n";
	  throw e;
	}
      catch (RVLException & e) 
	{
	  e<<"\ncalled from RealSim1Jet::getLin\n";
	  throw e;
	}
    }
    
    Sim<State> & getAdj() const 
    {
      try
	{
	  if (!adjinit) adjinit =  new ListAlg(fwdinitstep, adjinitstep);
	  if (!adjpre)	adjpre	= new TargetNextTimeStep<State, RefState>(fwd, adjstep);
	  if (!adjlist) adjlist = new TimeStepList<State> (adjstep, *adjpre, true);
	  if (!adj)	adj =  new StdSim<State>(*adjlist, adjterm, *adjinit);

	  if (adj) return *adj;

	  RVLException e; 
	  e<<"Error: Sim::getAdj - adj Sim member not constructed\n";
	  throw e;
	}

      catch(RVLException &e) 
	{
	  e<<"\ncalled from RealSim1Jet::getAdj\n";
	  throw e;
	}

    }
    
  };





  /** Jet for adaptive time stepping classes, but uses the class 
      TargetCurrTimeStep to dictate behavior of when to stop the
      reference simulator during the adjoint evolution
  */
  template<typename State, typename RefState>
  class RealSim1Jet_AC: public Sim1Jet<State> {

  private:

    Sim<State> & fwd;
    Algorithm &fwdinitstep;

    TimeStep<State> & linstep;    
    StateAlg<State> & lininitstep;
    TimeTerm & linterm;
    mutable ListAlg *lininit;
    mutable TargetCurrTimeStep<State, RefState> * linpost;
    mutable TimeStepList<State> * linlist;
    mutable StdSim<State> * lin;
   
    mutable ListAlg *adjinit;
    TimeStep<State> & adjstep;
    StateAlg<State> & adjinitstep;
    TimeTerm & adjterm;
    mutable TargetCurrTimeStep<State, RefState> * adjpre;
    mutable TimeStepList<State> * adjlist;
    mutable StdSim<State> * adj;

    RealSim1Jet_AC();
  
  public:

    /** copy constructor for jet; needed for external use of this class 
	with the alg::umin package */
    RealSim1Jet_AC(RealSim1Jet_AC<State, RefState> const & _j): fwd(_j.fwd), fwdinitstep(_j.fwdinitstep),
								linstep(_j.linstep), lininitstep(_j.lininitstep),
								linterm(_j.linterm), lininit(_j.lininit),
								linpost(_j.linpost), linlist(_j.linlist),
								lin(_j.lin), adjinit(_j.adjinit), adjstep(_j.adjstep),
								adjinitstep(_j.adjinitstep), adjterm(_j.adjterm),
								adjpre(_j.adjpre), adjlist(_j.adjlist), adj(_j.adj) {}

    

    /** main constructor: takes forward simulator and raw ingredients
	for linearized, adjoint simulators. */
    RealSim1Jet_AC(Sim<State> & _fwd,

		   TimeStep<State> & _linstep,
		   StateAlg<State> &_lininitstep,
		   TimeTerm & _linterm,	

		   TimeStep<State> & _adjstep,
		   StateAlg<State> &_adjinitstep,
		   TimeTerm & _adjterm
	
		   ) 
      : fwd(_fwd),
	fwdinitstep(_fwd.getInitStep()),
	linstep(_linstep),
	lininitstep(_lininitstep),
	linterm(_linterm),
	adjstep(_adjstep),
	adjinitstep(_adjinitstep),
	adjterm(_adjterm),

	lin(NULL),linpost(NULL),linlist(NULL), lininit(NULL),
	adj(NULL),adjpre(NULL), adjlist(NULL), adjinit(NULL)
    {}

    ~RealSim1Jet_AC() {

      if (lin) delete lin;
      if (linlist) delete linlist;
      if (linpost) delete linpost;
      if (lininit) delete lininit;

      if (adj) delete adj;
      if (adjpre) delete adjpre;
      if (adjlist) delete adjlist;
      if (adjinit) delete adjinit;
    }

    Sim<State> & getFwd() const { return fwd; }
    
    Sim<State> & getLin() const 
    {

      try  
	{
	  if (!lininit) lininit = new ListAlg(fwdinitstep, lininitstep);
	  if (!linpost) linpost = new TargetCurrTimeStep<State,RefState>(fwd, linstep);
	  if (!linlist) linlist = new TimeStepList<State>(linstep, *linpost, false);
	  if (!lin)	lin =  new StdSim<State>(*linlist, linterm, *lininit);
	  if (lin)  return *lin; 
	  RVLException e; 
	  e<<"Error: Sim::getLin - linearized Sim member not constructed\n";
	  throw e;
	}
      catch (RVLException & e) 
	{
	  e<<"\ncalled from RealSim1Jet::getLin\n";
	  throw e;
	}
    }
    
    Sim<State> & getAdj() const 
    {
      try
	{
	  if (!adjinit) adjinit =  new ListAlg(fwdinitstep, adjinitstep);
	  if (!adjpre)	adjpre	= new TargetCurrTimeStep<State, RefState>(fwd, adjstep);
	  if (!adjlist) adjlist = new TimeStepList<State> (adjstep, *adjpre, true);
	  if (!adj)	adj =  new StdSim<State>(*adjlist, adjterm, *adjinit);

	  if (adj) return *adj;

	  RVLException e; 
	  e<<"Error: Sim::getAdj - adj Sim member not constructed\n";
	  throw e;
	}

      catch(RVLException &e) 
	{
	  e<<"\ncalled from RealSim1Jet::getAdj\n";
	  throw e;
	}

    }
    
  };





}

#endif

