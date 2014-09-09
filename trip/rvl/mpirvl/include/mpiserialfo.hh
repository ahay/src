#ifndef __RVL_MPIFO
#define __RVL_MPIFO

#ifdef IWAVE_USE_MPI
#include "mpidatatransfer.hh"
#endif
#include "local.hh"

namespace RVL {

  /** Mixin asserting ability to broadcast internal data.
   */

  class MPISynchRoot {
  public:
    /** initialization of internal data */
    virtual void set() = 0;

    /** synchronizes internal data in all instances of 
	object, across processes
    */
    virtual void synch() = 0;
      virtual ~MPISynchRoot() {};
  };

  /** MPI-enabled function object for serial evaluation. Simplest
      parallel evaluation - not parallel at all! This implementation
      deals only with the world comm, and evaluates a standard FO on
      rk=0. Constructed to be a pass-through outside of MPI
      environment. Concrete class, confined to LDCs: depends on being
      able to cast input to LocalEvaluation. */

  template<class DataType>
  class MPISerialFunctionObject: public LocalFunctionObject<DataType> {

  private: 

    LocalEvaluation<DataType> * fptr;
    string fname;
    MPISerialFunctionObject();

  public:
  
    MPISerialFunctionObject(FunctionObject & f): fptr(NULL) {
      if (!(fptr=dynamic_cast<LocalEvaluation<DataType> *>(&f))) {
	RVLException e;
	e<<"Error: MPISerialFunctionObject constructor\n";
	e<<"input function, named "<<f.getName()<<" not LocalEvaluation\n";
	throw e;
      }
      fname=f.getName();
    }
    MPISerialFunctionObject(MPISerialFunctionObject<DataType> const & f)
      : fptr(f.fptr), fname(f.fname) {}
    ~MPISerialFunctionObject() {}

    void operator()(LocalDataContainer<DataType> & target,
		    vector<LocalDataContainer<DataType> const *> & sources) {
      try {
      
#ifdef IWAVE_USE_MPI
	int rk=0;
	MPI_Comm_rank(MPI_COMM_WORLD,&rk);
	if (rk==0) 
#endif
	  fptr->operator()(target,sources);
      }
      catch (RVLException & e) {
	e<<"\ncalled from MPISerialFunctionObject::operator()\n";
	throw e;
      }
    }

    string getName() const {
      string str="MPISerialFO wrapper around ";
      str+=fname;
      return str;
    }
  };

  /** MPI-enabled serial reduction of data container. Simplest
      parallel evaluation - not parallel at all! This implementation
      deals only with the world comm, and evaluates a standard FOR on
      rk=0. Initialization via set() on ALL PROCESSES. Broadcasts
      result to other processes via synch(). Thus all invocations must
      call set() before eval(), then synch(), to correctly compute
      reduction results and place on all processes. Constructed to be
      a pass-through outside of MPI environment. Concrete class,
      confined to LDCs: depends on being able to cast input to
      LocalConstEval, and to MPISynch. */

  template<typename DataType, typename ValType>
  class MPISerialFunctionObjectRedn: 
    public FunctionObjectScalarRedn<ValType>, 
    public LocalConstEval<DataType>,
    public MPISynchRoot {

  private: 

    LocalConstEval<DataType> * lrptr;          // reduction parent - data op
    FunctionObjectScalarRedn<ValType> & f;     // FO parent - redn value
    string fname;                              // buffer for function name
    int root;                                  // root process rank
    int rk;                                    // this process rank

#ifdef IWAVE_USE_MPI
    MPI_Comm comm;
#endif

    MPISerialFunctionObjectRedn();             

  public:
  
    MPISerialFunctionObjectRedn(FunctionObjectScalarRedn<ValType> & _f,
				int _root=0
#ifdef IWAVE_USE_MPI  				
				, MPI_Comm _comm=MPI_COMM_WORLD
#endif
				)
      : FunctionObjectScalarRedn<ValType>(_f.getValue()),
	lrptr(NULL), f(_f), root(_root), rk(_root)
#ifdef IWAVE_USE_MPI
      , comm(_comm)
#endif
    {
#ifdef IWAVE_USE_MPI
	MPI_Comm_rank(comm,&rk);
#endif
      if (!(lrptr=dynamic_cast<LocalConstEval<DataType> *>(&f))) {
	RVLException e;
	e<<"Error: MPISerialFunctionObjectRedn constructor\n";
	e<<"input function, named "<<f.getName()<<" not LocalReduction\n";
	throw e;
      }
      fname=f.getName();
    }

    MPISerialFunctionObjectRedn(MPISerialFunctionObjectRedn<DataType,ValType> const & f)
      : lrptr(f.lrptr), f(f.f), root(f.root), rk(f.root),
#ifdef IWAVE_USE_MPI
	comm(f.comm), 
#endif
	fname(f.fname) {
#ifdef IWAVE_USE_MPI
	MPI_Comm_rank(comm,&rk);
#endif
    }

    ~MPISerialFunctionObjectRedn() {}

    void operator()(vector<LocalDataContainer<DataType> const *> & sources) {
      try {
	//	cerr<<"MPISFOR: operator() on "<<foptr->getName()<<endl;
	if (rk==root) {
	  //	  cerr<<"  before: val = "<<foptr->getValue()<<endl;
	  lrptr->operator()(sources);
	  //	  cerr<<"  after: val = "<<foptr->getValue()<<endl;
	}
	// thus far, the FOSR data member has been evaluated, which
	// presumably updates its internal data Scalar, on root. Next
	// step: still on root, transfer to the internal data Scalar
	// of this.
	ScalarRedn<ValType>::setValue(f.getValue());
      }
      catch (RVLException & e) {
	e<<"\ncalled from MPISerialFunctionObjectRedn::operator()\n";
	throw e;
      }
    }
    
    void synch() {
      try {
	//	cerr<<"rk="<<rk<<" MPISFOR::synch - FOR="<<foptr->getName()<<endl;
	//	cerr<<"rk="<<rk<<" MPISFOR::synch - extract FOR value - legit on root, else default"<<endl;
	ValType a = ScalarRedn<ValType>::getValue();
#ifdef IWAVE_USE_MPI
	//	cerr<<"rk="<<rk<<" MPISFOR::synch - broadcast val="<<a<<" from root"<<endl;
	MPI_Broadcaster<ValType> bc(root,comm);
	bc(a);
#endif
	//	cerr<<"rk="<<rk<<" MPISFOR::synch - set val="<<a<<" on all processes"<<endl;
	ScalarRedn<ValType>::setValue(a);
	// for safety's sake do same for FOSR data member
	f.setValue(a);
	//	cerr<<"rk="<<rk<<" MPISFOR::synch - return"<<endl;
      }
      catch (RVLException & e) {
	e<<"\ncalled from MPISerialFunctionObjectRedn::synch\n";
	throw e;
      }
    }

    /** implement this reinitialization function by using the default value
	from the FOSR data member */
    void setValue() {
      f.setValue();
      ScalarRedn<ValType>::setValue(f.getValue());
    }

    /** implement "blind" reinitialization from parent class */
    void set() {
      this->setValue();
    }

    string getName() const {
      string str="MPISerialFOR wrapper around ";
      str+=fname;
      return str;
    }
  };

}
#endif
