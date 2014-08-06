#include "mpiserialdc.hh"

namespace RVL {

  void MPISerialDC::eval(FunctionObject & f, 
			 std::vector<DataContainer const * > &x) {
    try {
      if (rk==0) {
	// to eval, must extract DC data members - these will 
	// be type-checked against this DC in many cases.
	std::vector<DataContainer const *> myx(x.size());
	MPISerialDC const * testdc;
	for (size_t i=0;i<x.size();i++) {
	  testdc=NULL;
	  if (!(testdc=dynamic_cast<MPISerialDC const *>(x[i]))) {
	    RVLException e;
	    e<<"Error: MPISerialDC::eval(FO)\n";
	    e<<"component "<<i<<" of input vector<DC *> not MPISerialDC\n";
	    throw e;
	  }
	  myx[i]=testdc->mydc;
	}
	mydc->eval(f,myx);
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from MPISerialDC::eval(for), for = "<<f.getName()<<"\n";
      throw e;
    }
  }
    
  void MPISerialDC::eval(FunctionObjectConstEval & f, 
			 vector<DataContainer const *> & x) const {
    try {
      MPISynchRoot & s = dynamic_cast<MPISynchRoot &>(f);
      // initialize internal data of FOR on all processes
      s.set();
      if (rk==0) {
	// to eval, must extract DC data members - these will 
	// be type-checked against this DC in many cases.
	std::vector<DataContainer const *> myx(x.size());
	MPISerialDC const * testdc;
	for (size_t i=0;i<x.size();i++) {
	  testdc=NULL;
	  if (!(testdc=dynamic_cast<MPISerialDC const *>(x[i]))) {
	    RVLException e;
	    e<<"Error: MPISerialDC::eval(FOR)\n";
	    e<<"component "<<i<<" of input vector<DC *> not MPISerialDC\n";
	    throw e;
	  }
	  myx[i]=testdc->mydc;
	}
	// evaluate FOR
	mydc->eval(f,myx);
      }
      // cerr<<"MPISerialDC: call synch on "<<f.getName()<<endl;
      s.synch();
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: MPISerialDC::eval called on FOR object not of MPISynchRoot type\n";
      e<<"does not know how to synchronize "<<f.getName()<<" result across comm\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from MPISerialDC::eval(for), for = "<<f.getName()<<"\n";
      throw e;
    }
  }

  ostream & MPISerialDC::write(ostream & str) const {
    str<<"MPISerialDataContainer *** beg rk = "<<rk<<" **********************\n";
    if (rk==0) mydc->write(str);
    else str<<"(no content)\n";
    str<<"************************** end rk = "<<rk<<" **********************\n";
    return str;
  }
  
  MPISerialDC * MPISerialDCF::buildMPISerialDC() const {
    try {
      return new MPISerialDC(f);
    }
    catch (RVLException & e) {
      e<<"\ncalled from MPISerialDCF::buildPC\n";
      throw e;
    }
  }

  bool MPISerialDCF::compare(DataContainerFactory const & dcf) const {
    try {
      MPISerialDCF const & c = dynamic_cast<MPISerialDCF const &>(dcf);
      return this->f.compare(c.f);
    }
    catch (bad_cast) {
      return false;
    }
    catch (RVLException & e) {
      e<<"\ncalled from MPISerialDCF::compare\n";
      throw e;
    }
  }

  bool MPISerialDCF::isCompatible(DataContainer const & dc) const {
    try {
      MPISerialDC const & c = dynamic_cast<MPISerialDC const &>(dc);
      int ires=0;
      if (rk==0) ires = this->f.isCompatible(*(c.mydc));
#ifdef IWAVE_USE_MPI
      MPI_Bcast(&ires,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
      if (ires) return true;
      return false;
    }
    catch (bad_cast) {
      return false;
    }
    catch (RVLException & e) {
      e<<"\ncalled from MPISerialDC::isCompatible\n";
      throw e;
    }
  }

  ostream & MPISerialDCF::write(ostream & str) const {
    str<<"MPISerialDCFactory ******* beg rk = "<<rk<<" **********************\n";
    if (rk==0) f.write(str);
    else str<<"(no content)\n";
    str<<"************************** end rk = "<<rk<<" **********************\n";
    return str;
  }
  
}

