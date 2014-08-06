#ifndef __TSOPT_MPI_SEGYPP
#define __TSOPT_MPI_SEGYPP

#include "segypp.hh"
#include "mpiserialdc.hh"
#include "mpiseriallap.hh"

namespace TSOpt {

  using RVL::Space;
  using RVL::StdSpace;
  using RVL::MPISerialDC;
  using RVL::MPISerialDCF;
  using RVL::MPISerialLAP;
  using RVL::LinearAlgebraPackage;
  using RVL::DataContainerFactory;
  using RVL::ScalarFieldTraits;

  class MPISEGYSpace: public StdSpace<float,float>,
		      public ConstContainer<STRING_PAIR> { 

  private:
    
    SEGYDCF f;
    MPISerialDCF mpif; 
    MPISerialLAP<float> l;
    ostream & outfile;

  public:

    MPISEGYSpace(string hdr, string key
#ifdef IWAVE_USE_MPI 
		 , MPI_Comm _comm=retrieveGlobalComm()
#endif
		 , ostream & _outfile = cerr)
      : StdSpace<float,float>(), 
	ConstContainer<STRING_PAIR>(STRING_PAIR(key,hdr)),
	f(hdr,_outfile),
	mpif(f),
	l(ScalarFieldTraits<float>::One()
	  , 0
#ifdef IWAVE_USE_MPI 
	  , _comm
#endif
	  ),
	outfile(_outfile) {
      int rk = 0;
#ifdef IWAVE_USE_MPI 
      MPI_Comm_rank(_comm, &rk);
#endif
      if (rk==0) l.setScale(f.getDt());
    }
    
    MPISEGYSpace(MPISEGYSpace const & sp) 
      : StdSpace<float,float>(),
	ConstContainer<STRING_PAIR>(sp.get()),
	f(sp.get().val,sp.outfile),
	mpif(f),
	l(sp.l),
	outfile(sp.outfile){ }
    
    ~MPISEGYSpace() {}

    LinearAlgebraPackage<float> const & getLAP() const { return l; }

    DataContainerFactory const & getDCF() const { return mpif; }

    /** return time step for scaling l2 inner product */
    float getDt() const { 
      try { 
	return f.getDt(); }
      catch (RVLException & e) {
	e<<"\ncalled from SEGYSpace::getDt\n";
	throw e;
      }
    }
    /** return number of time samples */
    int getNt() const { 
      try {
	return f.getNt(); }
      catch  (RVLException & e) {
	e<<"\ncalled from SEGYSpace::getNt\n";
	throw e;
      }
    }

    ostream & write(ostream & str) const {
      str<<"MPISEGYSpace: MPI-enabled StdSpace based on SEGYDC\n";
      str<<"trace data container class\n";
      str<<" - header file = "<<this->get().val<<"\n";
      str<<" - data type   = "<<this->get().key<<"\n";
      return str;
    }
  };
}	      
      
#endif
