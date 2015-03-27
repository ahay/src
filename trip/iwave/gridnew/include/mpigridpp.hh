#ifndef __TSOPT_MPI_GRIDPP
#define __TSOPT_MPI_GRIDPP

#include "usempi.h"
#include "gridpp.hh"
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
  using RVL::STRING_PAIR;
  using RVL::ConstContainer;

  class MPIGridSpace: public StdSpace<ireal,ireal>,
		      public ConstContainer<STRING_PAIR> {

  private:
    
    GridDCF f;
    MPISerialDCF mpif; 
    MPISerialLAP<ireal> l;
    ostream & outfile;                 // verbose output unit

  public:

    MPIGridSpace(string hdr,
		 string dtype="notype",
		 bool _incore=false
#ifdef IWAVE_USE_MPI 
		 //      		 , MPI_Comm _comm=MPI_COMM_WORLD
       		 , MPI_Comm _comm=retrieveGlobalComm()
#endif
		 , ostream & _outfile = cerr)
      : StdSpace<ireal,ireal>(),
	ConstContainer<STRING_PAIR>(STRING_PAIR(dtype,hdr)),
	f(hdr,_incore,_outfile),
	mpif(f),
	l(f.getCellVol()
	  , 0
#ifdef IWAVE_USE_MPI 
	  , _comm
#endif
	  ), 
	outfile(_outfile) {}

    MPIGridSpace(MPIGridSpace const & sp) 
      : StdSpace<ireal,ireal>(), 
	ConstContainer<STRING_PAIR>(sp.get()),
	f(sp.f),
	mpif(f),
	l(sp.l),
        outfile(sp.outfile) {}

    ~MPIGridSpace() {}

    LinearAlgebraPackage<ireal> const & getLAP() const { return l; }

    DataContainerFactory const & getDCF() const { return mpif; }

    grid const & getGrid() const { 
      if (retrieveGlobalRank()) {
	RVLException e;
	e<<"Error: MPIGridSpace::getGrid\n";
	e<<"should not be called on ranks other than 0\n";
	throw e;
      }
      return f.getGrid(); 
    }

    bool isIncore() const { return f.isIncore(); }

    ostream & write(ostream & str) const {
      str<<"MPIGridSpace: MPI-enabled StdSpace based on RSF grid data\n";
      str<<" - header file = "<<f.getFilename()<<"\n";
      return str;
    }
  };
}	      
      
#endif

  
