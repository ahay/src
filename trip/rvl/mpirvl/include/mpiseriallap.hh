#ifndef __RVL_MPI_SERLAP
#define __RVL_MPI_SERLAP

#include "locallinalg.hh"
#include "mpiserialfo.hh"

namespace RVL {

  template<typename Scalar>
  class MPISerialLAP: public LinearAlgebraPackage<Scalar> {
    
  private:
    
    mutable RVLAssignConst<Scalar> this_zero;
    mutable RVLL2innerProd<Scalar> this_ser_inner;
    mutable RVLLinCombObject<Scalar> this_lco;
    
    int root;
#ifdef IWAVE_USE_MPI
    MPI_Comm comm;
#endif

    mutable MPISerialFunctionObjectRedn<Scalar,Scalar> this_mpi_inner;
    
  public:

    MPISerialLAP(Scalar ipscale = ScalarFieldTraits<Scalar>::One(),
		 int _root = 0
#ifdef IWAVE_USE_MPI 
		 , MPI_Comm _comm=MPI_COMM_WORLD
#endif
		 )
      : this_zero(ScalarFieldTraits<Scalar>::Zero()), 
	this_ser_inner(abs(ipscale)), this_lco(), 
	root(_root), 
#ifdef IWAVE_USE_MPI
	comm(_comm),
#endif
	this_mpi_inner(this_ser_inner,root
#ifdef IWAVE_USE_MPI
		       , comm
#endif
		       ) {}
    MPISerialLAP(const MPISerialLAP<Scalar> & p) 
      : this_zero(ScalarFieldTraits<Scalar>::Zero()),
	this_ser_inner(p.this_ser_inner.getScale()),this_lco(),
	root(p.root), 
#ifdef IWAVE_USE_MPI
	comm(p.comm),
#endif
	this_mpi_inner(this_ser_inner,root
#ifdef IWAVE_USE_MPI
		       , comm
#endif
		       ) {}
    ~MPISerialLAP() {}

    FunctionObjectScalarRedn<Scalar> & inner() const {
      return this_mpi_inner; 
    }
    FunctionObject & zero() const {
      return this_zero;
    }
    LinCombObject<Scalar> & linComb() const {
      return this_lco; 
    }

    virtual bool compare(LinearAlgebraPackage<Scalar> const & lap) const {
      MPISerialLAP<Scalar> const * tmp = NULL;
      if ((tmp=dynamic_cast<MPISerialLAP<Scalar> const *>(&lap))) return true;
      return false;
    }    

    /** added to spparate instantiation from initialization */
    void setScale(Scalar newscale) { this_ser_inner.setScale(newscale); }

    ostream & write(ostream & str) const {
      str<<"MPISerialLAP\n";
      return str;
    }
  };

}

#endif


