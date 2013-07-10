#ifndef __RVL_MPI_SERDC
#define __RVL_MPI_SERDC

#define FRUITCAKE

#include "mpiserialfo.hh"
#include "data.hh"

namespace RVL {

  class MPISerialDCF;

  class MPISerialDC: public DataContainer {
    
    friend class MPISerialDCF;

  private:
    
    MPISerialDC();
    MPISerialDC(MPISerialDC const &);

  protected:

    /** rank in multiprocess - =0 for default serial case */
    int rk;

    /** internal DC member */
    DataContainer * mydc;

  public:

    /** main constructor records reference to factory,
	allocates CP<T,M> buffer, and assigns a temp 
	filename */
    MPISerialDC(DataContainerFactory const & F):  rk(0), mydc(NULL) {
      /* set rank - all processes*/
#ifdef IWAVE_USE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD,&rk);
#endif
      /* construct internal PC */
      if (rk==0) mydc=F.build();
    }

    ~MPISerialDC() { if (mydc) delete mydc; }

    /** FO Eval defined by delegation to DC data member on rk=0
    */
    void eval(FunctionObject & f, 
	      std::vector<DataContainer const * > &x);
    
    /** FOCE Eval defined by delegation to DC data member on rk=0
    */
    void eval(FunctionObjectConstEval & f, 
	      vector<DataContainer const *> & x) const;

    ostream & write(ostream & str) const;
  };

  /** Factory class for MPISerialDCs, used both in MPISerialDC construction
      and initialization and in corresponding Space class. */
  class MPISerialDCF: public DataContainerFactory {

  private:

    int rk;

    /** captive DCF */
    DataContainerFactory const & f;

    MPISerialDC * buildMPISerialDC() const;

  public:

    MPISerialDCF(DataContainerFactory const & _f)
    : rk(0), f(_f) {
#ifdef IWAVE_USE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD,&rk);
#endif
    }
      
    MPISerialDCF(MPISerialDCF const & fact)
    : rk(0), f(fact.f) {
#ifdef IWAVE_USE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD,&rk);
#endif
    }
    ~MPISerialDCF() {}
    DataContainer * build() const { return buildMPISerialDC(); }
    bool compare(DataContainerFactory const & dcf) const;
    /** since MPISerialDCs only allocate data on rk=0, compatibility must
	be tested there and result broadcast. */
    bool isCompatible(DataContainer const & dc) const;
    ostream & write(ostream & str) const;
  };

}

#endif
