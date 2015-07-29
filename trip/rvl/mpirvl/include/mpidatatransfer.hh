// mpidatatransfer.H
// created by ADP 3/21/05
// modified by WWS spring 06

// This file defines a templated function for selecting the proper
// enumerated type to match a given basic type, and two functions for
// sending and receiving data to/from a LDC.

#ifndef __PARAMS_MPI_DATA_TRANSFER_
#define __PARAMS_MPI_DATA_TRANSFER_

#include "mpi.h"
#include "std_cpp_includes.hh"
#include "mpiutils.hh"
#include "except.hh"

namespace RVL {

  /** base template should always throw exception - only specializations
      are meaningful */
  template<typename Scalar>
  inline static MPI_Datatype findMPIDatatype() {
    RVLException e;
    e<<"Error: call to base findMPIDatatype template - undefined \n";
    throw e;
  }

  template<>
  inline MPI_Datatype findMPIDatatype<char>() {return MPI_CHAR;}

  template<>
  inline MPI_Datatype findMPIDatatype<short>() {return MPI_SHORT;}

  template<>
  inline MPI_Datatype findMPIDatatype<int>() {return MPI_INT;}

  template<>
  inline MPI_Datatype findMPIDatatype<long>() {return MPI_LONG;}

  template<>
  inline MPI_Datatype findMPIDatatype<unsigned char>() {return MPI_UNSIGNED_CHAR;}

  template<>
  inline MPI_Datatype findMPIDatatype<unsigned short>() {return MPI_UNSIGNED_SHORT;}

  template<>
  inline MPI_Datatype findMPIDatatype<unsigned int>() {return MPI_UNSIGNED;}

  template<>
  inline MPI_Datatype findMPIDatatype<unsigned long>() {return MPI_UNSIGNED_LONG;}

  template<>
  inline MPI_Datatype findMPIDatatype<float>() {return MPI_FLOAT;}

  template<>
  inline MPI_Datatype findMPIDatatype<double>() {return MPI_DOUBLE;}

  template<>
  inline MPI_Datatype findMPIDatatype<long double>() {return MPI_LONG_DOUBLE;}

  /** default implementations for built-in MPI_Datatypes */

  template<typename T>
  void Build_MPI_Datatype(MPI_Datatype * mpitype) {
    try {
      *mpitype=findMPIDatatype<T>();
    }
    catch (RVLException & e) {
      e<<"\ncalled from BuildMPI_Datatype\n";
      throw e;
    }
  }
      
  /** Suitable for use with any type T with fixed size, i.e.
      for which instantiation and allocation are the same. Should
      be specialized appropriately for types with dynamic size. 
  */

  template<typename T>
  class MPI_Sender {

  private:

    int dest;
    MPI_Comm comm;

  public:

    MPI_Sender(int _dest=0, MPI_Comm _comm = MPI_COMM_WORLD) 
      : dest(_dest), comm(_comm) {}
    MPI_Sender(MPI_Sender<T> const & s) 
      : dest(s.dest), comm(s.comm) {}
    ~MPI_Sender() {}

    bool setDestination( int d ){ return MPIRVL_SetRank(dest,d,comm); }

    bool operator()(T & x) const {
      try {
	int tag=0;
	MPI_Datatype dt;
	Build_MPI_Datatype<T>(&dt);
	if (MPI_SUCCESS == MPI_Send(&x,1,dt,dest,tag,comm)) return true;
	return false;
      }
      catch (RVLException & e) {
	e<<"\ncalled from MPI_Sender::operator()\n";
	throw e;
      }
    }
  };

  template<typename T>
  class MPI_Receiver {

    MPI_Status * status;
    int src;
    MPI_Comm comm;

  public:

    MPI_Receiver(MPI_Status * _status,
		 int _src=0, 
		 MPI_Comm _comm = MPI_COMM_WORLD) 
      : status(_status), src(_src), comm(_comm) {}
    MPI_Receiver(MPI_Receiver<T> const & s) 
      : status(s.status),src(s.src), comm(s.comm) {}
    ~MPI_Receiver() {}

    bool setSource( int s ){ return MPIRVL_SetRank(src,s,comm); }

    bool operator()(T & x) const {
      try {
	int tag=0;
	MPI_Datatype dt;
	Build_MPI_Datatype<T>(&dt);
	if (MPI_SUCCESS == MPI_Recv(&x,1,dt,src,tag,comm,status)) return true;
	return false;
      }
      catch (RVLException & e) {
	e<<"\ncalled from MPI_Receiver::operator()\n";
	throw e;
      }
    }
  };

  template<typename T>
  class MPI_Broadcaster {

    int root;
    MPI_Comm comm;

  public:

    MPI_Broadcaster(int _root=0, MPI_Comm _comm = MPI_COMM_WORLD)
      : root(_root), comm(_comm) {
      int rt;
      if (!MPIRVL_SetRank(rt,root,comm)) {
	RVLException e;
	e<<"Error: MPI_Broadcaster constructor\n";
	e<<"root outside range of ranks for assigned MPI_Comm\n";
	throw e;
      }
    }
    MPI_Broadcaster(MPI_Broadcaster<T> const & b)
      : root(b.root), comm(b.comm) {}
    ~MPI_Broadcaster() {}

    bool operator()(T & x) const {
      try {
	MPI_Datatype dt;
	Build_MPI_Datatype<T>(&dt);
	if (MPI_SUCCESS == MPI_Bcast(&x,1,dt,root,comm)) return true;
	return false;
      }
      catch (RVLException & e) {
	e<<"\ncalled from MPI_Broadcaster::operator()\n";
	throw e;
      }
    }
  };

  template<typename T>
  class MPI_Reducer {

    MPI_Op op;
    int root;
    MPI_Comm comm;

  public:

    MPI_Reducer(MPI_Op _op= MPI_SUM,
		int _root=0, 
		MPI_Comm _comm = MPI_COMM_WORLD)
      : op(_op), root(_root), comm(_comm) {
      int rt;
      if (!MPIRVL_SetRank(rt,root,comm)) {
	RVLException e;
	e<<"Error: MPI_Reducer constructor\n";
	e<<"root outside range of ranks for assigned MPI_Comm\n";
	throw e;
      }
    }
    MPI_Reducer(MPI_Reducer<T> const & b)
      : op(b.op), root(b.root), comm(b.comm) {}
    ~MPI_Reducer() {}

    bool operator()(T & xout, T const & xin) const {
      try {
	MPI_Datatype dt;
	Build_MPI_Datatype<T>(&dt);
	if (MPI_SUCCESS == MPI_Reduce(&xin,&xout,1,dt,op,root,comm)) return true;
	return false;
      }
      catch (RVLException & e) {
	e<<"\ncalled from MPI_Reducer::operator()\n";
	throw e;
      }
    }
  };

}

#endif // __PARAMS_MPI_DATA_TRANSFER_
