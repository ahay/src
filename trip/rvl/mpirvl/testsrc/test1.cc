#include "mpiserialdc.hh"
#include "functions.hh"
#include "rn.hh"

using namespace RVL;

int main(int argc,char ** argv) {

  int rk=0;

  try {

#ifdef IWAVE_USE_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);
#endif 

    PlantSeeds(19490615);

    int n=11;

    RnDataContainerFactory<float> f(n);
    MPISerialDCF mpif(f);

    DataContainer * d = mpif.build();
    vector<DataContainer const *> x(0);

    RVLRandomize<float> ac;
    
    d->eval(ac,x);
    
    RVLMax<float> mx;
    MPISerialFunctionObjectRedn<float,float> mpimx(mx);
    d->eval(mpimx,x);

    RVLMin<float> mn;
    MPISerialFunctionObjectRedn<float,float> mpimn(mn);
    d->eval(mpimn,x);
      
    cout<<"rk="<<rk<<" test1 result: max="<<mpimx.getValue()<<" min="<<mpimn.getValue()<<endl;
    
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
    return 0;
  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif    
    exit(1);
  }
}

