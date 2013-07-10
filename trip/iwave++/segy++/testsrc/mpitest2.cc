#include "rkstream.hh"
#include "mpisegypp.hh"

using namespace RVL;
using namespace TSOpt;

char ** xargv;

int main(int argc, char ** argv) {

  int rk=0;


#ifdef IWAVE_USE_MPI
  MPI_Init(&argc,&argv);
  storeComm(MPI_COMM_WORLD);
  MPI_Comm cml = retrieveComm();
  MPI_Comm_rank(cml, &rk);
#endif 

  try {

    ofstream str;
    makeRankStream(str,rk,"testsrc/mpitest2");

    string fname="testsrc/mpitest2/hdr.su";

    if (rk==0) {

      cout<<"MPI SEGYPP Unit Test 2"<<endl;
      cout<<"create SEGYSpace, construct 2 Vectors in space"<<endl;
      cout<<"assign const 1, 2 to these, execute axpy with a=1"<<endl;
      cout<<"first vector should contain all samples = 3.0"<<endl;
      cout<<"compute max and min"<<endl;

      system("sunull nt=101 ntr=11 dt=0.002|sushw key=sx a=1000|sushw key=gx a=2000 b=-100 > testsrc/mpitest2/hdr.su");
    }
    
#ifdef IWAVE_USE_MPI
    MPISEGYSpace sp(fname,cml,str);
#else
    SEGYSpace sp(fname,str);
#endif
    Vector<float> v1(sp);
    Vector<float> v2(sp);

    RVLAssignConst<float> ac1(1.0);
    RVLAssignConst<float> ac2(2.0);
    
    v1.eval(ac1);
    v2.eval(ac2);
    
    v1.linComb(1.0,v2);
    
    RVLMax<float> mx;
    MPISerialFunctionObjectRedn<float,float> mpimx(mx);
    v1.eval(mpimx);
      
    RVLMin<float> mn;
    MPISerialFunctionObjectRedn<float,float> mpimn(mn);
    v1.eval(mpimn);
    
    if (rk==0) cout<<"mpitest2 result: max="<<mpimx.getValue()<<" min="<<mpimn.getValue()<<endl;
    
    str.close();

#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif

    if (rk==0) iwave_fdestroy();

    return(0);

  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(cml,0);
#endif
    exit(1);
  }
}
