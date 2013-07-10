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
    RVL::makeRankStream(str,rk,"testsrc/mpitest1");
    
    if (rk==0) {

      system("sunull nt=101 ntr=11 dt=0.002|sushw key=sx a=1000|sushw key=gx a=2000 b=-100 > testsrc/mpitest1/hdr.su");
      system("suop op=exp <testsrc/mpitest1/hdr.su | sugain scale=12.0 > testsrc/mpitest1/mpitest1.su");

      cout<<"MPISEGYPP Unit Test 1"<<endl;
      cout<<"create MPISEGYSpace, construct Vector in space"<<endl;
      cout<<"assign to perm old file=testsrc/mpitest1/mpitest1.su"<<endl;
      cout<<"all data samples in file = 12.0, compute max and min"<<endl;
      
    }
    
    string fname="testsrc/mpitest1/hdr.su";
    string dname="testsrc/mpitest1/mpitest1.su";

#ifdef IWAVE_USE_MPI
    MPISEGYSpace sp(fname,cml,str);
#else
    SEGYSpace sp(fname,str);
#endif

    Vector<float> v(sp);

    AssignFilename af(dname);
    v.eval(af);

    RVLMax<float> mx;
    MPISerialFunctionObjectRedn<float,float> mpimx(mx);
    v.eval(mpimx);
      
    RVLMin<float> mn;
    MPISerialFunctionObjectRedn<float,float> mpimn(mn);
    v.eval(mpimn);
    
    if (rk==0) cout<<"mpitest1 result: max="<<mpimx.getValue()<<" min="<<mpimn.getValue()<<endl;
    
    str.close();

#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif

    if (rk==0)  iwave_fdestroy();
    
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
