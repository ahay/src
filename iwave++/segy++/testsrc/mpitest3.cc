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
    makeRankStream(str,rk,"testsrc/mpitest3");

    string fname="testsrc/mpitest3/hdr.su";

    if (rk==0) {

      cout<<"MPI SEGYPP Unit Test 3"<<endl;
      cout<<"Create SEGYSpace, two vectors in it, assign 1.0 to all samples,"<<endl;
      cout<<"print out inner product - should be "<<101*11*2.0<<endl;

      system("sunull nt=101 ntr=11 dt=0.002|sushw key=sx a=1000|sushw key=gx a=2000 b=-100 > testsrc/mpitest3/hdr.su");
    }

#ifdef IWAVE_USE_MPI
    MPISEGYSpace sp(fname,cml,str);
#else
    SEGYSpace sp(fname,str);
#endif
    
    Vector<float> x(sp);
    Vector<float> y(sp);
      
    RVLAssignConst<float>  ac(1.0);
    x.eval(ac);
    y.eval(ac);
      
    if (rk==0) cout<<"mpitest3 result: inner product = "<<y.inner(x)<<endl;
    
    str.close();

#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif

    if(rk == 0) iwave_fdestroy();

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
