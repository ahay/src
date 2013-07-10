#include "usempi.h"
#include "rkstream.hh"
#include "gridpp_top.hh"

extern "C" {
#include "iwave_fopen.h"
}

#include "create_hfile.hh"

using namespace RVL;
using namespace TSOpt;

int main(int argc, char ** argv) {

  int rk=0;

#ifdef IWAVE_USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rk);
#endif 

  try {

    string fname="testsrc/test13/dtestgrid.rsf";

    if (rk==0) {

      ofstream str;
      makeRankStream(str,rk,"testsrc/test13");
      FILE * fout = fopen("testsrc/test13/cout0.txt","r+");
      fseek(fout,0L,SEEK_END);

      cout<<"GRIDPP Unit Test 13"<<endl;
      cout<<"create GridSpace, construct 2 Vectors in space"<<endl;
      cout<<"assign const 1, 2 to these, execute axpy with a=1"<<endl;
      cout<<"first vector should contain all samples = 3.0"<<endl;
      cout<<"compute max and min"<<endl;

      create_hfile(fname,1);
 
      GridSpace<double> sp(fname,"notype",str);

      fprintf(fout,"file system before creation of tmp vectors, eval of \n");
      fprintf(fout,"AssignConst FOs and linComb\n");
      iwave_fprintall(fout);

      Vector<double> v1(sp);
      Vector<double> v2(sp);

      RVLAssignConst<double> ac1(1.0);
      RVLAssignConst<double> ac2(2.0);

      v1.eval(ac1);
      v2.eval(ac2);


      v1.linComb(1.0,v2);

      fprintf(fout,"file system before creation of tmp vectors, eval of \n");
      fprintf(fout,"AssignConst FOs and linComb\n");
      iwave_fprintall(fout);

      RVLMax<double> mx;
      v1.eval(mx);
      
      RVLMin<double> mn;
      v1.eval(mn);

      cout<<"test12 result: max="<<mx.getValue()<<" min="<<mn.getValue()<<endl;

      str.close();
      fclose(fout);
      iwave_fdestroy();

    }
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif

    return(0);

  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
}
