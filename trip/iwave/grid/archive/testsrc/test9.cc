#include "usempi.h"
#include "rkstream.hh"
#include "gridpp_top.hh"
extern "C" {
#include "iwave_fopen.h"
}
#include "create_hfile.hh"

int main(int argc, char ** argv) {

  using namespace RVL;
  using namespace TSOpt;
  using RVL::makeRankStream;

  int rk=0;

#ifdef IWAVE_USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rk);
#endif 

  try {

    if (rk==0) {

      ofstream str;
      makeRankStream(str,rk,"testsrc/test9");
      FILE * fout = fopen("testsrc/test9/cout0.txt","r+");
      fseek(fout,0L,SEEK_END);

      cout<<"GRIDPP Unit Test 9"<<endl;
      cout<<"create GridDC<double>, perm old file=dtest9.rsf"<<endl;
      cout<<"all data samples in file = 9.0, compute max and min"<<endl;

      string hfile="testsrc/test9/dtestgrid.rsf";
      string ffile="testsrc/test9/dtest9.rsf";

      create_hfile(hfile,1);
      create_hfile(ffile,1);
      
      vector<DataContainer const *> x(0);
      
      GridDC<double> d(hfile,str);

      AssignFilename af(ffile);
      d.eval(af,x);

      fprintf(fout,"file system before assignconst eval\n");
      iwave_fprintall(fout);

      RVLAssignConst<double> ac(9.0);
      d.eval(ac,x);

      fprintf(fout,"file system after assignconst eval\n");
      iwave_fprintall(fout);

      RVLMax<double> mx;
      d.eval(mx,x);
      
      RVLMin<double> mn;
      d.eval(mn,x);

      cout<<"test9 result: max="<<mx.getValue()<<" min="<<mn.getValue()<<endl;

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
