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

    if (rk==0) {
      
      ofstream str;
      makeRankStream(str,rk,"testsrc/test8");
      FILE * fout = fopen("testsrc/test8/cout0.txt","r+");
      fseek(fout,0L,SEEK_END);

      cout<<"GRIDPP Unit Test 8"<<endl;
      cout<<"create GridDCF, GridDC, perm old file=testgrid.rsf, assign"<<endl;
      cout<<"8.0 to all data samples, double case, compute max and min"<<endl;

      string fname="testsrc/test8/dtestgrid.rsf";
    
      create_hfile(fname,1);

      GridDC<double> d(fname,str);

      vector<DataContainer const *> x(0);
      
      AssignFilename af(fname);
      d.eval(af,x);

      fprintf(fout,"file system before assignconst eval\n");
      iwave_fprintall(fout);

      RVLAssignConst<double>  ac(8.0);
      d.eval(ac,x);
      
      fprintf(fout,"file system after assignconst eval\n");
      iwave_fprintall(fout);

      RVLMax<double> mx;
      d.eval(mx,x);
      
      RVLMin<double> mn;
      d.eval(mn,x);
      
      cout<<"test8 result: max="<<mx.getValue()<<" min="<<mn.getValue()<<endl;

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
