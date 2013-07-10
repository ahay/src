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

  try {

#ifdef IWAVE_USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rk);
#endif 

    if (rk==0) {

      ofstream str;
      makeRankStream(str,rk,"testsrc/test7");
      FILE * fout = fopen("testsrc/test7/cout0.txt","r+");
      fseek(fout,0L,SEEK_END);

      cout<<"GRIDPP Unit Test 7"<<endl;
      cout<<"create GridDCF, GridDC, perm new file=test7.rsf, assign"<<endl;
      cout<<"7.0 to all data samples, compute max and min"<<endl;

      /* write header file */
      string hname="testsrc/test7/testgrid.rsf";
      create_hfile(hname,0);
      
      GridDC<float> d(hname,str);

      vector<DataContainer const *> x(0);

      //    cerr<<"rk="<<rk<<" eval cc"<<endl;
      string dname="testsrc/test7/test7.rsf";
      AssignFilename af(dname);

      d.eval(af,x);

      fprintf(fout,"file system before assignconst eval\n");
      iwave_fprintall(fout);

      //    cerr<<"rk="<<rk<<" eval ac"<<endl;
      RVLAssignConst<float>  ac(7.0);
      d.eval(ac,x);

      fprintf(fout,"file system after assignconst eval\n");
      iwave_fprintall(fout);

      //    cerr<<"rk="<<rk<<" eval mx"<<endl;
      RVLMax<float> mx;
      d.eval(mx,x);
      //    cerr<<"rk="<<rk<<" eval mn"<<endl;
      RVLMin<float> mn;
      d.eval(mn,x);
      
      cout<<"test7 result: max="<<mx.getValue()<<" min="<<mn.getValue()<<endl;

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
