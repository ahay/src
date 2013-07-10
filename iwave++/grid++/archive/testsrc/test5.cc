#include "usempi.h"
#include "rkstream.hh"
#include "gridpp_top.hh"
#include "create_hfile.hh"

using namespace RVL;
using namespace TSOpt;
using namespace std;

int main(int argc, char ** argv) {

  int rk=0;

#ifdef IWAVE_USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rk);
#endif 

  try {

    if (rk==0) {

      ofstream str;
      makeRankStream(str,rk,"testsrc/test5");
      FILE * fout = fopen("testsrc/test5/cout0.txt","r+");
      fseek(fout,0L,SEEK_END);

      cout<<"GRIDPP Unit Test 5"<<endl;
      cout<<"create GridDCF, GridDC, tmp file not unlinked, assign"<<endl;
      cout<<"5.0 to all data samples, compute max and min"<<endl;

      string hname="testsrc/test5/testgrid.rsf";

      create_hfile(hname,0);

      // dummy vector of dc's
      vector<DataContainer const *> x(0);

      GridDC<float> d(hname,str);

      fprintf(fout,"file status before eval of FO:\n");
      iwave_fprintall(fout);

      RVLAssignConst<float>  ac(5.0);
      d.eval(ac,x);

      fprintf(fout,"file status after eval of FO:\n");
      iwave_fprintall(fout);

      RVLMax<float> mx;
      d.eval(mx,x);

      RVLMin<float> mn;
      d.eval(mn,x);
  
      cout<<"test5 result: max="<<mx.getValue()<<" min="<<mn.getValue()<<endl;
      cout<<"temp files not unlinked\n";
      cout<<endl;

      str.close();
      fclose(fout);
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
