#include "gridpp_top.hh"
#include "mpigridpp.hh"
#include "gridops.hh"
#include "griddiffops.hh"

#include "adjtest.hh"

extern "C" {
#include "utils.h"
#include "parser.h"
}

using RVL::Space;
using RVL::Vector;
using RVL::StdProductSpace;
using RVL::Components;
using RVL::LinearOp;
using RVL::LinearOpFO;
using RVL::AssignFilename;
using RVL::RVLException;
using RVL::ScalarFieldTraits;
using RVL::ContentPackage;
using RVL::OperatorEvaluation;
using RVL::RVLRandomize;


#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
#else
using TSOpt::GridSpace;
#endif
using TSOpt::GridWindowFO;
using TSOpt::Grid;

using TSOpt::GridReader;
using TSOpt::GridWriter;
//using TSOpt::GridProductDiffOp;
using TSOpt::GridDiffOp;

int main(int argc, char ** argv) {

  int rk=0;

#ifdef IWAVE_USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rk);
  storeGlobalComm(MPI_COMM_WORLD);
#endif 

  try {
    
    PARARRAY par;
    ps_createargs(&par,argc,argv);
    
    char * tmp;
    
    if (ps_ffcstring(par,"in1",&tmp)) {
      RVLException e;
      e<<"Error: rwtest - failed to read key=in1 from pararray\n";
      throw e;
    }
    string ing1=tmp;
    free(tmp);

    if (ps_ffcstring(par,"in2",&tmp)) {
      RVLException e;
      e<<"Error: rwtest - failed to read key=in2 from pararray\n";
      throw e;
    }
    string ing2=tmp;
    free(tmp);
   

//     if (ps_ffcstring(par,"out1",&tmp)) {
//       RVLException e;
//       e<<"Error: rwtest - failed to read key=out1 from pararray\n";
//       throw e;
//     }
//     string outg1=tmp;
//     free(tmp);

//     if (ps_ffcstring(par,"out2",&tmp)) {
//       RVLException e;
//       e<<"Error: rwtest - failed to read key=out2 from pararray\n";
//       throw e;
//     }
//     string outg2=tmp;
//     free(tmp);

//     if (ps_ffcstring(par,"aout1",&tmp)) {
//       RVLException e;
//       e<<"Error: rwtest - failed to read key=aout1 from pararray\n";
//       throw e;
//     }
//     string aoutg1=tmp;
//     free(tmp);
 
//     if (ps_ffcstring(par,"aout2",&tmp)) {
//       RVLException e;
//       e<<"Error: rwtest - failed to read key=aout2 from pararray\n";
//       throw e;
//     }
//     string aoutg2=tmp;
//     free(tmp);
    
#ifdef IWAVE_USE_MPI
    MPIGridSpace<float> gsp1(ing1);
    MPIGridSpace<float> gsp2(ing2);
#else
    GridSpace<float> gsp1(ing1);
    GridSpace<float> gsp2(ing2);
#endif
    StdProductSpace<float> gsp(gsp1,gsp2);
//     Vector<float> xin(gsp);
//     Vector<float> xout(gsp);
//     Vector<float> xaout(gsp);

//     Components<float> cmin(xin);
//     Components<float> cmout(xout);
//     Components<float> cmaout(xaout);
 
//     /* assign files */
//     AssignFilename afin1(ing1);
//     AssignFilename afin2(ing2);
//     AssignFilename afout1(outg1);
//     AssignFilename afout2(outg2);
//     AssignFilename afaout1(aoutg1);
//     AssignFilename afaout2(aoutg2);

//     cmin[0].eval(afin1);
//     cmin[1].eval(afin2);
//     cmout[0].eval(afout1);
//     cmout[1].eval(afout2);
//     cmaout[0].eval(afaout1);
//     cmaout[1].eval(afaout2);
 
//     xout.zero();
//     xaout.zero();
  
//    GridProductDiffOp<float> gpdcdiff(gsp,2);

    GridDiffOp<float> gpdcdiff(gsp,2);

//     gpdcdiff.applyOp(xin,xout);
       
//     gpdcdiff.applyAdjOp(xin,xaout);
   
    srand(getpid()); 
    RVLRandomize<float> adc(getpid(),-1.0,1.0);
    AdjointTest(gpdcdiff,adc,cout);

    if (rk==0) iwave_fdestroy();    
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
    /* clean up */
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


    
