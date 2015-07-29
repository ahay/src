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
    
    if (ps_ffcstring(par,"in",&tmp)) {
      RVLException e;
      e<<"Error: rwtest - failed to read key=in from pararray\n";
      throw e;
    }
    string ing=tmp;
    free(tmp);
      
//     if (ps_ffcstring(par,"out",&tmp)) {
//       RVLException e;
//       e<<"Error: rwtest - failed to read key=out from pararray\n";
//       throw e;
//     }
//     string outg=tmp;
//     free(tmp);

//     if (ps_ffcstring(par,"out2",&tmp)) {
//       RVLException e;
//       e<<"Error: rwtest - failed to read key=out from pararray\n";
//       throw e;
//     }
//     string outg2=tmp;
//     free(tmp);
    
#ifdef IWAVE_USE_MPI
    MPIGridSpace<float> gsp(ing);
#else
    GridSpace<float> gsp(ing);
#endif
    
    // Vector<float> xin(gsp);
    //    Vector<float> xout(gsp);
    //    Vector<float> xout2(gsp);    

    // AssignFilename afin(ing);
    // AssignFilename afout(outg);
    // AssignFilename afout2(outg2);   
 
    // xin.eval(afin);
    // xout.eval(afout);
    // xout2.eval(afout2);
   
    // xout.zero();
    // xout2.zero();
  
    GridDiffOp<float> gdiff(gsp,2);

    //gdiff.applyOp(xin,xout);
    //gdiff.applyAdjOp(xout,xout2);    

    srand(getpid()); 
    RVLRandomize<float> adc(getpid(),-1.0,1.0);
    AdjointTest(gdiff,adc,cout);


    if (rk==0) iwave_fdestroy();

#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif

    return(0);
  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
    exit(1);
  }
  
}


    
