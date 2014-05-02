#include "acd_defn.hh"
#include "grid.h"
#include "istate.hh"
#include "iwop.hh"
#include "functions.hh"
#include "cgnealg.hh"
#include "adjtest.hh"
#ifdef IWAVE_USE_MPI
#include "mpisegypp.hh"
#include "mpigridpp.hh"
#else
#include "segypp.hh"
#include "gridpp.hh"
#endif
#include "segyops.hh"
#include "gridops.hh"


//#define GTEST_VERBOSE
IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq",    0, true,  true },
  {"data",   1, false, true },
  {"source", 1, true,  false},
  {"",       0, false, false}
};

using RVL::parse;
using RVL::RVLException;
using RVL::Vector;
using RVL::Components;
using RVL::Operator;
using RVL::OperatorEvaluation;
using RVL::LinearOp;
using RVL::LinearOpFO;
using RVL::OpComp;
using RVL::SymmetricBilinearOp;
using RVL::AssignFilename;
using RVL::AssignParams;
using RVL::RVLRandomize;
using RVL::AdjointTest;
using TSOpt::IWaveEnvironment;
using TSOpt::IWaveTree;
using TSOpt::IWaveSampler;
using TSOpt::IWaveSim;
using TSOpt::TASK_RELN;
using TSOpt::IOTask;
using TSOpt::IWaveOp;
using TSOpt::SEGYLinMute;
using TSOpt::GridWindowOp;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
using TSOpt::MPISEGYSpace;
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
#endif

using RVLUmin::CGNEAlg;

int xargc;
char **xargv;

int main(int argc, char ** argv) {

  try {

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif

    PARARRAY * pars = NULL;
    FILE * stream = NULL;
    IWaveEnvironment(argc, argv, 0, &pars, &stream);
    
    // the Op
    IWaveOp iwop(*pars,stream);
      
    // assign window widths - default = 0;
    RPNT wind;
    RASN(wind,RPNT_0);
    wind[0]=valparse<float>(*pars,"windw1",0.0f);
    wind[1]=valparse<float>(*pars,"windw2",0.0f);
    wind[2]=valparse<float>(*pars,"windw3",0.0f);
//
//    SEGYLinMute mute(valparse<float>(*pars,"mute_slope",0.0f),
//                     valparse<float>(*pars,"mute_zotime",0.0f),
//                     valparse<float>(*pars,"mute_width",0.0f));
//      
//    LinearOpFO<float> muteop(iwop.getRange(),iwop.getRange(),mute,mute);
//
//    Vector<ireal> ddin(muteop.getRange());
//    Vector<ireal> ddout(muteop.getRange());
    Vector<ireal> min(iwop.getDomain());
    Vector<ireal> mout(iwop.getDomain());

    AssignFilename minfn(valparse<std::string>(*pars,"csqin"));
    Components<ireal> cmin(min);
    cmin[0].eval(minfn);
      
    AssignFilename moutfn(valparse<std::string>(*pars,"csqout"));
    Components<ireal> cmout(mout);
    cmout[0].eval(moutfn);

    GridWindowOp wop(iwop.getDomain(),min,wind);

    RVLRandomize<float> rnd(getpid(),-1.0,1.0);
    
    OperatorEvaluation<ireal> opeval(wop,min);
    AdjointTest<float>(opeval.getDeriv(),rnd,cerr);
    
    opeval.getDeriv().applyOp(min,mout);

    
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
  
}
