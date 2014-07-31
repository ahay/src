#include "parser.h"
#ifdef IWAVE_USE_MPI
#include "mpisegypp.hh"
#include "mpiserialfo.hh"
#else
#include "segypp.hh"
#endif
#include "segyops.hh"
#include "op.hh"

using RVL::valparse;
using RVL::RVLException;
using RVL::Vector;
using RVL::Components;
using RVL::LinearOpFO;
using RVL::AssignFilename;

using TSOpt::SEGYLinMute;

#ifdef IWAVE_USE_MPI
using RVL::MPISerialFunctionObject;
using TSOpt::MPISEGYSpace;
typedef TSOpt::MPISEGYSpace tsp;
#else
using TSOpt::SEGYSpace;
typedef TSOpt::SEGYSpace tsp;
#endif

int xargc;
char **xargv;

int main(int argc, char ** argv) {

  try {

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif

    PARARRAY * pars = ps_new();
    
    if (ps_createargs(pars,argc-1,&(argv[1]))) {
      RVLException e;
      e<<"ERROR: mute from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }
    // since the product of grid spaces is not really an 
    // out-of-core structure, this driver operates on single
    // grid spaces
    string inp = valparse<string>(*pars,"input");
    string outp = valparse<string>(*pars,"output");

    tsp dom(inp,"data");
        
    Vector<ireal> ddin(dom);
    Vector<ireal> ddout(dom);
        
    AssignFilename ddinfn(inp);
    Components<ireal> cddin(ddin);
    cddin[0].eval(ddinfn);
        
    AssignFilename ddoutfn(outp);
    Components<ireal> cddout(ddout);
    cddout[0].eval(ddoutfn);

    SEGYLinMute mute(valparse<float>(*pars,"mute_slope",0.0f),
		     valparse<float>(*pars,"mute_zotime",0.0f),
		     valparse<float>(*pars,"mute_width",0.0f));
#ifdef IWAVE_USE_MPI
    MPISerialFunctionObject<float> mpimute(mute);
    LinearOpFO<float> muteop(dom,dom,mpimute,mpimute);
    muteop.applyOp(ddin,ddout);
    MPI_Finalize();
#else
    LinearOpFO<float> muteop(dom,dom,mute,mute);
    muteop.applyOp(ddin,ddout);
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
