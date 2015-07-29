#include "parser.h"
#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#else
#include "gridpp.hh"
#endif
#include "helmfftw.hh"

using RVL::valparse;
using RVL::RVLException;
using RVL::Vector;
using RVL::LinearOp;
using RVL::AssignFilename;

using TSOpt::GridHelmFFTWOp;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
typedef TSOpt::MPIGridSpace gsp;
#else
using TSOpt::GridSpace;
typedef TSOpt::GridSpace gsp;
#endif

int xargc;
char **xargv;

int main(int argc, char ** argv) {

  try {
#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);   
    storeGlobalComm(MPI_COMM_WORLD);
#endif

    PARARRAY * pars = ps_new();

    if (ps_createargs(pars,argc-1,&(argv[1]))) {
      RVLException e;
      e<<"ERROR: GridDerivOp from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }
    // since the product of grid spaces is not really an 
    // out-of-core structure, this driver operates on single
    // grid spaces
    string inp = valparse<string>(*pars,"in");
    string outp = valparse<string>(*pars,"out");
    float power = valparse<float>(*pars,"power");
    float datum = valparse<float>(*pars,"datum", 0.0f);
    RPNT weights;
      IPNT sbc;
      IPNT ebc;
    for (int i=0;i<RARR_MAX_NDIM;i++) {
      stringstream msg;
      msg<<"weight"<<i;
      weights[i]=valparse<float>(*pars,msg.str(),1.0f);
      stringstream smsg;
      smsg<<"sbc"<<i;
      sbc[i]=valparse<int>(*pars,smsg.str(),0);
      stringstream emsg;
      emsg<<"ebc"<<i;
      ebc[i]=valparse<int>(*pars,emsg.str(),0);
        //cerr << "sbc[" << i << "] = " << sbc[i] << endl;
        //cerr << "ebc[" << i << "] = " << ebc[i] << endl;
    }
    gsp sp(inp, "notype", true
#ifdef IWAVE_USE_MPI 
	   , retrieveGlobalComm()
#endif
	   );
    GridHelmFFTWOp op(sp,weights,sbc,ebc,power,datum);
    Vector<float> invec(sp);
    Vector<float> outvec(sp);
    AssignFilename afin(inp);
    AssignFilename afout(outp);
    invec.eval(afin);
    outvec.eval(afout);
    if (valparse<bool>(*pars,"adjoint",false)) {
      op.applyAdjOp(invec,outvec);
    }
    else {
      op.applyOp(invec,outvec);
    }
    ps_delete(&pars);
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


