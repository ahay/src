#include "parser.h"
#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#else
#include "gridpp.hh"
#endif
#include "helmfftw.hh"
#include "adjtest.hh"

using RVL::valparse;
using RVL::RVLException;
using RVL::Vector;
using RVL::LinearOp;
using RVL::AssignFilename;
using RVL::AdjointTest;


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
    string inpx = valparse<string>(*pars,"inx");
    string outpx = valparse<string>(*pars,"outx");
    string inpy = valparse<string>(*pars,"iny");
    string outpy = valparse<string>(*pars,"outy");
      
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
      sbc[i]=valparse<int>(*pars,smsg.str(),1);
      stringstream emsg;
      emsg<<"ebc"<<i;
      ebc[i]=valparse<int>(*pars,emsg.str(),1);
        //cerr << "sbc[" << i << "] = " << sbc[i] << endl;
        //cerr << "ebc[" << i << "] = " << ebc[i] << endl;
    }
    gsp sp(inpx, "notype", false
#ifdef IWAVE_USE_MPI 
	   , retrieveGlobalComm()
#endif
	   );
    GridHelmFFTWOp op(sp,weights,sbc,ebc,power,datum);
    Vector<float> inx(sp);
    Vector<float> outx(sp);
    Vector<float> iny(sp);
    Vector<float> outy(sp);

    AssignFilename afinx(inpx);
    AssignFilename afoutx(outpx);
    inx.eval(afinx);
    outx.eval(afoutx);
    AssignFilename afiny(inpy);
    AssignFilename afouty(outpy);
    iny.eval(afiny);
    outy.eval(afouty);
    
    op.applyOp(inx,outx);
    op.applyAdjOp(iny,outy);
      
    float axnm = outx.norm();
    float ynm  = iny.norm();
    float axy  = iny.inner(outx);
    float xaty = inx.inner(outy);
      
    cerr << "<Ax,    y> = " << axy << endl;
    cerr << "< x, A^Ty> = " << xaty << endl;
    cerr << " |Ax| * |y| = " << axnm * ynm << endl;
    cerr << " adjvalue = " << abs(axy - xaty)/(axnm*ynm) << endl;
      
    GridHelmFFTWOp op1(sp,weights,sbc,ebc,power,datum);
    int seed = 10001; //getpid();
    RVL::RVLRandomize<float> rnd(seed,-1.0f,1.0f);
    AdjointTest(op1,rnd,cerr);
      

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


