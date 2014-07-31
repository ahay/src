#include "parser.h"
#include "gridpp.hh"
#include "gridops.hh"

using RVL::valparse;
using RVL::RVLException;
using RVL::Vector;
using RVL::LinearOp;
using RVL::AssignFilename;

using TSOpt::GridDerivOp;
using TSOpt::GridSpace;
typedef TSOpt::GridSpace gsp;

int xargc;
char **xargv;

int main(int argc, char ** argv) {

  try {

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
    string inp = valparse<string>(*pars,"input");
    string outp = valparse<string>(*pars,"output");
    int axis = valparse<int>(*pars,"axis");

    ps_delete(&pars);

    gsp sp(inp,"notype",true);
    
    GridDerivOp op(sp,axis);

    Vector<float> invec(sp);
    Vector<float> outvec(sp);
    AssignFilename afin(inp);
    AssignFilename afout(outp);
    invec.eval(afin);
    outvec.eval(afout);

    op.applyOp(invec,outvec);

  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
  
}


