#include "gridpp_top.hh"
#include "mpigridpp.hh"
#include "gridops.hh"
#include "gridrwfos.hh"
#include "mpiserialfo.hh"

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
using RVL::MPISerialFunctionObjectRedn;

#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
#else
using TSOpt::GridSpace;
#endif
using TSOpt::GridWindowFO;
using TSOpt::Grid;

using TSOpt::GridReader;
using TSOpt::GridWriter;


int main(int argc, char ** argv) {
  try {

    int rk=0;

#ifdef IWAVE_USE_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);
    storeGlobalComm(MPI_COMM_WORLD);
    {
#endif 
      printf("rk = %d \n",rk);
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
      
      if (ps_ffcstring(par,"out",&tmp)) {
	RVLException e;
	e<<"Error: rwtest - failed to read key=out from pararray\n";
	throw e;
      }
      string outg=tmp;
      free(tmp);
      
#ifdef IWAVE_USE_MPI
      MPIGridSpace<float> gsp(ing);
#else
      GridSpace<float> gsp(ing);
#endif
      
      Vector<float> xin(gsp);
      Vector<float> xout(gsp);
      
      AssignFilename afin(ing);
      AssignFilename afout(outg);
      
      xin.eval(afin);
      xout.eval(afout);
      xout.zero();
      
      ContentPackage<float, Grid<float> > cpbuf;
      Grid<float> g(gsp.getGrid()); 
      
      // set up a subgrid
      g.getData()[0].n=100;
      g.getData()[1].n=100;
      g.getExdAxis().n=2;

      g.fprint(stderr);

      cpbuf.initialize(g);
      
      GridReader<float> r(cpbuf);
      xin.eval(r);
      
      GridWriter<float> w(cpbuf);
      xout.eval(w);
      
      if (rk==0) iwave_fdestroy();
      
#ifdef IWAVE_USE_MPI
    }
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


    
