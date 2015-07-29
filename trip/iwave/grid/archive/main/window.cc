#include "gridpp_top.hh"
#include "mpigridpp.hh"
#include "gridops.hh"

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

#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
#else
using TSOpt::GridSpace;
#endif
using TSOpt::GridWindowFO;
using TSOpt::Grid;

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
    
    ps_printall(par,stdout);

    char * tmp;
    
    if (ps_ffcstring(par,"in",&tmp)) {
      RVLException e;
      e<<"Error: window - failed to read key=in from pararray\n";
      throw e;
    }
    string ing=tmp;
    free(tmp);
      
    if (ps_ffcstring(par,"out",&tmp)) {
      RVLException e;
      e<<"Error: window - failed to read key=out from pararray\n";
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

    float w=ScalarFieldTraits<float>::Zero();
    ps_fffloat(par,"width",&w);

    Grid<float> const & g = gsp.getGrid();
    std::vector<float> a(g.getDim());
    std::vector<float> b(g.getDim());
    
    for (int i=0;i<RARR_MAX_NDIM;i++) {
      a[i]=g.getData()[i].o;
      b[i]=g.getData()[i].o + g.getData()[i].d*(g.getData()[i].n-1);
      stringstream l;
      stringstream r;
      l<<"l"<<i+1;
      r<<"r"<<i+1;
      ps_fffloat(par,l.str().c_str(),&(a[i]));
      ps_fffloat(par,r.str().c_str(),&(b[i]));
    }

    GridWindowFO<float> f(a,b,w);

    LinearOpFO<float> op(gsp,gsp,f,f);

    xout.eval(f,xin);

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


    
