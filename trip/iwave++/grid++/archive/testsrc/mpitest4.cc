#include "rkstream.hh"
#include "gridpp_top.hh"
#include "mpigridpp.hh"
#include "op.hh"
#include "ls.hh"
#include "functions.hh"
#include "LBFGSBT.hh"
extern "C" {
#include "iwave_fopen.h"
}
#include "create_hfile.hh"

using namespace RVL;
using namespace TSOpt;
using namespace RVLAlg;
using namespace RVLUmin;

class tfun: public BinaryLocalFunctionObject<float> {

public:

  using RVL::BinaryLocalEvaluation<float>::operator();
  void operator()(LocalDataContainer<float> & y,
		  LocalDataContainer<float> const & x) {
    int n=min(x.getSize(),y.getSize());
    for (int i=0;i<n;i++) {
      y.getData()[i] = (1.0+(float)n)*x.getData()[i];
    }
  }
  string getName() const { string tmp="tfun"; return tmp; }
};

class dtfun: public TernaryLocalFunctionObject<float> {

public:

  using RVL::TernaryLocalEvaluation<float>::operator();
  void operator()(LocalDataContainer<float> & dy,
		  LocalDataContainer<float> const & x,
		  LocalDataContainer<float> const & dx) {
    int n=min(x.getSize(),dy.getSize());
    n=min(n,dx.getSize());
    for (int i=0;i<n;i++) 
      dy.getData()[i] = (1.0+(float)n)*dx.getData()[i];
  }
  string getName() const { string tmp="dtfun"; return tmp; }
};

int main(int argc, char ** argv) {

  int rk=0;

#ifdef IWAVE_USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rk);
#endif 

  try {

    ofstream str;
    makeRankStream(str,rk,"testsrc/mpitest4");
      
    string fname="testsrc/mpitest4/testgrid.rsf";

    if (rk==0) {

      cout<<"MPI GRIDPP Unit Test 4"<<endl;
      cout<<"Least squares solution of square system"<<endl;
      cout<<"  y_i = i * x_i"<<endl;
      cout<<"Operator implemented as OpFO, using FOs at top of file"<<endl;
      cout<<"Target solution - x_i=1, i=1,...n"<<endl;
      cout<<"LBFGS iteration with backtracking line search"<<endl;
      cout<<"float version, solution in temp file so not saved"<<endl<<endl;
      
      create_hfile(fname,0);
    }

#ifdef IWAVE_USE_MPI    
    MPIGridSpace<float> sp(fname,MPI_COMM_WORLD,str);
#else
    GridSpace<float> sp(fname,"notype",str);
#endif
    
    Vector<float> x(sp);
    Vector<float> y(sp);
      
    str<<"assign const to ref\n";

    RVLAssignConst<float>  ac(1.0);
    x.eval(ac);
    
    str<<"declare FOs, create OpFO, opeval\n";
    tfun tf;
    dtfun dtf;
    OpFO<float> f(sp,sp,tf,dtf,dtf);
    
    // in separate scope so that OpEval is destroyed
    {
      OperatorEvaluation<float> feval(f,x);
      
      str<<"copy output to y"<<endl;
      y.copy(feval.getValue());
    }
    
    str<<"create least squares function"<<endl;
    StdLeastSquaresFcnlGN<float> j(f,y);
    
    str<<"create UMinTable"<<endl;
    Table par;
    par.putValue("DispFlag",2);
    par.putValue("AbsGradTol",0.0);
    par.putValue("RelGradTol",0.01);
    par.putValue("MaxItn",20);
    par.putValue("LS_MaxSample",10);
    par.putValue("LS_FirstStep",0.1);
    par.putValue("LS_MinStepTol",100*numeric_limits<float>::epsilon());
    par.putValue("LS_FractionOfMaxStep",0.9);
    par.putValue("MinDecrease",0.01);
    par.putValue("GoodDecrease",0.8);
    par.putValue("StepDecrFactor",0.5);
    par.putValue("StepIncrFactor",1.8);
    par.putValue("BFGS_InvHessianScale",1.0);
    par.putValue("BFGS_MaxUpdates",5);
    
    str<<"zero solution vector for initial guess"<<endl;
    x.zero();
    
    str<<"create UMinMethod (optimization algorithm object)"<<endl;
    ostream * outstr = &str;
    if (rk==0) outstr = &cout;
    LBFGSBT<float> umin(j,x,par,*outstr);

    str<<"run optimization"<<endl;
    umin.run();
    
    str<<"compute, print max, min - should be close to 1"<<endl;
    
    RVLMax<float> mx;
    MPISerialFunctionObjectRedn<float,float> mpimx(mx);
    x.eval(mpimx);
      
    RVLMin<float> mn;
    MPISerialFunctionObjectRedn<float,float> mpimn(mn);
    x.eval(mpimn);
    
    if (rk==0) cout<<"mpitest2 result: max="<<mpimx.getValue()<<" min="<<mpimn.getValue()<<endl;

    str.close();
    if (rk==0) iwave_fdestroy();

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
