#include "usempi.h"
#include "rkstream.hh"
#include "gridpp_top.hh"
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
      y.getData()[i] = 2.0*(1.0+(float)n)*x.getData()[i];
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
      dy.getData()[i] = 2.0*(1.0+(float)n)*dx.getData()[i];
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

    if (rk==0) {

      ofstream str;
      makeRankStream(str,rk,"testsrc/test16");
      FILE * fout = fopen("testsrc/test16/stdout0.txt","w");
      
      string fname="testsrc/test16/testgrid.rsf";
      
      cout<<"GRIDPP Unit Test 16"<<endl;
      cout<<"Least squares solution of square linear system"<<endl;
      cout<<"  y_i = i * x_i"<<endl;
      cout<<"Operator implemented as OpFO, using FOs at top of file"<<endl;
      cout<<"Target solution - x_i=1, i=1,...n"<<endl;
      cout<<"LBFGS iteration with backtracking line search"<<endl;
      cout<<"float version, solution in temp file so not saved"<<endl<<endl;
      
      create_hfile(fname,0);

      GridSpace<float> sp(fname,"notype",str);
    
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
      LBFGSBT<float> umin(j,x,par,cout);
    
      str<<"run optimization"<<endl;
      umin.run();
    
      str<<"compute, print max, min - should be close to 1"<<endl;
    
      RVLMax<float> mx;
      x.eval(mx);
    
      RVLMin<float> mn;
      x.eval(mn);
    
      cout<<"test16 result: max="<<mx.getValue()<<" min="<<mn.getValue()<<endl;

      fprintf(fout,"file system after LBFGS run\n");
      iwave_fprintall(fout);

      str.close();
      fclose(fout);

      iwave_fdestroy();
    }
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
