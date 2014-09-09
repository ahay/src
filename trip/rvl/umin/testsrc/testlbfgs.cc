#include "op.hh"
#include "rnop.hh"
#include "ls.hh"
#include "functions.hh"
#include "LBFGSBT.hh"
#include "rnspace.hh"
#include "expls.hh"

using namespace RVL;
using namespace RVLAlg;
using namespace RVLUmin;



int main(int argc, char ** argv) {

  try {

    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *            BEGIN UMin UNIT TEST 3             *"<<endl;
    cout<<" * Purpose: test basic lbfgs alg with            *"<<endl;
    cout<<" * with backtracking line search (float)         *"<<endl;
    cout<<" *************************************************/"<<endl; 
    cout<<endl;
    cout<<"Least squares solution of 10-diml diagonal system"<<endl;
    cout<<"  y_i = z - 0.1*z^2 + 0.05*z^4, z = (i * x_i)-"<<endl;
    cout<<endl;
    cout<<"Target solution - x_i=1, i=1,...10"<<endl;
    cout<<"LBFGS iteration with backtracking line search"<<endl;
    cout<<"declare victory when gradient decreases by 0.01\n";

      cout<<"create parameter Table - same for all subtests"<<endl;
      Table par;

      par.putValue("DispFlag",2);
      par.putValue("AbsGradTol",0.0);
      par.putValue("RelGradTol",0.01);
      par.putValue("MaxItn",200);
      par.putValue("LS_MaxSample",10);
      par.putValue("LS_FirstStep",0.001);
      par.putValue("LS_MinStepTol",100*numeric_limits<float>::epsilon());
      par.putValue("LS_FractionOfMaxStep",0.9);
      par.putValue("MinDecrease",0.01);
      par.putValue("GoodDecrease",0.8);
      par.putValue("StepDecrFactor",0.5);
      par.putValue("StepIncrFactor",1.8);
      par.putValue("BFGS_InvHessianScale",1.0);
      par.putValue("BFGS_MaxUpdates",5);
    
    cout<<endl<<"Test 3.1 : dim = 2" <<endl;
    cout<<"tabulation in lbfgs1.rpt"<<endl<<endl;

    {
    
      int dim=2;
      RnSpace<float> sp(dim);
    
      Vector<float> x(sp);
      Vector<float> y(sp);
    
      cout<<"assign const to ref\n";
    
      RVLAssignConst<float>  ac(1.0);
      x.eval(ac);
    
      cout<<"create GenOp\n";
      GenOp<float, diagquart<float>, diagquartjac<float> > f(sp,sp);
    
      // in separate scope so that OpEval is destroyed
      {
	cout<<"  evaluate at initial data\n";
	OperatorEvaluation<float> feval(f,x);
      
	cout<<"  copy output to y = RHS"<<endl;
	y.copy(feval.getValue());
      }
    
      cout<<"create least squares function"<<endl;
      StdLeastSquaresFcnlGN<float> j(f,y);
    
      cout<<"zero solution vector for initial guess"<<endl;
      x.zero();
      //    x.write(cout);

      cout<<"create LBFGSBT (optimization algorithm object)"<<endl;
      ofstream str("./lbfgs1.rpt");
      //    UMinMethod<float> umin(j,x,par,str);
      LBFGSBT<float> umin(j,x,par,str);
      cout << "\n before run x.norm() = " << x.norm() << endl;
      cout<<"run optimization"<<endl;
      umin.run();
      cout << "\n end run    x.norm() = " << x.norm() << endl;
    
      cout<<"iteration count = "<<umin.getCount()<<endl;

    }

    cout << endl << "Test 3.2: dim = 5" <<endl;
    cout<<"tabulation in lbfgs2.rpt"<<endl<<endl;

    par.putValue("LS_MaxSample",20);

    {
    

      int dim=5;
      RnSpace<float> sp(dim);
    
      Vector<float> x(sp);
      Vector<float> y(sp);
    
      cout<<"assign const to ref\n";
    
      RVLAssignConst<float>  ac(1.0);
      x.eval(ac);
    
      cout<<"create GenOp\n";
      GenOp<float, diagquart<float>, diagquartjac<float> > f(sp,sp);
    
      // in separate scope so that OpEval is destroyed
      {
	cout<<"  evaluate at initial data\n";
	OperatorEvaluation<float> feval(f,x);
      
	cout<<"  copy output to y = RHS"<<endl;
	y.copy(feval.getValue());
      }
    
      cout<<"create least squares function"<<endl;
      StdLeastSquaresFcnlGN<float> j(f,y);
    
      cout<<"zero solution vector for initial guess"<<endl;
      x.zero();
      //    x.write(cout);

      cout<<"create LBFGSBT (optimization algorithm object)"<<endl;
      ofstream str("./lbfgs2.rpt");
      LBFGSBT<float> umin(j,x,par,str);
    
      cout<<"run optimization"<<endl;
      umin.run();
    
      cout<<"iteration count = "<<umin.getCount()<<endl;

    }

    cout << endl << "Test 3.3: dim = 10" <<endl;
    cout<<"tabulation in lbfgs3.rpt"<<endl<<endl;

    {
    
      int dim=10;
      RnSpace<float> sp(dim);
    
      Vector<float> x(sp);
      Vector<float> y(sp);
    
      cout<<"assign const to ref\n";
    
      RVLAssignConst<float>  ac(1.0);
      x.eval(ac);
    
      cout<<"create GenOp\n";
      GenOp<float, diagquart<float>, diagquartjac<float> > f(sp,sp);
    
      // in separate scope so that OpEval is destroyed
      {
	cout<<"  evaluate at initial data\n";
	OperatorEvaluation<float> feval(f,x);
      
	cout<<"  copy output to y = RHS"<<endl;
	y.copy(feval.getValue());
      }
    
      cout<<"create least squares function"<<endl;
      StdLeastSquaresFcnlGN<float> j(f,y);
    
      cout<<"zero solution vector for initial guess"<<endl;
      x.zero();
      //    x.write(cout);

      cout<<"create LBFGSBT (optimization algorithm object)"<<endl;
      ofstream str("./lbfgs3.rpt");
      //    UMinMethod<float> umin(j,x,par,str);
      LBFGSBT<float> umin(j,x,par,str);
    
      cout<<"run optimization"<<endl;
      umin.run();
    
      cout<<"iteration count = "<<umin.getCount()<<endl;

    }

    cout << endl << "Test 3.4: dim = 10, double precision" <<endl;
    cout<<"tabulation in lbfgs4.rpt"<<endl<<endl;

    par.putValue("LS_MinStepTol",100*numeric_limits<double>::epsilon());

    {
    
      int dim=10;
      RnSpace<double> sp(dim);
    
      Vector<double> x(sp);
      Vector<double> y(sp);
    
      cout<<"assign const to ref\n";
    
      RVLAssignConst<double>  ac(1.0);
      x.eval(ac);
    
      cout<<"create GenOp\n";
      GenOp<double, diagquart<double>, diagquartjac<double> > f(sp,sp);
    
      // in separate scope so that OpEval is destroyed
      {
	cout<<"  evaluate at initial data\n";
	OperatorEvaluation<double> feval(f,x);
      
	cout<<"  copy output to y = RHS"<<endl;
	y.copy(feval.getValue());
      }
    
      cout<<"create least squares function"<<endl;
      StdLeastSquaresFcnlGN<double> j(f,y);
    
      cout<<"zero solution vector for initial guess"<<endl;
      x.zero();
      //    x.write(cout);

      cout<<"create LBFGSBT (optimization algorithm object)"<<endl;
      ofstream str("./lbfgs4.rpt");
      //    UMinMethod<double> umin(j,x,par,str);
      LBFGSBT<double> umin(j,x,par,str);
    
      cout<<"run optimization"<<endl;
      umin.run();
    
      cout<<"iteration count = "<<umin.getCount()<<endl;

    }

    cout << endl << "Test 3.5: Rosenbrock function, double precision" <<endl;
    cout<<"tabulation in lbfgs5.rpt"<<endl<<endl;
    par.putValue("RelGradTol",0.0001);
    {
    
      int dim=2;
      RnSpace<double> sp(dim);
    
      LocalVector<double> x(sp);
      cout<<"assign initial vector per TOMS"<<endl;
      x.getData()[0]=-1.2;
      x.getData()[1]=1.0;
    
      cout<<"create GenOp\n";
      GenOp<double, rosie<double>, rosiejac<double> > f(sp,sp);
    
      cout<<"create least squares function"<<endl;
      LeastSquaresFcnlGN<double> j(f);
    
      cout<<"create LBFGSBT (optimization algorithm object)"<<endl;
      ofstream str("./lbfgs5.rpt");
      LBFGSBT<double> umin(j,x,par,str);
    
      cout<<"run optimization"<<endl;
      umin.run();
    
      cout<<"iteration count = "<<umin.getCount()<<endl;
      str<<"\n\nFinal Iterate: \n";
      x.write(str);
    }
   
    cout << endl << "Test 3.6: Rosenbrock function, single precision" <<endl;
    cout<<"tabulation in lbfgs6.rpt"<<endl<<endl;
    par.putValue("RelGradTol",0.0001);
    {
    
      int dim=2;
      RnSpace<float> sp(dim);
    
      LocalVector<float> x(sp);
      cout<<"assign initial vector per TOMS"<<endl;
      x.getData()[0]=-1.2;
      x.getData()[1]=1.0;
    
      cout<<"create GenOp\n";
      GenOp<float, rosie<float>, rosiejac<float> > f(sp,sp);
    
      cout<<"create least squares function"<<endl;
      LeastSquaresFcnlGN<float> j(f);
    
      cout<<"create LBFGSBT (optimization algorithm object)"<<endl;
      ofstream str("./lbfgs6.rpt");
      LBFGSBT<float> umin(j,x,par,str);
    
      cout<<"run optimization"<<endl;
      umin.run();
    
      cout<<"iteration count = "<<umin.getCount()<<endl;
      str<<"\n\nFinal Iterate: \n";
      x.write(str);
    }

    cout << endl << "Test 3.7: Rosenbrock function, double precision, rank 2 updates" <<endl;
    cout<<"tabulation in lbfgs7.rpt"<<endl<<endl;
    par.putValue("BFGS_MaxUpdates",2);
    {
    
      int dim=2;
      RnSpace<double> sp(dim);
    
      LocalVector<double> x(sp);
      cout<<"assign initial vector per TOMS"<<endl;
      x.getData()[0]=-1.2;
      x.getData()[1]=1.0;
    
      cout<<"create GenOp\n";
      GenOp<double, rosie<double>, rosiejac<double> > f(sp,sp);
    
      cout<<"create least squares function"<<endl;
      LeastSquaresFcnlGN<double> j(f);
    
      cout<<"create LBFGSBT (optimization algorithm object)"<<endl;
      ofstream str("./lbfgs7.rpt");
      LBFGSBT<double> umin(j,x,par,str);
    
      cout<<"run optimization"<<endl;
      umin.run();
    
      cout<<"iteration count = "<<umin.getCount()<<endl;
      str<<"\n\nFinal Iterate: \n";
      x.write(str);
    }
   
    cout << endl << "Test 3.8: Rosenbrock function, double precision, steepest descent" <<endl;
    cout<<"tabulation in lbfgs8.rpt"<<endl<<endl;
    par.putValue("BFGS_MaxUpdates",0);
    par.putValue("MaxItn",3000);
    {
    
      int dim=2;
      RnSpace<double> sp(dim);
    
      LocalVector<double> x(sp);
      cout<<"assign initial vector per TOMS"<<endl;
      x.getData()[0]=-1.2;
      x.getData()[1]=1.0;
    
      cout<<"create GenOp\n";
      GenOp<double, rosie<double>, rosiejac<double> > f(sp,sp);
    
      cout<<"create least squares function"<<endl;
      LeastSquaresFcnlGN<double> j(f);
    
      cout<<"create LBFGSBT (optimization algorithm object)"<<endl;
      ofstream str("./lbfgs8.rpt");
      LBFGSBT<double> umin(j,x,par,str);
    
      cout<<"run optimization"<<endl;
      umin.run();
    
      cout<<"iteration count = "<<umin.getCount()<<endl;
      str<<"\n\nFinal Iterate: \n";
      x.write(str);
    }
   
    return(0);
   
 
  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
