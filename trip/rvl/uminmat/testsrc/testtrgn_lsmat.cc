#include "op.hh"
#include "rnop.hh"
#include "functions.hh"
#include "TRGNAlg.hh"
#include "rnspace.hh"
#include "expls.hh"
#include "TRLSMat.hh"

using namespace RVL;
using namespace RVLAlg;
using namespace RVLUmin;

int main(int argc, char ** argv) {

  try {

    cout<<"create parameter Table - same for all subtests"<<endl;
    Table par;
    
    par.putValue("ResidualTol",0.0);          // jtol
    par.putValue("DispFlag",2);
    par.putValue("RelGradTol",0.01);          // gtol - relativize before call
    par.putValue("MaxItn",40);                // maxcount
    par.putValue("MinDecrease",0.01);         // eta1
    par.putValue("GoodDecrease",0.8);         // eta2
    par.putValue("StepDecrFactor",0.5);       // gamma1
    par.putValue("StepIncrFactor",1.8);       // gamma2
    par.putValue("TR_Delta",1.0);             // TR initial trust radius
    par.putValue("CGNE_ResTol",0.0);          // CGNE residual tolerance
    par.putValue("CGNE_GradTol",0.001);       // CGNE normal residual tolerance
    par.putValue("CGNE_MaxItn",10);           // CGNE CG step limit

    
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *            BEGIN UMin UNIT TEST 4             *"<<endl;
    cout<<" * Purpose: test TRGNAlg with LAPACK             *"<<endl;
    cout<<" * Rosenbrock function - standard problem        *"<<endl;
    cout<<" *************************************************/"<<endl; 
    cout<<endl;

    cout << endl << "Test 4.1: Rosenbrock function, double precision" <<endl;
    cout<<"tabulation in TRLAP1.rpt"<<endl<<endl;
    {
    
      int dim=2;
      RnSpace<double> sp(dim);
    
      LocalVector<double> x(sp);
      cout<<"assign initial vector per TOMS"<<endl;
      x.getData()[0]=-1.2;
      x.getData()[1]=1.0;
    
      cout<<"create GenOp\n";
      GenOp<double, rosie<double>, rosiejac<double> > f(sp,sp);
    
      cout<<"create TRGNAlg (optimization algorithm object)"<<endl;
      ofstream str("TRLAP1.rpt");
      TRGNAlg<double, LSMatPolicy<double> > umin(f,x,par,str);
 
      // ordinarily would have to assign policy parameters here, but
      // LSMatPolicy has none

      cout<<"run optimization"<<endl;
      umin.run();
    
      cout<<"iteration count = "<<umin.getCount()<<endl;
      str<<"\n\nFinal Iterate: \n";
      x.write(str);
    }

    return 0;
  }
  catch (RVLException & e) {
    e.write(cerr);
    return 1;
  }
}
