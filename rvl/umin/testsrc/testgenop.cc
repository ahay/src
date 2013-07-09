#include "rnop.hh"

using namespace RVL;

// Rosenbrock function
template<typename T>
int rosie(int m, int n, T const * x, T * y) {
  if ((m!=2) || (n!=2)) return 1;
  y[0]=10*(x[1]-x[0]*x[0]);
  y[1]=ScalarFieldTraits<T>::One()-x[0];
  return 0;
}

// Rosenbrock jacobian
template<typename T>
int rosiejac(int m, int n, T const * x, T * y) {
  if ((m!=2) || (n!=2)) return 1;
  // df[0]/dx[0]
  y[0]=-20*x[0];
  // df[1]/dx[0]
  y[1]=-ScalarFieldTraits<T>::One();
  // df[0]/dx[1]
  y[2]=10*ScalarFieldTraits<T>::One();
  // df[1]/dx[1]
  y[3]=ScalarFieldTraits<T>::Zero();
  return 0;
}

int main() {

  cout<<"******************************************************"<<endl;
  cout<<"***       Test of GenOp Rn operator class          ***"<<endl;
  cout<<"*** Rosenbrock function (Mor\'{e} et all TOMS 566) ***"<<endl;
  cout<<"***                                                ***"<<endl;
  cout<<"*** m=2 n=2                                        ***"<<endl;
  cout<<"***                                                ***"<<endl;
  cout<<"*** initial point from TOMS:                       ***"<<endl;
  cout<<"******************************************************"<<endl;
  // space, operator
  int m=2;
  RnSpace<float> dom(m);
  GenOp<float, rosie<float>, rosiejac<float> > op(dom,dom);
  
  // vector workspace
  LocalVector<float> x(dom);
  LocalVector<float> y(dom);

  // assign values to x - increments version, makes it "fresh"
  x.getData()[0]=-1.2f;
  x.getData()[1]=1.0f;

  // evaluate operator at x
  OperatorEvaluation<float> opeval(op,x);

  // extract the image of x under op (which causes it to be computed)
  y.copy(opeval.getValue());

  // report
  cout<<"\ninput vector:\n";
  x.write(cout);
  cout<<"\noutput vector:\n";
  y.write(cout);
  cout<<"\nJacobian:\n";
 
  // a bit complicated to pull the Jacobian data out, but it's
  // accessible. First, grab the Operator copy being used with
  // Operator::getOp. It's a GenOp, so cast it, making the additional
  // GenOp attributes available. One of these is the internal GenMat
  // used to store the Jacobian. The Jacobian has been initialized
  // because the operator has been applied to a vector, so you can
  // access the GenMat, and then its internal float * data.

  GenOp<float, rosie, rosiejac> const & currop 
    = dynamic_cast<GenOp<float, rosie, rosiejac> const &>(opeval.getOp());
  float const * j = currop.getGenMat().getData();
  cout<<"j[0,0]="<<j[0]<<endl;
  cout<<"j[1,0]="<<j[1]<<endl;
  cout<<"j[0,1]="<<j[2]<<endl;
  cout<<"j[1,1]="<<j[3]<<endl<<endl;  

  cout<<"******************************************************"<<endl;
  cout<<"*** global solution for LS problem from TOMS:      ***"<<endl;
  cout<<"******************************************************"<<endl;

  // non-const access to data causes version to increment, hence all opeval 
  // attributes depending on the eval point (x) will be recomputed.
  x.getData()[0]=1.0f;
  x.getData()[1]=1.0f;

  // invoking getValue recomputes everything
  y.copy(opeval.getValue());

  // report again
  cout<<"\ninput vector:\n";
  x.write(cout);
  cout<<"\noutput vector:\n";
  y.write(cout);
  cout<<"\nJacobian:\n";
  cout<<"j[0,0]="<<j[0]<<endl;
  cout<<"j[1,0]="<<j[1]<<endl;
  cout<<"j[0,1]="<<j[2]<<endl;
  cout<<"j[1,1]="<<j[3]<<endl;
}

  
