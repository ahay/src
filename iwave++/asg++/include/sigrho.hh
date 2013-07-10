/* Rami here is a possible solution, for the
   (sigma,rho)-to-(kappa,buoy) case. The idea is to treat it like
   InjectOp (which is completely generic, hence part of RVL) but make
   the method implementations particular to this case. 
*/

#include "fcnfo.hh"
//#include "sigrho_to_kappabuoy.hh"
#include "op.hh"

#define VELRHO
//#define SIGRHO
//#define VELONLY

#define rho_const 1000.0

#ifdef VELONLY
float kappa(float sig,float rho){
  return sig*sig*rho_const;
}

float buoy(float sig,float rho){
  return 1/rho_const;
}

float dkappadsig(float sig, float rho, float deltasig){
  return 2*sig*deltasig*rho_const;
}

float dkappadrho(float sig,float rho, float deltarho){
  return 0;
}

float dbuoydrho(float sig,float rho, float deltarho){
  return 0;
}

float dbuoydsig(float sig,float rho, float deltasig){
  return 0;
}
#else
#ifdef VELRHO
float kappa(float sig,float rho){
  return sig*sig*rho;
}

float buoy(float sig,float rho){
  return 1/rho;
}

float dkappadsig(float sig, float rho, float deltasig){
  return 2*sig*deltasig*rho;
}

float dkappadrho(float sig,float rho, float deltarho){
  return sig*sig*deltarho;
}

float dbuoydrho(float sig,float rho, float deltarho){
  return -deltarho/(rho*rho);
}

float dbuoydsig(float sig,float rho, float deltasig){
  return 0;
}
#else

float kappa(float sig,float rho){
  return sig*sig/rho;
}

float buoy(float sig,float rho){
  return 1/rho;
}

float dkappadsig(float sig, float rho, float deltasig){
  return 2*sig*deltasig/rho;
}

float dkappadrho(float sig,float rho, float deltarho){
  return -sig*sig*deltarho/(rho*rho);
}

float dbuoydrho(float sig,float rho, float deltarho){
  return -deltarho/(rho*rho);
}

float dbuoydsig(float sig,float rho, float deltasig){
  return 0;
}

#endif
#endif

using namespace RVL;

namespace ASG {

  /* lets just go right ahead and write the op we want from scratch,
     more or less - I started by copying InjectOp from blockop.hh and
     changing the names, adding all those ScalarFOs from asginv_rami
     as private members - might as well be private, will never expose
     them outside. Also it's only going to be used for floats so
     changed that too.
  */
  class SigRhoOp: public Operator<float> {
    
  private:

    Space<float> const & dom;
    Space<float> const & rng;

    // these FOs must be declared mutable because in principle
    // eval of an FO may change its internal state (does not happen
    // in this case, of course)
    mutable ScalarFO2<float,kappa> kappa_sigrho;
    mutable ScalarFO2<float,buoy> buoy_sigrho;

    mutable ScalarFO3<float,dkappadsig> dkappa_dsig;
    mutable ScalarFO3<float,dkappadrho> dkappa_drho;
    mutable ScalarFO3<float,dbuoydrho> dbuoy_drho;
    mutable ScalarFO3<float,dbuoydsig> dbuoy_dsig;
    
  protected:
    
    void apply(Vector<float> const & x,
	       Vector<float> & y) const;

    void applyDeriv(Vector<float> const & x,
		    Vector<float> const & dx,
		    Vector<float> & dy) const;

    void applyAdjDeriv(Vector<float> const & x,
		       Vector<float> const & dy,
		       Vector<float> & dx) const;

    Operator<float> * clone() const { return new SigRhoOp(*this); }

  public:

    // the constructor is considerably different. Basically you give
    // it domain and range spaces, rather than construct them. The
    // spaces contain information about the types of data
    // (ie. essentially unit information), so you can't interchange
    // (sig,rho) with (kappa,buoy) - have to treat them as different.

    SigRhoOp(Space<float> const & _dom,
	     Space<float> const & _rng);

    SigRhoOp(SigRhoOp const & f)
      : dom(f.dom), rng(f.rng) {}

    ~SigRhoOp() {}

    Space<float> const & getDomain() const { return dom; }
    Space<float> const & getRange() const { return rng; }

    ostream & write(ostream & str) const {
      str<<"SigRhoOp object\n";
      // you should probably invent some more interesting stuff to say
      // about this op
      return str;
    }
  };
}
  

