#ifndef __RVL_GRIDMUL
#define __RVL_GRIDMUL

#include "gridwindow.hh"
#include "griddata.hh"
#include "cubicfo.hh"

namespace RVLGrid {

  using namespace RVL;

  class mfit: public UnaryLocalFunctionObjectScalarRedn<float> {
    AxisSegs ax;
    GridData const & din;
    GridData const & dout;
    float s;
    float lam;
  public: 
    mfit(AxisSegs & _ax, 
	 GridData const & _din,
	 GridData const & _dout,
	 float _s, float _lam)
      : ax(_ax), din(_din), dout(_dout), s(_s), lam(_lam) {}
    mfit(mfit const & m)
      : ax(m.ax), din(m.din), dout(m.dout), s(m.s), lam(m.lam) {}
    ~mfit() {}

    void operator()(LocalDataContainer<float> const & in);

    std::string getName() const { std::string tmp="mfit"; return tmp; }
  };

  class gfit: public BinaryLocalFunctionObject<float> {
    AxisSegs ax;
    GridData const & din;
    GridData const & dout;
    float s;
    float lam;
  public: 
    gfit(AxisSegs & _ax, 
	 GridData const & _din,
	 GridData const & _dout,
	 float _s, float _lam)
      : ax(_ax), din(_din), dout(_dout), s(_s), lam(_lam) {}
    gfit(gfit const & m)
      : ax(m.ax), din(m.din), dout(m.dout), s(m.s), lam(m.lam) {}
    ~gfit() {}

    void operator()(LocalDataContainer<float> & grad,
		    LocalDataContainer<float> const & in);

    std::string getName() const { std::string tmp="gfit"; return tmp; }
  };

  class hfit: public TernaryLocalFunctionObject<float> {
  public: 
    hfit() {}
    hfit(hfit const &) {}
    ~hfit() {}

    void operator()(LocalDataContainer<float> & out,
		    LocalDataContainer<float> const & refin,
		    LocalDataContainer<float> const & pertin) {
      RVLException e;
      e<<"Error: hfit::operator() not implemented yet\n";
      throw e;
    }

    std::string getName() const { std::string tmp="hfit"; return tmp; }
  };


  template<typename Scalar>
  class mul: public TernaryLocalFunctionObject<Scalar> {
  private:
    Scalar s;
  public: 
    mul(Scalar _s = ScalarFieldTraits<Scalar>::One()):s(_s) {}
    mul(mul<Scalar> const & m):s(m.s) {}
    ~mul() {}

    void operator()(LocalDataContainer<float> & out,
		    LocalDataContainer<float> const & in1,
		    LocalDataContainer<float> const & in2) {
      int n = min(out.getSize(),in1.getSize());
      n = min(n,in2.getSize());
      float * p = out.getData();
      float const * q1 = in1.getData();
      float const * q2 = in2.getData();
      for (int i=0;i<n;i++) p[i]=s*q1[i]*q2[i];
    }

    std::string getName() const { std::string tmp="mul"; return tmp; }
  };

  template<typename Scalar>
  class expo: public UnaryLocalFunctionObject<Scalar> {
  private:
    Scalar s,t;
  public:
    expo(Scalar _s = ScalarFieldTraits<Scalar>::One(),
	 Scalar _t = ScalarFieldTraits<Scalar>::One()):s(_s),t(_t) {}
    expo(expo<Scalar> const & e):s(e.s),t(e.t) {}
    ~expo() {}

    void operator()(LocalDataContainer<Scalar> & x) {
      int n=x.getSize();
      Scalar * p = x.getData();
      for (int i=0;i<n;i++) {
	p[i]=s*exp(t*p[i]);
      }
    }

    std::string getName() const { std::string tmp="expo"; return tmp; }
  };
}

#endif
