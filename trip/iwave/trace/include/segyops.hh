#ifndef __TSOPT__SEGY__OPS__
#define __TSOPT__SEGY__OPS__

#include "segypp.hh"
#include "local.hh"

namespace TSOpt {

  using RVL::BinaryLocalFunctionObject;
  using RVL::LocalDataContainer;

  class SEGYLinMute: public BinaryLocalFunctionObject<float> {
    
  private:
    float s;   // type 0: slope of mute (dt/dx)
               // type 1: minimum value of gx where taper starts (amp=0.)

    float tm;  // type 0: mute onset at zero offset, time AFTER first sample
               // type 1: maximum value of gx where taper starts (amp=0.)    

    float w;   // width of mute zone
   
    int mute_type; // muting type: 0, conic-like (point source); 1, rectangle (plane-wave src)
    
    float gxbeg;  
    float gxend; 
  public:
    
    SEGYLinMute(float _s=0.0f, float _tm=0.0f, float _w=0.0f, int _type = 0) 
      : s(_s),tm(_tm),w(_w),mute_type(_type) {}
    
    SEGYLinMute(SEGYLinMute const & m) 
      : s(m.s),tm(m.tm),w(m.w),mute_type(m.mute_type),gxbeg(m.gxbeg),gxend(m.gxend) {}
    void set(float _s, float _tm, float _w, int _type = 0){
      mute_type = _type; 
      if (!mute_type) {
	s = _s;
	tm = _tm;
	w = _w;
      }
      else {
	s = _s; //gxbeg
	tm = _tm; //gxend
	w = _w;
      }
    }
    
    ~SEGYLinMute() {}
    
    using RVL::LocalEvaluation<float>::operator();    
    void operator()(LocalDataContainer<float> &,
		    LocalDataContainer<float> const &);
    
    string getName() const { string tmp = "SEGYLinMute"; return tmp; }
    
  };

  class SEGYFwdInt: public BinaryLocalFunctionObject<float> {
    
  private:
    
    int nint;  // number of trace integrations
    
  public:
    
    SEGYFwdInt(int _nint): nint(_nint) {}
    
    SEGYFwdInt(SEGYFwdInt const & m): nint(m.nint) {}
    
    ~SEGYFwdInt() {}

    using RVL::LocalEvaluation<float>::operator();    
    void operator()(LocalDataContainer<float> &,
		    LocalDataContainer<float> const &);
    
    string getName() const { string tmp = "SEGYFwdInt"; return tmp; }
    
  };

  class SEGYAdjInt: public BinaryLocalFunctionObject<float> {
    
  private:
    
    int nint;  // number of trace integrations
    
  public:
    
    SEGYAdjInt(int _nint): nint(_nint) {}
    
    SEGYAdjInt(SEGYAdjInt const & m): nint(m.nint) {}
    
    ~SEGYAdjInt() {}
    
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> &,
		    LocalDataContainer<float> const &);
    
    string getName() const { string tmp = "SEGYAdjInt"; return tmp; }
    
  };

}

#endif
