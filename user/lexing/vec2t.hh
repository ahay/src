#ifndef  _VEC2T_HH_
#define  _VEC2T_HH_

#include "commoninc.hh"

using std::istream;
using std::ostream;
using std::min;
using std::max;
using std::abs;

//-----------------------------------------------------------------------------------------------------
///Common VECtor Template
template <class F>
class Vec2T {
private:
  F _v[2];
public:
  enum{ X=0, Y=1 };
  //------------CONSTRUCTOR AND DESTRUCTOR 
  Vec2T()              { _v[0]=F(0);    _v[1]=F(0); }  //Vec2T(F f)           { _v[0]=f;       _v[1]=f; }
  Vec2T(const F* f)    { _v[0]=f[0];    _v[1]=f[1]; }
  Vec2T(F a,F b)       { _v[0]=a;       _v[1]=b; }
  Vec2T(const Vec2T& c){ _v[0]=c._v[0]; _v[1]=c._v[1]; }
  ~Vec2T() {}
  //------------POINTER and ACCESS
  operator F*()             { return &_v[0]; }
  operator const F*() const { return &_v[0]; }
  F* array()                { return &_v[0]; }  //access array
  F& operator()(int i)             { assert(i<2); return _v[i]; }
  const F& operator()(int i) const { assert(i<2); return _v[i]; }
  F& operator[](int i)             { assert(i<2); return _v[i]; }
  const F& operator[](int i) const { assert(i<2); return _v[i]; }
  //------------ASSIGN
  Vec2T& operator= ( const Vec2T& c ) { _v[0] =c._v[0]; _v[1] =c._v[1]; return *this; }
  Vec2T& operator+=( const Vec2T& c ) { _v[0]+=c._v[0]; _v[1]+=c._v[1]; return *this; }
  Vec2T& operator-=( const Vec2T& c ) { _v[0]-=c._v[0]; _v[1]-=c._v[1]; return *this; }
  Vec2T& operator*=( const F& s )     { _v[0]*=s;       _v[1]*=s;       return *this; }
  Vec2T& operator/=( const F& s )     { _v[0]/=s;       _v[1]/=s;       return *this; }
  //-----------LENGTH...
  F l1( void )     const  { F sum=F(0); for(int i=0; i<2; i++) sum=sum+abs(_v[i]); return sum; }
  F linfty( void ) const  { F cur=F(0); for(int i=0; i<2; i++) cur=max(cur,abs(_v[i])); return cur; }
  F l2( void )     const  { F sum=F(0); for(int i=0; i<2; i++) sum=sum+_v[i]*_v[i]; return sqrt(sum); }
  //F length( void ) const  { return l2(); }
  //Vec2T dir( void )    const  { F a=l2(); return (*this)/a; }
};

//-----------LEX COMPARE
template <class F> inline bool operator==(const Vec2T<F>& a, const Vec2T<F>& b) {
  return (a[0]==b[0] && a[1]==b[1]);
}
template <class F> inline bool operator!=(const Vec2T<F>& a, const Vec2T<F>& b) {
  return !(a==b);
}
template <class F> inline bool operator> (const Vec2T<F>& a, const Vec2T<F>& b) {
  for(int i=0; i<2; i++) {
	if(     a[i]>b[i])	  return true;
	else if(a[i]<b[i])	  return false;
  }
  return false;
}
template <class F> inline bool operator< (const Vec2T<F>& a, const Vec2T<F>& b) {
  for(int i=0; i<2; i++) {
	if(     a[i]<b[i])	  return true;
	else if(a[i]>b[i])	  return false;
  }
  return false;
}
template <class F> inline bool operator>=(const Vec2T<F>& a, const Vec2T<F>& b) {
  for(int i=0; i<2; i++) {
	if(     a[i]>b[i])	  return true;
	else if(a[i]<b[i])	  return false;
  }
  return true;
}
template <class F> inline bool operator<=(const Vec2T<F>& a, const Vec2T<F>& b) {
  for(int i=0; i<2; i++) {
	if(     a[i]<b[i])	  return true;
	else if(a[i]>b[i])	  return false;
  }
  return true;
}

//-----------NUMERICAL OPS
template <class F> inline Vec2T<F> operator- (const Vec2T<F>& a) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = -a[i]; return r;
}
template <class F> inline Vec2T<F> operator+ (const Vec2T<F>& a, const Vec2T<F>& b) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = a[i]+b[i]; return r; 
}
template <class F> inline Vec2T<F> operator- (const Vec2T<F>& a, const Vec2T<F>& b) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = a[i]-b[i]; return r;
}
template <class F> inline Vec2T<F> operator* (F scl, const Vec2T<F>& a) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = scl*a[i];  return r;
}
template <class F> inline Vec2T<F> operator* (const Vec2T<F>& a, F scl) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = scl*a[i];  return r;
}
template <class F> inline Vec2T<F> operator/ (const Vec2T<F>& a, F scl) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = a[i]/scl;  return r;
}
template <class F> inline F operator* (const Vec2T<F>& a, const Vec2T<F>& b) {
  F sum=F(0); for(int i=0; i<2; i++) sum=sum+a(i)*b(i); return sum;
}
template <class F> inline F dot       (const Vec2T<F>& a, const Vec2T<F>& b) {
  return a*b;
}

//-------------ew NUMERICAL OPS
template <class F> inline Vec2T<F> ewmin(const Vec2T<F>& a, const Vec2T<F>& b) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = min(a[i], b[i]); return r;
}
template <class F> inline Vec2T<F> ewmax(const Vec2T<F>& a, const Vec2T<F>& b) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = max(a[i], b[i]); return r;
}
template <class F> inline Vec2T<F> ewabs(const Vec2T<F>& a) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = abs(a[i]); return r;
}
template <class F> inline Vec2T<F> ewmul(const Vec2T<F>&a, const Vec2T<F>& b) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = a[i]*b[i]; return r;
}
template <class F> inline Vec2T<F> ewdiv(const Vec2T<F>&a, const Vec2T<F>& b) { 
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = a[i]/b[i]; return r;
}
template <class F> inline Vec2T<F> ewrnd(const Vec2T<F>&a) { //round
  Vec2T<F> r;  for(int i=0; i<2; i++)	r[i] = round(a[i]);  return r;
}

//---------------INOUT
template <class F> istream& operator>>(istream& is, Vec2T<F>& a) {
  for(int i=0; i<2; i++) is>>a[i]; return is;
}
template <class F> ostream& operator<<(ostream& os, const Vec2T<F>& a) { 
  for(int i=0; i<2; i++) os<<a[i]<<" "; return os;
}

//---------------------------------------------------------
/// MOST COMMONLY USED
typedef Vec2T<float> Point2;
typedef Vec2T<int>    Index2;



#endif
