#ifndef  _VEC3T_HH_
#define  _VEC3T_HH_

#include "commoninc.hh"

using std::istream;
using std::ostream;
using std::min;
using std::max;
using std::abs;

///Common VECtor Template
template <class F>
class Vec3T {
private:
  F _v[3];
public:
  enum{ X=0, Y=1, Z=2 };
  //------------CONSTRUCTOR AND DESTRUCTOR 
  Vec3T()              { _v[0]=F(0);    _v[1]=F(0);    _v[2]=F(0); }
  Vec3T(const F* f)    { _v[0]=f[0];    _v[1]=f[1];    _v[2]=f[2]; }
  Vec3T(F a,F b,F c)   { _v[0]=a;       _v[1]=b;       _v[2]=c; }
  Vec3T(const Vec3T& c){ _v[0]=c._v[0]; _v[1]=c._v[1]; _v[2]=c._v[2]; }
  ~Vec3T() {}
  //------------POINTER and ACCESS
  operator F*()             { return &_v[0]; }
  operator const F*() const { return &_v[0]; }
  F* array()                { return &_v[0]; }  //access array
  F& operator()(int i)             { assert(i<3); return _v[i]; }
  const F& operator()(int i) const { assert(i<3); return _v[i]; }
  F& operator[](int i)             { assert(i<3); return _v[i]; }
  const F& operator[](int i) const { assert(i<3); return _v[i]; }
  //------------ASSIGN
  Vec3T& operator= ( const Vec3T& c ) { _v[0] =c._v[0]; _v[1] =c._v[1]; _v[2] =c._v[2]; return *this; }
  Vec3T& operator+=( const Vec3T& c ) { _v[0]+=c._v[0]; _v[1]+=c._v[1]; _v[2]+=c._v[2]; return *this; }
  Vec3T& operator-=( const Vec3T& c ) { _v[0]-=c._v[0]; _v[1]-=c._v[1]; _v[2]-=c._v[2]; return *this; }
  Vec3T& operator*=( const F& s )     { _v[0]*=s;       _v[1]*=s;       _v[2]*=s;       return *this; }
  Vec3T& operator/=( const F& s )     { _v[0]/=s;       _v[1]/=s;       _v[2]/=s;       return *this; }
  //-----------LENGTH...
  F l1( void )     const  { F sum=F(0); for(int i=0; i<3; i++) sum=sum+abs(_v[i]); return sum; }
  F linfty( void ) const  { F cur=F(0); for(int i=0; i<3; i++) cur=max(cur,abs(_v[i])); return cur; }
  F l2( void )     const  { F sum=F(0); for(int i=0; i<3; i++) sum=sum+_v[i]*_v[i]; return sqrt(sum); }
};

//-----------LEX COMPARE
template <class F> inline bool operator==(const Vec3T<F>& a, const Vec3T<F>& b) {
  return (a[0]==b[0] && a[1]==b[1] && a[2]==b[2]);
}
template <class F> inline bool operator!=(const Vec3T<F>& a, const Vec3T<F>& b) {
  return !(a==b);
}
template <class F> inline bool operator> (const Vec3T<F>& a, const Vec3T<F>& b) {
  for(int i=0; i<3; i++) {
	if(     a[i]>b[i])	  return true;
	else if(a[i]<b[i])	  return false;
  }
  return false;
}
template <class F> inline bool operator< (const Vec3T<F>& a, const Vec3T<F>& b) {
  for(int i=0; i<3; i++) {
	if(     a[i]<b[i])	  return true;
	else if(a[i]>b[i])	  return false;
  }
  return false;
}
template <class F> inline bool operator>=(const Vec3T<F>& a, const Vec3T<F>& b) {
  for(int i=0; i<3; i++) {
	if(     a[i]>b[i])	  return true;
	else if(a[i]<b[i])	  return false;
  }
  return true;
}
template <class F> inline bool operator<=(const Vec3T<F>& a, const Vec3T<F>& b) {
  for(int i=0; i<3; i++) {
	if(     a[i]<b[i])	  return true;
	else if(a[i]>b[i])	  return false;
  }
  return true;
}

//-----------NUMERICAL OPS
template <class F> inline Vec3T<F> operator- (const Vec3T<F>& a) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = -a[i]; return r;
}
template <class F> inline Vec3T<F> operator+ (const Vec3T<F>& a, const Vec3T<F>& b) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = a[i]+b[i]; return r; 
}
template <class F> inline Vec3T<F> operator- (const Vec3T<F>& a, const Vec3T<F>& b) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = a[i]-b[i]; return r;
}
template <class F> inline Vec3T<F> operator* (F scl, const Vec3T<F>& a) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = scl*a[i];  return r;
}
template <class F> inline Vec3T<F> operator* (const Vec3T<F>& a, F scl) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = scl*a[i];  return r;
}
template <class F> inline Vec3T<F> operator/ (const Vec3T<F>& a, F scl) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = a[i]/scl;  return r;
}
template <class F> inline F operator* (const Vec3T<F>& a, const Vec3T<F>& b) {
  F sum=F(0); for(int i=0; i<3; i++) sum=sum+a(i)*b(i); return sum;
}
template <class F> inline F dot       (const Vec3T<F>& a, const Vec3T<F>& b) {
  return a*b;
}
template <class F> inline Vec3T<F> operator^ (const Vec3T<F>& a, const Vec3T<F>& b) {
  return Vec3T<F>(a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0));
}
template <class F> inline Vec3T<F> cross     (const Vec3T<F>& a, const Vec3T<F>& b) { 
  return a^b; 
}

//-------------ew NUMERICAL OPS
template <class F> inline Vec3T<F> ewmin(const Vec3T<F>& a, const Vec3T<F>& b) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = min(a[i], b[i]); return r;
}
template <class F> inline Vec3T<F> ewmax(const Vec3T<F>& a, const Vec3T<F>& b) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = max(a[i], b[i]); return r;
}
template <class F> inline Vec3T<F> ewabs(const Vec3T<F>& a) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = abs(a[i]); return r;
}
template <class F> inline Vec3T<F> ewmul(const Vec3T<F>&a, const Vec3T<F>& b) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = a[i]*b[i]; return r;
}
template <class F> inline Vec3T<F> ewdiv(const Vec3T<F>&a, const Vec3T<F>& b) { 
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = a[i]/b[i]; return r;
}
template <class F> inline Vec3T<F> ewrnd(const Vec3T<F>&a) { //round
  Vec3T<F> r;  for(int i=0; i<3; i++)	r[i] = round(a[i]);  return r;
}

/*
//-----------accumulative BOOLEAN OPS
template <class F> inline bool allequ(const Vec3T<F>& a, const Vec3T<F>& b) {
  bool res = true;  for(int i=0; i<3; i++)   res = res && (a(i)==b(i));  return res;
}
template <class F> inline bool allneq(const Vec3T<F>& a, const Vec3T<F>& b) {
  return !(a==b);
}
template <class F> inline bool allgtt(const Vec3T<F>& a, const Vec3T<F>& b) {
  bool res = true;  for(int i=0; i<3; i++)   res = res && (a(i)> b(i));  return res; 
}
template <class F> inline bool alllst(const Vec3T<F>& a, const Vec3T<F>& b) {
  bool res = true;  for(int i=0; i<3; i++)   res = res && (a(i)< b(i));  return res; 
}
template <class F> inline bool allgoe(const Vec3T<F>& a, const Vec3T<F>& b) {
  bool res = true;  for(int i=0; i<3; i++)	res = res && (a(i)>=b(i));  return res; 
}
template <class F> inline bool allloe(const Vec3T<F>& a, const Vec3T<F>& b) {
  bool res = true;  for(int i=0; i<3; i++)   res = res && (a(i)<=b(i));  return res; 
}
*/

//---------------INOUT
template <class F> istream& operator>>(istream& is, Vec3T<F>& a) {
  for(int i=0; i<3; i++) is>>a[i]; return is;
}
template <class F> ostream& operator<<(ostream& os, const Vec3T<F>& a) { 
  for(int i=0; i<3; i++) os<<a[i]<<" "; return os;
}

//---------------------------------------------------------
/// MOST COMMONLY USED
typedef Vec3T<float> Point3;
typedef Vec3T<int>    Index3;

#endif

