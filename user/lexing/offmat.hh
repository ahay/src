#ifndef _OFFMAT_HH_
#define _OFFMAT_HH_

#include "offvec.hh"

template <class F>
class OffMat
{
public:
  int _m, _n;
  int _s, _t;
  bool _owndata;
  F* _data;
public:
  OffMat(int m=0, int n=0, int s=0, int t=0): _m(m), _n(n), _s(s), _t(t), _owndata(true) {
	if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
  }
  OffMat(int m, int n, int s, int t, bool owndata, F* data): _m(m), _n(n), _s(s), _t(t), _owndata(owndata) {
	if(_owndata) {
	  if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
	  if(_m>0 && _n>0) { for(int i=0; i<_m*_n; i++) _data[i] = data[i]; }
	} else {
	  _data = data;
	}
  }
  OffMat(const OffMat& C): _m(C._m), _n(C._n), _s(C._s), _t(C._t), _owndata(C._owndata) {
	if(_owndata) {
	  if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
	  if(_m>0 && _n>0) { for(int i=0; i<_m*_n; i++) _data[i] = C._data[i]; }
	 } else {
		_data = C._data;
	 }
  }
  ~OffMat() { 
	 if(_owndata) { 
		if(_m>0 && _n>0) { delete[] _data; _data = NULL; } 
	 }
  }
  OffMat& operator=(const OffMat& C) {
	 if(_owndata) { 
		if(_m>0 && _n>0) { delete[] _data; _data = NULL; } 
	 }
	 _m = C._m; _n=C._n; _s = C._s; _t=C._t; _owndata=C._owndata;
	 if(_owndata) {
		if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
		if(_m>0 && _n>0) { for(int i=0; i<_m*_n; i++) _data[i] = C._data[i]; }
	 } else {
		_data = C._data;
	 }
	 return *this;
  }
  void resize(int m, int n, int s, int t)  {
	assert( _owndata==true );
	if(_m!=m || _n!=n) {
	  if(_m>0 && _n>0) { delete[] _data; _data = NULL; }
	  _m = m; _n = n; _s = s; _t = t;
	  if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
	}
  }
  const F& operator()(int i, int j) const  {
	assert( i>=_s && i<_m+_s && j>=_t && j<_n+_t );
	return _data[(i-_s) + (j-_t)*_m];
  }
  F& operator()(int i, int j)  {
	assert( i>=_s && i<_m+_s && j>=_t && j<_n+_t );
	return _data[(i-_s) + (j-_t)*_m];
  }
  
  F* data() const { return _data; }
  F* clmdata(int j) { return &(_data[j*_m]); }
  int m() const { return _m; }
  int n() const { return _n; }
  int s() const { return _s; }
  int t() const { return _t; }
};

template <class F> inline ostream& operator<<( ostream& os, const OffMat<F>& mat)
{
  os<<mat.m()<<" "<<mat.n()<<" "<<mat.s()<<" "<<mat.t()<<endl;
  os.setf(ios_base::scientific, ios_base::floatfield);
  for(int i=mat.s(); i<mat.s()+mat.m(); i++) {
	for(int j=mat.t(); j<mat.t()+mat.n(); j++)
	  os<<" "<<mat(i,j);
	os<<endl;
  }
  return os;
}

template <class F> inline void setvalue(OffMat<F>& M, F val)
{
  for(int i=M.s(); i<M.s()+M.m(); i++)
	for(int j=M.t(); j<M.t()+M.n(); j++)
	  M(i,j) = val;
}

template <class F> inline float energy(OffMat<F>& M)
{
  float sum = 0;
  for(int i=M.s(); i<M.s()+M.m(); i++)
	for(int j=M.t(); j<M.t()+M.n(); j++)
	  sum += abs(M(i,j)*M(i,j));
  return sum;
}


typedef OffMat<bool>   BolOffMat;
typedef OffMat<int>    IntOffMat;
typedef OffMat<double> DblOffMat;
typedef OffMat<float>  FltOffMat;
typedef OffMat<cpx>    CpxOffMat;

#endif




