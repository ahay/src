#ifndef _NUMTNS_HH_
#define _NUMTNS_HH_

#include "nummat.hh"

template <class F>
class NumTns
{
public:
  int _m, _n, _p;
  bool _owndata;
  F* _data;
public:
  NumTns(int m=0, int n=0, int p=0): _m(m), _n(n), _p(p), _owndata(true) {
	if(_m>0 && _n>0 && _p>0) { _data = new F[_m*_n*_p]; assert( _data!=NULL ); } else _data=NULL;
  }
  NumTns(int m, int n, int p, bool owndata, F* data): _m(m), _n(n), _p(p), _owndata(owndata) {
	if(_owndata) {
	  if(_m>0 && _n>0 && _p>0) { _data = new F[_m*_n*_p]; assert( _data!=NULL ); } else _data=NULL;
	  if(_m>0 && _n>0 && _p>0) { for(int i=0; i<_m*_n*_p; i++) _data[i] = data[i]; }
	} else {
	  _data = data;
	}
  }
  NumTns(const NumTns& C): _m(C._m), _n(C._n), _p(C._p), _owndata(C._owndata) {
	if(_owndata) {
	  if(_m>0 && _n>0 && _p>0) { _data = new F[_m*_n*_p]; assert( _data!=NULL ); } else _data=NULL;
	  if(_m>0 && _n>0 && _p>0) { for(int i=0; i<_m*_n*_p; i++) _data[i] = C._data[i]; }
	} else {
	  _data = C._data;
	}
  }
  ~NumTns() { 
	if(_owndata) { 
	  if(_m>0 && _n>0) { delete[] _data; _data = NULL; } 
	}
  }
  NumTns& operator=(const NumTns& C) {
	if(_owndata) { 
	  if(_m>0 && _n>0 && _p>0) { delete[] _data; _data = NULL; } 
	}
	_m = C._m; _n=C._n; _p=C._p; _owndata=C._owndata;
	if(_owndata) {
	  if(_m>0 && _n>0 && _p>0) { _data = new F[_m*_n*_p]; assert( _data!=NULL ); } else _data=NULL;
	  if(_m>0 && _n>0 && _p>0) { for(int i=0; i<_m*_n*_p; i++) _data[i] = C._data[i]; }
	} else {
	  _data = C._data;
	}
	return *this;
  }
  void resize(int m, int n, int p)  {
	assert( _owndata==true );
	if(_m!=m || _n!=n || _p!=p) {
	  if(_m>0 && _n>0 && _p>0) { delete[] _data; _data = NULL; } 
	  _m = m; _n = n; _p=p;
	  if(_m>0 && _n>0 && _p>0) { _data = new F[_m*_n*_p]; assert( _data!=NULL ); } else _data=NULL;
	}
  }
  const F& operator()(int i, int j, int k) const  {
	assert( i>=0 && i<_m && j>=0 && j<_n && k>=0 && k<_p);
	return _data[i+j*_m+k*_m*_n];
  }
  F& operator()(int i, int j, int k)  {
	assert( i>=0 && i<_m && j>=0 && j<_n && k>=0 && k<_p);
	return _data[i+j*_m+k*_m*_n];
  }
  
  F* data() const { return _data; }
  int m() const { return _m; }
  int n() const { return _n; }
  int p() const { return _p; }
};

template <class F> inline ostream& operator<<( ostream& os, const NumTns<F>& tns)
{
  os<<tns.m()<<" "<<tns.n()<<" "<<tns.p()<<endl;
  os.setf(ios_base::scientific, ios_base::floatfield);
  for(int i=0; i<tns.m(); i++) {
	 for(int j=0; j<tns.n(); j++) {
		for(int k=0; k<tns.p(); k++) {
		  os<<" "<<tns(i,j,k);
		}
		os<<endl;
	 }
	 os<<endl;
  }
  return os;
}
template <class F> inline void setvalue(NumTns<F>& T, F val)
{
  for(int i=0; i<T.m(); i++)	for(int j=0; j<T.n(); j++)	  for(int k=0; k<T.p(); k++)		T(i,j,k) = val;
  return;
}
template <class F> inline float energy(NumTns<F>& T)
{
  float sum = 0;
  for(int i=0; i<T.m(); i++)	for(int j=0; j<T.n(); j++)	  for(int k=0; k<T.p(); k++)		sum += abs(T(i,j,k)*T(i,j,k));
  return sum;
}

typedef NumTns<bool>   BolNumTns;
typedef NumTns<int>    IntNumTns;
typedef NumTns<double> DblNumTns;
typedef NumTns<float>  FltNumTns;
typedef NumTns<cpx>    CpxNumTns;

#endif




