#ifndef _NUMMAT_HH_
#define _NUMMAT_HH_

#include "numvec.hh"

template <class F>
class NumMat
{
public:
    int _m, _n;
    bool _owndata;
    F* _data;
public:
    NumMat(int m=0, int n=0): _m(m), _n(n), _owndata(true) {
	if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
    }
    NumMat(int m, int n, bool owndata, F* data): _m(m), _n(n), _owndata(owndata) {
	if(_owndata) {
	    if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
	    if(_m>0 && _n>0) { for(int i=0; i<_m*_n; i++) _data[i] = data[i]; }
	} else {
	    _data = data;
	}
    }
    NumMat(const NumMat& C): _m(C._m), _n(C._n), _owndata(C._owndata) {
	if(_owndata) {
	    if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
	    if(_m>0 && _n>0) { for(int i=0; i<_m*_n; i++) _data[i] = C._data[i]; }
	} else {
	    _data = C._data;
	}
    }
    ~NumMat() {
	if(_owndata) {
	    if(_m>0 && _n>0) { delete[] _data; _data = NULL; }
	}
    }
    NumMat& operator=(const NumMat& C) {
	if(_owndata) {
	    if(_m>0 && _n>0) { delete[] _data; _data = NULL; }
	}
	_m = C._m; _n=C._n; _owndata=C._owndata;
	if(_owndata) {
	    if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
	    if(_m>0 && _n>0) { for(int i=0; i<_m*_n; i++) _data[i] = C._data[i]; }
	} else {
	    _data = C._data;
	}
	return *this;
    }
    void resize(int m, int n)  {
	assert( _owndata==true );
	if(_m!=m || _n!=n) {
	    if(_m>0 && _n>0) { delete[] _data; _data = NULL; }
	    _m = m; _n = n;
	    if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
	}
    }
    const F& operator()(int i, int j) const  { 
	assert( i>=0 && i<_m && j>=0 && j<_n );
	return _data[i+j*_m];
    }
    F& operator()(int i, int j)  { 
	assert( i>=0 && i<_m && j>=0 && j<_n );
	return _data[i+j*_m];
    }
  
    F* data() const { return _data; }
    F* clmdata(int j) { return &(_data[j*_m]); }
    int m() const { return _m; }
    int n() const { return _n; }
};

template <class F> inline ostream& operator<<( ostream& os, const NumMat<F>& mat)
{
    os<<mat.m()<<" "<<mat.n()<<endl;
    os.setf(ios_base::scientific, ios_base::floatfield);
    for(int i=0; i<mat.m(); i++) {
	for(int j=0; j<mat.n(); j++)
	    os<<" "<<mat(i,j);
	os<<endl;
    }
    return os;
}
template <class F> inline void setvalue(NumMat<F>& M, F val)
{
    for(int i=0; i<M.m(); i++)
	for(int j=0; j<M.n(); j++)
	    M(i,j) = val;
}
template <class F> inline float energy(NumMat<F>& M)
{
    float sum = 0;
    for(int i=0; i<M.m(); i++)	for(int j=0; j<M.n(); j++)	  sum += abs(M(i,j)*M(i,j));
    return sum;
}

typedef NumMat<bool>   BolNumMat;
typedef NumMat<int>    IntNumMat;
typedef NumMat<double> DblNumMat;
typedef NumMat<float>  FltNumMat;
typedef NumMat<cpx>    CpxNumMat;
typedef NumMat<zpx>    ZpxNumMat;

#endif




