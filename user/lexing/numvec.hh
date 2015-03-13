#ifndef _NUMVEC_HH_
#define _NUMVEC_HH_

#include "commoninc.hh"

template <class F>
class NumVec
{
public:
  int  _m;
  bool _owndata;
  F* _data;
public:
  NumVec(int m=0): _m(m), _owndata(true)  {
	if(_m>0) { _data = new F[_m]; assert(_data!=NULL); } else _data=NULL;
  }
  NumVec(int m, bool owndata, F* data): _m(m), _owndata(owndata) {
	if(_owndata) {
	  if(_m>0) { _data = new F[_m]; assert(_data!=NULL); } else _data=NULL;
	  if(_m>0) { for(int i=0; i<_m; i++) _data[i] = data[i]; }
	} else {
	  _data = data;
	}
  }
  NumVec(const NumVec& C): _m(C._m), _owndata(C._owndata)  {
	if(_owndata) {
	  if(_m>0) { _data = new F[_m]; assert(_data!=NULL); } else _data=NULL;
	  if(_m>0) { for(int i=0; i<_m; i++) _data[i] = C._data[i]; }
	} else {
	  _data = C._data;
	}
  }
  ~NumVec() {
	if(_owndata) {
	  if(_m>0) { delete[] _data; _data = NULL; }
	}
  }
  NumVec& operator=(const NumVec& C)  {
	if(_owndata) {
	  if(_m>0) { delete[] _data; _data = NULL; }
	}
	_m = C._m; _owndata=C._owndata;
	if(_owndata) {
	  if(_m>0) { _data = new F[_m]; assert(_data!=NULL); } else _data=NULL;
	  if(_m>0) { for(int i=0; i<_m; i++) _data[i] = C._data[i]; }
	} else {
	  _data =C._data;
	}
	return *this;
  }
  void resize(int m)  {
	assert(_owndata==true);
	if(m !=_m) {
	  if(_m>0) { delete[] _data; _data = NULL; }
	  _m = m;
	  if(_m>0) { _data = new F[_m]; assert(_data!=NULL); } else _data=NULL;
	}
  }
  const F& operator()(int i) const  {
	assert(i>=0 && i<_m);
	return _data[i]; 
  }
  F& operator()(int i)  {
	assert(i>=0 && i<_m);
	return _data[i]; 
  }
  
  F* data() const { return _data; }
  int m () const { return _m; }
};

template <class F> inline ostream& operator<<( ostream& os, const NumVec<F>& vec)
{
  os<<vec.m()<<endl;
  os.setf(ios_base::scientific, ios_base::floatfield);
  for(int i=0; i<vec.m(); i++)	 os<<" "<<vec(i);
  os<<endl;
  return os;
}

template <class F> inline void setvalue(NumVec<F>& vec, F val)
{
  for(int i=0; i<vec.m(); i++)
    vec(i) = val;
}
template <class F> inline float energy(NumVec<F>& vec)
{
  float sum = 0;
  for(int i=0; i<vec.m(); i++)    sum += abs(vec(i)*vec(i));
  return sum;
}  
template <class F> inline float energy(const NumVec<F>& vec)
{
  float sum = 0;
  for(int i=0; i<vec.m(); i++)    sum += abs(vec(i)*vec(i));
  return sum;
}  


typedef NumVec<bool>   BolNumVec;
typedef NumVec<int>    IntNumVec;
typedef NumVec<double> DblNumVec;
typedef NumVec<float>  FltNumVec;
typedef NumVec<cpx>    CpxNumVec;
typedef NumVec<zpx>    ZpxNumVec;

#endif


