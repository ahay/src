#ifndef _OFFVEC_HH_
#define _OFFVEC_HH_

#include "commoninc.hh"

template <class F>
class OffVec
{
public:
  int  _m;
  int  _s;
  bool _owndata;
  F* _data;
public:
  OffVec(int m=0, int s=0): _m(m), _s(s), _owndata(true)  {
	if(_m>0) { _data = new F[_m]; assert(_data!=NULL); } else _data=NULL;
  }
  OffVec(int m, int s, bool owndata, F* data): _m(m), _s(s), _owndata(owndata) {
	if(_owndata) {
	  if(_m>0) { _data = new F[_m]; assert(_data!=NULL); } else _data=NULL;
	  if(_m>0) { for(int i=0; i<_m; i++) _data[i] = data[i]; }
	} else {
	  _data = data;
	}
  }
  OffVec(const OffVec& C): _m(C._m), _s(C._s), _owndata(C._owndata)  {
	if(_owndata) {
	  if(_m>0) { _data = new F[_m]; assert(_data!=NULL); } else _data=NULL;
	  if(_m>0) { for(int i=0; i<_m; i++) _data[i] = C._data[i]; }
	} else {
	  _data = C._data;
	}
  }
  ~OffVec() {
	if(_owndata) {
	  if(_m>0) { delete[] _data; _data = NULL; }
	}
  }
  OffVec& operator=(const OffVec& C)  {
	if(_owndata) {
	  if(_m>0) { delete[] _data; _data = NULL; }
	}
	_m=C._m; _s=C._s; _owndata=C._owndata;
	if(_owndata) {
	  if(_m>0) { _data = new F[_m]; assert(_data!=NULL); } else _data=NULL;
	  if(_m>0) { for(int i=0; i<_m; i++) _data[i] = C._data[i]; }
	 } else {
		_data =C._data;
	 }
	 return *this;
  }
  void resize(int m, int s)  {
	assert(_owndata==true);
	if(m !=_m) {
	  if(_m>0) { delete[] _data; _data = NULL; }
	  _m = m;		_s = s;
	  if(_m>0) { _data = new F[_m]; assert(_data!=NULL); } else _data=NULL;
	}
  }
  const F& operator()(int i) const  {
	assert(i>=_s && i<_m+_s);
	return _data[i-_s]; 
  }
  F& operator()(int i)  {
	assert(i>=_s && i<_m+_s);
	return _data[i-_s]; 
  }
  
  F* data() const { return _data; }
  int m() const { return _m; }
  int s() const { return _s; }
};

template <class F> inline ostream& operator<<( ostream& os, const OffVec<F>& vec)
{
  os<<vec.m()<<" "<<vec.s()<<endl;
  os.setf(ios_base::scientific, ios_base::floatfield);
  for(int i=vec.s(); i<vec.m()+vec.s(); i++)	 os<<" "<<vec(i);
  os<<endl;
  return os;
}
template <class F> inline void setvalue(OffVec<F>& vec, F val)
{
  for(int i=vec.s(); i<vec.s()+vec.m(); i++)	 vec(i) = val;
}
template <class F> inline float energy(OffVec<F>& vec)
{
  float sum = 0;
  for(int i=vec.s(); i<vec.s()+vec.m(); i++)
	sum += abs(vec(i)*vec(i));
  return sum;
}

typedef OffVec<bool>   BolOffVec;
typedef OffVec<int>    IntOffVec;
typedef OffVec<double> DblOffVec;
typedef OffVec<float>  FltOffVec;
typedef OffVec<cpx>    CpxOffVec;
typedef OffVec<cpx16>  ZpxOffVec;

#endif


