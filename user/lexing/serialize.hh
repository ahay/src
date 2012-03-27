#ifndef _SERIALIZE_HH_
#define _SERIALIZE_HH_

#include "commoninc.hh"
#include "vec2t.hh"
#include "vec3t.hh"
#include "numtns.hh"

using std::vector;
using std::set;
using std::map;
using std::pair;
using std::istream;
using std::ostream;
using std::istringstream;
using std::ostringstream;
using std::string;
using std::cerr;

template<class T, class S> int serialize(const pair<T,S>& val, ostream& os, const vector<int>& mask);
template<class T, class S> int deserialize(pair<T,S>& val, istream& is, const vector<int>& mask);

//-------------------
//char
inline int serialize(const char& val, ostream& os, const vector<int>& mask)
{
  os.write((char*)&val, sizeof(char));
  return 0;
}

inline int deserialize(char& val, istream& is, const vector<int>& mask)
{
  is.read((char*)&val, sizeof(char));
  return 0;
}

//-------------------
//int
inline int serialize(const int& val, ostream& os, const vector<int>& mask)
{
  os.write((char*)&val, sizeof(int));
  return 0;
}

inline int deserialize(int& val, istream& is, const vector<int>& mask)
{
  is.read((char*)&val, sizeof(int));
  return 0;
}

//-------------------
//float
inline int serialize(const float& val, ostream& os, const vector<int>& mask)
{
  os.write((char*)&val, sizeof(float));
  return 0;
}

inline int deserialize(float& val, istream& is, const vector<int>& mask)
{
  is.read((char*)&val, sizeof(float));
  return 0;
}

//-------------------
//cpx
inline int serialize(const cpx& val, ostream& os, const vector<int>& mask)
{
  os.write((char*)&val, sizeof(cpx));
  return 0;
}

inline int deserialize(cpx& val, istream& is, const vector<int>& mask)
{
  is.read((char*)&val, sizeof(cpx));
  return 0;
}

//-------------------
//Index2
inline int serialize(const Index2& val, ostream& os, const vector<int>& mask)
{
  os.write((char*)&(val[0]), 2*sizeof(int));
  return 0;
}

inline int deserialize(Index2& val, istream& is, const vector<int>& mask)
{
  is.read((char*)&(val[0]), 2*sizeof(int));
  return 0;
}

//-------------------
//Point2
inline int serialize(const Point2& val, ostream& os, const vector<int>& mask)
{
  os.write((char*)&(val[0]), 2*sizeof(float));
  return 0;
}

inline int deserialize(Point2& val, istream& is, const vector<int>& mask)
{
  is.read((char*)&(val[0]), 2*sizeof(float));
  return 0;
}

//-------------------
//Index3
inline int serialize(const Index3& val, ostream& os, const vector<int>& mask)
{
  os.write((char*)&(val[0]), 3*sizeof(int));
  return 0;
}

inline int deserialize(Index3& val, istream& is, const vector<int>& mask)
{
  is.read((char*)&(val[0]), 3*sizeof(int));
  return 0;
}

//-------------------
//Point3
inline int serialize(const Point3& val, ostream& os, const vector<int>& mask)
{
  os.write((char*)&(val[0]), 3*sizeof(float));
  return 0;
}

inline int deserialize(Point3& val, istream& is, const vector<int>& mask)
{
  is.read((char*)&(val[0]), 3*sizeof(float));
  return 0;
}

//-------------------
//vector
template<class T>
int serialize(const vector<T>& val, ostream& os, const vector<int>& mask)
{
  int sz = val.size();
  os.write((char*)&sz, sizeof(int));
  for(int k=0; k<sz; k++)
	serialize(val[k], os, mask);
  return 0;
}

template<class T>
int deserialize(vector<T>& val, istream& is, const vector<int>& mask)
{
  int sz;
  is.read((char*)&sz, sizeof(int));
  val.resize(sz);
  for(int k=0; k<sz; k++)
    deserialize(val[k], is, mask);
  return 0;
}

//-------------------
//set
template<class T>
int serialize(const set<T>& val, ostream& os, const vector<int>& mask)
{
  int sz = val.size();
  os.write((char*)&sz, sizeof(int));
  for(typename set<T>::const_iterator mi=val.begin(); mi!=val.end(); mi++) 
	serialize((*mi), os, mask);
  return 0;
}

template<class T>
int deserialize(set<T>& val, istream& is, const vector<int>& mask)
{
  val.clear();
  int sz;
  is.read((char*)&sz, sizeof(int));
  for(int k=0; k<sz; k++) {
	T t; deserialize(t, is, mask);
	val.insert(t);
  }
  return 0;
}

//-------------------
//map
template<class T, class S>
int serialize(const map<T,S>& val, ostream& os, const vector<int>& mask)
{
  int sz = val.size();
  os.write((char*)&sz, sizeof(int));
  for(typename map<T,S>::const_iterator mi=val.begin(); mi!=val.end(); mi++) {
	serialize((*mi).first, os, mask);
	serialize((*mi).second, os, mask);
  }
  return 0;
}

template<class T, class S>
int deserialize(map<T,S>& val, istream& is, const vector<int>& mask)
{
  val.clear();
  int sz;
  is.read((char*)&sz, sizeof(int));  //cerr<<"DES sz "<<sz<<endl;
  for(int k=0; k<sz; k++) {    //cerr<<"  k stt "<<k<<" "<<sz<<endl;
    T t;	deserialize(t, is, mask);
    S s;	deserialize(s, is, mask);
    val[t] = s;    //cerr<<"  k end "<<k<<" "<<sz<<endl;
  }
  return 0;
}

//-------------------
//pair
template<class T, class S>
int serialize(const pair<T,S>& val, ostream& os, const vector<int>& mask)
{
  serialize(val.first, os, mask);
  serialize(val.second, os, mask);
  return 0;
}

template<class T, class S>
int deserialize(pair<T,S>& val, istream& is, const vector<int>& mask)
{
  deserialize(val.first, is, mask);
  deserialize(val.second, is, mask);
  return 0;
}


//-------------------
//BolNumVec
inline int serialize(const BolNumVec& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();
  os.write((char*)&m, sizeof(int));
  os.write((char*)(val.data()), m*sizeof(bool));
  return 0;
}

inline int deserialize(BolNumVec& val, istream& is, const vector<int>& mask)
{
  int m;
  is.read((char*)&m, sizeof(int));
  val.resize(m);
  is.read((char*)(val.data()), m*sizeof(bool));
  return 0;
}

//-------------------
//BolNumMat
inline int serialize(const BolNumMat& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();
  int n = val.n();
  os.write((char*)&m, sizeof(int));
  os.write((char*)&n, sizeof(int));
  os.write((char*)(val.data()), m*n*sizeof(bool));
  return 0;
}

inline int deserialize(BolNumMat& val, istream& is, const vector<int>& mask)
{
  int m;
  int n;
  is.read((char*)&m, sizeof(int));
  is.read((char*)&n, sizeof(int));
  val.resize(m,n);
  is.read((char*)(val.data()), m*n*sizeof(bool));
  return 0;
}

//-------------------
//BolNumTns
inline int serialize(const BolNumTns& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();  int n = val.n();  int p = val.p();
  os.write((char*)&m, sizeof(int));
  os.write((char*)&n, sizeof(int));
  os.write((char*)&p, sizeof(int));
  os.write((char*)(val.data()), m*n*p*sizeof(bool));
  return 0;
}

inline int deserialize(BolNumTns& val, istream& is, const vector<int>& mask)
{
  int m,n,p;
  is.read((char*)&m, sizeof(int));
  is.read((char*)&n, sizeof(int));
  is.read((char*)&p, sizeof(int));
  val.resize(m,n,p);
  is.read((char*)(val.data()), m*n*p*sizeof(bool));
  return 0;
}


//-------------------
//IntNumVec
inline int serialize(const IntNumVec& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();
  os.write((char*)&m, sizeof(int));
  os.write((char*)(val.data()), m*sizeof(int));
  return 0;
}

inline int deserialize(IntNumVec& val, istream& is, const vector<int>& mask)
{
  int m;
  is.read((char*)&m, sizeof(int));
  val.resize(m);
  is.read((char*)(val.data()), m*sizeof(int));
  return 0;
}

//-------------------
//IntNumMat
inline int serialize(const IntNumMat& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();
  int n = val.n();
  os.write((char*)&m, sizeof(int));
  os.write((char*)&n, sizeof(int));
  os.write((char*)(val.data()), m*n*sizeof(int));
  return 0;
}

inline int deserialize(IntNumMat& val, istream& is, const vector<int>& mask)
{
  int m;
  int n;
  is.read((char*)&m, sizeof(int));
  is.read((char*)&n, sizeof(int));
  val.resize(m,n);
  is.read((char*)(val.data()), m*n*sizeof(int));
  return 0;
}

//-------------------
//IntNumTns
inline int serialize(const IntNumTns& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();  int n = val.n();  int p = val.p();
  os.write((char*)&m, sizeof(int));
  os.write((char*)&n, sizeof(int));
  os.write((char*)&p, sizeof(int));
  os.write((char*)(val.data()), m*n*p*sizeof(int));
  return 0;
}

inline int deserialize(IntNumTns& val, istream& is, const vector<int>& mask)
{
  int m,n,p;
  is.read((char*)&m, sizeof(int));
  is.read((char*)&n, sizeof(int));
  is.read((char*)&p, sizeof(int));
  val.resize(m,n,p);
  is.read((char*)(val.data()), m*n*p*sizeof(int));
  return 0;
}


//-------------------
//FltNumVec
inline int serialize(const FltNumVec& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();
  os.write((char*)&m, sizeof(int));
  os.write((char*)(val.data()), m*sizeof(float));
  return 0;
}

inline int deserialize(FltNumVec& val, istream& is, const vector<int>& mask)
{
  int m;
  is.read((char*)&m, sizeof(int));
  val.resize(m);
  is.read((char*)(val.data()), m*sizeof(float));
  return 0;
}

//-------------------
//FltNumMat
inline int serialize(const FltNumMat& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();
  int n = val.n();
  os.write((char*)&m, sizeof(int));
  os.write((char*)&n, sizeof(int));
  os.write((char*)(val.data()), m*n*sizeof(float));
  return 0;
}

inline int deserialize(FltNumMat& val, istream& is, const vector<int>& mask)
{
  int m;
  int n;
  is.read((char*)&m, sizeof(int));
  is.read((char*)&n, sizeof(int));
  val.resize(m,n);
  is.read((char*)(val.data()), m*n*sizeof(float));
  return 0;
}

//-------------------
//FltNumTns
inline int serialize(const FltNumTns& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();  int n = val.n();  int p = val.p();
  os.write((char*)&m, sizeof(int));
  os.write((char*)&n, sizeof(int));
  os.write((char*)&p, sizeof(int));
  os.write((char*)(val.data()), m*n*p*sizeof(float));
  return 0;
}

inline int deserialize(FltNumTns& val, istream& is, const vector<int>& mask)
{
  int m,n,p;
  is.read((char*)&m, sizeof(int));
  is.read((char*)&n, sizeof(int));
  is.read((char*)&p, sizeof(int));
  val.resize(m,n,p);
  is.read((char*)(val.data()), m*n*p*sizeof(float));
  return 0;
}

//-------------------
//CpxNumVec
inline int serialize(const CpxNumVec& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();
  os.write((char*)&m, sizeof(int));
  os.write((char*)(val.data()), m*sizeof(cpx));
  return 0;
}

inline int deserialize(CpxNumVec& val, istream& is, const vector<int>& mask)
{
  int m;
  is.read((char*)&m, sizeof(int));
  val.resize(m);
  is.read((char*)(val.data()), m*sizeof(cpx));
  return 0;
}

//-------------------
//CpxNumMat
inline int serialize(const CpxNumMat& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();
  int n = val.n();
  os.write((char*)&m, sizeof(int));
  os.write((char*)&n, sizeof(int));
  os.write((char*)(val.data()), m*n*sizeof(cpx));
  return 0;
}

inline int deserialize(CpxNumMat& val, istream& is, const vector<int>& mask)
{
  int m;
  int n;
  is.read((char*)&m, sizeof(int));
  is.read((char*)&n, sizeof(int));
  val.resize(m,n);
  is.read((char*)(val.data()), m*n*sizeof(cpx));
  return 0;
}

//-------------------
//CpxNumTns
inline int serialize(const CpxNumTns& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();  int n = val.n();  int p = val.p();
  os.write((char*)&m, sizeof(int));
  os.write((char*)&n, sizeof(int));
  os.write((char*)&p, sizeof(int));
  os.write((char*)(val.data()), m*n*p*sizeof(cpx));
  return 0;
}

inline int deserialize(CpxNumTns& val, istream& is, const vector<int>& mask)
{
  int m,n,p;
  is.read((char*)&m, sizeof(int));
  is.read((char*)&n, sizeof(int));
  is.read((char*)&p, sizeof(int));
  val.resize(m,n,p);
  is.read((char*)(val.data()), m*n*p*sizeof(cpx));
  return 0;
}

//-------------------
//NumVec
template<class T>
int inline serialize(const NumVec<T>& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();
  os.write((char*)&m, sizeof(int));
  for(int i=0; i<m; i++)
	serialize(val(i), os, mask);
  return 0;
}
template<class T>
int inline deserialize(NumVec<T>& val, istream& is, const vector<int>& mask)
{
  int m;
  is.read((char*)&m, sizeof(int));
  val.resize(m);
  for(int i=0; i<m; i++)
	deserialize(val(i), is, mask);
  return 0;
}

//-------------------
//NumMat
template<class T>
int inline serialize(const NumMat<T>& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();
  int n = val.n();
  os.write((char*)&m, sizeof(int));
  os.write((char*)&n, sizeof(int));
  for(int j=0; j<n; j++)
	for(int i=0; i<m; i++)
	  serialize(val(i,j), os, mask);
  return 0;
}
template<class T>
int inline deserialize(NumMat<T>& val, istream& is, const vector<int>& mask)
{
  int m;
  int n;
  is.read((char*)&m, sizeof(int));
  is.read((char*)&n, sizeof(int));
  val.resize(m,n);
  for(int j=0; j<n; j++)
	for(int i=0; i<m; i++)
	  deserialize(val(i,j), is, mask);
  return 0;
}

//-------------------
//NumTns
template<class T>
int inline serialize(const NumTns<T>& val, ostream& os, const vector<int>& mask)
{
  int m = val.m();
  int n = val.n();
  int p = val.p();
  os.write((char*)&m, sizeof(int));
  os.write((char*)&n, sizeof(int));
  os.write((char*)&p, sizeof(int));
  for(int k=0; k<p; k++)
	for(int j=0; j<n; j++)
	  for(int i=0; i<m; i++)
		serialize(val(i,j,k), os, mask);
  return 0;
}
template<class T>
int inline deserialize(NumTns<T>& val, istream& is, const vector<int>& mask)
{
  int m;
  int n;
  int p;
  is.read((char*)&m, sizeof(int));
  is.read((char*)&n, sizeof(int));
  is.read((char*)&p, sizeof(int));
  val.resize(m,n,p);
  for(int k=0; k<p; k++)
	for(int j=0; j<n; j++)
	  for(int i=0; i<m; i++)
		deserialize(val(i,j,k), is, mask);
  return 0;
}


#endif



