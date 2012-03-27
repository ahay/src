#ifndef _BFIO_HH_
#define _BFIO_HH_

#include <rsf.hh>

#include "comobject.hh"
#include "numtns.hh"
#include "offtns.hh"
#include "vec2t.hh"
#include "vecmatop.hh"

using std::vector;
using std::pair;
using std::map;

class Entry
{
public:
  FltNumVec _grid;
  NumVec<CpxNumMat> _mats;
public:
  Entry() {;}
  ~Entry() {;}
  FltNumVec& grid() { return _grid; }
  NumVec<CpxNumMat>& mats() { return _mats; }
};

int serialize(const Entry&, ostream&, const vector<int>&);
int deserialize(Entry&, istream&, const vector<int>&);

//-----------------------------
class BFIO: public ComObject
{
public:
  int _EPSx1;
  int _EPSx2;
  int _EPSp1;
  int _EPSp2;
  int _fi;
  int _EL;
  map<int, Entry> _e2dmap;
  float tmin, tmax, pmin, pmax;
  float wmin, wmax, zmin, zmax; 
public:
  BFIO(const string& p): ComObject(p) {;}
  ~BFIO() {;}
  int& EPSx1() { return _EPSx1; }
  int& EPSx2() { return _EPSx2; }
  int& EPSp1() { return _EPSp1; }
  int& EPSp2() { return _EPSp2; }
  int& fi() { return _fi; }
  int& EL() { return _EL; }
  map<int, Entry>& e2dmap() { return _e2dmap; }
  //
  int setup(iRSF& par, iRSF& inp);
  int setup32(iRSF& par, iRSF& inp);
  int setup23(iRSF& par, iRSF& inp);
  int eval(int N, const CpxNumMat& f, const FltNumVec& w, const FltNumVec& z, CpxNumMat& u, const FltNumVec& t, const FltNumVec& p);
  int kernel(int N, vector<Point2>& trg, vector<Point2>& src, CpxNumMat& res);
  int check(int N, const CpxNumMat& f, const FltNumVec& w, const FltNumVec& z, const CpxNumMat& u, const FltNumVec& t, const FltNumVec& p, int NC, float& relerr);
  //
  int prep_aux(FltNumVec& grid, vector<float>& ts, CpxNumMat& tmp);
};

#endif



