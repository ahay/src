#ifndef _BFIO_HH_
#define _BFIO_HH_

#include <rsf.hh>

#include "comobject.hh"
#include "numtns.hh"
#include "offtns.hh"
#include "vec2t.hh"
#include "vec3t.hh"

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
  int _EPSx3;
  int _EPSk1;
  int _EPSk2;
  int _EPSk3;
  int _fi;
  int _EL;
  map<int, Entry> _e2dmap;
  float taumin, taumax, pmin, pmax, qmin, qmax;
  float wmin, wmax, xmin, xmax, ymin, ymax;
public:
  BFIO(const string& p): ComObject(p) {;}
  ~BFIO() {;}
  int& EPSx1() { return _EPSx1; }
  int& EPSx2() { return _EPSx2; }
  int& EPSx3() { return _EPSx3; }
  int& EPSk1() { return _EPSk1; }
  int& EPSk2() { return _EPSk2; }
  int& EPSk3() { return _EPSk3; }
  int& fi() { return _fi; }
  int& EL() { return _EL; }
  map<int, Entry>& e2dmap() { return _e2dmap; }
  //
  int setup2(iRSF& par, iRSF& inp);
  int setup32(iRSF& par, iRSF& inp);
  int setup3(iRSF& par, iRSF& inp);
  //
  int kernel2(int N, vector<Point2>& trg, vector<Point2>& src, CpxNumMat& res);
  int apkernel2(int N, vector<Point2>& trg, vector<Point2>& src, CpxNumMat& res, const float xx);
  int kernel3(int N, vector<Point3>& trg, vector<Point3>& src, CpxNumMat& res);
  int dikernel3(const int fi, const float tau, const float p, const float q, const float x, const float y, float& t);
  int kernel34(int N, vector<Point3>& trg, vector<Point3>& src, CpxNumMat& res, const float xx);
  //
  int check2(int N, const CpxNumMat& f, const FltNumVec& w, const FltNumVec& x, const CpxNumMat& u, const FltNumVec& tau, const FltNumVec& p, int NC, float& relerr);
  int apcheck2(int N, const CpxNumMat& f, const FltNumVec& w, const FltNumVec& x, const CpxNumMat& u, const FltNumVec& tau, const FltNumVec& p, const float xx, int NC, float& relerr);
  int check3(int N, const CpxNumTns& f, const FltNumVec& w, const FltNumVec& x, const FltNumVec& y, const CpxNumTns& u, const FltNumVec& tau, const FltNumVec& p, const FltNumVec& q, int NC, float& relerr);
  int check34(int N, const CpxNumTns& f, const FltNumVec& w, const FltNumVec& x, const FltNumVec& y, const CpxNumTns& u, const FltNumVec& tau, const FltNumVec& p, const FltNumVec& q, const float xx, int NC, float& relerr);
  //
  int prep_aux(FltNumVec& grid, vector<float>& ts, CpxNumMat& tmp);
  int eval_addaux(const CpxNumTns& ext, CpxNumTns& all, CpxNumMat& m1, CpxNumMat& m2, CpxNumMat& m3);
  //
  int eval2(int N, const CpxNumMat& f, const FltNumVec& w, const FltNumVec& x, CpxNumMat& u, const FltNumVec& tau, const FltNumVec& p);
  int apeval2(int N, const CpxNumMat& f, const FltNumVec& w, const FltNumVec& x, CpxNumMat& u, const FltNumVec& tau, const FltNumVec& p, const float xx);
  int eval3(int N, const CpxNumTns& f, const FltNumVec& w, const FltNumVec& x, const FltNumVec& y, CpxNumTns& u, const FltNumVec& tau, const FltNumVec& p, const FltNumVec& q);
  int eval34(int N, const CpxNumTns& f, const FltNumVec& w, const FltNumVec& x, const FltNumVec& y, CpxNumTns& u, const FltNumVec& tau, const FltNumVec& p, const FltNumVec& q, const float xx);
  //
};

#endif



