#ifndef _BFIO_HPP_
#define _BFIO_HPP_

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
    DblNumVec _grid;
    NumVec<CpxNumMat> _mats;
    NumVec<CpxNumMat> _dirc;
public:
    Entry() {;}
    ~Entry() {;}
    DblNumVec& grid() { return _grid; }
    NumVec<CpxNumMat>& mats() { return _mats; }
    NumVec<CpxNumMat>& dirc() { return _dirc; }
};

int serialize(const Entry&, ostream&, const vector<int>&);
int deserialize(Entry&, istream&, const vector<int>&);

//-----------------------------
class BFIO: public ComObject
{
public:
    int _EPSx1;
    int _EPSx2;
    int _EPSk1;
    int _EPSk2;
    int _fi;
    map<int, Entry> _e2dmap;
    double tmin, tmax, pmin, pmax;
    double wmin, wmax, zmin, zmax; 
public:
    BFIO(const string& p): ComObject(p) {;}
    ~BFIO() {;}
    int& EPSx1() { return _EPSx1; }
    int& EPSx2() { return _EPSx2; }
    int& EPSk1() { return _EPSk1; }
    int& EPSk2() { return _EPSk2; }
    int& fi() { return _fi; }
    map<int, Entry>& e2dmap() { return _e2dmap; }
    //
    int setup(iRSF &par, iRSF &inp);
    int eval(int N, const CpxNumMat& f, CpxNumMat& u);
    int kernel(int N, vector<Point2>& trg, vector<Point2>& src, CpxNumMat& res);
    int check(int N, const CpxNumMat& f, const CpxNumMat& u, int NC, double& relerr);
};

#endif



