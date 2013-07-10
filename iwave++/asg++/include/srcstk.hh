#ifndef __TSOPT_SRCSTACK__
#define __TSOPT_SRCSTACK__

#include "seamx_headers.hh"

#include "rnspace.hh"
#include "linop.hh"

#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#include "mpisegypp.hh"
#else
#include "gridpp.hh"
#include "segypp.hh"
#endif

namespace TSOpt {

  using RVL::RnSpace;
  using RVL::Vector;
  using RVL::LocalVector;
  using RVL::LinearOp;
  using RVL::RVLException;
  using RVL::AssignFilename;
  using RVL::AssignTag;
  using RVL::AssignParams;
  
  /* this class collects source gather data from a SEGY file
     containing headers for a survey, and checks that every source
     gather has the same number of traces. It creates SU files for a
     common source gather and a common receiver gather; in the latter,
     the source coordinates are transposed to receiver coordinates, so
     that the CRG can be used as a header file for an array source in
     IWAVE. For the same reason, an optional source for time data permits
     the CRG to be defined with different time sampling.
  */
  
  class SrcData {

  private:
    int nt;
    int ntr;
    int nsrc;
    int kmute;
    bool isfix;
    string line;
    string csg;
    string crg;
    string time;

    FILE * fpl;
    FILE * fpg;
    FILE * fpr;
    FILE * fpt;

    unsigned short nsmod;
    unsigned short dtmod;
    short t0mod;    

    SrcData();

  public:

    SrcData(string line, 
	    string csg, 
	    string crg, 
	    string time, 
	    float tmute=0.0);
    SrcData(SrcData const & sd)
      : nt(sd.nt),
	ntr(sd.ntr),
	nsrc(sd.nsrc),
	kmute(sd.kmute),
	isfix(sd.isfix),
	csg(sd.csg),
	crg(sd.crg),
	time(sd.time),
        fpl(sd.fpl),
        fpg(sd.fpg),
        fpr(sd.fpr),
        fpt(sd.fpt),
	nsmod(sd.nsmod),
	dtmod(sd.dtmod),
	t0mod(sd.t0mod) {}
    ~SrcData();
    int get_nt() const { return nt; }
    int get_ntr() const { return ntr; }
    int get_nsrc() const { return nsrc; }
    int get_kmute() const { return kmute; }
    bool isFixedSpread() const { return isfix; }

    string getLine() const { return line; }
    string getCSG() const { return csg; }
    string getCRG() const { return crg; }
    string getSrc() const { return time; }

    unsigned short get_ntsrc() const { return nsmod; }
    unsigned short get_dtsrc() const { return dtmod; }
    short get_t0src() const { return t0mod; }

    ostream & write(ostream & str) const;
  };

  class SrcOp: public LinearOp<float> {

  private:

    FILE * fpl;
    SrcData const & sd;
    SEGYSpace rng;
    RnSpace<float> dom;
    float * bufg;
    SrcOp();

  protected:

    LinearOp<float> * clone() const {
      return new SrcOp(*this);
    }

    void apply(Vector<float> const & x,
	       Vector<float> & y) const;
    void applyAdj(Vector<float> const & x,
		  Vector<float> & y) const;

  public:
    
    SrcOp(SrcData const & _sd)
      : sd(_sd), 
	rng(sd.getCSG()), 
	dom(sd.get_nsrc()),
	bufg(new float[sd.get_ntr()*sd.get_nt()]) {
      if (!(fpl=iwave_const_fopen(sd.getLine().c_str(),"r",NULL,stderr))) {
	RVLException e;
	e<<"Error: SrcOp constructor\n";
	e<<"failed to open file for line data - file = "<<sd.getLine()<<"\n";
	throw e;
      }
      if (!sd.isFixedSpread()) {
	RVLException e;
	e<<"Error: SrcOp constructor\n";
	e<<"line data is not fixed spread - file = "<<sd.getLine()<<"\n";
	throw e;
      }
    }

    SrcOp(SrcOp const & a)
      : fpl(a.fpl),
	sd(a.sd),
	dom(a.dom),
	rng(a.rng),
	bufg(new float[sd.get_ntr()*sd.get_nt()]) 
    {}
    
    ~SrcOp() { iwave_fclose(fpl); delete [] bufg; }

    Space<float> const & getDomain() const { return dom; }
    Space<float> const & getRange() const { return rng; }
    
    ostream & write(ostream & str) const;
  };
}
    
#endif
