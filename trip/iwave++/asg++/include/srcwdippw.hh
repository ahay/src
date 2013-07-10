#ifndef __TSOPT_SRCWPW__
#define __TSOPT_SRCWPW__

#include "rnspace.hh"
#include "linop.hh"

#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#include "mpisegypp.hh"
#else
#include "gridpp.hh"
#include "segypp.hh"
#endif

#include "srcstk.hh"

namespace TSOpt {

  using RVL::RnSpace;
  using RVL::Vector;
  using RVL::LocalVector;
  using RVL::LinearOp;
  using RVL::RVLException;
  using RVL::AssignFilename;
  using RVL::AssignTag;
  using RVL::AssignParams;
  
  /* this class defines a rule for mapping an array of weights to a dipping array 
     source, each trace of which is a reference pulse scaled by a weight and with 
     different starting time computed from dipping factor. 
     So this operator generates planewave with different incidence angle.
     
     The fundamental relation between incidence angle and dipping factor used 
     in the following computation:
     let \xi be slowness (1/km), v be velocity at source points (km/s), 
         dip be time-dip of plane-wave (ms/trace), \theta be the incidence angle,
	 dx be the horizontal space interval between two nearby source array/trace (m),
     then 
         \xi = sin(\theta)/v , dip = dx * \xi
 */
  
  class SrcWDipPwOp: public LinearOp<float> {

  private:

    SrcData const & sd;       /* trace data object, creates gather files */
#if IWAVE_USE_MPI
    MPISEGYSpace rng;
#else    
    SEGYSpace rng;            /* space defined by CRG */
#endif
    RnSpace<float> dom;       /* space of weight vectors */
    segy w;                   /* wavelet */

    float dip;	              /* time-dip of plane-wave (ms/trace)	*/
    //  float max_offset;         /* dist(src_start, src_end)*/
    // int ct1,cx1;	      /* center of plane (sample and trace)	*/
    // float xi;                 /* slowness */
    // float dx;                 /* space interval between nearby source arrays*/
   
    SrcWDipPwOp();

  protected:

    LinearOp<float> * clone() const {
      return new SrcWDipPwOp(*this);
    }

    void apply(Vector<float> const & x,
	       Vector<float> & y) const;
    void applyAdj(Vector<float> const & x,
		  Vector<float> & y) const;

  public:
    
    SrcWDipPwOp(string wavelet,
		 SrcData const & sd, float _dip=0);

    SrcWDipPwOp(SrcWDipPwOp const & a)
      : sd(a.sd), rng(a.rng), dom(a.dom), dip(a.dip){}
    
    ~SrcWDipPwOp() {}

    Space<float> const & getDomain() const { return dom; }
    Space<float> const & getRange() const { return rng; }
    
    ostream & write(ostream & str) const;
  };
}
    
#endif
