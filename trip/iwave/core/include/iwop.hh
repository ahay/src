#ifndef __IWAVE_OP
#define __IWAVE_OP

#define DEFAULT_SNAPS 10

//#include "alg.hh"
#include "op.hh"
#include "ocdc.hh"
#include "gridpp.hh"
#include "segypp.hh"
#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#include "mpisegypp.hh"
#endif
#include "istate.hh"

namespace TSOpt {

  //using namespace RVLAlg;
  using RVLAlg::ListAlg;
  using RVL::DataContainer;
  using RVL::ProductDataContainer;
  using RVL::StdProductDataContainer;
  using RVL::Space;
  using RVL::SpaceDCF;
  using RVL::ProductSpace;
  using RVL::ConstContainer;
  using RVL::STRING_PAIR;
  using RVL::Vector;
  using RVL::Components;
  using RVL::FunctionObject;
  using RVL::Operator;
  using RVL::Writeable;
  using RVL::AssignParams;
  
  class IWaveSpace: public ProductSpace<ireal> {

  private:

    /** vector of const pointers to Space */
    std::vector< Space<ireal> const * > _s;
    std::vector< std::string > _keys;
    
  public:

    IWaveSpace(PARARRAY const & par, 
	       IWaveInfo const & ic,
	       bool input,
	       ostream & outfile = cerr);

    IWaveSpace(IWaveSpace const & sp);
    
    ~IWaveSpace();

    /** implements virtual DataContainer constructor via
	StdProductDataContainer class. */
    DataContainer * buildDataContainer() const;

    size_t getSize() const;
    Space<ireal> const & operator[](size_t i) const;
    std::vector<std::string> getKeys() const;
  };

  class IWaveOp: public Operator<ireal>  {
      
  private:

    IWaveInfo ic;
    IWaveSpace dom;
    IWaveSpace rng;
    mutable FILE * stream;              /* output stream            */
    PARARRAY * pars;            /* parameter array ref copy */

    // verbosity control
    int dump_steps;
    int dump_pars;
    int dump_term;

    // other verbosity control handled within iwave code
    
    // dry run option
    bool dryrun;
    ostream & drystr;
    
    // verbose announcements
    ostream & announce;

    // convenience filename transfer - weak sanity check, presume 
    // that membership is already established
    void param_set(RVL::Vector<ireal> const & x, 
		   PARARRAY & pars, 
		   IWaveSpace const & sp,
		   std::string const & suf,
		   FILE * stream=stderr) const;

    IWaveOp();
      
  protected:
      
    void apply(const Vector<ireal> & x, 
	       Vector<ireal> & y) const;

    void applyDeriv(const Vector<ireal> & x, 
		    const Vector<ireal> & dx,
		    Vector<ireal> & dy) const;
      
    void applyAdjDeriv(const Vector<ireal> & x, 
		       const Vector<ireal> & dy,
		       Vector<ireal> & dx) const;
      
    void applyDeriv2(const Vector<ireal> & x, 
		     const Vector<ireal> & dx0,
		     const Vector<ireal> & dx1,
		     Vector<ireal> & dy) const;

    void applyAdjDeriv2(const Vector<ireal> & x, 
			const Vector<ireal> & dx0,
			const Vector<ireal> & dy,
			Vector<ireal> & dx1) const;

    Operator<ireal> * clone() const;
    
  public:
    
    IWaveOp(PARARRAY _pars, FILE * _stream,
	    bool _dryrun=false,
	    ostream & _drystr=cerr,
	    ostream & _announce=cerr);
    
    IWaveOp(IWaveOp const & x);

    ~IWaveOp();

    const Space<ireal> & getDomain() const;
    const Space<ireal> & getRange() const;

    // added 23.06.10 to facilitate using source as variable
    // without admitting that it's part of domain
    PARARRAY & getPar();
    PARARRAY const & getPar() const;

    ostream & write(ostream & str) const;
  };
}

#endif
