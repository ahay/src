#ifndef __TSOPT_SEGYPP
#define __TSOPT_SEGYPP

#include "space.hh"
#include "locallinalg.hh"
#include "ocdc.hh"
#include "su.h"
#include "header.h"
#include "segy.h"
#include "parserdecl.hh"
#include "iwave_fopen.h"
#include "utils.h"

// uncomment for verbose output
//#define FRUITCAKE

/* forward declarations for CP<segy> functions, and specialization
   of pair extractor
*/

namespace RVL {

  template<>
  size_t getDataSize<segy>(segy const & md);
  template<>
  float * newData<float,segy>(segy & md);
  template<>
  void deleteData<float,segy>(float ** d,segy ** md);
  template<>
  ostream & writeMeta<segy>(segy const & md, ostream & s);
  /*
  template<>
  class PairExtractor< segy > {
  public:
    PairExtractor() {}
    PairExtractor(PairExtractor< segy > const &) {}
    ~PairExtractor() {}
    
    static STRING_PAIR extract(string in) {
      STRING_PAIR pr;
      pr.key="datafile";
      pr.val=in;
      return pr;
    }
  };
  */
}

// standalone global namespace function
int count_unique_flthdr(std::string fpath, std::string key, float tol=1.e-3);

namespace TSOpt {

  using RVL::getDataSize;
  using RVL::newData;
  using RVL::deleteData;
  using RVL::ContentPackage;
  using RVL::PackageContainer;
  using RVL::PackageContainerFactory;
  using RVL::OCDC;
  using RVL::STRING_PAIR;
  //  using RVL::PairExtractor;
  using RVL::ConstContainer;
  using RVL::Space;
  using RVL::StdSpace;
  using RVL::Vector;
  using RVL::DataContainer;
  using RVL::DataContainerFactory;  
  using RVL::FunctionObject;
  using RVL::FunctionObjectScalarRedn;
  using RVL::LinearAlgebraPackage;
  using RVL::RVLLinearAlgebraPackage;
  using RVL::RVLException;
  using RVL::Writeable;

  /** identify segy trace as ContentPackage */
  typedef ContentPackage<float,segy> segytrace;

  /* forward declaration */
  class segygen;

  /** OCDC specialization. Provides access to segy traces stored in
      files, as PackageContainer instance. A prototype filename is
      stored as object data, and all files used as out-of-core data
      for SEGYDC objects must match this prototype file in
      structure. File regarded as temporary if no name supplied by
      eval of AssignFilename prior to any other FO or FOR eval. In
      that case, file has w+ access and may be unlinked by
      iwave_fdestroy if this function is called in driver. If filename
      is supplied, then file may either exist or not. Constructor
      first attempts to open with r+ access. Total length is checked
      on open by iwave_fopen. If the file fails to open with r+ access
      (as will occur if the file does not exist, for example), or if
      total length does not match that of the prototype file, then
      regarded as a new file and opened with w+ permission instead -
      note that a side effect is that any existing data is truncated
      in this case! If open with r+ access succeeds, then number of
      samples per trace and time step are checked against prototype
      (for first trace - all traces assumed to be uniform in these
      headers).

      PackageContainer access functions (both flavors of get, also put
      and reset) implemented via calls to SU i/o functions. Uses
      ContentPackage constructed out of segy trace struct from SU,
      with data part of trace identified as data part of
      ContentPackage - that is, in this case the data is part of the
      metadata!

      To enable correct non-const FO evals, each trace is flagged as
      it is read - the flag signifies "read but not yet written". A
      write operation begins by seeking to the beginning of the trace,
      the offset of which is cached, if the flag is set (but not
      otherwise). The write unsets the flag. This design supports the
      purely sequential structure of PackageContainer function
      evaluation. A more flexible design with some degree of random
      access would require greater exposure and interpretation of file
      position information. For RTOp-type ops the present design is
      adequate. Operators implementing more complex interactions with
      trace data should probably be implemented as standalone
      out-of-core functions, on the IWAVE model.
  */
  class SEGYDC: public OCDC<float,segy> {

  private:

    string hdr;            // name of prototype SU file
    mutable segy tr;       // workspace for reads, proto for segytrace
    segytrace * buf;       // read buffer, get returns ref
    mutable FILE * fp;     // stream for data i/o
    mutable bool rd;       // set if last op is read, unset on write
    mutable off_t off_cur; // offset of current trace in bytes
    ostream & outfile;     // verbose output file

    FILE * fph;            // stream for prototype i/o
    off_t off_eof;         // prototype end-of-file offset in bytes
    int nb;                // prototype number of bytes per trace
    int ns;                // prototype number of samples per trace
    int dt;                // prototype sample rate

    SEGYDC();
    SEGYDC(SEGYDC const &);

  protected:

    void open_p() const;
    segytrace & get(bool & more);
    segytrace const & get(bool & more) const;
    void put() const;
    void reset() const;

  public:

    /** only legal constructor takes prototype filename and
	flag for unlinking tmp files */
    SEGYDC(string _hdr, ostream & _outfile = cerr);
    /** destructor */
    ~SEGYDC();

    string getHdr() const { return hdr; }
    ostream & write(ostream & str) const;

  };

  /** DC Factory methods for OC classes must depend only on
      initialization of the prototype filename string. All
      other class methods should be limited to rank 0 in 
      parallel execution, and DCF methods should not depend
      on them. This feature should probably be abstracted in
      a superclass.
  */
  class SEGYDCF: public PackageContainerFactory<float,segy> {

  private:

    /* initialized on construction */
    string hdr;        // name of prototype file
    ostream & outfile; // verbose output stream

    /* post-construction initialization */
    mutable int nb;    // number of bytes for all traces
    mutable int nt;    // number of time samples
    mutable float dt;  // time step (scale in L2 norm)
    mutable bool iflag;
    void init() const;

    SEGYDCF();

  protected:

    /** must initialize before invoking SEGYDC constr -
        registers hdr with file mgmt */
    SEGYDC * buildSEGY() const {
      init();
      return new SEGYDC(hdr,outfile);
    }
    PackageContainer<float,segy> * buildPC() const {
      return buildSEGY();
    }

  public: 

    /** bare constructor - does not initialize */
    SEGYDCF(string _hdr, ostream & _outfile=cerr);

    /** bare copy constructor - does not initialize */
    SEGYDCF(SEGYDCF const & f) 
      : hdr(f.hdr), outfile(f.outfile), iflag(false) {}

    ~SEGYDCF() {}

    /** bare clone method - does not initialize */
    SEGYDCF * cloneSEGY() const { 
      return new SEGYDCF(*this); 
    }

    PackageContainerFactory<float,segy> * clone() const { 
      return cloneSEGY();
    } 

    /** return header file name */
    string getFilename() const { return hdr; }
    
    /** return time step for scaling l2 inner product - initialization
        required */
    float getDt() const { init(); return dt; }

    /** return number of time samples - initialization required */
    int getNt() const { init(); return nt; }

    /** return number of bytes per trace - initialization required */
    int getNB() const { init(); return nb; }

    /** bare compare method - needs only hdr filename */
    bool compare( DataContainerFactory const & dcf) const;
    /** bare compatibility method - needs only hdr filename */
    bool isCompatible( DataContainer const & dc ) const;
    /** bare write method - needs only hdr filename */
    ostream & write(ostream & str) const;
  };

  /** Space class encapsulating SEGY data */
  class SEGYSpace: public StdSpace<float,float>,
		   public ConstContainer<STRING_PAIR> {
    
  private:
    
    SEGYDCF f;
    RVLLinearAlgebraPackage<float> l;

  public:

    SEGYSpace(string hdr, string key, ostream & outfile = cerr)
      : StdSpace<float,float>(),
	//	ConstContainer<STRING_PAIR>(PairExtractor< segy >::extract(hdr)),
	ConstContainer<STRING_PAIR>(STRING_PAIR(key,hdr)),
	f(hdr,outfile),
	l(f.getDt()) {}
    SEGYSpace(SEGYSpace const & sp) 
      : StdSpace<float,float>(sp),
	ConstContainer<STRING_PAIR>(sp.get()),
	f(sp.f), l(sp.l) {}
    ~SEGYSpace() {}

    // can only be called by other spaces, since base class operator new is protected
    //    SEGYSpace * clone() { return new SEGYSpace(*this); }

    DataContainerFactory const & getDCF() const { return f; }
    LinearAlgebraPackage<float> const & getLAP() const { return l; }

    /** return header file name */
    string getPrototypeFilename() const { return f.getFilename(); }
    
    /** return time step for scaling l2 inner product */
    float getDt() const { 
      try { 
	return f.getDt(); }
      catch (RVLException & e) {
	e<<"\ncalled from SEGYSpace::getDt\n";
	throw e;
      }
    }
    /** return number of time samples */
    int getNt() const { 
      try {
	return f.getNt(); }
      catch  (RVLException & e) {
	e<<"\ncalled from SEGYSpace::getNt\n";
	throw e;
      }
    }
	
    ostream & write(ostream & str) const;
  };    
    
}

#endif
