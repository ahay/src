#ifndef __TSOPT_GRIDPP
#define __TSOPT_GRIDPP

//#define FRUITCAKE
//#define VERBOSE

#define MAX_NAMELEN 128

#include "space.hh"
#include "locallinalg.hh"
#include "ocdc.hh"
#include "parserdecl.hh"
#include "iwave_fopen.h"
#include "grid.h"
#include "gridio.h"
#include "rarray.h"
#include "write.hh"

namespace RVL {

  template<>
  size_t getDataSize<RARR>(RARR const & a);

  template<>
  ireal * newData<ireal,RARR>(RARR & md);

  template<>
  void deleteData<ireal,RARR>(ireal ** d,RARR ** md);

  template<>
  ostream & writeMeta<RARR>(RARR const & md, ostream & s);

}

namespace TSOpt {

  using RVL::OCDC;
  using RVL::ConstContainer;
  using RVL::STRING_PAIR;
  using RVL::RVLException;
  using RVL::ScalarFieldTraits;
  using RVL::DataContainer;
  using RVL::ContentPackage;
  using RVL::PackageContainer;
  using RVL::DataContainerFactory;
  using RVL::PackageContainerFactory;
  using RVL::getDataSize;
  using RVL::Writeable;
  using RVL::StdSpace;
  using RVL::LinearAlgebraPackage;
  using RVL::RVLLinearAlgebraPackage;
  using RVL::parse;

  class GridDC: public OCDC<ireal, RARR> {

  private:

    /* constructor argument copies */
    string const & protohdr;  // name of prototype RSF header file
    string const & protodata; // name of prototype RSF data file
    grid const & protog;      // created from prototype RSF header file
    string data_format;       // read from prototype RSF header file
    ostream & outfile;        // verbose output stream
    
    /* internally generated in open_p - must be mutable */
    mutable FILE * fph;       // header file stream
    mutable FILE * fp;        // data file stream
    mutable string datafile;  // guess

    // read buffer
    mutable ContentPackage<ireal,RARR> buf;
    mutable int panelindex;   // chunk index
    int panelnum;             // number of chunks
    mutable bool rd;          // set if last op is read, unset on write

    GridDC();
    GridDC(GridDC const &);

  protected:

    /** private file access method - used in reset */

    void open_p() const;
    ContentPackage<ireal,RARR> & get(bool & more);
    ContentPackage<ireal,RARR> const & get(bool & more) const;
    void put() const;
    void reset() const;

  public:
    
    /** standalone constructor takes prototype RSF header
	filename. NOTE: in this version, prototype RSF data file must
	also exist and have correct length. Since all internal data
        depend on this file, consistency is guaranteed. 

	If incore flag is set, then g.dim = g.gdim and entire
	data volume is read/written, rather than panel-by-panel.
    */
    GridDC(string const & _protohdr,
	   string const & _protodata,
	   grid const & _protog,
	   string const & _data_format,
	   string data_type,
	   bool incore = false,
	   ostream & _outfile = cerr);
    /** destructor */
    ~GridDC();
    bool isInCore() const;
    string getProtohdr() const;
    ostream & write(ostream & str) const;

  };
  
  /** Factory class for GridDC. */
  
  class GridDCF: public PackageContainerFactory<ireal, RARR> {
    
  private:
    
    string protohdr;        // name of prototype RSF metadata file
    string protodata;       // name of prototype RSF data file
    string data_format;     // data format
    string data_type;       // data type tag
    ostream & outfile;      // externally supplied verbose output unit
    grid protog;            // grid determined by metadata
    ireal scfac;            // scale factor (from units)
    bool incore;            // RVL::Vector ops done incore
    ireal vol;              // cell volume
    
    GridDCF();

  protected:

    GridDC * buildGridDC() const;
    PackageContainer<ireal, RARR> * buildPC() const;

  public:

    /** main constructor. Stores name of prototype header 
	file (first arg), also extracts filename information 
	and ensures that prototype files register in the file
	manager database. Proto data extracted here is passed
	to GridDC constructor. */

    GridDCF(string _protohdr, bool _incore=false, ostream & _outfile = cerr);

    /** copy constructor. Since some copy of arg must have been built
	with main constructor, no need to worry about file manager. */
    GridDCF(GridDCF const & f);
    ~GridDCF();
    
    PackageContainerFactory<ireal, RARR> * clone() const;
    grid const & getGrid() const;
    ireal getCellVol() const;
    ireal getScaleFactor() const;
    string getFilename() const;
    bool compare( DataContainerFactory const & dcf) const;

    // compatible if either DFC incore flag is unset (whether data is incore or not)
    // or DFC incore flag is set and DC data is incore for any reason (i.e. internal
    // grid has dim=gdim).
    bool isCompatible( DataContainer const & dc ) const;
    bool isIncore() const;
    ostream & write(ostream & str) const;

  };

  /** Space class for Grid data */
  class GridSpace: public StdSpace<ireal,ireal>,
		   public ConstContainer<STRING_PAIR> {
    
  private:
    
    GridDCF f; 
    RVLLinearAlgebraPackage<ireal> l;
    
    /** function to build access STRING_PAIR (tag=hdr) from basic data
	note that the thdr string is passed by non-const reference, so is
	possibly altered on return - this is for use in the constructor,
	which requires the header file name.
    */
    STRING_PAIR getPTblEntry(grid const & g,
			     std::string tag,
			     std::string & thdr,
			     std::string fmt);
  public:

    /** normal constructor - assumes hdr is name of RSF header file 

	data type serves as key in IWAVE parameter tables. will eventually
	be deduced from units specified in metadata.
 
	optional specification of incore vector arithmetic - causes
	grid physical dimension to be overwritten in DC class by grid
	global dimension, so that entire data volume is read into memory
	rather than physical panels.

	@param[in] hdr      - path to prototype RSF metadata file
	@param[in] dtype    - data type (density, bulkmod, velocity,...)
	@param[in] incore   - implement vector arithmetic in core if true
	@param[in] outfile  - verbose output stream

    */
    GridSpace(string hdr, 
	      string dtype = "notype", 
	      bool incore = false,
	      ostream & outfile = cerr);

    /** constructor from grid data. useful when many intermediate spaces
	with same or related grids must be constructed.
	@param[in] g   - prototype grid
	@param[in] tag - data type (standin for units)
	@param[in] thdr - optional header file name - default is tmp name
	@param[in] fmt - optional data format - default is native_ireal
	@param[in] outfile - verbose output stream
    */
    GridSpace(grid const & g,
	      std::string tag,
	      std::string thdr="",
	      std::string fmt="native_ireal",
	      bool incore=false,
	      ostream & outfile = cerr);
    GridSpace(GridSpace const & sp);
    ~GridSpace();

    LinearAlgebraPackage<ireal> const & getLAP() const;
    DataContainerFactory const & getDCF() const;
    grid const & getGrid() const;
    bool isIncore() const;
    ostream & write(ostream & str) const;
  };
}	      
      
#endif
