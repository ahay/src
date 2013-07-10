#ifndef __TSOPT_GRIDPP
#define __TSOPT_GRIDPP

//#define FRUITCAKE
//#define VERBOSE

#define MAX_NAMELEN 128

#include "space.hh"
#include "locallinalg.hh"
#include "ocdc.hh"
#include "axispp.hh"
#include "axisdec.hh"
#include "parserdecl.hh"
extern "C" {
#include "iwave_fopen.h"
}

namespace TSOpt {

  using RVL::writeMeta;
  using RVL::ContentPackage;
  using RVL::ScalarFieldTraits;
  using RVL::RVLException;
  using RVL::Writeable;
  using RVL::parse;

  /** as a data structure, a Grid is simply a ContentPackage with
      int metadata and Axis data. The specialization includes 
      a specific initialization-on-construction. Since Grid is
      going to function as metadata for an o/c PackageContainer,
      do not bother to equip it with the standard incore CP functions.
  */

  template<typename Scalar>
  class Grid: public ContentPackage<Axis<Scalar>,int> {

  private:

    Axis<Scalar> exd;

  public:
    Grid()
      : ContentPackage<Axis<Scalar>,int>() {
      this->initialize(RARR_MAX_NDIM);      
      for (int i=0;i<RARR_MAX_NDIM;i++) this->getData()[i].id=i+1;
      exd.id=RARR_MAX_NDIM+1;
    }
    Grid(Grid<Scalar> const & g)
      : ContentPackage<Axis<Scalar>,int>(g), exd(g.exd) {
      if (!g.isInitialized())
	this->initialize(RARR_MAX_NDIM);
      for (int i=0;i<RARR_MAX_NDIM;i++)	this->getData()[i].id
	=g.getData()[i].id;
    }
    ~Grid() {}

    /** bounds-checked axis access */
    Axis<Scalar> const & getAxis(int i) const {
      if (i<0 || i > RARR_MAX_NDIM-1) {
	RVLException e;
	e<<"Error: Axis::getAxis\n";
	e<<"requested index = "<<i<<" out of bounds [0,"<<RARR_MAX_NDIM<<"\n";
	throw e;
      }
      return this->getData()[i];
    }

    Axis<Scalar> const & getExdAxis() const {
      return exd;
    }
    Axis<Scalar> & getExdAxis() {
      return exd;
    }

    /** dimension */
    int getDim() const {
      for (int i=RARR_MAX_NDIM-1;i>-1;i--) 
	if (this->getData()[i].n > 1) return i+1;
      return 0;
    }

    int getExdDim() const {
      if (exd.n > 1) return getDim()+1;
      return getDim();
    }

    /** size of corresponding data array - should probably be size_t */
    int getDataSize() const {
      int j=1;
      for (int i=0;i<RARR_MAX_NDIM;i++)
	j *= this->getData()[i].n;
      j *= exd.n;
      return j;
    }

    /** cell volume */
    Scalar getCellVol() const {
      try {
	Scalar vol=ScalarFieldTraits<Scalar>::One();
	int dim=this->getDim();
	for (int i=0;i<dim;i++)
	  vol *= this->getData()[i].d;
	// comment out at this time, for no extension scale in segy-space now
	// if (exd.n > 1) vol *= this->getExdAxis().d;
	return vol;
      }
      catch (RVLException & e) {
	e<<"\ncalled from Grid::getCellVol\n";
	throw e;
      }
    }
	
    /** comparison */
    bool operator==(Grid<Scalar> const & g) const {
      int dim=this->getDim();
      for (int i=0;i<dim;i++) 
	if (this->getAxis(i)!=g.getAxis(i)) return false;
      if (exd != g.exd) return false;
      return true;
    }
    bool operator!=(Grid<Scalar> const & a) const { return !operator==(a); }

    /** print method to reproduce RSF file */
    void fprint(FILE * fp) const {
      if (!fp) {
	RVLException e;
	e<<"Error: Axis::fprint\n";
	e<<"file pointer is null\n";
	throw e;
      }
      for (int i=0;i<RARR_MAX_NDIM;i++) {
	fprintf(fp,"n%d=%d\n",i+1,this->getData()[i].n);
	fprintf(fp,"d%d=%e\n",i+1,this->getData()[i].d);
	fprintf(fp,"o%d=%e\n",i+1,this->getData()[i].o);	
	fprintf(fp,"id%d=%d\n",i+1,this->getData()[i].id);	
      }
      if (exd.n > 1) {
	fprintf(fp,"en=%d\n",this->exd.n);
	fprintf(fp,"ed=%e\n",this->exd.d);
	fprintf(fp,"eo=%e\n",this->exd.o);	
	fprintf(fp,"eid=%d\n",this->exd.id);	      
      }
    }

    /** overload of CP::write */
    ostream & write(ostream & str) const {
      if (this->isInitialized()) {
	str<<"Grid axis list (non-defaults only)\n";
	for (int i=0;i<RARR_MAX_NDIM;i++) {
	  if (this->getData()[i].n > 1) 
	    writeMeta<Axis<Scalar> >(this->getData()[i],str);
	}
	if (exd.n > 1) 
	  writeMeta<Axis<Scalar> >(exd,str);	    
      }
      else {
	str<<"Grid - not initialized\n";
      }
      return str;
    }

    void readPar(PARARRAY par) {
      string key;
      for (int i=0;i<RARR_MAX_NDIM;i++) {
	std::ostringstream str;
	str<<i+1;
	key="n"+str.str();
	/* if you can read n, then try for the rest */
	if (parse<int>(par,key,this->getData()[i].n)) {
	  key="d"+str.str();
	  parse<Scalar>(par,key,this->getData()[i].d);
	  key="o"+str.str();
	  parse<Scalar>(par,key,this->getData()[i].o);
	  key="id"+str.str();
	  parse<int>(par,key,this->getData()[i].id); 
	}
      }
      key="en";
      if (parse<int>(par,key,this->exd.n)) {
	key="ed";
	parse<Scalar>(par,key,this->exd.d);
	key="eo";
	parse<Scalar>(par,key,this->exd.o);
	key="eid";
	parse<int>(par,key,this->exd.id); 
      }
    }

    void readFile(string filename) {
      PARARRAY * par = ps_new();
      if (ps_createfile(par,filename.c_str())) {
	RVLException e;
	e<<"Error: Grid::readFile\n";
	e<<"failed to create param array from file "<<filename<<"\n";
	throw e;
      }
      readPar(*par);
      ps_delete(&par);
    }
  };

}

#include "griddec.hh"

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

  template<typename Scalar>
  class GridDC: public OCDC<Scalar, Grid<Scalar> > {

  private:

    /* inherited from OCDC mixin: filename */
    string hdr;               // name of prototype RSF header file
    string protodata;         // name of prototype RSF data file
    PARARRAY * par;           // created from prototype RSF header file
    Grid<Scalar> g;           // created from prototype RSF header file
    int ns;                   // length of first axis (chunk)
    off_t off_eof;            // offset of proto data eof, in bytes
    // read buffer, refs returned by get
    ContentPackage<Scalar, Grid<Scalar> > * buf; 
    mutable IPNT ioff;        // integer offset tuple
    mutable int exdoff;       // extended axis offset index
    mutable off_t off_cur;    // nominal offset of current axis 0 copy
    mutable FILE * fph;       // header file pointer
    mutable FILE * fp;        // stream for data i/o
    mutable string datafile;  // data file name
    mutable bool rd;          // set if last op is read, unset on write
    ostream & outfile;        // verbose output stream

    GridDC();
    GridDC(GridDC<Scalar> const &);

    bool update_ioff() const {
      int carry=1;
      for (int i=1;i<g.getDim();i++) {
	ioff[i]+=carry;
	if (ioff[i]<g.getData()[i].n) carry=0;
	else {
	  carry=1;
	  ioff[i]=0;
	}
      }
      exdoff+=carry;
      if (exdoff<g.getExdAxis().n) carry=0;

      if (carry) return false;
      return true;
    }
	
  protected:

    /** private file access method - used in reset */

    void open_p() const {

      // always an error if data file is already open
      if (fp) {
	RVLException e;
	e<<"Error: GridDC::open_p\n";
	e<<"error state on call to open_p: file pointer set\n";
	throw e;
      }

      char * name = NULL;

      // temp case
      
      if (this->getFilename() == "") {

        // First, deal with header file 
#ifdef FRUITCAKE
	outfile<<"GridDC::open_p: attempt to open new temp file for hdr\n";
	outfile<<"  proto = "<<hdr<<endl;
#endif
	fph=iwave_fopen(&name,"w+",hdr.c_str(),stderr);
	
	if (!fph) {
	  RVLException e;
	  e<<"Error: GridDC::open_p\n";
	  e<<"temp file case: failed to open header file\n";
	  throw e;
	}
	if (!name) {
	  RVLException e;
	  e<<"Error: GridDC::open_p\n";
	  e<<"temp file case: failed to assign name\n";
	  throw e;
	}

	// postcondition - file is open with w+ access,
	// contents = copy of hdr
	// set filename
	string sname = name;
	this->setFilename(sname);
	userfree_(name);

#ifdef FRUITCAKE
	outfile<<"open_p: temp hdr file "<<this->getFilename()<<" opened on file ptr "<<fph<<"\n";
#endif
	
	// data file name will be temp name returned by iwave_fopen
	name=NULL;

#ifdef FRUITCAKE
	outfile<<"GridDC::open_p: attempt to open new temp file for data\n";
	outfile<<"  proto = "<<protodata<<endl;
#endif
	// open data file 
	fp=iwave_fopen(&name,"w+",protodata.c_str(),stderr);
	if (!fp) {
	  RVLException e;
	  e<<"Error: GridDC::open_p\n";
	  e<<"temp file case: failed to open data file = "<<datafile<<"\n";
	  throw e;
	}

	// record data name, both in string buffer and in header file
	datafile=name;
	string inkey="in";

#ifdef FRUITCAKE
	outfile<<"open_p: temp data file "<<datafile<<" opened on file ptr "<<fp<<"\n";
#endif

	/* add new in=datafile pair at bottom */
	if (ps_slcstring(*par,inkey.c_str(),datafile.c_str())) {
	  RVLException e;
	  e<<"Error GridDC::open_p\n";
	  e<<"failed to add in=... line to header file\n";
	  e<<"new archival file case\n";
	  throw e;
	}

	// overwrite header file - adding on to end - repeats lots 
	// of stuff but so what
	fseek(fph,0L,SEEK_END);
	fprintf(fph,"GRID++ new file branch\n\n");
	//	ps_printall(*par,fph);
	fprintf(fph,"in=%s\n",datafile.c_str());
	fflush(fph);
	
	// this is dangerous - relies on all access to header file
	// reading only last instance of key in.

	// even more dangerous - must flush! else will be buffered, and
	// file may be used before write - eg. with different file ptr
	fflush(fph);

	userfree_(name);

      }

      else {
	
	// archival file - new or old
	
	name=(char *)usermalloc_((this->getFilename().size()+1)*sizeof(char));
	strcpy(name,this->getFilename().c_str());

	// first try old-file open
	fph=iwave_fopen(&name,"r+",hdr.c_str(),stderr);
	
	if (!fph) {
#ifdef FRUITCAKE
	  outfile<<"open_p: new archival file branch, file = "
		 <<this->getFilename()<<endl;
#endif

	  /* create header file with new "in=" line at top
	     start with data file name and key */
	  datafile=this->getFilename()+"@";
	  string inkey="in";

	  /* new in=datafile pair at bottom */
	  if (ps_slcstring(*par,inkey.c_str(),datafile.c_str())) {
	    RVLException e;
	    e<<"Error GridDC::open_p\n";
	    e<<"failed to add in=... line to header file\n";
	    e<<"new archival file case\n";
	    throw e;
	  }

	  /* open header file for write */
	  if (!(fph=iwave_fopen(&name,"w",hdr.c_str(),stderr))) {      
	    RVLException e;
	    e<<"Error: GridDC::open_p\n";
	    e<<"failed to open file "<<this->getFilename()
	     <<" for read/write\n";
	    throw e;
	  }

	  /* overwrite with new parray */
	  fseek(fph,0L,SEEK_END);
	  fprintf(fph,"GRID++ new file branch\n\n");
	  //	  ps_printall(*par,fph);
	  fprintf(fph,"in=%s\n",datafile.c_str());
	  fflush(fph);

	  /* open data file for write + read */

	  userfree_(name);
	  name=(char *)usermalloc_((datafile.size()+1)*sizeof(char));
	  strcpy(name,datafile.c_str());

	  // open data file 
	  fp=iwave_fopen(&name,"w+",protodata.c_str(),stderr);
	  if (!fp) {
	    RVLException e;
	    e<<"Error: GridDC::open_p\n";
	    e<<"temp file case: failed to open data file = "<<datafile<<"\n";
	    throw e;
	  }
	  userfree_(name);
	  
	}
      
	/* old file option - check for grid compatibility, extract
	   datafile name, test for readability. Note that any additional
	   compatibility tests with header must be inserted by hand. */
	
	else {

#ifdef FRUITCAKE
	  outfile<<"open_p: old file branch, file = "<<this->getFilename()<<endl;
#endif
	  /* compatibility check inescapably requires Grid comparison.
	     NOTE that other attributes (data type, units, ...) are NOT
	     CHECKED.- only attributes of Grid type. This could be changed 
	     later. 
	  */
	  // 10.04.10: no longer needed, as ref grid is member data
	  //	  Grid<Scalar> g;
	  //	  g.readPar(par);
	  PARARRAY * partest  = ps_new();
	  if (ps_createfile(partest,this->getFilename().c_str())) {
	    RVLException e;
	    e<<"Error: GridDC::open_p\n";
	    e<<"failed to read file "<<this->getFilename()<<" into PARARRAY\n";
	    throw e;
	  }
	  Grid<Scalar> gtest;
	  gtest.readPar(*partest);
	  if (gtest != g) {
	    RVLException e;
	    e<<"Error: GridDC::open_p\n";
	    e<<"input RSF file "<<this->getFilename()<<" incompatible with\n";
	    e<<"prototype grid:\n";
	    g.RVL::Writeable::write(e);
	    e<<"\ngrid from "<<this->getFilename()<<":\n";
	    gtest.RVL::Writeable::write(e);
	    throw e;
	  }
	
	  /* data file must exist in this branch - parse it out and
	     test with fopen - fails if not found with correct prototype
	     and r+ access
	  */
	  string inkey="in";
	  if (!(parse<string>(*partest,inkey,datafile))) {
	    RVLException e;
	    e<<"Error: GridDC::open_p\n";
	    e<<"input RSF file "<<this->getFilename()<<" failed parse key=in\n";
	    throw e;
	  }

	  // done with partest
	  ps_delete(&partest);

	  // open data file

	  userfree_(name);
	  name=(char *)usermalloc_((datafile.size()+1)*sizeof(char));
	  strcpy(name,datafile.c_str());

	  if (!(fp=iwave_fopen(&name,"r+",protodata.c_str(),stderr))) {      
	    RVLException e;
	    e<<"Error: GridDC::open_p\n";
	    e<<"failed to open file "<<datafile<<" for read/write\n";
	    e<<"identified as input file in header file"<<this->getFilename()<<"\n";
	    throw e;
	  }	
	  userfree_(name);
  
	  // this call initializes off_cur, off_eof properly
	  if (fseeko(fp,0L,SEEK_SET)) {
	    RVLException e;
	    e<<"Error: GridDC::open_p\n";
	    e<<"seek to start of file "<<datafile<<" failed\n";
	    throw e;
	  }
	  off_cur = ftello(fp);

	  // check data file length against proto
	  if (fseeko(fp,0L,SEEK_END)) {
	    RVLException e;
	    e<<"Error: GridDC::open_p\n";
	    e<<"seek to end of file "<<datafile<<" failed\n";
	    throw e;
	  }
	  // temp fix to accomodate DS
	  //	  if (off_eof != ftello(fp)) {
	  if (off_eof > ftello(fp)) {
	    RVLException e;
	    e<<"Error: GridDC::open_p\n";
	    e<<"ftello indicates wrong size of file "
	     <<datafile<<"\n";
	    e<<"is "<<(long)ftello(fp)<<" should be "<<(long)off_eof<<"\n";
	    throw e;
	  }

	  // seek back to file begin
	  if (fseeko(fp,off_cur,SEEK_SET)) {
	    RVLException e;
	    e<<"Error: GridDC::open_p\n";
	    e<<"seek to current position in file "
	     <<datafile<<" failed\n";
	    throw e;
	  }
	}
      }
    }

    ContentPackage<Scalar, Grid<Scalar> > & get(bool & more) {
      try {
	// record current file position
	off_cur=ftello(fp);

	// record current physical position
	for (int i=1;i<g.getDim();i++) 
	  (buf->getMetadata()).getData()[i].o 
	    = ioff[i]*g.getData()[i].d + g.getData()[i].o;
	buf->getMetadata().getExdAxis().o =
	  exdoff*g.getExdAxis().d + g.getExdAxis().o;

#ifdef FRUITCAKE 
	outfile<<"gridpp::get - read 1D buffer at location ";
	for (int i=0;i<g.getDim();i++) 
	  outfile<<"ioff["<<i<<"]="<<ioff[i]<<" o["<<i<<"]="
		 <<(buf->getMetadata()).getData()[i].o<<"\n";
#endif
	// read data into buffer
	int nb = fread(buf->getData(),sizeof(Scalar),ns,fp);

	// this should not happen
	if (nb != ns) {
	  RVLException e;
	  e<<"Error: GridDC::get_p\n";
	  e<<"did not read prototype number of samples\n";
	  e<<"proto = "<<ns<<" this = "<<nb<<"\n";
	  throw e;
	}

	// set rd flag for successful read
	rd=true;
	
	// determine whether this is last trace
	more=false;
	off_t off_ptr=ftello(fp);
	if ((off_ptr<off_eof) &&
	    update_ioff()) more=true;
#ifdef FRUITCAKE
	outfile<<"gridpp::get - end of read rk=0 off_eof="<<off_eof
	       <<" off_cur="<<off_cur<<" off_ptr="<<off_ptr<<" more="
	       <<more<<endl;
#endif
	return *buf;
      }
      catch (RVLException & e) {
	e<<"\ncalled from GridDC::get (mutable)\n";
	throw e;
      }
    }

    ContentPackage<Scalar, Grid<Scalar> > const & get(bool & more) const {
      try {
	// record current file position
	off_cur=ftello(fp);

	// record current physical position
	for (int i=1;i<g.getDim();i++) 
	  (buf->getMetadata()).getData()[i].o 
	    = ioff[i]*g.getData()[i].d + g.getData()[i].o;
	buf->getMetadata().getExdAxis().o =
	  exdoff*g.getExdAxis().d + g.getExdAxis().o;

#ifdef FRUITCAKE 
	outfile<<"gridpp::get - read 1D buffer at location ";
	for (int i=0;i<g.getDim();i++) 
	  outfile<<"ioff["<<i<<"]="<<ioff[i]<<" o["<<i<<"]="
		 <<(buf->getMetadata()).getData()[i].o<<"\n";
#endif
	// read data into buffer
	int nb = fread(buf->getData(),sizeof(Scalar),ns,fp);

	// this should not happen
	if (nb != ns) {
	  RVLException e;
	  e<<"Error: GridDC::get_p\n";
	  e<<"did not read prototype number of samples\n";
	  e<<"proto = "<<ns<<" this = "<<nb<<"\n";
	  throw e;
	}

	// set rd flag for successful read
	rd=true;
	
	// determine whether this is last trace
	more=false;
	off_t off_ptr=ftello(fp);
	if ((off_ptr<off_eof) &&
	    update_ioff()) more=true;
#ifdef FRUITCAKE
	outfile<<"gridpp::get - end of read rk=0 off_eof="<<off_eof
	       <<" off_cur="<<off_cur<<" off_ptr="<<off_ptr<<" more="
	       <<more<<endl;
#endif
	return *buf;
      }
      catch (RVLException & e) {
	e<<"\ncalled from GridDC::get (const)\n";
	throw e;
      }
    }

    void put() const {
      try {
	//if (!fp) open_p();

	// note that put does not require any update of buffer origin,
	// since buffer is left with same data - does not change column
	// stored in buffer. Only a read can change the data in the 
	// buffer.

	// if read flag is set, overwrite most recently read segment
	if (rd) {
	  if (fseeko(fp,off_cur,SEEK_SET)) {
	    RVLException e;
	    e<<"Error: GridDC::put_p\n";
	    e<<"seek failed on file "<<this->getFilename()<<"\n";
	    throw e;
	  }	
	}
	int nb=fwrite(buf->getData(),sizeof(Scalar),ns,fp);
	if (nb != ns) {
	  RVLException e;
	  e<<"Error: GridDC::put_p\n";
	  e<<"did not write prototype number of samples\n";
	  e<<"proto = "<<ns<<" this = "<<nb<<"\n";
	  throw e;
	}
	fflush(fp);
	off_cur=ftello(fp);
#ifdef FRUITCAKE
	outfile<<"GridDC::put: ns="<<ns<<" nb="<<nb<<" off_cur="<<off_cur<<"\n";
#endif
	// unset read flag
	rd=false;
      }
      catch (RVLException & e) {
	e<<"\ncalled from GridDC::put\n";
	throw e;
      }
    }    
    
    void reset() const { 
      try {
#ifdef FRUITCAKE
	outfile<<"reset: filename="<<this->getFilename()<<endl;
#endif
	if (!fp) open_p();
	if (fseeko(fp,0L,SEEK_SET)) {
	  RVLException e;
	  e<<"Error: GridDC::reset_p\n";
	  e<<"seek failed on file "<<datafile<<"\n";
	  throw e;
	}
	// reset origin in CP buffer metadata Grid to reference
	for (int i=0;i<g.getDim();i++) {
	  ioff[i]=0;
	  buf->getMetadata().getData()[i].o=g.getData()[i].o;
	}
	exdoff=0;
	buf->getMetadata().getExdAxis().o = g.getExdAxis().o;
      }
      catch (RVLException & e) {
	e<<"\ncalled from GridDC::reset\n";
	throw e;
      }
    }

  public:
    
    /** standalone constructor takes prototype RSF header
	filename. NOTE: in this version, prototype RSF data file must
	also exist and have correct length. Since all internal data
        depend on this file, consistency is guaranteed. */
    GridDC(string _hdr, ostream & _outfile = cerr)
      : OCDC<Scalar, Grid<Scalar> >(_outfile), 
	hdr(_hdr),
	datafile(""),
	fp(NULL),
	fph(NULL),
	rd(false),
	par(ps_new()),
	outfile(_outfile) {

#ifdef FRUITCAKE
      outfile<<"GridDC constructor, header file = "<<hdr<<endl;
#endif
      if (ps_createfile(par,hdr.c_str())) { 
	RVLException e;
	e<<"Error: GridDC constructor\n";
	e<<"failed to generate parameter array from proto header file "<<hdr<<"\n";
	throw e;
      }
#ifdef FRUITCAKE
      cerr<<"par file from "<<hdr<<endl;
      ps_printall(*par,stderr);
#endif
      // extract tag from header, insert in OCDC database
      // OK if none is present
      string tagstr;
      if (parse(*par,"data_type",tagstr)) this->setTag(tagstr);

#ifdef FRUITCAKE
      outfile<<"GridDC constr read grid\n";
#endif
      g.readPar(*par);
      ns=g.getAxis(0).n;
#ifdef FRUITCAKE
      outfile<<"GridDC constr buf\n";
#endif

      // set up prototype grid 'gtmp' for single column buffer
      Grid<Scalar> gtmp(g);           
      for (int i=1;i<g.getDim();i++) gtmp.getData()[i].n=1;
      gtmp.getExdAxis().n = 1;

      // set up CP for read of single column
      buf=new ContentPackage<Scalar, Grid<Scalar> >;

      // origin vector left initialized for first col
      buf->initialize(gtmp);

      // initialize buffer index array
      IASN(ioff,IPNT_0);
      exdoff=0;

      // compute eof from grid info
      off_eof=1;
      for (int i=0;i<RARR_MAX_NDIM;i++) 
	off_eof *= g.getAxis(i).n;
      off_eof *= g.getExdAxis().n;
      off_eof *= sizeof(Scalar);
      off_cur = 0L;

      // extract data prototype from hdr PARARRAY
      string inkey="in";
      if (!(parse<string>(*par,inkey,protodata))) {
	RVLException e;
	e<<"Error: GridDC constructor\n";
	e<<"input RSF file "<<hdr<<" failed parse key=in\n";
	throw e;
      }
      
#ifdef FRUITCAKE
      outfile<<"exit GridDC constr off_eof="<<off_eof<<endl;
#endif
    }
    
    /** destructor */
    ~GridDC() {
      delete buf;
      ps_delete(&par);
      if (fph) iwave_fclose(fph);
      if (fp) iwave_fclose(fp); 
    }

    string getHdr() const { return hdr; }

    ostream & write(ostream & str) const {
      str<<"Grid Data Container object \n";
      str<<"based on prototype header file "<<hdr<<"\n";
      str<<"header file = "<<this->getFilename()<<" pointer = "<<fph<<endl;
      str<<"data file = "<<datafile<<" pointer = "<<fp<<endl;
#ifdef VERBOSE
      bool more=true;
      this->reset();
      while (more) {
	ContentPackage<Scalar, Grid<Scalar> > const & ref = this->get(more);
	ref.write(str);
      }
#endif
      return str;
    }
  };

  /** Factory class for GridDC. */
  
  template<typename Scalar>
  class GridDCF: public PackageContainerFactory<Scalar, Grid<Scalar> > {

  private:

    string hdr;             // name of prototype RSF file
    ostream & outfile;      // externally supplied verbose output unit
    Grid<Scalar> g;         // grid determined by proto file
    Scalar vol;             // cell volume

    GridDCF();

  protected:

    GridDC<Scalar> * buildGridDC() const {
      return new GridDC<Scalar>(hdr,outfile);
    }

    PackageContainer<Scalar, Grid<Scalar> > * buildPC() const {
      return buildGridDC();
    }

  public:

    /** main constructor. Stores name of prototype header 
	file (first arg), also extracts filename information 
	and ensures that prototype files register in the file
	manager database. Involves repetition of some steps needed
	in construction of GridDC, so some data could be stored
	in this class and simply transferred to GridDCs constructed
	by build method. However the amount is small enough, and the
	opportunity for incoherency strictly controlled, so for now
	permit these classes to extract the same info from the 
	prototype header file independently. */

    GridDCF(string _hdr, ostream & _outfile = cerr)
      : hdr(_hdr), outfile(_outfile), vol(ScalarFieldTraits<Scalar>::Zero()) {

#ifdef IWAVE_USE_MPI
      if (retrieveGlobalRank() == 0) {
#endif
      cerr<<"rank = "<<retrieveGlobalRank()<<" GridDCF("<<hdr<<")\n";
      // sanity check
      char * buf = new char[MAX_NAMELEN];
      if (hdr.size()>MAX_NAMELEN-10) {
	RVLException e;
	e<<"Error: GridDCF constructor\n";
	e<<"proto header filename = "<<hdr<<" too large - > "<<MAX_NAMELEN-10<<"\n";
	throw e;
      }

      // assure presence in file mgr
      strcpy(buf,hdr.c_str());
      FILE * fp=NULL;
      if (!(fp=iwave_fopen(&buf,"r",NULL,stderr))) {
	RVLException e;
	e<<"Error: GridDCF constructor\n";
	e<<"failed to open proto header file "<<hdr<<" in file mgr\n";
	e<<"rank = "<<retrieveGlobalRank()<<"\n";
	throw e;
      }      
      iwave_fclose(fp);

      // extract parameter array
      PARARRAY * par = ps_new();	
      if (ps_createfile(par,hdr.c_str())) { 
	RVLException e;
	e<<"Error: GridDCF constructor\n";
	e<<"failed to generate parameter array from proto header file "<<hdr<<"\n";
	throw e;
      }

      // store grid, cell volume
      g.readPar(*par);
      vol=g.getCellVol();

      // done with pararray

      // extract prototype data file name
      string protodata="";
      string inkey="in";
      if (!(parse<string>(*par,inkey,protodata))) {
	RVLException e;
	e<<"Error: GridDCF constructor\n";
	e<<"proto RSF file "<<hdr<<" failed to parse key=in\n";
	e<<"rank = "<<retrieveGlobalRank()<<"\n";
	throw e;
      }      

      // assure presence in file mgr
      strcpy(buf,protodata.c_str());
      FILE * fpd=NULL;
      if (!(fpd=iwave_fopen(&buf,"r",NULL,stderr))) {
	RVLException e;
	e<<"Error: GridDCF constructor\n";
	e<<"failed to open proto data file "<<protodata<<" in file mgr\n";
	e<<"rank = "<<retrieveGlobalRank()<<"\n";
	throw e;
      }      
      iwave_fclose(fpd);

      // clean up
      delete [] buf;
      ps_delete(&par);

#ifdef IWAVE_USE_MPI
      }
      // cast to double, Bcast, assign
      double tvol = vol;
      MPI_Bcast(&tvol,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      vol = tvol;
#endif

    }

    /** copy constructor. Since some copy of arg must have been built
	with main constructor, no need to worry about file manager. */
    GridDCF(GridDCF<Scalar> const & f) 
      : hdr(f.hdr),
	outfile(f.outfile),
        g(f.g),
        vol(f.vol) {}
    
    ~GridDCF() {}

    PackageContainerFactory<Scalar, Grid<Scalar> > * clone() const {
      return new GridDCF<Scalar>(*this);
    }

    Grid<Scalar> const & getGrid() const { 
#ifdef IWAVE_USE_MPI
      if (retrieveGlobalRank()) {
	RVLException e;
	e<<"Error: GridDCF::getGrid rank="<<retrieveGlobalRank()<<"\n";
	e<<"should not be called except on rank 0\n";
	throw e;
      }
#endif
      return g;
    }

    Scalar getCellVol() const { return vol; }

    string getFilename() const { return hdr; }

    bool compare( DataContainerFactory const & dcf) const {
      GridDCF<Scalar> const * f = NULL;
      f = dynamic_cast<GridDCF<Scalar> const *>(&dcf);
      if ((f) && (this->hdr == f->hdr)) return true;
      return false;
    }

    bool isCompatible( DataContainer const & dc ) const {
      GridDC<Scalar> const * gc = NULL;
      gc = dynamic_cast<GridDC<Scalar> const *>(&dc);
      if (gc && gc->getHdr()==hdr) return true; 
      return false;
    }

    ostream & write(ostream & str) const {
      str<<"GridDCF: factory class for GridDC\n";
      str<<"precision: decimal digits = "<<std::numeric_limits<Scalar>::digits10<<"\n";
      str<<"Prototype RSF filename = "<<hdr<<"\n";
      return str;
    }
  };

  /** Space class for Grid data */
  template<typename Scalar>
  class GridSpace: public StdSpace<Scalar,Scalar>,
		   public ConstContainer<STRING_PAIR> {

  private:
    
    GridDCF<Scalar> f; 
    RVLLinearAlgebraPackage<Scalar> l;

    /** function to build access STRING_PAIR (tag=hdr) from basic data
     note that the thdr string is passed by non-const reference, so is
     possibly altered on return - this is for use in the constructor,
     which requires the header file name.
    */
    STRING_PAIR getPTblEntry(Grid<Scalar> const & g,
			     std::string tag,
			     std::string & thdr,
			     std::string fmt) {
      /** open file for construction */
      char * name = NULL;
      if (thdr != "") {
	name=(char *)usermalloc_(thdr.size()+1);
	strcpy(name,thdr.c_str());
      }
      FILE * fph = iwave_fopen(&name,"w",NULL,stderr);
      /* now have header name */
      if (fph) thdr=name;
      else {
	RVLException e;
	e<<"Error: GridSpace constructor from grid \n";
	e<<"failed to open header file = "<<name<<"\n";
	throw e;
      }	

      /* write grid to file */
      g.fprint(fph);
      /* write other stuff to file */
      fprintf(fph,"data_format = %s\n",fmt.c_str());
      fprintf(fph,"data_type   = %s\n",tag.c_str());

      /* the next block is only necessary if we keep the 
	 requirement that proto header files have to have
	 proto data files to go with them 
      */
      userfree_(name);
      name=NULL;
      FILE * fp = iwave_fopen(&name,"w",NULL,stderr);
      /* now have data name */
      fprintf(fph,"in          = %s\n",name);
      int n1 = g.getAxis(0).n;
      int nt = g.getDataSize();
      char * x = (char *)usermalloc_(n1*sizeof(char)*sizeof(float));
      memset(x,0,n1*sizeof(float));
      while (nt>0) {
	fwrite(x,sizeof(char),n1*sizeof(float),fp);
	nt-=n1;
      }
      userfree_(name);
      userfree_(x);
      if (ftello(fp) != g.getDataSize()*sizeof(float)) {
	RVLException e;
	e<<"Error: GridSpace constructor from grid \n";
	e<<"failed to initialize data file\n";
	throw e;
      }
      fseek(fp,0L,SEEK_SET);
      
      STRING_PAIR pr;
      pr.key=tag;
      pr.val=thdr;
      
      return pr;
    }

  public:

    /** normal constructor - assumes hdr is name of RSF header file */
    GridSpace(string hdr, string dtype = "notype", ostream & outfile = cerr)
      : StdSpace<Scalar,Scalar>(),
	//	ConstContainer<STRING_PAIR>(PairExtractor< Grid<Scalar> >::extract(hdr)),
	//	ConstContainer<STRING_PAIR>(extract_pair(hdr,"data_type")),
	ConstContainer<STRING_PAIR>(STRING_PAIR(dtype,hdr)),
	f(hdr,outfile), 
	l(f.getCellVol()) {
      //      cerr<<"GridSpace constructor: key="<<this->get().key<<" val="<<this->get().val<<"\n";
    }

    /** constructor from grid data. useful when many intermediate spaces
	with same or related grids must be constructed.
	@param[in] g   - prototype grid
	@param[in] tag - data type (standin for units)
	@param[in] thdr - optional header file name - default is tmp name
	@param[in] fmt - optional data format - default is native_float
	@param[in] outfile - verbose output stream
    */
    GridSpace(Grid<Scalar> const & g,
	      std::string tag,
	      std::string thdr="",
	      std::string fmt="native_float",
	      ostream & outfile = cerr)
      : StdSpace<Scalar,Scalar>(),
	ConstContainer<STRING_PAIR>(getPTblEntry(g,tag,thdr,fmt)),
	f(thdr,outfile), 
	l(f.getCellVol()) {}

    GridSpace(GridSpace<Scalar> const & sp) 
      : StdSpace<Scalar,Scalar>(), 
	ConstContainer<STRING_PAIR>(sp.get()),
	f(sp.get().val,sp.outfile), 
	l(f.getCellVol()) {}

    ~GridSpace() {}

    LinearAlgebraPackage<Scalar> const & getLAP() const { return l; }

    DataContainerFactory const & getDCF() const { return f; }

    Grid<Scalar> const & getGrid() const { return f.getGrid(); }

    ostream & write(ostream & str) const {
      str<<"GridSpace: StdSpace based on RSF grid data\n";
      str<<" - data type   = "<<this->get().key<<"\n";
      str<<" - header file = "<<this->get().val<<"\n";
      return str;
    }

  };
}	      
      
#endif
