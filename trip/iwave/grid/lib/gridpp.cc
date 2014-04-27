#define MAX_NAMELEN 128

#include "gridpp.hh"
#include "gridfun.hh"

/**#define DUMP_RARR*/

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

  void GridDC::open_p() const {

    /* first task: file registration, including
       creation of temp files */

    char * name = NULL;

    // temporary file case: signalled by metadata filename
    // not initialized  
    if (this->getFilename() == "") {

      // First, deal with header file 
#ifdef FRUITCAKE
      outfile<<"GridDC::open_p: open new temp hdr file \n";
      outfile<<"  proto = "<<protohdr<<endl;
#endif
      // open metadata file
      fph=iwave_fopen(&name,"w+",protohdr.c_str(),stderr);
	
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
      // contents = copy of protohdr
      // set filename
      string sname = name;
      this->setFilename(sname);
      userfree_(name);

#ifdef FRUITCAKE
      outfile<<"  open_p: temp hdr file "<<this->getFilename()<<" opened on file ptr "<<fph<<"\n";
#endif
	
      // data file name will be temp name returned by iwave_fopen
      name=NULL;

#ifdef FRUITCAKE
      outfile<<"GridDC::open_p: open new temp data file\n";
      outfile<<"  proto = "<<protodata<<endl;
#endif
      // open data file 
      // NOTE: 28.05.13 iwave_fopen automatically adds DATAPATH to tmp file
      // names
      fp=iwave_fopen(&name,"w+",protodata.c_str(),stderr);
      if (!fp) {
	RVLException e;
	e<<"Error: GridDC::open_p\n";
	e<<"temp file case: failed to open data file\n";
        e<<"with prototype data "<<protodata<<"\n";
	throw e;
      }

      // claim real estate
      ireal * tmpbuf = new ireal[get_datasize_grid(protog)];
      memset(tmpbuf,0,get_datasize_grid(protog)*sizeof(ireal));
      for (int ip=0; ip<get_panelnum_grid(protog); ip++) {
	if (get_datasize_grid(protog) != 
	    (int)fwrite(tmpbuf,sizeof(ireal),get_datasize_grid(protog),fp)) {
	  RVLException e;
	  e<<"Error: GridDC::open_p\n";
	  e<<"  failed to write panel "<<ip<<" of "<<get_datasize_grid(protog)<<" zeroes to data file\n";
	  e<<"  "<<name<<"  with prototype data "<<protodata<<"\n";
	  throw e;
	}
      }
      delete [] tmpbuf;
      fseeko(fp,0L,SEEK_SET);

      // record datafile for future use
      datafile = name;
	
#ifdef FRUITCAKE
      outfile<<"  open_p: temp data file "<<datafile<<" opened on file ptr "<<fp<<"\n";
      fseeko(fp,0L,SEEK_SET);
      off_t fbeg = ftello(fp);
      fseeko(fp,0L,SEEK_END);
      outfile<<"  file length = "<<ftello(fp)-fbeg<<" bytes\n";
      fseeko(fp,0L,SEEK_SET);
#endif
      /*************************************************************
      //      if (retrieveGlobalRank()==0) {
	//	cerr<<"  open_p: temp hdr file "<<this->getFilename()<<
	//  " opened on file ptr "<<fph<<"\n";
	//	cerr<<"  open_p: temp data file "<<datafile<<
	//  " opened on file ptr "<<fp<<"\n";
	//	iwave_fprintall(stderr);
      //      }
      ************************************************************/

      // add new datafile name to end of header file
      fseek(fph,0L,SEEK_END);
      fprintf(fph,"\nGridDC::open_p - new temp file branch\n\n");
      fprintf(fph,"in=%s\n",name);
      fflush(fph);

      userfree_(name);

    }

    // existing file - new or old
	
    else {
	
      name=(char *)usermalloc_((this->getFilename().size()+1)*sizeof(char));
      strcpy(name,this->getFilename().c_str());

      // first try old-file open
      fph=iwave_fopen(&name,"r+",protohdr.c_str(),stderr);
	
      // new archival file - does not exist

      if (!fph) {
#ifdef FRUITCAKE
	outfile<<"GridDC::open_p: open new archival hdr file \n";
	outfile<<"  proto = "<<protohdr<<endl;
#endif

	/* open metadata file for write - note that this action
	   automatically copies prototype metadata file */
	if (!(fph=iwave_fopen(&name,"w+",protohdr.c_str(),stderr))) {      
	  RVLException e;
	  e<<"Error: GridDC::open_p\n";
	  e<<"failed to open file "<<this->getFilename()
	   <<" for write with prototype "<<protohdr<<"\n";
	  throw e;
	}

#ifdef FRUITCAKE
	outfile<<"  open_p: archival hdr file "<<this->getFilename()<<" opened on file ptr "<<fph<<"\n";
#endif

	/* open data file for write + read */
	/* 28.05.13: prepend DATAPATH - the following code implements "split", 
           essentially - it strips off the base path from a fully qualified path,
           if given, otherwise uses the given filename locally
	*/
	userfree_(name);
	char * c_dpath = getenv("DATAPATH"); // standard says: return value shall not be modified!
	//	cerr<<"datapath = "<<c_dpath<<endl;
	//        cerr<<"filename = "<<this->getFilename()<<endl;
	// split up filename to find root
	char * fname = (char *)usermalloc_((this->getFilename().size()+1)*sizeof(char));
	strcpy(fname,this->getFilename().c_str());
	//	fprintf(stderr,"fname = %s\n",fname);
	char * c_dname = NULL;
	char * jnk = strtok(fname, "/"); // no memory allocated here!
	while (jnk) {
	  //	  fprintf(stderr,"token = %s\n",jnk);
	  jnk = strtok(NULL, "/");
          // if more than one segment found, use last
	  if (jnk) c_dname = jnk;
          // else use filename
	}
        if (!c_dname) c_dname=fname;
	//	fprintf(stderr,"c_dname = %s\n",c_dname);
	//        fprintf(stderr,"last non-NULL token = %s\n",c_dname);
        string dpath(c_dpath);
	string dname(c_dname);
	userfree_(fname);
	dname = dpath + dname;
	name=(char *)usermalloc_((dname.size()+2)*sizeof(char));
	strcpy(name,dname.c_str());
	strcat(name,"@");
	/////
	//	cerr<<"name = "<<name<<endl;

#ifdef FRUITCAKE
	outfile<<"GridDC::open_p: open new archival file for data\n";
	outfile<<"  proto = "<<protodata<<endl;
#endif
	// open data file 
	fp=iwave_fopen(&name,"w+",protodata.c_str(),stderr);
	if (!fp) {
	  RVLException e;
	  e<<"Error: GridDC::open_p\n";
	  e<<"new archival file case: failed to open data file = "<<name<<"\n";
	  e<<"with prototype "<<protodata<<" mode=w+\n";
	  throw e;
	}

	// record data file name
	datafile = name;

	/* add in = datafile at end */
	fseek(fph,0L,SEEK_END);
	fprintf(fph,"\nGridDC::open_p - new file branch\n\n");
	fprintf(fph,"in=%s\n",datafile.c_str());
	fflush(fph);

	userfree_(name);
	  
#ifdef FRUITCAKE
	outfile<<"  open_p: archival data file "<<datafile<<" opened on file ptr "<<fp<<"\n";
#endif

      }
      
      /* old file option - check for grid compatibility, extract
	 datafile name, test for readability. Note that any additional
	 compatibility tests with header must be inserted by hand. */
	
      else {

#ifdef FRUITCAKE
	outfile<<"GridDC::open_p: open existing archival hdr file \n";
	outfile<<"  proto = "<<protohdr<<endl;
#endif
#ifdef FRUITCAKE
	outfile<<"  open_p: archival hdr file "<<this->getFilename()<<" opened on file ptr "<<fph<<"\n";
#endif
	/* compatibility check inescapably requires Grid comparison.
	   NOTE that other attributes (data format, data type, ...) are NOT
	   CHECKED.- only attributes of Grid type. This could be changed 
	   later. 
	*/
	PARARRAY * par  = ps_new();
	if (ps_createfp(par,fph)) {
	  RVLException e;
	  e<<"Error: GridDC::open_p\n";
	  e<<"failed to read file "<<this->getFilename()<<" into PARARRAY\n";
	  throw e;
	}
	grid gtest;
	init_grid(&gtest,protog.dim,protog.gdim);
	par_grid(&gtest,*par,stderr);
	if (compare_grid(protog,gtest)) {
	  RVLException e;
	  e<<"Error: GridDC::open_p\n";
	  e<<"input RSF file "<<this->getFilename()<<"incompatible with\n";
	  e<<"prototype grid defined in "<<protohdr<<"\n";
	  throw e;
	}
	  
	/* data file must exist in this branch - parse it out and
	   test with fopen - fails if not found with correct prototype
	   and r+ access
	*/
	string inkey="in";
	if (!(parse(*par,inkey,datafile))) {
	  RVLException e;
	  e<<"Error: GridDC::open_p\n";
	  e<<"metadata file "<<this->getFilename()<<" failed parse key=in\n";
	  throw e;
	}

	// done with partest
	ps_delete(&par);

	// open data file - remains open for lifetime of object,
	// other accesses to data file return same FILE*

#ifdef FRUITCAKE
	outfile<<"GridDC::open_p: open existing archival file for data\n";
	outfile<<"  proto = "<<protodata<<endl;
#endif
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

#ifdef FRUITCAKE
	outfile<<"  open_p: archival data file "<<datafile<<" opened on file ptr "<<fp<<"\n";
#endif
  
      }
    }
    // in effect set file pointer to start
    panelindex=0;
    rd=false;
  }

  ContentPackage<ireal,RARR> & GridDC::get(bool & more) {

    try {

#ifdef FRUITCAKE 
      outfile<<"GridDC::get (mutable)\n";
      outfile<<"  hdr = "<<this->getFilename()<<" more = "<<more<<" rd = "<<rd<<" panelindex = "<<panelindex<<endl;
#endif

      // if more = false on call you are a bad, bad person
      /*
	if (!more) {
	RVLException e;
	e<<"Error: GridDC::get (mutable)\n";
	e<<"  called with more = false - this should not happen!\n";
	throw e;
	}
      */
      // if panelindex too big, you are a bad, bad person
      if (panelindex < 0 || panelindex >= panelnum) {
	RVLException e;
	e<<"Error: GridDC::get (mutable)\n";
	e<<"reading file = "<<this->getFilename()<<"\n";
	e<<"panelindex = "<<panelindex<<" out of range [0, "
	 <<panelnum<<"\n";
	throw e;
      }

      // read data into buffer
      size_t ngrid;
      ra_a_datasize(&(buf.getMetadata()),&ngrid);
      off_t cur_pos = panelindex * ngrid * sizeof(float);
      if (fseeko(fp,cur_pos,SEEK_SET)) {
	RVLException e;
	e<<"Error: GridDC::get (mutable)\n";
	e<<"  from fseeko, metafile = "<<this->getFilename()<<"panelindex="<<panelindex<<"\n";
	throw e;
      }
      if (ngrid != fread(buf.getData(),sizeof(float),ngrid,fp)) {
	RVLException e;
	e<<"Error: GridDC::get (mutable)\n";
	e<<"  from fread, metafile = "<<this->getFilename()<<" panelindex="<<panelindex<<"\n";
	e<<"  panelnum="<<panelnum<<"\n";
	throw e;
      }
      // set rd flag for successful read
      rd=true;
	
      // determine whether this is last trace
      if (panelindex < panelnum) {
	more=true;
	panelindex++;
      }
      if (panelindex >= panelnum) more=false;

#ifdef FRUITCAKE
      outfile<<"GridDC::get (mutable) - finish\n";
      outfile<<"  hdr = "<<this->getFilename()<<" more = "<<more<<" rd = "<<rd<<" panelindex = "<<panelindex<<endl;
#endif
      return buf;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridDC::get (mutable)\n";
      throw e;
    }
  }

  ContentPackage<ireal, RARR> const & GridDC::get(bool & more) const {
    try {
#ifdef FRUITCAKE 
      outfile<<"GridDC::get (const)\n";
      outfile<<"  hdr = "<<this->getFilename()<<" more = "<<more<<" rd = "<<rd<<" panelindex = "<<panelindex<<endl;
#endif

      // if more = false on call you are a bad, bad person
      /*
	if (!more) {
	RVLException e;
	e<<"Error: GridDC::get (const)\n";
	e<<"  called with more = false - this should not happen!\n";
	throw e;
	}
      */
      // if panelindex too big, you are a bad, bad person
      if (panelindex < 0 || panelindex >= panelnum) {
	RVLException e;
	e<<"Error: GridDC::get (const) \n";
	e<<"reading file = "<<this->getFilename()<<"\n";
	e<<"panelindex = "<<panelindex<<" out of range [0, "
	 <<panelnum<<"\n";
	throw e;
      }

      // read data into buffer
      size_t ngrid;
      ra_a_datasize(&(buf.getMetadata()),&ngrid);
      off_t cur_pos = panelindex * ngrid * sizeof(float);
      if (fseeko(fp,cur_pos,SEEK_SET)) {
	RVLException e;
	e<<"Error: GridDC::get (const)\n";
	e<<"  from fseeko, metafile = "<<this->getFilename()<<"panelindex="<<panelindex<<"\n";
	throw e;
      }
      if (ngrid != fread(buf.getData(),sizeof(float),ngrid,fp)) {
	RVLException e;
	e<<"Error: GridDC::get (const)\n";
	e<<"  from fread, metafile = "<<this->getFilename()<<"panelindex="<<panelindex<<"\n";
	throw e;
      }
      // set rd flag for successful read
      rd=true;
	
      // determine whether this is last trace
      // determine whether this is last trace
      if (panelindex < panelnum) {
	more=true;
	panelindex++;
      }
      if (panelindex >= panelnum) more=false;
      

#ifdef FRUITCAKE
      outfile<<"GridDC::get (const) - finish\n";
      outfile<<"  hdr = "<<this->getFilename()<<" more = "<<more<<" rd = "<<rd<<" panelindex = "<<panelindex<<endl;
#endif
      return buf;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridDC::get (const)\n";
      throw e;
    }
  }

  void GridDC::put() const {
    try {
#ifdef FRUITCAKE 
      outfile<<"GridDC::put\n";
      outfile<<"  hdr = "<<this->getFilename()<<" rd = "<<rd<<" panelindex = "<<panelindex<<endl;
#endif

      // note that put does not require any update of buffer origin,
      // since buffer is left with same data - does not change column
      // stored in buffer. Only a read can change the data in the 
      // buffer.

      // if read flag is set, overwrite most recently read segment
      if (rd) panelindex--;

      // if panelindex too big, you are a bad, bad person
      if (panelindex < 0 || panelindex >= panelnum) {
	RVLException e;
	e<<"Error: GridDC::put\n";
	e<<"writing file = "<<this->getFilename()<<"\n";
	e<<"panelindex = "<<panelindex<<" out of range [0, "
	 <<panelnum<<"\n";
	throw e;
      }

#ifdef FRUITCAKE
      outfile<<"GridDC::put - call fwrite, panelindex = "<<panelindex<<endl;
      outfile<<"  hdr = "<<this->getFilename()<<" data = "<<datafile<<endl;
      
#endif
      size_t ngrid;
      ra_a_datasize(&(buf.getMetadata()),&ngrid);
      off_t cur_pos = panelindex * ngrid * sizeof(float);
      if (fseeko(fp,cur_pos,SEEK_SET)) {
	RVLException e;
	e<<"Error: GridDC::put\n";
	e<<"  from fseeko, metafile = "<<this->getFilename()<<"panelindex="<<panelindex<<"\n";
	throw e;
      }
      if (ngrid != fwrite(buf.getData(),sizeof(float),ngrid,fp)) {
	RVLException e;
	e<<"Error: GridDC::put\n";
	e<<"  from fwrite, metafile = "<<this->getFilename()<<"panelindex="<<panelindex<<"\n";
	throw e;
      }
      
      fflush(fp);
#ifdef FRUITCAKE
      iwave_fprintall(stderr);
      outfile<<"  write "<<ngrid<<" floats at offset "<<cur_pos<<" FILE* = "<<fp<<endl;
      outfile<<"  from meta: getSize = "<<getDataSize<RARR>(buf.getMetadata())<<endl;
      //      for (int i=0;i<ngrid;i++) outfile<<"  buf["<<i<<"] = "<<buf.getData()[i]<<endl;
      //      for (int i=0;i<ngrid;i++) outfile<<"  s0["<<i<<"] = "<<buf.getMetadata()._s0[i]<<endl;
#endif
      // unset read flag
      rd=false;

      // advance panelindex - ready to write next panel
      panelindex++;
#ifdef FRUITCAKE
      outfile<<"GridDC::put - finish\n";
      outfile<<"  hdr = "<<this->getFilename()<<" rd = "<<rd<<" panelindex = "<<panelindex<<endl;
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridDC::put\n";
      throw e;
    }
  }    
    
  void GridDC::reset() const { 
    try {
#ifdef FRUITCAKE
      outfile<<"GridDC::reset"<<endl;
      outfile<<"  hdr = "<<this->getFilename()<<endl;
#endif
      if (!fp) open_p();
      panelindex=0;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridDC::reset\n";
      throw e;
    }
  }

  GridDC::GridDC(string const & _protohdr,
		 string const & _protodata,
		 grid const & _protog,
		 string const & _data_format,
		 string data_type,
		 bool incore,
		 ostream & _outfile)
    : OCDC<ireal, RARR>(_outfile), 
      protohdr(_protohdr),
      protodata(_protodata),
      protog(_protog),
      data_format(_data_format),
      outfile(_outfile),
      fph(NULL),
      fp(NULL),
      panelindex(0),
      panelnum(1),
      rd(false)
       {

    // NB: no sanity tests here, compatibility presumed
    // inherited from DCF.
#ifdef FRUITCAKE
    outfile<<"GridDC constructor:\n";
    outfile<<"metadata = "<<protohdr<<endl;
    outfile<<"data     = "<<protodata<<endl;
    outfile<<"format   = "<<data_format<<endl;
    outfile<<"datatype = "<<data_type<<endl;
#endif

    // pass data_type token
    if (data_type.size()>0) this->setTag(data_type);

    // extract diml arrays from grid
    IPNT ng;
    IPNT sg;
    get_n(ng, protog);
    get_gs(sg, protog);

    // incore branch: set dim = gdim in copy of grid, 
    // compute panelnum
    grid g;
    copy_grid(&g,&protog);
    if (incore) g.dim=g.gdim;
    panelnum = get_panelnum_grid(g);

    // how this works: ra_declare assigns the IPNT data members of
    // the RARR workspace tmp. CP::initialize uses operator new and
    // tmp to dynamically allocate its internal RARR, which does a
    // deep copy of the IPNTs and copies the pointer to ireal data,
    // nothing more. Then the newData template specialization calls
    // ra_allocate to finish the construction of the internal RARR,
    // by allocating memory to the ireal pointer. The RARR workspace
    // tmp never has its ireal pointer allocated, and goes out of
    // scope at the end of the constructor, so this construction is
    // memory-efficient. In effect, as far as metadata role goes,
    // the ra_declare step fully initializes the RARR tmp. The data
    // buffer (internal RARR _s0) gets allocated in the initialize
    // step - in contrast to segy case, where data buffer has fixed
    // length and is deep-copied.
    RARR tmp;
    ra_declare_s(&tmp,g.dim,sg,ng);
    buf.initialize(tmp);

#ifdef DUMP_RARR
    ra_dump(&(buf.getMetadata()),stderr);
    fprint_grid(stderr,g);
#endif
  }
    
  /** destructor */
  GridDC::~GridDC() {
    iwave_fclose(fp);
    iwave_fclose(fph);
  }

  bool GridDC::isInCore() const { 
    if (panelnum==1) return true;
    return false;
  }

  string GridDC::getProtohdr() const { return protohdr; }

  ostream & GridDC::write(ostream & str) const {
    str<<"Grid Data Container object \n";
    str<<"based on prototype metadata "<<protohdr<<"\n";
    str<<"metadata file = "<<this->getFilename()<<" pointer = "<<fph<<endl;
    str<<"data file = "<<this->datafile<<" pointer = "<<fp<<endl;
#ifdef VERBOSE
    bool more=true;
    this->reset();
    while (more) {
      ContentPackage<ireal, grid> const & ref = this->get(more);
      ref.write(str);
    }
#endif
    return str;
  }
  
  GridDC * GridDCF::buildGridDC() const {
    return new GridDC(protohdr,
		      protodata,
		      protog,
		      data_format,
		      data_type,
		      incore,
		      outfile);
  }

  PackageContainer<ireal, RARR> * GridDCF::buildPC() const {
    return buildGridDC();
  }

  GridDCF::GridDCF(string _protohdr, bool _incore, ostream & _outfile)
    : protohdr(_protohdr),
      protodata(""),
      data_format("native_float"),
      data_type(""),
      outfile(_outfile), 
      scfac(REAL_ONE),
      incore(_incore),
      vol(REAL_ONE) {

    // sanity check
    char * buf = new char[MAX_NAMELEN];
    if (protohdr.size()>MAX_NAMELEN-10) {
      RVLException e;
      e<<"Error: GridDCF constructor\n";
      e<<"proto header filename = "<<protohdr<<" too large - > "<<MAX_NAMELEN-10<<"\n";
      throw e;
    }

    // assure presence in file mgr
    //    cerr<<"protohdr = "<<protohdr<<endl;
    strcpy(buf,protohdr.c_str());
    FILE * fp=NULL;
    if (!(fp=iwave_fopen(&buf,"r+",NULL,stderr))) {
      RVLException e;
      e<<"Error: GridDCF constructor\n";
      e<<"failed to open proto header file "<<protohdr<<" in file mgr\n";
      e<<"rank = "<<retrieveGlobalRank()<<"\n";
      throw e;
    }      

    // extract parameter array
    PARARRAY * par = ps_new();	
    if (ps_createfp(par,fp)) { 
      RVLException e;
      e<<"Error: GridDCF constructor\n";
      e<<"failed to generate parameter array from proto metadata file "<<protohdr<<"\n";
      throw e;
    }
      
    iwave_fclose(fp);
      
    // store grid, cell volume - should be global
    par_grid(&protog,*par,stderr);
    vol=get_global_cellvol_grid(protog);

    
    //    cerr<<"in GridDCF:\n";
    //    cerr<<"par\n";
    //    ps_printall(*par,stderr);
    //    cerr<<"grid\n";
    //    fprint_grid(stderr,protog);
    

    string inkey;

    //    cerr<<"extract prototype data file name\n";
    inkey="in";
    if (!(parse(*par,inkey,protodata))) {
      RVLException e;
      e<<"Error: GridDCF constructor\n";
      e<<"proto metadata file "<<protohdr<<" failed to parse key=in\n";
      e<<"rank = "<<retrieveGlobalRank()<<"\n";
      throw e;
    }      

    //   cerr<<"assure presence in file mgr\n";
    strcpy(buf,protodata.c_str());
    FILE * fpd=NULL;
    if (!(fpd=iwave_fopen(&buf,"r+",NULL,stderr))) {
      RVLException e;
      e<<"Error: GridDCF constructor\n";
      e<<"failed to open proto data file "<<protodata<<" in file mgr\n";
      e<<"rank = "<<retrieveGlobalRank()<<"\n";
      throw e;
    }      
    iwave_fclose(fpd);

    //   cerr<<"extract other metadata - defaults OK\n";
    inkey="data_format";
    parse(*par,inkey,data_format);
    inkey="scale";
    int scale = 0;
    parse(*par,inkey,scale);
    inkey="data_type";
    parse(*par,inkey,data_type);

    //   cerr<<"scale factor\n";
    for (int i=0;i<abs(scale);i++) {
      if (scale>0) scfac *= 10;
      if (scale<0) scfac *= 0.1;
    }

    //   cerr<<"clean up\n";
    delete [] buf;
    ps_delete(&par);

  }

  GridDCF::GridDCF(GridDCF const & f) 
    : protohdr(f.protohdr),
      protodata(f.protodata),
      outfile(f.outfile),
      protog(f.protog),
      scfac(f.scfac),
      incore(f.incore),
      vol(f.vol)
       {}
    
  GridDCF::~GridDCF() {}
    
  PackageContainerFactory<ireal, RARR> * GridDCF::clone() const {
    return new GridDCF(*this);
  }

  grid const & GridDCF::getGrid() const { 
    return protog;
  }
    
  ireal GridDCF::getCellVol() const { return vol; }
  ireal GridDCF::getScaleFactor() const { return scfac; }
    
  string GridDCF::getFilename() const { return protohdr; }
    
  bool GridDCF::compare( DataContainerFactory const & dcf) const {
    GridDCF const * f = NULL;
    f = dynamic_cast< GridDCF const *>(&dcf);
    if ((f) && (this->protohdr == f->protohdr) &&
	(this->incore == f->incore)) return true;
    return false;
  }

  // compatible if either DFC incore flag is unset (whether data is incore or not)
  // or DFC incore flag is set and DC data is incore for any reason (i.e. internal
  // grid has dim=gdim).
  bool GridDCF::isCompatible( DataContainer const & dc ) const {
    GridDC const * gc = NULL;
    gc = dynamic_cast<GridDC const *>(&dc);
    if (gc && (gc->getProtohdr()==this->protohdr)) 
      if ((this->incore && gc->isInCore()) || !(this->incore)) return true;
    return false;
  }

  bool GridDCF::isIncore() const { return incore; }
    
  ostream & GridDCF::write(ostream & str) const {
    str<<"GridDCF: factory class for GridDC\n";
    str<<"Prototype RSF filename = "<<protohdr<<"\n";
    str<<"Incore vector arithmetic = "<<incore<<"\n";
    return str;
  }

  STRING_PAIR GridSpace::getPTblEntry(grid const & g,
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
    fprint_grid(fph,g);
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

    int n1 = get_datasize_grid(g);
    int nt = get_global_datasize_grid(g);
    char * x = (char *)usermalloc_(n1*sizeof(char)*sizeof(ireal));
    memset(x,0,n1*sizeof(ireal));
    while (nt>0) {
      fwrite(x,sizeof(char),n1*sizeof(ireal),fp);
      nt-=n1;
    }
    userfree_(name);
    userfree_(x);
    if (ftello(fp) != (off_t) (get_global_datasize_grid(g)*sizeof(ireal))) {
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

  GridSpace::GridSpace(string hdr, 
		       string dtype, 
		       bool incore,
		       ostream & outfile)
    : StdSpace<ireal,ireal>(),
      ConstContainer<STRING_PAIR>(STRING_PAIR(dtype,hdr)),
      f(hdr,incore,outfile), 
      l(f.getCellVol()*f.getScaleFactor()*f.getScaleFactor()) {
  }

  GridSpace::GridSpace(grid const & g,
		       std::string tag,
		       std::string thdr,
		       std::string fmt,
		       bool incore,
		       ostream & outfile)
    : StdSpace<ireal,ireal>(),
      ConstContainer<STRING_PAIR>(getPTblEntry(g,tag,thdr,fmt)),
      f(thdr,incore,outfile), 
      l(f.getCellVol()*f.getScaleFactor()*f.getScaleFactor()) {}

  GridSpace::GridSpace(GridSpace const & sp) 
    : StdSpace<ireal,ireal>(), 
      ConstContainer<STRING_PAIR>(sp.get()),
      //	f(sp.get().val,sp.outfile), 
      f(sp.f),
      l(f.getCellVol()*f.getScaleFactor()*f.getScaleFactor()) {}

  GridSpace::~GridSpace() {}

  LinearAlgebraPackage<ireal> const & GridSpace::getLAP() const { return l; }

  DataContainerFactory const & GridSpace::getDCF() const { return f; }

  grid const & GridSpace::getGrid() const { return f.getGrid(); }

  bool GridSpace::isIncore() const { return f.isIncore(); }

  ostream & GridSpace::write(ostream & str) const {
    str<<"GridSpace: StdSpace based on RSF grid data\n";
    str<<" - data type   = "<<this->get().key<<"\n";
    str<<" - header file = "<<this->get().val<<"\n";
    return str;
  }

}	      
      

