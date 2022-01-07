#include "segypp.hh"
#include "segyfun.hh"

namespace TSOpt {

  using RVL::RVLException;
  using RVL::ContentPackage;
  using RVL::PackageContainer;
  using RVL::PackageContainerFactory;

  SEGYDC::SEGYDC(string _hdr, ostream & _outfile)
    : OCDC<float,segy>(_outfile), 
      hdr(_hdr),
      fp(NULL),
      rd(false),
      off_cur(0L),
      outfile(_outfile),
      off_eof(0L) {
    buf=new ContentPackage<float,segy>;
    segy tr;
    buf->initialize(tr);
    // 18.11.09: for the nonce, the universal tag is "datafile"
    string tagstr="datafile";
    this->setTag(tagstr);
    // 27.03.10: iwave_fopen based implementation presumes that 
    // datafile is initialized

    // open header file - should ALWAYS be opened with these permissions
    // note that if it's already open, this call simply returns the 
    // associated pointer
    FILE * fph = iwave_const_fopen(hdr.c_str(),"r",NULL,stderr);
    if (!fph) {
      RVLException e;
      e<<"Error: SEGYCFC constructor\n";
      e<<"failed to access proto file ptr, filename="<<hdr<<"\n";
      throw e;
    }

    if (fseeko(fph,0L,SEEK_SET)) {
      RVLException e;
      e<<"Error: SEGYDC constructor\n";
      e<<"existing file branch\n";
      e<<"seek to begin of prototype file "<<hdr<<" failed\n";
      throw e;
    }

    // set nb, ns and dt
    if (!(nb=fgettr(fph,&tr))) {
      RVLException e;
      e<<"Error: SEGYDC constructor\n";
      e<<"existing file branch\n";
      e<<"failed to read first trace from prototype file "<<hdr<<"\n";
      throw e;
    }
    ns = tr.ns; 
    dt = tr.dt;

    // set off_eof
    if (fseeko(fph,0L,SEEK_END)) {
      RVLException e;
      e<<"Error: SEGYDC constructor\n";
      e<<"seek to end of proto file "<<hdr<<" failed\n";
      throw e;
    }    
    off_eof=ftello(fph);
    if (fseeko(fph,0L,SEEK_SET)) {
      RVLException e;
      e<<"Error: SEGYDC constructor\n";
      e<<"seek to begin of proto file "<<hdr<<" failed\n";
      throw e;
    }    

    // unset in-use flag to permit reuse of proto file
    iwave_fclose(fph);

#ifdef FRUITCAKE
    outfile <<"SEGYDC constructor: hdr="<<hdr<<" nb="<<nb<<" ns="<<ns<<" dt="<<dt<<" off_eof="<<off_eof<<"\n";
#endif
    ///////////////
    //    cerr <<"SEGYDC constructor: hdr="<<hdr<<" nb="<<nb<<" ns="<<ns<<" dt="<<dt<<" off_eof="<<off_eof<<"\n";
    ///////////////
  }

  SEGYDC::~SEGYDC() {
    if (fp) iwave_fclose(fp); 
    fp=NULL;
    if (fph) iwave_fclose(fph);
    if (buf) delete buf;
    ///////////////
    //    cerr <<"SEGYDC::~SEGYDC: close file "<<this->getFilename()<<"\n";
    //    iwave_fprintall(stderr);
    ///////////////
  }
    
  void SEGYDC::open_p() const {

    // no try block as no function calls which might throw exception

    if (fp) {
      RVLException e;
      e<<"Error: SEGYDC::open_p\n";
      e<<"error state on call to open_p: file pointer set\n";
      throw e;
    }

#ifdef FRUITCAKE
    outfile<<"SEGYDC::open_p hdr="<<hdr<<" file="<<this->getFilename()<<"\n";
#endif    
      
    /* set up inputs to iwave_fopen
     * name - equal to filename c string if nonzero len, else NULL
     * proto - header file name
     * mode - w+ if temp file (name not set on call), else r+ if possible,
     *   else w+.
     */

    bool old = false;

    if (this->getFilename() == "") {
      char * name = NULL;

      fp = iwave_fopen(&name,"w+",hdr.c_str(),stderr);
      if (!fp) {
	RVLException e;
	e<<"Error: SEGYDC::open_p\n";
	e<<"temp file case: failed to open file\n";
	throw e;
      }
      if (!name) {
	RVLException e;
	e<<"Error: SEGYDC::open_p\n";
	e<<"temp file case: failed to set name\n";
	throw e;
      }

      //////////////
      //      cerr<<"SEGYDC->open_p->iwave_fopen name="<<name<<" hdr="<<hdr<<endl;
      //      iwave_fprintall(stderr);
      //////////////

      string sname = name;
      this->setFilename(sname);
      if (name) userfree_(name);

    }
    else {
      if (!(fp = iwave_const_fopen(this->getFilename().c_str(),"r+",hdr.c_str(),stderr))) {
	fp = iwave_const_fopen(this->getFilename().c_str(),"w+",hdr.c_str(),stderr);
	if (!fp) {
	  RVLException e;
	  e<<"Error: SEGYDC::open_p";
	  e<<"perm file case: failed to open file\n";
	  throw e;
	}
	//////////////
	//	cerr<<"SEGYDC->open_p->iwave_const_fopen name = "<<this->getFilename()<<" hdr = "<<hdr<<endl;
	//	iwave_fprintall(stderr);
	//////////////
      }
      else {
	old=true;
	//////////////
	//	cerr<<"SEGYDC->open_p old file name = "<<this->getFilename()<<" hdr = "<<hdr<<endl;
	//////////////
      }
    }

    /* now check that ns and dt are same as prototype, for existing
       files - then H-space structures are same. iwave_fopen has
       already checked that file lengths are same, so this check
       guarantees that traces are same length as well.
    */
    if (old) {

      if (fseeko(fp,0L,SEEK_SET)) {
	RVLException e;
	e<<"Error: SEGYDC::open_p\n";
	e<<"existing file branch\n";
	e<<"seek to begin of data file "<<this->getFilename()<<" failed\n";
	throw e;
      }
      
      int nbtmp = fgettr(fp,&tr);
      if (nb!=nbtmp) {
	RVLException e;
	e<<"Error: SEGYDC::open_p\n";
	e<<"existing file branch, file = "<<this->getFilename()<<" proto="<<hdr<<"\n";
	e<<"first trace from prototype file has incorrect number of bytes"<<hdr<<"\n";
	e<<"proto = "<<(int)nb<<" this file = "<<(int)nbtmp<<"\n";
	throw e;
      }
      if ((ns != tr.ns) || (dt != tr.dt)) {
	RVLException e;
	e<<"Error: SEGYDC::open_p\n";
	e<<"existing file branch\n";
	e<<"data file "<<this->getFilename()<<" differs in structure from prototype\n";
	e<<"file "<<hdr<<": either number of bytes per trace, number of samples per\n";
	e<<"trace, or time sample rate are different\n";
	throw e;
      }

      // reposition file at begin for subsequent use
      if (fseeko(fp,0L,SEEK_SET)) {
	RVLException e;
	e<<"Error: SEGYDC::open_p\n";
	e<<"existing file branch\n";
	e<<"seek to begin of data file "<<this->getFilename()<<" failed\n";
	throw e;
      }
      // end existing file check
    }

    // clean up char buffers

  }

  segytrace & SEGYDC::get(bool & more) {
    try {
#ifdef FRUITCAKE
      outfile<<"SEGYDC::get (mutable) -> open_p file="<<this->getFilename()<<" fileptr="<<fp<<"\n";
#endif
      //      if (!fp) open_p();

      // record current position
      off_cur=ftello(fp);
#ifdef FRUITCAKE
      outfile<<"SEGYDC::get_p->fgettr off_cur = "<<off_cur<<"\n";
#endif
      if (!(fgettr(fp,&(buf->getMetadata())))) {
	outfile<<"SEGYDC::get (mutable) - failure\n";
	RVLException e;
	e<<"Error: SEGYDC::get\n";
	e<<"failed to read trace on file "<<this->getFilename()<<"\n";
	throw e;
      }
      // set rd flag for successful read
      rd=true;
      // determine whether this is last trace
      more=false;
#ifdef FRUITCAKE
      outfile<<"SEGYDC::get_p->ftello\n";
#endif
      if (ftello(fp)<off_eof) more=true;
#ifdef FRUITCAKE
      outfile<<"SEGYDC::get_p - returning with more="<<more<<endl;
      outfile<<"SEGYDC::get_p - ftello = "<<ftello(fp)<<" off_eof = "<<off_eof<<endl;
#endif
      /*
      int i;
      fprintf(stderr,"at end of SEGYDC::get (mutable)\n");
      scanf("%d",&i);
      */
      return *buf;
    }
    catch (RVLException & e) {
      e<<"\ncalled from SEGYDC::get (nonconst)\n";
      throw e;
    }
  }

  segytrace const & SEGYDC::get(bool & more) const {
    try {
#ifdef FRUITCAKE
      outfile<<"SEGYDC::get (const) -> open_p file="<<this->getFilename()<<" fileptr="<<fp<<"\n";
#endif
      //      if (!fp) open_p();

      // record current position
      off_cur=ftello(fp);
#ifdef FRUITCAKE
      outfile<<"SEGYDC::get_p->fgettr off_cur = "<<off_cur<<"\n";
#endif
      if (!(fgettr(fp,&(buf->getMetadata())))) {
	outfile<<"SEGYDC::get (const) - failure\n";
	RVLException e;
	e<<"Error: SEGYDC::get\n";
	e<<"failed to read trace on file "<<this->getFilename()<<"\n";
	throw e;
      }
      // set rd flag for successful read
      rd=true;
      // determine whether this is last trace
      more=false;
#ifdef FRUITCAKE
      outfile<<"SEGYDC::get_p->ftello\n";
#endif
      if (ftello(fp)<off_eof) more=true;
#ifdef FRUITCAKE
      outfile<<"SEGYDC::get_p - returning with more="<<more<<endl;
      outfile<<"SEGYDC::get_p - ftello = "<<ftello(fp)<<" off_eof = "<<off_eof<<endl;
#endif
      /*
      int i;
      fprintf(stderr,"at end of SEGYDC::get (const)\n");
      scanf("%d",&i);
      */
      return *buf;
    }
    catch (RVLException & e) {
      e<<"\ncalled from SEGYDC::get (const)\n";
      throw e;
    }
  }

  void SEGYDC::put() const {
    try {
#ifdef FRUITCAKE
      outfile<<"SEGYDC::put -> open_p file="<<this->getFilename()<<" fileptr="<<fp<<"\n";
#endif
      //      if (!fp) open_p(); 

      // if read flag is set, seek to recorded trace begin 
      if (rd) {
	if (fseeko(fp,off_cur,SEEK_SET)) {
	  RVLException e;
	  e<<"Error: SEGYDC::put_p\n";
	  e<<"seek failed to "<<off_cur<<"on file "<<this->getFilename()<<"\n";
	  throw e;
	}	
      }
      off_cur=ftello(fp);
#ifdef FRUITCAKE
      outfile<<"SEGYDC::put_p -> fputtr at "<<off_cur<<"\n";
#endif
      fputtr(fp,&(buf->getMetadata()));
      fflush(fp);
      // unset read flag
      rd=false;

      /*
      int i;
      fprintf(stderr,"at end of SEGYDC::put\n");
      scanf("%d",&i);
      */
#ifdef FRUITCAKE
      outfile<<"SEGYDC::put_p - return\n";
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from SEGYDC::put\n";
      throw e;
    }
  }    

  void SEGYDC::reset() const { 
    try {
#ifdef FRUITCAKE
      outfile<<"SEGYDC::reset -> open_p file="<<this->getFilename()<<" fileptr="<<fp<<"\n";
#endif
      if (!fp) open_p();
#ifdef FRUITCAKE
      outfile<<"SEGYDC::reset -> open_p file="<<this->getFilename()<<" fileptr="<<fp<<"\n";
      off_t i = ftello(fp);
      outfile<<"SEGYDC::reset at offset"<<i<<"\n";
#endif
      fseeko(fp,0L,SEEK_SET);
#ifdef FRUITCAKE
      i=ftello(fp);
      outfile<<"SEGYDC::reset after fseeko at offset"<<i<<"\n";
#endif
      /*
      fprintf(stderr,"at end of SEGYDC::reset\n");
      scanf("%d",&i);
      */
    }
    catch (RVLException & e) {
      e<<"\ncalled from SEGYDC::reset\n";
      throw e;
    }
  }

  ostream & SEGYDC::write(ostream & str) const {
    str<<"SEGYDC: Data Container class for SEGY trace data\n";
    str<<"Prototype SU filename = "<<hdr<<"\n";
    str<<"Data file             = "<<this->getFilename()<<"\n";
    return str;
  }    

  SEGYDCF::SEGYDCF(string _hdr, ostream & _outfile)
    : hdr(_hdr), outfile(_outfile), iflag(false) {
    //    cerr<<"**** built SEGYDCF hdr = "<<hdr<<endl; 
  }

  void SEGYDCF::init() const {
    
    if (retrieveGlobalRank() != 0) {
      RVLException e;
      e<<"Error: SEGYDCF::init() called on non-root process\n";
      e<<"no dependency on second-level initialization permitted\n";
      e<<"except on root\n";
      throw e;
    }
    
    if (iflag) return;

    /* ensure that proto file is opened */
    FILE * fp = NULL;
    //    cerr<<"SEGYDCF: opening proto file "<<hdr<<endl;
    if (!(fp=iwave_const_fopen(hdr.c_str(),"r",NULL,stderr))) {
      RVLException e;
      e<<"Error: SEGYDCF constructor\n";
      e<<"failed to open proto file "<<hdr<<"\n";
      throw e;
    }
    
    /* read first trace, set number of bytes for all traces 
       (assumed const!!!) */
    segy m;
    nb=fgettr(fp,&m);
    if (!nb) {
      RVLException e;
      e<<"Error: SEGYDCF::getDt\n";
      e<<"read zero bytes on first trace of proto file "<<hdr<<"\n";
      throw e;
    }
    
    // extract float sample interval ms
    dt = 0.001*((float)(m.dt));
    // cerr<<"m.dt="<<m.dt<<" dt="<<dt<<endl;
    nt = (int)m.ns;

    // unset use flag
    iwave_fclose(fp);

    iflag=true;
  }

  bool SEGYDCF::compare( DataContainerFactory const & dcf) const {
    SEGYDCF const * f = NULL;
    f = dynamic_cast<SEGYDCF const *>(&dcf);
    if ((f) && (this->hdr == f->hdr)) return true;
    return false;
  }

  bool SEGYDCF::isCompatible( DataContainer const & dc ) const {
    SEGYDC const * gc = NULL;
    gc = dynamic_cast<SEGYDC const *>(&dc);
    if (gc && gc->getHdr()==hdr) return true; 
    return false;
  }

  ostream & SEGYDCF::write(ostream & str) const {
    str<<"SEGYDCF: factory class for SEGYDC\n";
    str<<"Prototype SU filename = "<<hdr<<"\n";
    return str;
  }

  ostream & SEGYSpace::write(ostream & str) const {
    str<<"SEGYSpace: StdSpace based on SEGYDC trace data container class\n";
    str<<" - header file = "<<this->get().val<<"\n";
    return str;
  }
}
