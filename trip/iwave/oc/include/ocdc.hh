#ifndef __RVL_OCDC
#define __RVL_OCDC

//#define FRUITCAKE

#include "usempi.h"
#include "contentpackage.hh"
#include "mpiserialfo.hh"
#include "parser.h"

namespace RVL {

  /** FO wrapper around a const object. A DC may eval this FO and do
      anything it damn well pleases with the object data. For example:
      use ConstContainer<string> to assign filenames to out-of-core
      DataContainer objects. 
  */
  template<typename T>
  class ConstContainer: public FunctionObject {
  private:
    T t;
    ConstContainer();
  public:
    ConstContainer(T _t): t(_t) {}
    ConstContainer(ConstContainer const & af): t(af.t) {}
    virtual ~ConstContainer() {}
    T get() const { return t; }
    virtual string getName() const { string ret="ConstContainer"; return ret; }
  };

  /** Const FO (FOR) wrapper around a mutable object. A DC may eval
      this FOR and alter the data in some constructive way. For
      example, use Container<PARARRAY> to transfer key=value pairs to
      parameter tables (PARARRAY from iwave package, or any equivalent
      type). The return value is essentially a dummy. This may look
      like an abuse of the FOR concept, and in fact if we introduced a
      const FO type we could use it instead of FOR as a
      parent. Clearly the data type T should not be required to behave
      in any sense as small data, so T is not a RetType and the
      set/getResult interfaces are not appropriate. Hence the
      provision of an alternate get of a non-const reference.
  */
  template<typename T>
#ifdef IWAVE_USE_MPI
  class Container: public FunctionObjectConstEval,
		   public MPISynchRoot 
#else
  class Container: public FunctionObjectConstEval 
#endif
  {

  private:
    
    T & t;
    
    Container();
    Container(Container const &);
    
  public:
    /** main constructor: record T reference */
    Container(T & _t): t(_t) {}
    virtual ~Container() {}
    
    virtual string getName() const { string ret="Container"; return ret; }
    
    /** mutable access to T reference */
    T & get() { return t; }
    /** const access to T reference */
    T const & get() const { return t; }
    
#ifdef IWAVE_USE_MPI
    /** MPISynchRoot functions are no-ops. - T member data is not
	reduction target, has no a-priori defined reference value, and
	does not need to be broadcast. Any communication of the T data
	is the responsibility of the evaluator DC. */
    void set() {}
    virtual void synch() {}
#endif
  };
    
  /** Two special subclasses of ConstContainer have special meaning
      for OCDCs. The first sets the filename, and can only be invoked 
      before the data file is open.
  */
  class AssignFilename: public ConstContainer<string> {
  private:
    AssignFilename();
  public:
    AssignFilename(string s): ConstContainer<string>(s) {}
    AssignFilename(AssignFilename const & af): ConstContainer<string>(af) {}
    ~AssignFilename() {}
    string getName() const { string ret="AssignFilename"; return ret; }
  };

  /** The second assigns a tag */ 
  class AssignTag: public ConstContainer<string> {
  private:
    AssignTag();
  public:
    AssignTag(string s): ConstContainer<string>(s) {}
    AssignTag(AssignTag const & af): ConstContainer<string>(af) {}
    ~AssignTag() {}
    string getName() const { string ret="AssignTag"; return ret; }
  };

  /** FOR whose eval transfers tag=filename pair from OCDC to PARARRAY.
      permits override of tag (by data string), so that OCDC can be used
      in arbitrary role as indicated by key, also addition of pre- and
      post-fix strings to modify internal tag data of OCDC. */
  class AssignParams: public Container<PARARRAY> {
  private:
    string tag;           // optional tag (key) spec
    string pre;           // prefix to add to all tags
    string post;          // postfix to add to all tags

    mutable string skey;  // workspace for key synch
    mutable string sval;  // workspace for value synch

    FILE * stream;        // verbose output

    AssignParams();
    AssignParams(AssignParams const &);
  public:
    AssignParams(PARARRAY & par, FILE * _stream=stderr)
      : Container<PARARRAY>(par), tag(""), pre(""), post(""), 
	skey(""), sval(""), stream(_stream) {}
    AssignParams(PARARRAY & par, string _tag, FILE * _stream=stderr)
      : Container<PARARRAY>(par), tag(_tag), pre(""), post(""), 
	skey(""), sval(""), stream(_stream) {}
    AssignParams(PARARRAY & par, 
		 string _pre, string _post, FILE * _stream=stderr)
      : Container<PARARRAY>(par), tag(""), pre(_pre), post(_post), 
	skey(""), sval(""), stream(_stream) {}
    AssignParams(PARARRAY & par, string _tag, 
		 string _pre, string _post, FILE * _stream=stderr)
      : Container<PARARRAY>(par), tag(_tag), pre(_pre), post(_post), 
	skey(""), sval(""), stream(_stream) {}
    string getTag() const { return tag; }
    string getPrefix() const { return pre; }
    string getPostfix() const {return post; }

    void set_Bcast_key(string _key) { skey=_key; }
    void set_Bcast_val(string _val) { sval=_val; }

    void synch() {
#ifdef IWAVE_USE_MPI
      int nkey = skey.size()+1;
      int nval = sval.size()+1;
      MPI_Bcast(&nkey,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nval,1,MPI_INT,0,MPI_COMM_WORLD);
      if (nkey > 1 && nval > 1) {
	char * tmp = new char[nkey+nval];
	strcpy(tmp,skey.c_str());
	MPI_Bcast(tmp,nkey,MPI_CHAR,0,MPI_COMM_WORLD);
	skey=tmp;
	strcpy(tmp,sval.c_str());
	MPI_Bcast(tmp,nval,MPI_CHAR,0,MPI_COMM_WORLD);
	sval=tmp;
	delete [] tmp;

	// for verbose output
	int rk;
	MPI_Comm_rank(MPI_COMM_WORLD,&rk);
	fprintf(stream,"Assign_Params::synch: rk=%d skey=%s sval=%s\n",
		rk,skey.c_str(),sval.c_str());
	//

	// no need to add pair again on rk 0 - already done
	if (rk) {
	  if (ps_slcstring(this->get(),skey.c_str(),sval.c_str())) {
	    RVLException e;
	    e<<"AssignParams::synch\n";
	    e<<"failed to add key "<<skey<<" and/or value "<<sval<<"\n";
	    e<<"to PARARRAY data member\n";
	    throw e;
	  }
	}
      }
#endif
    }

    string getName() const { string ret="AssignParams"; return ret; }
  };

  /** OCDC - simply PackageContainer with file management
      attributes */
  template<typename T, typename M>
  class OCDC: public PackageContainer<T,M> {

  private:

    ostream & outfile;                  // verbose output stream
    mutable string filename;            // filename 
    mutable string tag;                 // use tag - indicates data type

    OCDC(OCDC<T,M> const &);

  protected:

    /** set filename: should be used to assign filename by 
	subclass method that opens file */
    void setFilename(string s) const { filename=s; }
    /** use to test prior assignment of filename */
    string getFilename() const { return filename; }
    void setTag(string s) const { tag=s; }
    string getTag() const { return tag; }

  public:

    OCDC(ostream & _outfile = cerr) 
      : PackageContainer<T,M>(),
	outfile(_outfile),
	filename(""), 
	tag("")
	 {}

    virtual ~OCDC() {}
    
    /** intercept eval of a string container, interpret as filename or
	tag assignment. Note how this works: evaluation of an FO
	BEFORE file is opened (which always assigns a nonempty
	filename) is an error, UNLESS the FO is either an
	AssignFilename or AssignTag instance. In those two cases,
	evaluation calls setFilename or setTag respectively. Otherwise
	an exception is thrown. Evaluation of an FO AFTER file is
	opened (signalled by a nonempty filename) is an error IF the
	FO is either an AssignFilename or AssignTag instance. Any
	other FO is delegated to PC::eval. Thus a calling unit can
	assign filenames or tags only BEFORE anything else happens
	that would cause a file to be opened.
    */
    void eval(FunctionObject & f, 
	      std::vector<DataContainer const * > &x) {
      try {

	// ConstContainer<string> should be eval'd on all 
	// processors in same way
	ConstContainer<string> * ff = NULL;
	if ((ff=dynamic_cast<ConstContainer<string> *>(&f))) {
	  AssignFilename * af = NULL;
	  AssignTag * at = NULL;
	  if ((af=dynamic_cast<AssignFilename *>(ff))) {
	    if (filename!="") {
	      RVLException e;
	      e<<"Error: OCDC::eval(FO)\n";
	      e<<"AssignFile called after file is opened\n";
	      e<<"current filename = "<<this->filename<<"\n";
	      e<<"attempted to assign "<<ff->get()<<"\n";
	      throw e;
	    }	    
#ifdef FRUITCAKE
	    outfile<<" ocdc::eval of AssignFile, assign fname="
		   <<ff->get()<<"\n";
#endif
	    this->setFilename(ff->get());
	    
#ifdef FRUITCAKE
	    outfile<<" new value of filename = "<<this->filename<<endl;
	    outfile.flush();
#endif
	    // a call to reset is mandatory here - if the filename is
	    // assigned, it still won't be registered with the file
	    // system until a call to reset. subclass-specific databases
	    // may also need to be updated to reflect the file assignment.
	    // any eval of an RVL::FO or FOR will automatically invoke
	    // reset, due to its use in PC::eval, but any other interaction
	    // with data via a non-FO app will not. so it has to happen here.
	  
	    this->reset();
	  }
	  else if ((at=dynamic_cast<AssignTag *>(ff))) {
	    this->setTag(ff->get());
	  }
	  else {
	    RVLException e;
	    e<<"Error: OCDC::eval\n";
	    e<<"eval only defined for special subtypes of ConstContainer<string>\n";
	    throw e;
	  }
	}
	
	else {
#ifdef FRUITCAKE
	  outfile<<" socdc::eval of generic FO = "<<f.getName()<<" on OCDC file = "<<this->filename<<"\n";
	  outfile.flush();
#endif
	  // rely on reset to correctly open, position file
	  PackageContainer<T,M>::eval(f,x);
	}
#ifdef FRUITCAKE
	outfile<<" socdc::eval FO return from "<<f.getName()<<" on OCDC file = "<<this->filename<<"\n";	 
	outfile.flush();
#endif
      }
      catch (RVLException & e) {
	e<<"\ncalled from SOCDC::eval(fo), filename="<<this->filename
	 <<" fo = "<<f.getName()<<"\n";
	throw e;
      }
    }

    // eval of FORs always delegated to PC::eval, EXCEPT for 
    // Container<PARARRAY> instances - in that case, add tag=filename

    virtual void eval(FunctionObjectConstEval & f, 
		      vector<DataContainer const *> & x) const {
      try {
	AssignParams * ff = NULL;
	if ((ff=dynamic_cast<AssignParams *>(&f))) {

	  // reset has to be called here in case the data is temp -
	  // then AssignFilename will not have been called prior to
	  // this point, and reset needed for same reason as explained
	  // in the source for AssignFilename, above.
	  this->reset();

	  //	  cerr<<" socdc::eval of AssignParams on OCDC file = "<<this->filename<<" tag = "<<this->getTag()<<"\n";
#ifdef FRUITCAKE
	  outfile<<" socdc::eval of AssignParams on OCDC file = "<<this->filename<<" tag = "<<this->getTag()<<"\n";
	  outfile.flush();
#endif	
	  if (this->getFilename()=="") {
	    RVLException e;
	    e<<"Error: OCDC::eval (AssignParams)\n";
	    e<<"eval attempted before filename assigned\n";
	    throw e;
	  }	    

	  /* test to see if tag override supplied - if so use, if not don't */
	  string tagstr = ff->getTag();
	  if (tagstr.size()==0) tagstr = this->getTag();
	  if (tagstr.size()==0) {
	    RVLException e;
	    e<<"Error: OCDC eval(FOR) - special AssignParams cast\n";
	    e<<"neither DC nor FOR supplies valid tag\n";
	    throw e;
	  }
	  string key = ff->getPrefix()+tagstr+ff->getPostfix();

	  if (ps_slcstring(ff->get(),key.c_str(),
			   (this->getFilename()).c_str())) {
	    RVLException e;
	    e<<"Error: OCDC eval(FOR) - special AssignParam case\n";
	    e<<"failed to add key "<<key<<" and/or value "
	     <<this->getFilename()<<" to PARARRAY data member of FOR\n";
	    throw e;
	  }	

	  ff->set_Bcast_key(key);
	  ff->set_Bcast_val(this->getFilename());
	}
	else {

#ifdef FRUITCAKE
	  outfile<<" socdc::eval of generic FOR = "<<f.getName()<<" on OCDC file = "<<this->filename<<"\n";
	  outfile.flush();
#endif
	  PackageContainer<T,M>::eval(f,x);
	}
#ifdef FRUITCAKE
	outfile<<" socdc::eval FO return from "<<f.getName()<<" on OCDC file = "<<this->filename<<"\n";	 
	outfile.flush();
#endif
      }
      catch (RVLException & e) {
	e<<"\ncalled from SOCDC::eval(fo), filename="<<this->filename
	 <<" fo = "<<f.getName()<<"\n";
	throw e;
      }
    }
  };

}

#endif
