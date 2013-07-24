/************************************************************************

Copyright Rice University, 2006
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/

#ifndef __RVL_TOOLS_NEW_CONTENT_PACKAGE_
#define __RVL_TOOLS_NEW_CONTENT_PACKAGE_

#include "local.hh"
#include "dejavu.hh"

namespace RVL {
  
  // rewrite of ADP's ContentPackage - now stores only a single 
  // metadata object, which is supposed to be collective
  
  /** Helper function template to extract the size of a data array
      implicit in a MetaType object. The base template will work for
      any intrinsic integer MetaType (size_t, int, long, uint,...).
      More complex MetaTypes will require template specialization. */

  template<class MetaType>
  size_t getDataSize(MetaType const & md) {
    size_t t = md;
    return t;
  }

  /** Helper template for implementation of standard write method. Will
      give sensible if minimal output. */
  template<class MetaType>
  ostream & writeMeta(MetaType const & md, ostream & e) { 
    e<<md<<"\n"; return e; }

  /** Helper template to allocate data memory - permits specialization
      in case eg. metatype object already allocates required memory */
  template<typename DataType,typename MetaType>
  DataType * newData(MetaType & md) {
    return new DataType[getDataSize<MetaType>(md)];
  }

  /** Helper template to delete data memory - paired with allocData.
      specialization permits other arrangements. */
  template<typename DataType,typename MetaType>
  void deleteData(DataType ** d, MetaType ** md) {
    delete [] *d; *d=NULL;
  }

  /** A normal pairing of an array of data along with Metadata
      object. 

      The only requirement on DataType: 

      * must have static size - and object of type DataType must have
      * a definite size known at compile time - and possess a default
      * constructor which allocates a memory block of constant,
      * definite size - sizeof(DataType) must make sense, and actually
      * give the size in bytes of a word of type DataType. Essentially
      * the only possibilities are an intrinsic type or a struct made
      * out of intrinsic types or a struct made out of such structs
      * or...

      The only requirements on MetaType: 

      * a deep copy constructor; 

      * correct behaviour of the getDataSize<MetaType>() helper
      * function, if necessary by template specialization, returning
      * the size of the DataType array specified by a MetaType object.

      These assumptions are implicit in the presumption that all
      information necessary to initialize and define the behaviour of
      an object of this type is present in the MetaType data member
      and the DataType type specification..

      This class follows the "initialization separate from
      instantiation" paradigm. The added flexibility is useful eg. in
      distributed computing - a recipient object can be constructed
      before the data it is to receive has been transmitted. Since the
      MetaType data member determines the size of the DataType array,
      the alternative paradigm would require some kind of auxiliary
      factory object.
  */
  
  template<class DataType, class MetaType >
  class ContentPackage: public LocalDataContainer<DataType> {
    
  private:
    
    mutable MetaType * md;
    mutable DataType * d;
    bool own;

  public:    

    /** copy constructor - deep copy, assumes const word length for
	DataType returned by sizeof. */
    ContentPackage(const ContentPackage<DataType,MetaType> & bin)
      : md(NULL), d(NULL), own(bin.own) {
      if (bin.md) {
	md = new MetaType(*(bin.md));
	if (own) {
	  //	  d = new DataType[getDataSize<MetaType>(*md)];
	  d=newData<DataType,MetaType>(*md);
	  memcpy(d,bin.d,getDataSize<MetaType>(*md)*sizeof(DataType));
	}
	else {
	  d = bin.d;
	}
      }
      //cerr<<"CP copy constr: md = "<<md<<" d = "<<d<<" own = "<<own<<" this = "<<this<<endl;
      //cerr<<"*** metadata:"<<endl;
      //      writeMeta<MetaType>(*md,cerr);
    } 

    /** main constructor */
    ContentPackage() 
      :md(NULL), d(NULL), own(true) {
      //cerr<<"CP main constr: md = "<<md<<" d = "<<d<<" this = "<<this<<endl;
    }

    virtual ~ContentPackage() {
      //cerr<<"CP destr begin, this = "<<this<<"\n";
      //cerr<<"   d = "<<d<<" own = "<<own<<endl;
      //if (d && own) delete [] d; d = NULL;
      if (d && own) deleteData<DataType,MetaType>(&d,&md);
      //cerr<<"  d = "<<d<<endl;
      //cerr<<"*** metadata: md="<<md<<endl;
      //writeMeta<MetaType>(*md,cerr);
      if (md) {
	//cerr<<"  delete md\n";
	delete md; md = NULL;
	//cerr<<"  md = "<<md<<endl;
      }
      //cerr<<"CP destr finish\n";
    }
        
    /** post-construction initialization - can only be done once! 
	Optional reference to external storage permits a CP to 
	act as a "view" of another hunk of memory. 
    */
    bool initialize(MetaType const & _md, 
		    LocalDataContainer<DataType> * p = NULL,
		    size_t start = 0) {
      if (md) return false;
      md = new MetaType(_md);
      //cerr<<"CP::initialize - md = "<<md<<" this = "<<this<<endl;
      //cerr<<"*** metadata:"<<endl;
      //writeMeta<MetaType>(*md,cerr);
      if (p) {
	if (getDataSize<MetaType>(*md) <= p->getSize()-start) {
	  d = &((p->getData())[start]);
	  own = false;
	  return true;
	}
	else {
	  //cerr<<"warning: ContentPackage::initialize\n";
	  //	  cerr<<"attempted to build view of external LDC\n";
	  //	  cerr<<"of length "<<getDataSize<MetaType>(*md)<<" starting at word "
	  //	   <<start<<" of LDC of length "<<p->getSize()<<"\n";
	  delete md;
	  md=NULL;
	  return false;
	}
      }
      else {
	//	d = new DataType[getDataSize<MetaType>(*md)];
	d=newData<DataType,MetaType>(*md);
	return true;
      }
    }

    /** initialization status query */
    bool isInitialized() const { if (md) return true; return false; }

    /** assignment */
    ContentPackage<DataType,MetaType> & 
    operator=(ContentPackage<DataType,MetaType> const & src) { 
      if (this != &src) {
	if (*md != *(src.md)) {
	  RVLException e;
	  e<<"Error: ContentPackage assignment op\n";
	  e<<"metadata of this differs from metadata of source\n";
	  throw e;
	}
	size_t n = getDataSize<MetaType>(*md);
	if (!d) d = new DataType[n];
	memcpy(d,src.d,n*sizeof(DataType));
      }
      return *this;
    }
    /** comparison */
    bool operator==(ContentPackage<DataType,MetaType> const & a) const {
      if (this == &a) return true;
      if (*md != *(a.md)) return false;
      if (getSize() != a.getSize()) return false;
      for (size_t i=0;i<getSize();i++) {
	if (getData()[i] != a.getData()[i]) return false;
      }
      return true;
    }
    bool operator!=(ContentPackage<DataType,MetaType> const & a) const { 
      return !operator==(a); 
    }

    /** access metadata, mutable */
    MetaType & getMetadata() { 
      if (!md) {
	RVLException e;
	e<<"Error: ContentPackage::getMetaData\n";
	e<<"metadata not initialized yet\n";
	throw e;
      }
      return *md;
    }

    /** access metadata, const */
    MetaType const & getMetadata() const { 
      if (!md) {
	RVLException e;
	e<<"Error: ContentPackage::getMetaData\n";
	e<<"metadata not initialized yet\n";
	throw e;
      }
      return *md;
    }

    /** access data array size */
    size_t getSize() const { 
      if (md) return getDataSize<MetaType>(*md); 
      return 0;
    }

    /** access data, mutable */
    DataType * getData() {
      if (!d) {
	RVLException e;
	e<<"Error: ContentPackage::getData\n";
	e<<"data not initialized yet\n";
	throw e;
      }
      return d;
    } 
    /** access data, const */
    DataType const * getData() const {
      if (!d) {
	RVLException e;
	e<<"Error: ContentPackage::getData\n";
	e<<"data not initialized yet\n";
	throw e;
      }
      return d;
    } 

    ostream & write(ostream & str) const {
      str<<"ContentPackage LDC, MetaData = \n";
      if (md) {
	writeMeta<MetaType>(*md,str);
      }
      else {
	str<<"(not initialized)\n";
      }
      str<<"  and data array\n";
      for (size_t i=0;i<this->getSize(); i++) {
	str<<"  "<<i<<". ";
	writeMeta<DataType>(d[i],str);
      }
      return str;
    }
  };


  /** Copy function for CPs. Default implementation may be specialized
      for more subtle copies, eg. of windows into multidimensional
      arrays or sectors of unstructured meshes or... */
  template<typename Source, typename Target>
  void copyCP(Source const & in, Target & out) {
    try {
      size_t n = min(in.getSize(),out.getSize());
      for (size_t i=0;i<n;i++) out.getData()[i]=in.getData()[i];
    }
    catch (RVLException & e) {
      e<<"\ncalled from copyCP\n";
      throw e;
    }
  }

  /** A abstract DC construct, implementing an implicit product
      structure with ContentPackage components. An explicit product
      structure could be overlain easily. Also exhibits
      RepeatableAccess interface and a private version of a server
      interface, which directly exposes internal data or stores
      it. This is more appropriate than the RecordServer interface, as
      a DC such as this one owns its data, and either stores it
      available for use or accesses windows into it over the
      internet. A standard RecordServer interface would imply needless
      copying.

      Only reset, get, and put - plus the usual write methods - need
      be implemented in child classes, as the eval methods of DC are
      implemented in the base class. 

      This class poses an interesting aliasing problem. The eval
      methods define data access patterns as for all
      DCs. RepeatableAccess models sequential access to a Cartesian
      product of DCs (as would happen wvia sequential file i/o). In
      evals with more than one input argument, several of these may be
      aliased. Reading a component DC from argument i will therefore
      read the wrong data if another argument, say j, aliases with i.

      The fix adopted here is to test each DC input argument for
      identity (by addresss) with previously accessed args. The CP in
      which PC::get() places the component data is dynamically
      allocated if the DC argument doesn't appear earlier on the
      list. Otherwise the pointer CP is set to the previously used
      pointer, thus explicitly aliasing the input data by address.

      The rationale for this design lies in the unique correspondence
      between RVL::Vector instances and their DCs. The only case in
      which Vectors share DCs comes in the Components construction for
      Vectors in ProductSpaces. Since evaluation happens at the level
      of components anyway, it is impossible for non-aliased Vector
      args to yield aliased DC args. In any case it is the DC args
      which pose the immediate aliasing issue.

      Of course, aliasing an input argument with the output argument
      violates RVL design principles. Since the eval methods are
      checking addresses in any case, violation of the anti-aliasing
      rule is also checked and throws an exception.

      Because FO evaluation is so central, and can accomodate so many
      different purposes, in RVL applications, we have left the eval
      methods virtual - child classes can override them to implement
      special cases, while delegating to the parent class methods for
      general evaluations.
  */

  template<typename DataType, typename MetaType>
  class PackageContainer: 
    public DataContainer {

  protected:

    /** Sequential access to component CPs. Should return more=false
	when there are no more CPs to be returned, else true. Thus a
	typical use is
	
	bool more = true;
	while (more) {
	  ContentPackage<...> const & cp = get(more);
          (do something)
        }

	Implementations should test for more=true and throw an
	exception if not, as an additional check.

	Since access is sequential, revisiting the contents of the 
	PC requires a call to reset().

	Child classes must store an internal CP
	buffer so that they can always return a ref to a CP.
    */
    virtual ContentPackage<DataType,MetaType> 
    & get(bool & more) = 0;
    virtual ContentPackage<DataType,MetaType> const 
    & get(bool & more) const = 0;

    /** store the internal CP buffer - to allow for out-of-order
	execution (eg. in parallel), child classes can index into
	global storage using metadata attributes of their internal CP
	buffer. */
    virtual void put() const = 0;

    /** restore PC to its original state. For out-of-core PCs this 
	should seek the file pointer to BOF. */
    virtual void reset() const = 0;

  public:

    PackageContainer() {}
    PackageContainer(PackageContainer<DataType,MetaType> const &) {}
    virtual ~PackageContainer() {}
    
    /** FO eval method - virtual to allow specialization in subclasses */
    virtual void eval(FunctionObject & f, 
                      std::vector<DataContainer const *> & x) {
      try {
	//	cerr<<"pc: before reset\n";
	//	this->write(cerr);
	this->reset();
	//	cerr<<"pc: after reset\n";
	//	this->write(cerr);
	std::vector<PackageContainer<DataType,MetaType> const *> pcx(x.size());
	for (size_t i=0;i<x.size();i++) {
	  // first check against aliasing
	  if (this == pcx[i]) {
	    RVLException e;
	    e<<"Error: PackageContainer::eval(FO)\n";
	    e<<"input argument "<<i<<" aliased with output\n";
	    throw e;
	  }
	  if (!(pcx[i] = dynamic_cast<PackageContainer<DataType,MetaType> 
		const *>(x[i]))) {
	    RVLException e;
	    e<<"Error: PackageContainer::eval(FO)\n";
	    e<<"at least one other arg is not a package container.\n";
	    for (size_t j=0;j<x.size();j++) {
	      e<<"*** arg "<<j<<":\n";
	      x[j]->write(e);
	    }
	    throw e;
	  }
	  pcx[i]->reset();
	}
	
	bool more = true;
	std::vector<DataContainer const *> cpx(x.size());

	//	cerr<<"cp: size="<<x.size()<<endl;
	//	cerr<<"cp: more="<<more<<endl;
	while (more) {
	  for (size_t i=0;i<x.size();i++) {
	    // check for alias
	    size_t j=i;
	    dejavu<PackageContainer<DataType,MetaType> const>(&j,pcx);
	    // if j has not changed value, then this is a non-aliased arg
	    if (j==i) { cpx[i]=&(pcx[i]->get(more)); }
	    // otherwise we've already seen it
	    else { cpx[i] = cpx[j]; }
	  }
	  //	  cerr<<"cp->get\n";
	  //	  this->write(cerr);
	  ContentPackage<DataType,MetaType> & cp = this->get(more);
	  //	  cerr<<"cp->eval\n";
	  cp.eval(f,cpx);
	  //	  cerr<<"this->put\n";
	  this->put();
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from PackageContainer::eval(FO)\n";
	throw e;
      }
    }
	          
    /** FOR eval method - virtual to allow specialization in subclasses */
    virtual void eval(FunctionObjectConstEval & f, 
                      vector<DataContainer const *> & x) const {
      try {
	this->reset();
	std::vector<PackageContainer<DataType,MetaType> const *> pcx(x.size());
	for (size_t i=0;i<x.size();i++) {
	  if (!(pcx[i] = 
		dynamic_cast<PackageContainer<DataType,MetaType> 
		const *>(x[i]))) {
	    RVLException e;
	    e<<"Error: PackageContainer::eval(FOR)\n";
	    e<<"at least one other arg is not a package container.\n";
	    throw e;
	  }
	  pcx[i]->reset();
	}
	bool more = true;
	std::vector<DataContainer const *> cpx(x.size());

	while (more) {
	  ContentPackage<DataType,MetaType> const & cp = this->get(more);
	  for (size_t i=0;i<x.size();i++) {
	    // check for alias with target
	    if (this == pcx[i]) {
	      cpx[i] = &cp;
	    }
	    else {
	      // check for alias
	      size_t j=i;
	      dejavu<PackageContainer<DataType,MetaType> const>(&j,pcx);
	      // if j has not changed value, then this is a non-aliased arg
	      if (j==i) { 
		cpx[i]=&(pcx[i]->get(more)); }
	      // otherwise we've already seen it
	      else { 
		cpx[i] = cpx[j]; 
	      }
	    }
	  }
	  cp.eval(f,cpx);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from PackageContainer::eval(FOR)\n";
	throw e;
      }
    }

  };

  /** An almost-concrete factory class for PackageContainers, lacking
      only one protected and one public method to be complete. Added
      clone method because copying of children in classes owning one
      of these will require virtual construction.

      The comparison methods may not be sufficiently discriminatory
      for various applications, as they merely typecheck, So they are
      made virtual and can be overridden in child classes.
  */
  template<typename DataType,typename MetaType>
  class PackageContainerFactory: public DataContainerFactory {
  protected:

    virtual PackageContainer<DataType,MetaType> * buildPC() const = 0;
    
  public:

    /** Typically implemnented by op new. */    
    virtual PackageContainerFactory<DataType,MetaType> * clone() const = 0;

    DataContainer * build() const { return buildPC(); }

    virtual bool compare( DataContainerFactory const & dcf) const {
      PackageContainerFactory<DataType,MetaType> const * ptr = NULL;
      if ((ptr = 
	   dynamic_cast< PackageContainerFactory<DataType,MetaType> const * >(&dcf)))
	return true;
      return false;
    }

    virtual bool isCompatible(DataContainer const & dc) const {
      PackageContainer<DataType,MetaType> const * ptr = NULL;
      if ((ptr = 
	   dynamic_cast< PackageContainer<DataType,MetaType> const * >(&dc)))
	return true;
      return false;
    }
      
    virtual ostream & write(ostream & str) const {
      str<<"PackageContainerFactory\n";
      return str;
    }    
  };
    
  /** SingleDataContainer is a wrapper for a single ContentPackage,
      and is a child rather than a template specialization of
      PackageContainer. The additional attributes are server
      functions (get and put). Any number of replications of the
      ContentPackage may be served up; when the specified number have
      been served, get returns false.q
  */
  
  template<typename Datatype, typename Metatype>
  class SingleDataContainer: public PackageContainer<Datatype,Metatype> {
    
  private:
    
    ContentPackage<Datatype,Metatype> gd;
    mutable int nrep;
    int nrep0;

  protected:

    ContentPackage<Datatype,Metatype> & get(bool & more) {
      nrep--;
      if (nrep<1) more = false; 
      return gd;
    }
    
    ContentPackage<Datatype,Metatype> const & get(bool & more) const {
      nrep--;
      if (nrep<1) more = false; 
      return gd;
    }

    void put() const {}

    void reset() const { nrep=nrep0; }

  public:
    
    SingleDataContainer(int _nrep = 1)
      : PackageContainer<Datatype,Metatype>(), 
	gd(), nrep(_nrep), nrep0(nrep) {}
    
    SingleDataContainer(SingleDataContainer<Datatype,Metatype> const & g)
      : PackageContainer<Datatype,Metatype>(g), gd(g.gd), 
	nrep(g.nrep), nrep0(g.nrep0) {}
    ~SingleDataContainer() {}
    
    void initialize(Metatype const & g) { gd.initialize(g); }
    Metatype const & getMetadata() const {
      if (gd.isInitialized()) return gd.getMetadata();
      RVLException e;
      e<<"Error: SingleDataContainer::getMetadata\n";
      e<<"cannot return reference to internal Metadata, as\n";
      e<<"this is not initialized yet.\n";
      throw e;
    }
    
    ostream & write(ostream & e) const {
      e<<"SingleDataContainer encapsulating ContentPackage object:\n";
      gd.write(e);
      return e;
    }
    
  };  

  /** Factory class for SingleDataContainer. */
  template<typename Datatype, typename Metatype>
  class SingleDataContainerFactory: 
    public PackageContainerFactory<Datatype,Metatype> {
  private:
    int nrep;
    Metatype * g;

  protected:
    PackageContainer<Datatype,Metatype> * buildPC() const {
      SingleDataContainer<Datatype,Metatype> * gdc 
	= new SingleDataContainer<Datatype,Metatype>(nrep);
      if (g) {
	gdc->initialize(*g);
	return gdc;
      }
      else {
	RVLException e;
	e<<"Error: SingleDataContainerFactory::buildPC\n";
	e<<"internal Metadata not initialized\n";
	throw e;
      }
    }
    
  public:
    SingleDataContainerFactory(int _nrep=1): nrep(_nrep), g(NULL) {}
    SingleDataContainerFactory
    (SingleDataContainerFactory<Datatype,Metatype> const & f)
      : nrep(f.nrep), g(NULL) { 
      if (f.g) g=new Metatype(*(f.g)); }
    ~SingleDataContainerFactory() { if (g) delete g; }
  
    PackageContainerFactory<Datatype,Metatype> * clone() const {
      if (g) {
	SingleDataContainerFactory<Datatype,Metatype> * f 
	  = new SingleDataContainerFactory<Datatype,Metatype>;
	f->initialize(*g);
	return f;
      }
      RVLException e;
      e<<"Error: SingleDataContainerFactory::clone\n";
      e<<"cannot construct copy because not initialized\n";
      throw e;
    }
    void initialize(Metatype const & _g) { g = new Metatype(_g); } 

    bool compare( DataContainerFactory const & dcf) const {
      try {
	if (!g) {
	  cerr<<"Warning: SingleDataContainerFactory::compare\n";
	  cerr<<"cannot compare uninitialized factory to others\n";
	  return false;
	}
	SingleDataContainerFactory<Datatype,Metatype> const * dcfptr = NULL;
	if ((dcfptr = dynamic_cast<SingleDataContainerFactory<Datatype,Metatype> const *>(dcfptr))) {
	  if (!(dcfptr->g)) {
	    cerr<<"Warning: GridDataContainerFactory::compare\n";
	    cerr<<"cannot compare this to uninitialized factory\n";
	    return false;
	  }
	  if ( *g == *(dcfptr->g)) return true;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from SingleDataContainerFactory::compare\n";
	throw e;
      }
    }	  

    bool isCompatible(DataContainer const & dc) const {
      try {
	if (!g) {
	  cerr<<"Warning: SingleDataContainerFactory::isCompatible\n";
	  cerr<<"cannot decide compatibility for uninitialized factory\n";
	  return false;
	}
	SingleDataContainer<Datatype,Metatype> const * dcptr = NULL;
	if ((dcptr = dynamic_cast<SingleDataContainer<Datatype,Metatype> const *>(&dc))) 
	if (*g == dcptr->getMetadata()) return true;
	return false;
      }
      catch (RVLException & e) {
	e<<"\ncalled from GridDataContainerFactory::isCompatible\n";
	throw e;
      }
    }   
    
    ostream & write(ostream & str) const {
      str<<"SingleDataContainerFactory";
      if (g) {
	str<<" using Metadata:\n";
	g->write(str);
      }
      else {
	str<<" (uninitialized)\n";
      }
      return str;
    }      
  
  };

}

#endif
