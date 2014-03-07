/*************************************************************************

Copyright Rice University, 2004, 2005, 2006.
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

#ifndef __RVL_VEC
#define __RVL_VEC

#include "data.hh"
#include "linalg.hh"
#include "write.hh"

namespace RVL {

  /** RVL abstract base class for Hilbert Spaces. A vector space over a
      field of scalars is set of objects, called vectors, together with
      two operations, scalar multiplication and vector addition, obeying
      certain axioms. Hilbert spaces add inner product op. 

      Access to vector objects is implemented through access to their
      underlying data structures (DataContainers). The build() method
      returns a pointer to a dynamically allocated DataContainer.

      The three linear algebra methods "teach the DataContainer
      how to behave like a vector" (S. Scott, Y2K). 
  */

  template<class Scalar> 
  class Space: public Writeable {

    /** convenient typedef for positive real type */
    typedef typename ScalarFieldTraits<Scalar>::AbsType NormRetType;

  protected:

    /** Spaces should always be allocated on the stack. Various RVL
	classes store Space references as data members, and RVL does
	not reference count. Therefore dynamic allocation of Spaces
	risks dangling references. However RVL does not virtually
	construct spaces: all actual space constructors are
	concrete. Therefore dynamic allocation of Spaces is seldom
	necessary in RVL applications, except possibly in child
	classes of Space (for instance to define a std::vector of
	subspaces). Since it is unnecessary and potentially leads to
	memory management errors, we forbid it to the "general
	public". 

	Version 1.0: for those who must, we will allow dynamically
	allocated spaces - caveat emptor!
    */

#ifndef RVL_OPERATOR_NEW_ENABLED
    void * operator new(size_t size) { 
      void * ptr;
      ptr = (void *) ::new unsigned char[size]; 
      return ptr;
    }
#endif

  public:

    Space() {}
    Space(const Space<Scalar> & sp) {}
    virtual ~Space(){}

    /** returns a dynamically allocated DataContainer. The type
	of this datacontainer encodes the data structure which realizes
	the storage appropriate for vectors in this space.
    */
    virtual DataContainer * buildDataContainer() const  = 0;

    /** compare two vector spaces - useful for sanity 
	checking in applications */
    virtual bool operator ==(const Space<Scalar> & sp) const = 0;
    bool operator !=(const Space<Scalar> & sp) const {
      try {
	return !operator==(sp);
      }
      catch (RVLException & e) {
	e<<"\ncalled from Space::operator!=\n";
	throw e;
      }
    }

    /** detect compatibility with foreign data container -
	also useful for sanity checking */
    virtual bool isCompatible(DataContainer const & dc) const = 0;

    /** core linear algebra methods */

    /** inner product */
    virtual Scalar inner(DataContainer const & x, 
			 DataContainer const & y) const = 0;

    /** assign to data of zero vector */
    virtual void zero(DataContainer & x) const = 0;

    /** linear combination y=ax+by. Output argument y may not
	be alised with input argument x. 
     */
    virtual void linComb(Scalar a, DataContainer const & x,
			 Scalar b, DataContainer & y) const = 0;

    /** virtual convenience methods, implemented in the base class
	using linComb. If linComb is implemented carefully these
	implementations are often adequate. In some cases it may be
	more efficient to override them. */

    /** Copy: target may not be aliased with source. */
    virtual void copy(DataContainer & tgt, DataContainer const & src) const {
      try {
	Scalar one = ScalarFieldTraits<Scalar>::One();
	Scalar zip = ScalarFieldTraits<Scalar>::Zero();
	linComb(one,src,zip,tgt);
      }
      catch (RVLException & e) {
	e<<"\ncalled fom Space::copy\n";
	throw e;
      }
    }

    /** negate, unary version. Can be efficient with careful
        implementation of DC classes and linComb, so that 
        unreferenced memory is not allocated. Else this implementation
	can be overridden. 
    */
    virtual void negate(DataContainer & tgt) const {
      try {
	DataContainer const * jnk = this->buildDataContainer();
	Scalar zip = ScalarFieldTraits<Scalar>::Zero();
	Scalar one = ScalarFieldTraits<Scalar>::One();
	linComb(zip,*jnk,-one,tgt);
	delete jnk;
      }
      catch (RVLException & e) {
	e<<"\ncalled fom Space::negate\n";
	throw e;
      }
    }

    /** negation, binary version. Target may not be aliased with source.*/
    virtual void negate(DataContainer & tgt,
			DataContainer const & src) const {
      try {
	Scalar one = ScalarFieldTraits<Scalar>::One();
	Scalar zip = ScalarFieldTraits<Scalar>::Zero();
	linComb(-one,src,zip,tgt);
      }
      catch (RVLException & e) {
	e<<"\ncalled fom Space::negate\n";
	throw e;
      }
    }

    /** scale, unary version. Can be efficient with careful
        implementation of DC classes and linComb, so that 
        unreferenced memory is not allocated. Else this implementation
	can be overridden.
    */
    virtual void scale(DataContainer & tgt,
		       Scalar c) const {
      try {
	DataContainer * jnk = this->buildDataContainer();
	zero(*jnk);
	Scalar zip = ScalarFieldTraits<Scalar>::Zero();
	linComb(zip,*jnk,c,tgt);
	delete jnk;
      }
      catch (RVLException & e) {
	e<<"\ncalled fom Space::scale\n";
	throw e;
      }
    }

    /** Scale, binary version. Target may not be aliased with source. */
    virtual void scale(DataContainer & tgt,
		       Scalar c,
		       DataContainer const & src) const {
      try {
	Scalar zip = ScalarFieldTraits<Scalar>::Zero();
	linComb(c,src,zip,tgt);
      }
      catch (RVLException & e) {
	e<<"\ncalled fom Space::scale\n";
	throw e;
      }
    }

    /** return Norm Squared.  The norm always returns a non-negative
	real number */
    NormRetType normsq(DataContainer const & x) const {
      try {
	return abs(inner(x,x));
      }
      catch (RVLException & e) {
	e<<"\n*** called from Space::normsq\n";
	throw e;
      }
    }

    /** return Norm. The norm always returns a non-negative real number */
    NormRetType norm(DataContainer const & x) const {
      try {
	NormRetType f=normsq(x);
	return (NormRetType) sqrt(f);
      }
      catch (RVLException & e) {
	e<<"\n*** called from Space::norm\n";
	throw e;
      }
    }

  };

  /** Standard modular RVL space class. This construction implements
      Space by pairing a DataContainerFactory (access to vectors) and a
      LinearAlgebraPackage (access to ops), both defined abstractly
      through methods returning appropriate references. All other
      methods (expressing vector space behaviour) implemented in this
      base class. 

      A typical implementation will provide a LinAlgPack and a 
      DataContainerFactory, either as arguments to the constructor
      or through some internal construction mechanism, and return
      references to these. Thus only two methods need be implemented.

      We could provide a concrete subclass which takes the LAP and
      DCF as constructor args. However this class would hardly ever
      be used, as LAPs and DCFs are not natural data in most apps.

      Note that the LinearAlgebraPackage needs a separate template
      parameter for the DataType in the LDCs.  In most cases, this is
      the same as the Scalar type.
  */

  template<class Scalar, class DataType = Scalar> 
  class StdSpace: public Space<Scalar> {

  public:

    StdSpace() {}
    StdSpace(const StdSpace<Scalar, DataType> & sp) {}
    virtual ~StdSpace(){}

    /** access to DataContainerFactory */
    virtual  DataContainerFactory const & getDCF() const = 0;
    /** access to LinearAlgebraPackage */
    virtual  LinearAlgebraPackage<Scalar> const & getLAP() const = 0;

    // default virtual data container constructor
    DataContainer * 
    buildDataContainer() const { return getDCF().build(); }

    /** comparison: first check for identity of addresses. If these
	are equal, a very low-cost comparison has succeeded. Otherwise,
	delegate the decision to the comparison methods of the DC factory
	and LA package. */
    bool operator ==(const Space<Scalar> & sp) const {
      try {
	if (this == &sp) return true;
	const StdSpace<Scalar, DataType> & stdsp = 
	  dynamic_cast<const StdSpace<Scalar, DataType> &>(sp);
	return (getDCF().compare(stdsp.getDCF()) && 
		getLAP().compare(stdsp.getLAP()));
      }
      catch (bad_cast) {
	return 0;
      }
    }

    // detect compatibility with data container
    bool isCompatible(DataContainer const & dc) const {
      return getDCF().isCompatible(dc);
    }

    // inner product
    Scalar inner(DataContainer const & x, 
		 DataContainer const & y) const {
      try {
	getLAP().inner().setValue();
	vector<DataContainer const *> vy(1);
	vy[0]=&y;
	x.eval(getLAP().inner(),vy);
	return getLAP().inner().getValue();
      }
      catch (RVLException & e) {
	e<<"\ncalled from Space::inner\n";
	throw e;
      }
    }

    // zero vector
    void zero(DataContainer & x) const {
      try {
	vector<DataContainer const *> vy(0);
	x.eval(getLAP().zero(),vy);
      }
      catch (RVLException & e) {
	e<<"\ncalled from StdSpace::zero\n";
	throw e;
      }
    }

    // linear combination
    void linComb(Scalar a, DataContainer const & x,
		 Scalar b, DataContainer & y) const {
      (getLAP().linComb()).setScalar(a,b);
      try {
	vector<DataContainer const *> vx(1);
	vx[0]=&x;
	y.eval(getLAP().linComb(),vx);
      }
      catch (RVLException & e) {
	e<<"\ncalled from StdSpace::linComb\n";
	throw e;
      }
    }

    ostream & write(ostream & str) const {
      str<<"StdSpace defined by DataContainerFactory\n";
      getDCF().write(str);
      str<<"and LinearAlgebraPackage\n";
      getLAP().write(str);
      return str;
    }
  };

  /** A converse to StdSpace: takes any space and makes manifest its
      inner DataContainerFactory. Note that a SpaceDCF object
      instantiated from a StdSpace is functionally equivalent to the
      latter's DCF. */
  template<class Scalar> 
  class SpaceDCF: public DataContainerFactory {

  private:

    Space<Scalar> const & sp;
    SpaceDCF();

  public:
    
    SpaceDCF(Space<Scalar> const & _sp): sp(_sp) {}
    SpaceDCF(SpaceDCF<Scalar> const & f): sp(f.sp) {}
    ~SpaceDCF() {}

    DataContainer * build() const { return sp.buildDataContainer(); }

    Space<Scalar> const & getSpace() const { return sp; }

    bool compare( DataContainerFactory const & dcf ) const {
      SpaceDCF<Scalar> const * p = NULL;
      p = dynamic_cast< SpaceDCF<Scalar> const * >(&dcf);
      if (p) return (this->getSpace()==p->getSpace());
      return false;
    }

    bool isCompatible(DataContainer const & dc) const {
      return sp.isCompatible(dc);
    }

    ostream & write(ostream & str) const {
      str<<"Space-derived DataContainerFactory based on Space:\n";
      sp.write(str);
      return str;
    }
  };
    
  template<class Scalar> 
  class OpComp;

  // forward declaration
  template<class Scalar>
  class Components;

  /** RVL Vector class. Concrete class with all aspects of a natural
      vector interface implemented through delegation to Space and
      DataContainer data members. Provision for evaluation of
      arbitrary FunctionObject and FunctionObjectConstEval instances
      implemented through delegation to DataContainer::eval.
  */

  template<class Scalar>
  class Vector: public Writeable {

    friend class Components<Scalar>;

  private: 

    // the space to which the vector belongs
    const Space<Scalar> & sp;
    // the data container for the vector data
    // mutable to allow delayed allocation in a const member function
    mutable DataContainer * d; 
    // flags:
    //   own:       true if allocation of DataContainer is internal, false 
    //              if it's externally allocated
    //   initizero: initialize with all zeros if true, leave uninitialized else
    bool own, initzero;
    // reference to ver
    unsigned int & verref;
    // version index - incremented at every assignment
    mutable unsigned int ver;

    // default construction disallowed
    Vector();

  protected:

    /** The following three functions expose the internal data of
	Vector to its friends and children classes (there are only two
	- Components and LocalVector) and enable construction of a
	Vector from its constituent parts.  General objects should not
	be able to either access the private data of a Vector or build
	a Vector from random pieces, so these functions are
	protected. */

    /** access to DataContainer */
    /*
    DataContainer const * getDataContainer() const { 
      if( !d) {
	d = sp.buildDataContainer();
	if(initzero) sp.zero(*d);
      }
      return d; 
    }
    */

    DataContainer * getDataContainer() const { 
      if( !d) {
	d = sp.buildDataContainer();
	if(initzero) sp.zero(*d);
      }
      return d; 
    }

    /** Protected vector copy constructor, shallow copy semantics.
	In default case, this object owns only references to external
	data.
    */
    Vector(const Space<Scalar> & _sp, 
	   DataContainer * _d, 
	   unsigned int & _verref, 
	   bool _own = false)
      :sp(_sp), d(_d), own(_own), verref(_verref) {
      if (!(sp.isCompatible(*getDataContainer()))) {
	RVLException e; e<<"*** Error: Vector constructor (sp,dc)\n";
	e<<"*** input data container not compatible with space\n";
	e<<"*** this space:\n";
	sp.write(e);
	e<<"*** data container:\n";
	d->write(e);
	throw e;
      }
      ver = verref;
    }

    /** by hiding the constructor, makes it available to child 
	classes */
    Vector<Scalar> * build_from_kit(const Space<Scalar> & _sp, 
				    DataContainer * _d, 
				    unsigned int & _verref, 
				    bool _own = false) {
      return new Vector<Scalar>(_sp,_d,_verref,_own);
    }

    /** Protected constructor which lets child classes wrap
	either space or datacontainer data members in other
	interfaces. */
    Vector(const Vector<Scalar> * v)
      : sp(v->sp), d(v->d), own(false), ver(0), verref(ver) {}

    /** The Evaluation classes have Vectors as data
	members. Accordingly, general RVL apps should avoid dynamic
	allocation of Vectors, just as they should avoid dynamic
	allocation of Spaces - to eliminate the risk of dangling
	references, and because it is unnecessary. However the
	Components and LocalVector classes do need to construct
	vectors dynamically, using the protected constructors
	described above. Therefore we make operator new protected. 

	Version 1.0: user control, for the bold. "There are old
	sailors, and there are bold sailors, but there are no old,
	bold sailors".
    */
#ifndef RVL_OPERATOR_NEW_ENABLED
    void * operator new(size_t size) { 
      void * ptr;
      ptr = (void *) ::new unsigned char[size]; 
      return ptr;
    }
#endif
    unsigned int & getVersionRef() const { return verref; }
    
  public:

    using RVL::Writeable::write;

    /** convenient typedef for positive real type */
    typedef typename ScalarFieldTraits<Scalar>::AbsType NormRetType;

    /** Copy constructor. Implements deep copy, as data is actually
	owned by this object. */
    Vector(const Vector<Scalar> & x)
      : sp(x.sp), 
	d(sp.buildDataContainer()), 
	own(true), 
	initzero(false), 
	verref(ver) {
      try {
	sp.copy(*getDataContainer(),*(x.getDataContainer()));
	ver=0;
      }
      catch (RVLException & e) {
	e<<"\ncalled from Vector copy constructor\n";
      }
    }

    /** Standard constructor. Constructs a member of a given
	vector space. Note that data is not initialized. For 
	large vector sizes, this constructor is therefore to
	be preferred to the copy constructor when a choice
	exists and data copying is not needed.

	The data container can be initialized to zero by setting
	the optional second parameter to true.
    */
    Vector(const Space<Scalar> & _sp, bool _initZero = false)
      : sp(_sp), d(NULL), own(true), 
	initzero(_initZero), verref(ver), ver(0) {}
  
    /** Destructor. Deallocates dynamically allocated
	DataContainer. */
    ~Vector() { if (own&&(d!=NULL)) delete d; }

    /** access to Space */
    const Space<Scalar> & getSpace() const { return sp; }

    /** returns nonzero if this vector is a member of this space. */
    bool inSpace(const Space<Scalar> & sp1) const { return (sp==sp1); }
    /** returns nonzero if this vector is member of same space as
	argument vector. */
    bool inSameSpace(const Vector<Scalar> & x) const { return inSpace(x.sp); }

    /** Version numbers are automatically incremented in the eval
	methods for which this is target. Else they are left the
	same. If a version number seems to be getting to high, make
	sure that the writesData() method of all FunctionObjects are
	returning correct information.
    */
    unsigned int getVersion() const { return ver; }
    void incrementVersion() { 
      if (ver < numeric_limits<unsigned int>::max()) { ver++; }
      else {
	RVLException e;
	e<<"Error: Vector::incrementVersion\n";
	e<<"run out of unsigned ints!!!\n";
	throw e;
      }
    }
    
    /** generic evaluation of a FunctionObject */
    void eval(FunctionObject & f,
	      vector<Vector<Scalar> const *> & x) {
      try {
	vector<DataContainer const *> dx(x.size());
	for (size_t i=0;i<x.size();i++) {
	  dx[i] = x[i]->getDataContainer();
	}
	(this->getDataContainer())->eval(f,dx);
	incrementVersion();
      }
      catch (RVLException & e) {
	e<<"\ncalled from Vector::eval (generic)\n";
	throw e;
      }
    }

    /** eval convenience interface: no additional vector args */
    void eval(FunctionObject & f) {
      try {
	vector<Vector<Scalar> const *> v(0);
	this->eval(f,v);
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::eval(fo &)\n";
	throw e;
      }
    }

    /** eval convenience interface: one additional vector arg */
    void eval(FunctionObject & f, 
	      const Vector<Scalar> & x) {
      try {
	vector<Vector<Scalar> const *> v(1);
	v[0]=&x;
	this->eval(f,v);
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::eval(fo &, vec &)\n";
	throw e;
      }
    }
  
    /** eval convenience interface: two additional vector args */
    void eval(FunctionObject & f, 
	      const Vector<Scalar> & x,
	      const Vector<Scalar> & y) {
      try {
	vector<Vector<Scalar> const *> v(2);
	v[0]=&x;
	v[1]=&y;
	this->eval(f,v);
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::eval(fo &, vec &, vec &)\n";
	throw e;
      }
    }

    /** eval convenience interface: three additional vector args */
    void eval(FunctionObject & f, 
	      const Vector<Scalar> & x,
	      const Vector<Scalar> & y,
	      const Vector<Scalar> & z) {
      try {
	vector<Vector<Scalar> const *> v(3);
	v[0]=&x;
	v[1]=&y;
	v[2]=&z;
	this->eval(f,v);
      }      
      catch (RVLException & e) {
	e<<"\n*** called from Vector::eval(fo &, vec &, vec &, vec &)\n";
	throw e;
      }
    }

    /** generic evaluation of a FunctionObjectConstEval */
    void eval(FunctionObjectConstEval & f,
	      vector<Vector<Scalar> const *> & x) const {
      try {
	vector<DataContainer const *> dx(x.size());
	for (int i=0;i<(int)x.size();i++) {
	  dx[i] = x[i]->getDataContainer();
	}
	(this->getDataContainer())->eval(f,dx);
      }
      catch (RVLException & e) {
	e<<"\ncalled from Vector::eval (FOR, generic)\n";
	throw e;
      }
    }

    /** eval convenience interface: no additional vector args */
    void eval(FunctionObjectConstEval & f) const {
      try {
	vector<Vector<Scalar> const *> v(0);
	this->eval(f,v);
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::eval(for &)\n";
	throw e;
      }
    }

    /** eval convenience interface: one additional vector arg */
    void eval(FunctionObjectConstEval & f, 
	      const Vector<Scalar> & x) const {
      try {
	vector<Vector<Scalar> const *> v(1);
	v[0]=&x;
	this->eval(f,v);
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::eval(for &, vec &)\n";
	throw e;
      }
    }
  
    /** eval convenience interface: two additional vector args */
    void eval(FunctionObjectConstEval & f, 
	      const Vector<Scalar> & x,
	      const Vector<Scalar> & y) const {
      try {
	vector<Vector<Scalar> const *> v(2);
	v[0]=&x;
	v[1]=&y;
	this->eval(f,v);
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::eval(for &, vec &, vec &)\n";
	throw e;
      }
    }

    /** eval convenience interface: three additional vector args */
    void eval(FunctionObjectConstEval & f, 
	      const Vector<Scalar> & x,
	      const Vector<Scalar> & y,
	      const Vector<Scalar> & z) const {
      try {
	vector<Vector<Scalar> const *> v(3);
	v[0]=&x;
	v[1]=&y;
	v[2]=&z;
	this->eval(f,v);
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::eval(for &, vec &, vec &, vec &)\n";
	throw e;
      }
    }

    /** Inner product. */
    Scalar inner(const Vector<Scalar> & y) const {
      try {
	return sp.inner(*getDataContainer(),*(y.getDataContainer()));
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::inner(Vector &)\n";
	throw e;
      }
    }

    /** Linear Combination. Default: this = a*x + this (i.e. axpy). */
    void linComb(Scalar a, const Vector<Scalar> & x,
		 Scalar b = ScalarFieldTraits<Scalar>::One()) {
      try {
	sp.linComb(a,*(x.getDataContainer()),
		   b,*(this->getDataContainer()));
	incrementVersion();
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::linComb (this = a x + b this)\n";
	throw e;
      }
    }

    /** Assignment to zero vector */
    void zero() {
      try {
	sp.zero(*(this->getDataContainer()));
	incrementVersion();
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::zero()\n";
	throw e;
      }
    }

    /** convenience methods, defined in terms of space convenience methods */

    /** Copy */
    void copy( const Vector<Scalar> & x) {
      try {
	sp.copy(*(this->getDataContainer()),*(x.getDataContainer()));
	incrementVersion();
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::copy()\n";
	throw e;
      }
    }

    /** Scale, unary version. If low-level DCs defer allocation
	and linComb is properly implemented, no actual data beyond a
	few pointers will be allocated or accessed for the workspace
	vector z.
    */
    void scale(Scalar c) {
      try {
	sp.scale(*(this->getDataContainer()),c);
	incrementVersion();
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::scale()\n";
	throw e;
      }
    }

    /** Scale, binary version */
    void scale(Scalar c, const Vector<Scalar> & x) {
      try {
	sp.scale(*(this->getDataContainer()),c,*(x.getDataContainer()));
	incrementVersion();
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::scale()\n";
	throw e;
      }
    }

    /** Negation, unary version. */
    void negate() {
      try {
	sp.negate(*(this->getDataContainer()));
	incrementVersion();
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::negate()\n";
	throw e;
      }
    }

    /** Negation, binary version */
    void negate(const Vector<Scalar> & x) {
      try {
	sp.negate(*(this->getDataContainer()),*(x.getDataContainer()));
	incrementVersion();
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::negate()\n";
	throw e;
      }
    }

    /** return Norm Squared.  The norm always returns a non-negative
	real number */
    NormRetType normsq() const {
      try {
	return sp.normsq(*(this->getDataContainer()));
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::normsq\n";
	throw e;
      }
    }

    /** return Norm. The norm always returns a non-negative real number */
    NormRetType norm() const {
      try {
	return sp.norm(*(this->getDataContainer()));
      }
      catch (RVLException & e) {
	e<<"\n*** called from Vector::norm\n";
	throw e;
      }
    }

    ostream & write(ostream & str) const {
      str<<"Vector Object"<<endl;
      str<<"member of space:"<<endl;
      sp.write(str);
      str<<"data container:"<<endl;
      (this->getDataContainer())->write(str);
      return str;
    }
  };
    
  /** This class references a vector and will store a version number.
      The update method returns true if an update was necessary (meaning
      the vector had been modified since the last update).
  */
  template<class Scalar>
  class WatchedVecRef {

  private:

    Vector<Scalar> & x;
    mutable unsigned int ver;
    WatchedVecRef() {}

  public:

    WatchedVecRef(const Vector<Scalar> & _x)
      : x(const_cast< Vector<Scalar> & >(_x)), 
	ver(x.getVersion()) {}
    WatchedVecRef(const WatchedVecRef<Scalar> & w)
      : x(w.x), ver(w.x.getVersion()) {}
    virtual ~WatchedVecRef() {}

    Vector<Scalar> & get() { return x; }
    Vector<Scalar> const & get() const { return x; }

    bool update() const {
      try {
	if (ver < x.getVersion()) {
	  ver = x.getVersion();
	  return true;
	}
	return false;
      }
      catch (RVLException & e) {
	e<<"\ncalled from WatchedVecRef::update\n";
	throw e;
      }
    }
  };

  /** space membership test - turns standard test into 
      one-liner 
  */
  template<typename Scalar> 
  void SpaceTest(Space<Scalar> const & sp, 
		 Vector<Scalar> const & v,
		 std::string msg) {
    if (!v.inSpace(sp)) {
      RVLException e; 
      e<<"Error: "<<msg<<"\n";
      e<<"vector not in space\n";
      e<<"vector:\n";
      v.write(e);
      e<<"space:\n";
      sp.write(e);
      throw e;
    }
  }

}

#endif



