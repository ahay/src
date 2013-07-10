#ifndef __TSOPT_GRIDDIFF_OPS__
#define __TSOPT_GRIDDIFF_OPS__

//#define VERBOSE

#include "mpigridpp.hh"
#include "mpiserialfo.hh"
#include "griddec.hh"
#include "local.hh"
#include "linop.hh"
#include "productspace.hh"
#include "gridrwfos.hh"

namespace TSOpt {
  using RVL::ScalarFieldTraits;
  using RVL::SpaceTest;
  using RVL::LinearOp;
  using RVL::Space;
  using RVL::ProductSpace;
  using RVL::Vector;
  using RVL::Components;
  //  using RVL::MPISerialFunctionObjectRedn;

  /** Differetiation Operator defined on GridSpace<Scalar> or its produts 
   * differentiating along given axis
   */
  template<typename Scalar> 
  class GridDiffOp: public LinearOp<Scalar> {
 
  private:
    // runtime error when dom not GridSpace<Scalar> or its products 
    Space<Scalar> const & dom;
    
    int const diffaxisid;
    int buf_capacity;   // units: number of Scalar type data 
    int sg_size;        // local subgrid size

    mutable IPNT ioff; 
    mutable int exdoff;
    
    GridDiffOp();

    bool update_ioff() const {
      /*
      Grid<Scalar> const & g = dom.getGrid();
      int carry=1;
      for (int i=1;i<g.getDim();i++) {
	ioff[i]+=carry;
	if (ioff[i]<g.getData()[i].n) carry=0;
	else {
	  carry=1;
	  if (i != diffaxisid - 1) ioff[i]=0;
	  else ioff[i] = g.getAxis(i).n - 1;
	}
      }
      exdoff+=carry;
      if (exdoff<g.getExdAxis().n) carry=0;

      if (carry) return false;
     
      return true;
      */
      return false;
    }

    void reset_ioff(Grid<Scalar> const & gref) const{
      IASN(ioff,IPNT_0);
      if (diffaxisid == gref.getExdAxis().id) 
	exdoff = gref.getExdAxis().n - 1;
      else {
	ioff[diffaxisid-1] = gref.getAxis(diffaxisid-1).n - 1;
	exdoff=0;
      }      
    }    

    void Dfun(ContentPackage< Scalar, Grid<Scalar> > & _cpbuf) const{
      Grid<Scalar> & _g = _cpbuf.getMetadata();
      Scalar * _buf = _cpbuf.getData();
      int cpbuf_size = _g.getDataSize();
      
      int nd;
      Scalar dd;
      int dnd = 1;
      if (diffaxisid == _g.getExdAxis().id){
	nd = _g.getExdAxis().n;
	dd = _g.getExdAxis().d;
	dnd = _g.getDataSize()/nd;
      }
      else {
	nd = _g.getAxis(diffaxisid - 1).n;
	dd = _g.getAxis(diffaxisid - 1).d;
	for (int kk = 0; kk< diffaxisid -1; kk++)
	  dnd *= _g.getAxis(kk).n;
      }

      /*-- previous working approach --*/
      /*
#ifdef VERBOSE
      cerr<<"----- GridDiffOp::Dfun(), debugging output --------"<<endl;
      cerr<<"... dnd ="<<dnd<<", size ="<<_g.getDataSize()<<", diffaxisid ="<<diffaxisid<<", nd ="<<nd<<endl;
      int lcount = 0;
      int lloop = 0;
#endif

      int indd = nd-1;
      for (int jj = cpbuf_size - 1; jj>=0; jj--) {
	indd = (jj/dnd)%nd;

#ifdef VERBOSE
	lcount++;
	if (lcount == dnd) {
	  lloop++;
	  cerr<<"... lloop = "<<lloop<<", indd = "<<indd<<endl;
	  lcount=0;
	}
#endif

	if (indd>0)
	  _buf[jj] = (_buf[jj] - _buf[jj - dnd])/dd;
	else
	  _buf[jj] = 0;
      }
      */

      /*-- new approach (vecterizable) --*/
      int dslicesz = dnd*nd;
      int nstrid = cpbuf_size/dslicesz;
      Scalar * _lbuf = NULL;
      for (int kk= 0; kk < nstrid; kk++) {
	_lbuf = _buf + kk * dslicesz + dslicesz - dnd;

	for (int jj= nd - 1; jj> 0; jj--) {
	  for (int ii=0; ii < dnd; ii++) {
	    _lbuf[ii] = (_lbuf[ii] - _lbuf[ii - dnd])/dd;
	  }
	  _lbuf -= dnd;
	}

	memset((void *) (_lbuf), 0, sizeof(Scalar) * dnd);
	// for (int ii= 0; ii < dnd; ii++) {
	//   _lbuf[ii] = 0;
	// }
      }
      
#ifdef VERBOSE
      cerr<<"----- GridDiffOp::Dfun(), exit --------"<<endl;  
#endif
    }
    
    void DTfun(ContentPackage< Scalar, Grid<Scalar> > & _cpbuf) const{
      Grid<Scalar> & _g = _cpbuf.getMetadata();
      Scalar * _buf = _cpbuf.getData();
      int cpbuf_size = _g.getDataSize();
 
      int nd;
      Scalar dd;
      int dnd = 1;
      if (diffaxisid == _g.getExdAxis().id){
	nd = _g.getExdAxis().n;
	dd = _g.getExdAxis().d;
	dnd = _g.getDataSize()/nd;
      }
      else {
	nd = _g.getAxis(diffaxisid - 1).n;
	dd = _g.getAxis(diffaxisid - 1).d;
	for (int kk = 0; kk< diffaxisid -1; kk++)
	  dnd *= _g.getAxis(kk).n;
      }

      /*-- previous working approach --*/
      /*
#ifdef VERBOSE
      cerr<<"----- GridDiffOp::DTfun(), debugging output --------"<<endl;
      cerr<<"...adj... dnd ="<<dnd<<", size ="<<_g.getDataSize()<<", diffaxisid ="<<diffaxisid<<", nd ="<<nd<<endl; 
      int lcount = 0;
      int lloop = 0;
#endif
     
      int indd = 0;
      if (nd >= 3) 
	for (int jj =0; jj< cpbuf_size; jj++) {
	  indd = (jj/dnd)%nd;

#ifdef VERBOSE
	  lcount++;
	  if (lcount == dnd) {
	    lloop++;
	    cerr<<"...adj... lloop = "<<lloop<<", ind = "<<indd<<endl;
	    lcount=0;
	  }
#endif

	  if (indd == 0)
	    _buf[jj] = -_buf[jj + dnd]/dd;
	  else if (indd < nd - 1)
	    _buf[jj] = (_buf[jj] - _buf[jj + dnd])/dd;
	  else
	    _buf[jj] /= dd;
	}
      else // nd == 2
	for (int jj =0; jj< cpbuf_size; jj++) {
	  indd = (jj/dnd)%nd;

#ifdef VERBOSE
	  lcount++;
	  if (lcount == dnd) {
	    lloop++;
	    cerr<<"...adj... lloop = "<<lloop<<", ind = "<<indd<<endl;
	    lcount=0;
	  }
#endif

	  if (indd == 0)
	    _buf[jj] = -_buf[jj + dnd]/dd;
	  else
	    _buf[jj] /= dd;
	}
      */

      /*-- new approach (vecterizable) --*/   
      int dslicesz = dnd*nd;
      int nstrid = cpbuf_size/dslicesz;
      Scalar * _lbuf = NULL;
      for (int kk=0; kk < nstrid; kk++) {
	_lbuf = _buf + kk * dslicesz;
	
	//-- jj= 0
	for (int ii=0; ii < dnd; ii++) 
	  _lbuf[ii] = -_lbuf[ii + dnd]/dd;
	
	_lbuf += dnd;	  
	//-- jj= 1 to nd-2
	for (int jj= 1; jj< nd-1; jj++) {
	  for (int ii=0; ii < dnd; ii++) {
	    _lbuf[ii] = (_lbuf[ii] - _lbuf[ii + dnd])/dd;
	  }
	  _lbuf += dnd;	  
	}
	// jj= nd - 1
	for (int ii=0; ii < dnd; ii++) 
	  _lbuf[ii] /= dd;
      }      

#ifdef VERBOSE
      cerr<<"----- GridDiffOp::DTfun(), exit --------"<<endl;  
#endif
      
    }
    
  protected:

    using LinearOp<Scalar>::apply;
    
    void apply(Vector<Scalar> const & x, 
	       Vector<Scalar> & y) const{
      try{
	SpaceTest(this->getDomain(),x,"TSOpt::GridDiffOp::apply (dom)");
	SpaceTest(this->getRange(),y,"TSOpt::GridDiffOp::apply (rng)");
	Components<Scalar> cx(x);
	Components<Scalar> cy(y);
	int csize = cx.getSize();
	
	for (int k=0; k< csize; k++) {
	  // extract grid for current component
	  try{
#ifdef IWAVE_USE_MPI
	    MPIGridSpace<Scalar> const & xsp = 
	      dynamic_cast<MPIGridSpace<Scalar> const &>(cx[k].getSpace());
#else
	    GridSpace<Scalar> const & xsp = 
	      dynamic_cast<GridSpace<Scalar> const &>(cx[k].getSpace());	  
#endif
	  
	    Grid<Scalar> const & gref = xsp.getGrid();
	    int g_size = gref.getDataSize();
	    int g_eid = gref.getExdAxis().id;
	    //	    int g_dim = gref.getDim();
	    //	    int g_edim = gref.getExdDim();
	    int g_nd = 1; 
	    int g_dnd = 1;
	    if (diffaxisid == g_eid) {
	      g_nd = gref.getExdAxis().n;
	      if(g_nd < 2) {
		RVLException e;
		e<<"Error: GridDiffOp::applyOp \n";
		e<<"the number of grid points along differetiating direction is less than 2 \n";
		throw e;
	      }
	      g_dnd = g_size/g_nd;
	    }
	    else {
	      g_nd = gref.getAxis(diffaxisid-1).n;
	      if(g_nd < 2) {
		RVLException e;
		e<<"Error: GridDiffOp::applyOp \n";
		e<<"the number of grid points along differetiating direction is less than 2 \n";
		throw e;
	      }
	      for (int kk = 0; kk< diffaxisid -1; kk++)
		g_dnd *= gref.getAxis(kk).n;
	    }      

	    // create working buffer
	    ContentPackage<Scalar, Grid<Scalar> > cpbuf;
	    Grid<Scalar> sg(gref);
	    /*
	      for(int i = RARR_MAX_NDIM; i > 1; i--){
	      if (i != diffaxisid) sg.getData()[i-1].n = 1;	  
	      }
	    */
	    int sgeid = sg.getExdAxis().id;
	    //if (sgeid != diffaxisid) sg.getExdAxis().n = 1;
	    
	    cpbuf.initialize(sg);
	    reset_ioff(gref);		
	    
	    do {
	      // record current origin position
	      int itmp;
	      if (sgeid == diffaxisid) itmp = 0;
	      else itmp = exdoff;
	      sg.getExdAxis().o =
		itmp * gref.getExdAxis().d + gref.getExdAxis().o;
	      
	      for (int i=1;i<gref.getDim();i++){
		if (i == diffaxisid-1) itmp = 0;
		else itmp = ioff[i];
		sg.getData()[i].o =
		  itmp * gref.getData()[i].d + gref.getData()[i].o;
	      }
	  
 	      GridReader<Scalar> r(cpbuf);
	      cx[k].eval(r);
	      Dfun(cpbuf);
	      GridWriter<Scalar> w(cpbuf);
	      cy[k].eval(w);
	    }
	    while(update_ioff());
	    
	  }
	  catch (std::bad_cast &bd){
	    RVLException e;
	    e<<bd.what()<<"\n";
	    e<< "GridDiffOp::apply, vector component "<<k<<" is not compatiable with GridSpace/MPIGridSpace\n";
	    throw e;
	  } // try-dyn-cast ends
	} // component branch ends
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSOpt::GridDiffOp::apply\n";
	throw e;
      }
    }

    using LinearOp<Scalar>::applyAdj;

    void applyAdj(Vector<Scalar> const & x, 
		  Vector<Scalar> & y) const{
      try{
	SpaceTest(this->getDomain(),x,"TSOpt::GridDiffOp::applyAdj (dom)");
	SpaceTest(this->getRange(),y,"TSOpt::GridDiffOp::applyAdj (rng)");

	Components<Scalar> cx(x);
	Components<Scalar> cy(y);
	int csize = cx.getSize();
	
	for (int k=0; k< csize; k++) {
	  // extract grid for current component
	  try{
#ifdef IWAVE_USE_MPI
	    MPIGridSpace<Scalar> const & xsp = 
	      dynamic_cast<MPIGridSpace<Scalar> const &>(cx[k].getSpace());
#else
	    GridSpace<Scalar> const & xsp = 
	      dynamic_cast<GridSpace<Scalar> const &>(cx[k].getSpace());	  
#endif
	    Grid<Scalar> const & gref = xsp.getGrid();

	    if (diffaxisid == gref.getExdAxis().id) {
	      if(gref.getExdAxis().n < 2) {
		RVLException e;
		e<<"Error: GridDiffOp::applyOp \n";
		e<<"the number of grid points along differetiating direction is less than 2 \n";
		throw e;
	      }
	    }
	    else {
	      if(gref.getAxis(diffaxisid-1).n < 2) {
		RVLException e;
		e<<"Error: GridDiffOp::applyOp \n";
		e<<"the number of grid points along differetiating direction is less than 2 \n";
		throw e;
	      }
	    }      

	    // create working buffer
	    ContentPackage<Scalar, Grid<Scalar> > cpbuf;
	    Grid<Scalar> sg(gref);
	    /*
	      for(int i = RARR_MAX_NDIM; i > 1; i--){
	      if (i != diffaxisid) sg.getData()[i-1].n = 1;	  
	      }
	    */
	    int sgeid = sg.getExdAxis().id;
	    // if (sgeid != diffaxisid) sg.getExdAxis().n = 1;
	    cpbuf.initialize(sg);
	    reset_ioff(gref);
	    
	    do {
	      // record current origin position
	      int itmp;
	      if (sgeid == diffaxisid) itmp = 0;
	      else itmp = exdoff;
	      sg.getExdAxis().o =
		itmp * gref.getExdAxis().d + gref.getExdAxis().o;
	      
	      for (int i=1;i<gref.getDim();i++){
		if (i == diffaxisid-1) itmp = 0;
		else itmp = ioff[i];
		sg.getData()[i].o =
		  itmp * gref.getData()[i].d + gref.getData()[i].o;
	      }
	    
 	      GridReader<Scalar> r(cpbuf);
	      cx[k].eval(r);
	      DTfun(cpbuf);
	      GridWriter<Scalar> w(cpbuf);
	      cy[k].eval(w);
	    }
	    while(update_ioff());
	  }
	  catch (std::bad_cast &bd){
	    RVLException e;
	    e<<bd.what()<<"\n";
	    e<< "GridDiffOp::applyAdj, vector component "<<k<<" is not compatiable with GridSpace/MPIGridSpace\n";
	    throw e;
	  } // try-dyn-cast ends
	} // component branch ends
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSOpt::GridDiffOp::applyAdj\n";
	throw e;
      }
    }

    LinearOp<Scalar> * clone() const {
      return new GridDiffOp(*this);
    }
    
  public:
    
    GridDiffOp(Space<Scalar> const & _dom
	       , int _difid
	       , int _capacity = 536870912):
      dom(_dom), diffaxisid(_difid), buf_capacity(_capacity) {
      IASN(ioff,IPNT_0);
    }

    GridDiffOp(GridDiffOp<Scalar> const & dop):
      dom(dop.dom), diffaxisid(dop.diffaxisid), buf_capacity(dop.buf_capacity) {
      IASN(ioff,dop.ioff);
      exdoff = dop.exdoff;      
    }

    ~GridDiffOp(){ }
    
    Space<Scalar> const & getDomain() const { return dom; }
    Space<Scalar> const & getRange() const { return dom; }   
    ostream & write(ostream & str) const {
      str<<"GridDiffOp\n";
      return str;
    }
  };

}




#endif
