#include "gridops.hh"

/* helmholtz power function */
extern "C" void helm_(int,         /* bc */
                      integer *,   /* n1 */
                      integer *,   /* n2 */
                      float *,     /* d1 */
                      float *,     /* d2 */
                      float *,     /* w1 */
                      float *,     /* w2 */
                      float *,     /* p */
                      float *,     /* datum */
                      float *,     /* data in */
                      float *,     /* data out */
                      float *,     /* workspace */
                      integer *,   /* length of workspace */
                      integer *    /* error flag */
		      );

namespace TSOpt {

  using RVL::ScalarFieldTraits;
  using RVL::SpaceTest;
  using RVL::Operator;
  using RVL::LinearOp;
  using RVL::Space;
  using RVL::ProductSpace;
  using RVL::Vector;
  using RVL::Components;
  using RVL::ProtectedDivision;
  using RVL::RnArray;
  using RVL::RVLScale;
  using RVL::BinaryLocalFunctionObject;
  using RVL::RVLException;
  using RVL::ContentPackage;
  using RVL::LocalDataContainer;
  using RVL::MPISerialFunctionObject;
  using RVL::MPISerialFunctionObjectRedn;

  void GridMaskFO::operator()(LocalDataContainer<ireal> & x,
			      LocalDataContainer<ireal> const & y) {
    try {
      // cerr<<"GridWindowFO::operator() begin\n";
      ContentPackage< ireal, RARR > const & gy =
	dynamic_cast<ContentPackage< ireal, RARR > const &>(y);
            
      ContentPackage< ireal, RARR > & gx =
	dynamic_cast<ContentPackage< ireal, RARR > &>(x);
            
      // precondition - metadata are same dimn
      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();
      int dimx; int dimy;
      ra_ndim(&rax,&dimx);
      ra_ndim(&ray,&dimy);
      if (dimx != dimy) {
	RVLException e;
	e<<"Error: GridMaskFO::operator()\n";
	e<<"arguments have different dims:\n";
	e<<"dimx="<<dimx<<" dimy="<<dimy<<"\n";
	throw e;
      }
            
      // compute grid params
      IPNT gsx; IPNT gex;
      IPNT gsy; IPNT gey;
      IPNT s; IPNT e;
      ra_a_gse(&rax,gsx,gex);
      ra_a_gse(&ray,gsy,gey);
      // calculate grid overlap
      for (int ii=0;ii<dimx;ii++)  {
	s[ii]=max(gsy[ii],gsx[ii]);
	e[ii]=min(gey[ii],gex[ii]);
      }
            
      IPNT i;
      RPNT fac;
      RASN(fac,RPNT_1);
#if RARR_MAX_NDIM > 0
      if (dimx==1) {
#pragma ivdep
	for (i[0]=siw[0];i[0]<=e[0]-eiw[0];i[0]++) {
	  if (bias) {
	    rax._s1[i[0]]+=ray._s1[i[0]];
	  }
	  else{
	    rax._s1[i[0]]=ray._s1[i[0]];
	  }
	}
      }
#endif
#if RARR_MAX_NDIM > 1
      if (dimx==2) {
	for (i[1]=siw[1];i[1]<=e[1]-eiw[1];i[1]++) {
#pragma ivdep
	  for (i[0]=siw[0];i[0]<=e[0]-eiw[0];i[0]++) {
	    if (bias) {
	      rax._s2[i[1]][i[0]]+=ray._s2[i[1]][i[0]];
	    }
	    else{
	      rax._s2[i[1]][i[0]]=ray._s2[i[1]][i[0]];
	    }
	  }
	}
      }
#endif
#if RARR_MAX_NDIM > 2
      if (dimx==3) {
	for (i[2]=siw[2];i[2]<=e[2]-eiw[2];i[2]++) {
	  for (i[1]=siw[1];i[1]<=e[1]-eiw[1];i[1]++) {
#pragma ivdep
	    for (i[0]=siw[0];i[0]<=e[0]-eiw[0];i[0]++) {
	      if (bias) {
		rax._s3[i[2]][i[1]][i[0]]+=ray._s3[i[2]][i[1]][i[0]];
	      }
	      else{
		rax._s3[i[2]][i[1]][i[0]]=ray._s3[i[2]][i[1]][i[0]];
	      }
	    }
	  }
	}
      }
#endif
      if (dimx<1 || dimx>3) {
	RVLException e;
	e<<"Error: GridMaskFO::operator()\n";
	e<<"dim = "<<dimx<<" outside of admissible set {1, 2, 3}\n";
	throw e;
      }
      // cerr<<"GridWindowFO::operator() end\n";
    }
    catch (bad_cast) {
      RVLException e;
      e<<"\nError: GridMaskFO::operator()\n";
      e<<"at least one arg is not ContentPackage<ireal,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridMaskFO::operator()\n";
      throw e;
    }
  }
    
   
#ifdef IWAVE_USE_MPI
  typedef MPIGridSpace myGridSpace;
#else
  typedef GridSpace myGridSpace;
#endif
    
  GridMaskOp::GridMaskOp(Space<ireal> const & _dom,
			 Vector<ireal> const & _bg,
			 RPNT const & sw, RPNT const & ew)
    : dom(_dom), bg(_bg) {
    try {
      // generic initialization of iw
      IASN(siw,IPNT_0);
      IASN(eiw,IPNT_0);
            
      // branch on product structure - unfortunately required
      ProductSpace<ireal> const * pdom = dynamic_cast<ProductSpace<ireal> const *>(&dom);
      ProductSpace<ireal> const * prng = dynamic_cast<ProductSpace<ireal> const *>(&(bg.getSpace()));
      if (pdom && prng) {
	if (pdom->getSize() != prng->getSize()) {
	  RVLException e;
	  e<<"Error GridMaskOp constructor\n";
	  e<<"  domain and range are product spaces with different numbers of components\n";
	  throw e;
	}
	// check that component gridspaces are pairwise compatible and compatible
	// with 0th domain component - then they are all compatible
                
	myGridSpace const & gref = dynamic_cast<myGridSpace const &>((*pdom)[0]);
	for (int j=0; j<(int)pdom->getSize(); j++) {
	  myGridSpace const & gdom = dynamic_cast<myGridSpace const &>((*pdom)[j]);
	  myGridSpace const & grng = dynamic_cast<myGridSpace const &>((*prng)[j]);
                    
	  if (retrieveGlobalRank()==0) {
	    if (compatible_grid(gdom.getGrid(),grng.getGrid()) ||
		compatible_grid(gref.getGrid(),grng.getGrid())) {
	      RVLException e;
	      e<<"Error: GridWindowOp constructor\n";
	      e<<"  domain, range defined on incompatible grids\n";
	      e<<"  product case, component = "<<j<<"\n";
	      e<<"  domain:\n";
	      for (int i=0;i<gdom.getGrid().gdim;i++)
		e<<"    axis "<<i<<" d="<<gdom.getGrid().axes[i].d<<" o="<<gdom.getGrid().axes[i].o<<"\n";
	      e<<"  range:\n";
	      for (int i=0;i<grng.getGrid().gdim;i++)
		e<<"    axis "<<i<<" d="<<grng.getGrid().axes[i].d<<" o="<<grng.getGrid().axes[i].o<<"\n";
	      throw e;
	    }
	  }
	}
	if (retrieveGlobalRank()==0) {
	  grid const & g = gref.getGrid();
	  for (int i=0; i< g.dim; i++) {
	    siw[i]=(int) (sw[i]/(g.axes[i].d) + 0.1);
	    eiw[i]=(int) (ew[i]/(g.axes[i].d) + 0.1);

	    siw[i]=iwave_max(siw[i],1);
	    eiw[i]=iwave_max(eiw[i],1);
	  }
	}
      }
      else {
	myGridSpace const & gdom = dynamic_cast<myGridSpace const &> (dom);
	myGridSpace const & grng = dynamic_cast<myGridSpace const &>(bg.getSpace());
	if (retrieveGlobalRank()==0) {
	  if (compatible_grid(gdom.getGrid(),grng.getGrid())) {
	    RVLException e;
	    e<<"Error: GridMaskOp constructor\n";
	    e<<"  domain, range defined on incompatible grids\n";
	    e<<"  domain:\n";
	    for (int i=0;i<gdom.getGrid().gdim;i++) 
	      e<<"    axis "<<i<<" d="<<gdom.getGrid().axes[i].d<<" o="<<gdom.getGrid().axes[i].o<<"\n";
	    e<<"  range:\n";
	    for (int i=0;i<grng.getGrid().gdim;i++) 
	      e<<"    axis "<<i<<" d="<<grng.getGrid().axes[i].d<<" o="<<grng.getGrid().axes[i].o<<"\n";
	    throw e;
	  }
	  grid const & g = gdom.getGrid();
	  for (int i=0; i< g.dim; i++) {
	    siw[i]=(int) (sw[i]/(g.axes[i].d) + 0.1);
	    eiw[i]=(int) (ew[i]/(g.axes[i].d) + 0.1);

	    siw[i]=iwave_max(siw[i],0);
	    eiw[i]=iwave_max(eiw[i],0);
	  }
	}
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridWindowOp constructor\n";
      e<<"  either domain or range is neither product nor a GridSpace,\n";
      e<<"  or some component is not a GridSpace\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridWindowOp constructor\n";
      throw e;
    }
  }
    
  void GridMaskOp::apply(Vector<ireal> const & x,
			 Vector<ireal> & y) const {
    try {
      GridMaskFO op(siw,eiw,true);
      MPISerialFunctionObject<ireal> mpiop(op);
      y.copy(bg);
      y.eval(mpiop,x);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridMaskOp::apply\n";
      throw e;
    }
  }
    
  void GridMaskOp::applyDeriv(Vector<ireal> const & x,
			      Vector<ireal> const & dx,
			      Vector<ireal> & dy) const {
    try {
      GridMaskFO op(siw,eiw,false);
      MPISerialFunctionObject<ireal> mpiop(op);
      dy.zero();
      dy.eval(mpiop,dx);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridMaskOp::applyDeriv\n";
      throw e;
    }
  }
    
  void GridMaskOp::applyAdjDeriv(Vector<ireal> const & x,
				 Vector<ireal> const & dy,
				 Vector<ireal> & dx) const {
    try {
      GridMaskFO op(siw,eiw,false);
      MPISerialFunctionObject<ireal> mpiop(op);
      dx.zero();
      dx.eval(mpiop,dy);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridMaskOp::applyAdjDeriv\n";
      throw e;
    }
        
  }
  ostream & GridMaskOp::write(ostream & str) const {
    if (!retrieveGlobalRank()) {
      str<<"GridMaskOp\n";
    }
    return str;
  }
    
  void GridWindowFO::operator()(LocalDataContainer<ireal> & x,
				LocalDataContainer<ireal> const & y) {
    try {
      // cerr<<"GridWindowFO::operator() begin\n";
      ContentPackage< ireal, RARR > const & gy = 
	dynamic_cast<ContentPackage< ireal, RARR > const &>(y);

      ContentPackage< ireal, RARR > & gx = 
	dynamic_cast<ContentPackage< ireal, RARR > &>(x);
      
      // precondition - metadata are same dimn
      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();
      int dimx; int dimy;
      ra_ndim(&rax,&dimx);
      ra_ndim(&ray,&dimy);
      if (dimx != dimy) {
	RVLException e;
	e<<"Error: GridWindow::operator()\n";
	e<<"arguments have different dims:\n";
	e<<"dimx="<<dimx<<" dimy="<<dimy<<"\n";
	throw e;
      }

      // compute grid params
      IPNT gsx; IPNT gex;
      IPNT gsy; IPNT gey;
      IPNT s; IPNT e; 
      ra_a_gse(&rax,gsx,gex);
      ra_a_gse(&ray,gsy,gey);
      // calculate grid overlap
      for (int ii=0;ii<dimx;ii++)  {
	s[ii]=max(gsy[ii],gsx[ii]);
	e[ii]=min(gey[ii],gex[ii]);
      }	

      IPNT i;
      RPNT fac;
      RASN(fac,RPNT_1);
#if RARR_MAX_NDIM > 0
      if (dimx==1) {
#pragma ivdep
	for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	  fac[0] = iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(ireal(i[0]-s[0]+1))/ireal(iw[0]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(ireal(e[0]+1-i[0]))/ireal(iw[0]))));
	  if (bias) {
	    rax._s1[i[0]]+=fac[0]*ray._s1[i[0]];
	  }
	  else {
	    rax._s1[i[0]] =fac[0]*ray._s1[i[0]];
	  }
	}
      }
#endif
#if RARR_MAX_NDIM > 1 
      if (dimx==2) {
	for (i[1]=s[1];i[1]<=e[1];i[1]++) {
	  fac[1] = iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(ireal(i[1]-s[1]+1))/ireal(iw[1]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(ireal(e[1]+1-i[1]))/ireal(iw[1]))));
#pragma ivdep
	  for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	    fac[0] = fac[1]*iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(ireal(i[0]-s[0]+1))/ireal(iw[0]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(ireal(e[0]+1-i[0]))/ireal(iw[0]))));
	    if (bias) {
	      rax._s2[i[1]][i[0]]+=fac[0]*ray._s2[i[1]][i[0]];	  
	    }
	    else {
	      rax._s2[i[1]][i[0]] =fac[0]*ray._s2[i[1]][i[0]];	  
	    }		
	  }
	}
      }
#endif
#if RARR_MAX_NDIM > 2 
      if (dimx==3) {
	for (i[2]=s[2];i[2]<=e[2];i[2]++) {
	  fac[2] = iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(ireal(i[2]-s[2]+1))/ireal(iw[2]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(ireal(e[2]+1-i[2]))/ireal(iw[2]))));
	  for (i[1]=s[1];i[1]<=e[1];i[1]++) {
	    fac[1] = fac[2]*iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(ireal(i[1]-s[1]+1))/ireal(iw[1]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(ireal(e[1]+1-i[1]))/ireal(iw[1]))));
#pragma ivdep
	    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	      fac[0] = fac[1]*iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(ireal(i[0]-s[0]+1))/ireal(iw[0]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(ireal(e[0]+1-i[0]))/ireal(iw[0]))));
	      if (bias) {
		rax._s3[i[2]][i[1]][i[0]]+=fac[0]*ray._s3[i[2]][i[1]][i[0]];	  
	      }
	      else {
		rax._s3[i[2]][i[1]][i[0]] =fac[0]*ray._s3[i[2]][i[1]][i[0]];	  
	      }
	    }
	  }
	}
      }
#endif
      if (dimx<1 || dimx>3) {
	RVLException e;
	e<<"Error: GridWindowFO::operator()\n";
	e<<"dim = "<<dimx<<" outside of admissible set {1, 2, 3}\n";
	throw e;
      }
      // cerr<<"GridWindowFO::operator() end\n";
    }
    catch (bad_cast) {
      RVLException e;
      e<<"\nError: GridWindowFO::operator()\n";
      e<<"at least one arg is not ContentPackage<ireal,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridWindowFO::operator()\n";
      throw e;
    }
  }


  GridWindowOp::GridWindowOp(Space<ireal> const & _dom,
			     Vector<ireal> const & _bg,
			     RPNT const & w) 
    : dom(_dom), bg(_bg) {
    try {
      // generic initialization of iw
      IASN(iw,IPNT_0);

      // branch on product structure - unfortunately required
      ProductSpace<ireal> const * pdom = dynamic_cast<ProductSpace<ireal> const *>(&dom);
      ProductSpace<ireal> const * prng = dynamic_cast<ProductSpace<ireal> const *>(&(bg.getSpace()));
      if (pdom && prng) {
	if (pdom->getSize() != prng->getSize()) {
	  RVLException e;
	  e<<"Error GridWindowOp constructor\n";
	  e<<"  domain and range are product spaces with different numbers of components\n";
	  throw e;
	}
	// check that component gridspaces are pairwise compatible and compatible
	// with 0th domain component - then they are all compatible

	myGridSpace const & gref = dynamic_cast<myGridSpace const &>((*pdom)[0]);
	for (int j=0; j<(int)pdom->getSize(); j++) {
	  myGridSpace const & gdom = dynamic_cast<myGridSpace const &>((*pdom)[j]);
	  myGridSpace const & grng = dynamic_cast<myGridSpace const &>((*prng)[j]);
	  
	  if (retrieveGlobalRank()==0) {
	    if (compatible_grid(gdom.getGrid(),grng.getGrid()) ||
		compatible_grid(gref.getGrid(),grng.getGrid())) {
	      RVLException e;
	      e<<"Error: GridWindowOp constructor\n";
	      e<<"  domain, range defined on incompatible grids\n";
	      e<<"  product case, component = "<<j<<"\n";
	      e<<"  domain:\n";
	      for (int i=0;i<gdom.getGrid().gdim;i++) 
		e<<"    axis "<<i<<" d="<<gdom.getGrid().axes[i].d<<" o="<<gdom.getGrid().axes[i].o<<"\n";
	      e<<"  range:\n";
	      for (int i=0;i<grng.getGrid().gdim;i++) 
		e<<"    axis "<<i<<" d="<<grng.getGrid().axes[i].d<<" o="<<grng.getGrid().axes[i].o<<"\n";
	      throw e;
	    }
	  }
	}
	if (retrieveGlobalRank()==0) {
	  grid const & g = gref.getGrid();
	  for (int i=0; i< g.dim; i++) {
	    iw[i]=(int) (w[i]/(g.axes[i].d) + 0.1);
	    iw[i]=iwave_max(iw[i],1);
	  }
	}
      }
      else {
	myGridSpace const & gdom = dynamic_cast<myGridSpace const &> (dom);
	myGridSpace const & grng = dynamic_cast<myGridSpace const &>(bg.getSpace());
	if (retrieveGlobalRank()==0) {
	  if (compatible_grid(gdom.getGrid(),grng.getGrid())) {
	    RVLException e;
	    e<<"Error: GridWindowOp constructor\n";
	    e<<"  domain, range defined on incompatible grids\n";
	    e<<"  domain:\n";
	    for (int i=0;i<gdom.getGrid().gdim;i++) 
	      e<<"    axis "<<i<<" d="<<gdom.getGrid().axes[i].d<<" o="<<gdom.getGrid().axes[i].o<<"\n";
	    e<<"  range:\n";
	    for (int i=0;i<grng.getGrid().gdim;i++) 
	      e<<"    axis "<<i<<" d="<<grng.getGrid().axes[i].d<<" o="<<grng.getGrid().axes[i].o<<"\n";
	    throw e;
	  }
	  grid const & g = gdom.getGrid();
	  for (int i=0; i< g.dim; i++) {
	    iw[i]=(int) (w[i]/(g.axes[i].d) + 0.1);
	    iw[i]=iwave_max(iw[i],0);
	  }
	}
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridWindowOp constructor\n";
      e<<"  either domain or range is neither product nor a GridSpace,\n";
      e<<"  or some component is not a GridSpace\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridWindowOp constructor\n";
      throw e;
    }
  }
    
  void GridWindowOp::apply(Vector<ireal> const & x,
			   Vector<ireal> & y) const {
    try {
      GridWindowFO op(iw,true);
      MPISerialFunctionObject<ireal> mpiop(op);
      y.copy(bg);
      y.eval(mpiop,x);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridWindowOp::apply\n";
      throw e;
    }
  }

  void GridWindowOp::applyDeriv(Vector<ireal> const & x,
				Vector<ireal> const & dx,
				Vector<ireal> & dy) const {
    try {
      GridWindowFO op(iw,false);
      MPISerialFunctionObject<ireal> mpiop(op);
      dy.zero();
      dy.eval(mpiop,dx);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridWindowOp::applyDeriv\n";
      throw e;
    }
  }

  void GridWindowOp::applyAdjDeriv(Vector<ireal> const & x,
				   Vector<ireal> const & dy,
				   Vector<ireal> & dx) const {
    try {
      GridWindowFO op(iw,false);
      MPISerialFunctionObject<ireal> mpiop(op);    
      dx.zero();
      dx.eval(mpiop,dy);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridWindowOp::applyAdjDeriv\n";
      throw e;
    }

  }

  ostream & GridWindowOp::write(ostream & str) const {
    if (!retrieveGlobalRank()) {
      str<<"GridWindowOp\n";
    }
    return str;
  }

  void GridFwdDerivFO::operator()(LocalDataContainer<ireal> & x,
				  LocalDataContainer<ireal> const & y) {
    try {
      // cerr<<"GridFwdDerivFO::operator() begin\n";
      ContentPackage< ireal, RARR > const & gy = 
	dynamic_cast<ContentPackage< ireal, RARR > const &>(y);

      ContentPackage< ireal, RARR > & gx = 
	dynamic_cast<ContentPackage< ireal, RARR > &>(x);

      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();	  

      // sanity test 
      if (ra_compare_meta(&rax,&ray)) {
	RVLException e;
	e<<"Error: GridFwdDerivFO::operator()\n";
	e<<"  incompatible input, output LDC\n";
	e<<"  input RARR:\n";
	throw e;
      }

      IPNT s0, e0, i;
      ra_a_gse(&rax,s0,e0);

#if RARR_MAX_NDIM > 0

      if (dir==0 && rax.ndim==1) {
	rax._s1[e0[0]] = REAL_ZERO;
#pragma ivdep
	for (i[0]=s0[0];i[0]<e0[0];i[0]++) 
	  rax._s1[i[0]] = (ray._s1[i[0]+1]-ray._s1[i[0]])*fac;
      }

#if RARR_MAX_NDIM > 1

      else if (dir==0 && rax.ndim==2) {
	for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	  rax._s2[i[1]][e0[0]] = REAL_ZERO;
#pragma ivdep
	  for (i[0]=s0[0];i[0]<e0[0];i[0]++) 
	    rax._s2[i[1]][i[0]] = (ray._s2[i[1]][i[0]+1]-ray._s2[i[1]][i[0]])*fac;
	}
      }

      else if (dir==1 && rax.ndim==2) {
	for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	  rax._s2[e0[1]][i[0]] = REAL_ZERO;
	  for (i[1]=s0[1];i[1]<e0[1];i[1]++) {
	    rax._s2[i[1]][i[0]] = (ray._s2[i[1]+1][i[0]]-ray._s2[i[1]][i[0]])*fac;
	  }
	}
      }

#if RARR_MAX_NDIM > 2

      else if (dir==0 && rax.ndim==3) {
	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    rax._s3[i[2]][i[1]][e0[0]] = REAL_ZERO;
#pragma ivdep
	    for (i[0]=s0[0];i[0]<e0[0];i[0]++) 
	      rax._s3[i[2]][i[1]][i[0]] = (ray._s3[i[2]][i[1]][i[0]+1]-ray._s3[i[2]][i[1]][i[0]])*fac;
	  }
	}
      }
      else if (dir==1 && rax.ndim==3) {
	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
	    rax._s3[i[2]][e0[1]][i[0]] = REAL_ZERO;
#pragma ivdep
	    for (i[1]=s0[1];i[1]<e0[1];i[1]++) {
	      rax._s3[i[2]][i[1]][i[0]] = (ray._s3[i[2]][i[1]+1][i[0]]-ray._s3[i[2]][i[1]][i[0]])*fac;
	    }	
	  }
	}
      }

      else if (dir==2 && rax.ndim==3) {
	for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	    rax._s3[e0[2]][i[1]][i[0]] = REAL_ZERO;
#pragma ivdep
	    for (i[2]=s0[2];i[2]<e0[2];i[2]++) {
	      rax._s3[i[2]][i[1]][i[0]] = (ray._s3[i[2]+1][i[1]][i[0]]-ray._s3[i[2]][i[1]][i[0]])*fac;
	    }
	  }
	}
      }

#if RARR_MAX_NDIM > 3
	
      else if (dir==0 && rax.ndim==4) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      rax._s4[i[3]][i[2]][i[1]][e0[0]] = REAL_ZERO;
#pragma ivdep
	      for (i[0]=s0[0];i[0]<e0[0];i[0]++) 
		rax._s4[i[3]][i[2]][i[1]][i[0]] = (ray._s4[i[3]][i[2]][i[1]][i[0]+1]-
						   ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;
	    }
	  }
	}
      }
      else if (dir==1 && rax.ndim==4) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
	      rax._s4[i[3]][i[2]][e0[1]][i[0]] = REAL_ZERO;
#pragma ivdep
	      for (i[1]=s0[1];i[1]<e0[1];i[1]++) {
		rax._s4[i[3]][i[2]][i[1]][i[0]] = (ray._s4[i[3]][i[2]][i[1]+1][i[0]]-
						   ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;
	      }	
	    }
	  }
	}
      }
      else if (dir==2 && rax.ndim==4) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
	      rax._s4[i[3]][e0[2]][i[1]][i[0]] = REAL_ZERO;
#pragma ivdep
	      for (i[2]=s0[2];i[2]<e0[2];i[2]++) {
		rax._s4[i[3]][i[2]][i[1]][i[0]] = (ray._s4[i[3]][i[2]+1][i[1]][i[0]]-
						   ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;
	      }	
	    }
	  }
	}
      }

      else if (dir==3 && rax.ndim==4) {
	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
	      rax._s4[e0[3]][i[2]][i[1]][i[0]] = REAL_ZERO;
#pragma ivdep
	      for (i[3]=s0[3];i[3]<e0[3];i[3]++) {
		rax._s4[i[3]][i[2]][i[1]][i[0]] = (ray._s4[i[3]+1][i[2]][i[1]][i[0]]-
						   ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;
	      }	
	    }
	  }
	}
      }

#if RARR_MAX_NDIM > 4

      else if (dir==0 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
		rax._s5[i[4]][i[3]][i[2]][i[1]][e0[0]] = REAL_ZERO;
#pragma ivdep
		for (i[0]=s0[0];i[0]<e0[0];i[0]++) 
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = (ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]+1]-
							   ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;
	      }
	    }
	  }
	}
      }
      else if (dir==1 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
		rax._s5[i[4]][i[3]][i[2]][e0[1]][i[0]] = REAL_ZERO;
#pragma ivdep
		for (i[1]=s0[1];i[1]<e0[1];i[1]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = (ray._s5[i[4]][i[3]][i[2]][i[1]+1][i[0]]-
							   ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;
		}	
	      }
	    }
	  }
	}
      }
      else if (dir==2 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
		rax._s5[i[4]][i[3]][e0[2]][i[1]][i[0]] = REAL_ZERO;
#pragma ivdep
		for (i[2]=s0[2];i[2]<e0[2];i[2]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = (ray._s5[i[4]][i[3]][i[2]+1][i[1]][i[0]]-
							   ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;
		}	
	      }
	    }
	  }
	}
      }
      else if (dir==3 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
		rax._s5[i[4]][e0[3]][i[2]][i[1]][i[0]] = REAL_ZERO;
#pragma ivdep
		for (i[3]=s0[3];i[3]<e0[3];i[3]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = (ray._s5[i[4]][i[3]+1][i[2]][i[1]][i[0]]-
							   ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;
		}	
	      }
	    }
	  }
	}
      }

      else if (dir==4 && rax.ndim==5) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
		rax._s5[e0[4]][i[3]][i[2]][i[1]][i[0]] = REAL_ZERO;
#pragma ivdep
		for (i[4]=s0[4];i[4]<e0[4];i[4]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = (ray._s5[i[4]+1][i[3]][i[2]][i[1]][i[0]]-
							   ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;
		}	
	      }
	    }
	  }
	}
      }

#endif // RARR_MAX_NDIM > 4
#endif // > 3
#endif // > 2
#endif // > 1

      else {
	RVLException e;
	e<<"Error: GridFwdDerivFO::operator()\n";
	e<<"  attempt to apply divided diff in direction "<<dir<<" on array of dim "<<rax.ndim<<"\n";
	throw e;
      }
#endif // > 0
      // cerr<<"GridFwdDerivFO::operator() end\n";
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridFwdDerivFO::operator()\n";
      e<<"input type error - not CP<ireal,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridFwdDerivFO::operator()\n";
      throw e;
    }
  }

  void GridAdjDerivFO::operator()(LocalDataContainer<ireal> & x,
				  LocalDataContainer<ireal> const & y) {
    try {
      // cerr<<"GridAdjDerivFO::operator() begin\n";
      ContentPackage< ireal, RARR > const & gy = 
	dynamic_cast<ContentPackage< ireal, RARR > const &>(y);

      ContentPackage< ireal, RARR > & gx = 
	dynamic_cast<ContentPackage< ireal, RARR > &>(x);

      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();	  

      // sanity test 
      if (ra_compare_meta(&rax,&ray)) {
	RVLException e;
	e<<"Error: GridAdjDerivFO::operator()\n";
	e<<"  incompatible input, output LDC\n";
	throw e;
      }

      IPNT s0, e0, i;
      ra_a_gse(&rax,s0,e0);
      
#if RARR_MAX_NDIM > 0

      if (dir==0 && rax.ndim==1) {
	rax._s1[s0[0]] = -ray._s1[s0[0]]*fac;
	rax._s1[e0[0]] = ray._s1[e0[0]-1]*fac;
#pragma ivdep
	for (i[0]=s0[0]+1;i[0]<e0[0];i[0]++) 
	  rax._s1[i[0]]=(ray._s1[i[0]-1]-ray._s1[i[0]])*fac;      }

#if RARR_MAX_NDIM > 1

      else if (dir==0 && rax.ndim==2) {
	for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	  rax._s2[i[1]][s0[0]] = -ray._s2[i[1]][s0[0]]*fac;
	  rax._s2[i[1]][e0[0]] = ray._s2[i[1]][e0[0]-1]*fac;
#pragma ivdep
	  for (i[0]=s0[0]+1;i[0]<e0[0];i[0]++) {
	    rax._s2[i[1]][i[0]]=(ray._s2[i[1]][i[0]-1]-ray._s2[i[1]][i[0]])*fac;      
	  }
	}
      }
      
      else if (dir==1 && rax.ndim==2) {
	for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	  rax._s2[s0[1]][i[0]] = -ray._s2[s0[1]][i[0]]*fac;
	  rax._s2[e0[1]][i[0]] = ray._s2[e0[1]-1][i[0]]*fac;
	  for (i[1]=s0[1]+1;i[1]<e0[1];i[1]++) {
	    rax._s2[i[1]][i[0]]=(ray._s2[i[1]-1][i[0]]-ray._s2[i[1]][i[0]])*fac;      
	  }
	}
      }

#if RARR_MAX_NDIM > 2

      else if (dir==0 && rax.ndim==3) {
	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    rax._s3[i[2]][i[1]][s0[0]] = -ray._s3[i[2]][i[1]][s0[0]]*fac;
	    rax._s3[i[2]][i[1]][e0[0]] = ray._s3[i[2]][i[1]][e0[0]-1]*fac;
#pragma ivdep
	    for (i[0]=s0[0]+1;i[0]<e0[0];i[0]++) {
	      rax._s3[i[2]][i[1]][i[0]]=(ray._s3[i[2]][i[1]][i[0]-1]-ray._s3[i[2]][i[1]][i[0]])*fac;      
	    }
	  }
	}
      }

      else if (dir==1 && rax.ndim==3) {
	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	    rax._s3[i[2]][s0[1]][i[0]] = -ray._s3[i[2]][s0[1]][i[0]]*fac;
	    rax._s3[i[2]][e0[1]][i[0]] = ray._s3[i[2]][e0[1]-1][i[0]]*fac;
#pragma ivdep
	    for (i[1]=s0[1]+1;i[1]<e0[1];i[1]++) {
	      rax._s3[i[2]][i[1]][i[0]]=(ray._s3[i[2]][i[1]-1][i[0]]-ray._s3[i[2]][i[1]][i[0]])*fac;      
	    }
	  }
	}
      }

      else if (dir==2 && rax.ndim == 3) {
	for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	    rax._s3[s0[2]][i[1]][i[0]] = -ray._s3[s0[2]][i[1]][i[0]]*fac;
	    rax._s3[e0[2]][i[1]][i[0]] = ray._s3[e0[2]-1][i[1]][i[0]]*fac;
#pragma ivdep
	    for (i[2]=s0[2]+1;i[2]<e0[2];i[2]++) {
	      rax._s3[i[2]][i[1]][i[0]]=(ray._s3[i[2]-1][i[1]][i[0]]-ray._s3[i[2]][i[1]][i[0]])*fac;      
	    }
	  }
	}
      }

#if RARR_MAX_NDIM > 3
	
      else if (dir==0 && rax.ndim==4) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      rax._s4[i[3]][i[2]][i[1]][s0[0]] = -ray._s4[i[3]][i[2]][i[1]][s0[0]]*fac;
	      rax._s4[i[3]][i[2]][i[1]][e0[0]] = ray._s4[i[3]][i[2]][i[1]][e0[0]-1]*fac;
#pragma ivdep
	      for (i[0]=s0[0]+1;i[0]<e0[0];i[0]++) {
		rax._s4[i[3]][i[2]][i[1]][i[0]]=(ray._s4[i[3]][i[2]][i[1]][i[0]-1]-
						 ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;      
	      }
	    }
	  }
	}
      }
      else if (dir==1 && rax.ndim==4) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s4[i[3]][i[2]][s0[1]][i[0]] = -ray._s4[i[3]][i[2]][s0[1]][i[0]]*fac;
	      rax._s4[i[3]][i[2]][e0[1]][i[0]] = ray._s4[i[3]][i[2]][e0[1]-1][i[0]]*fac;
#pragma ivdep
	      for (i[1]=s0[1]+1;i[1]<e0[1];i[1]++) {
		rax._s4[i[3]][i[2]][i[1]][i[0]]=(ray._s4[i[3]][i[2]][i[1]-1][i[0]]-
						 ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;      
	      }
	    }
	  }
	}
      }
      else if (dir==2 && rax.ndim==4) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s4[i[3]][s0[2]][i[1]][i[0]] = -ray._s4[i[3]][s0[2]][i[1]][i[0]]*fac;
	      rax._s4[i[3]][e0[2]][i[1]][i[0]] = ray._s4[i[3]][e0[2]-1][i[1]][i[0]]*fac;
#pragma ivdep
	      for (i[2]=s0[2]+1;i[2]<e0[2];i[2]++) {
		rax._s4[i[3]][i[2]][i[1]][i[0]]=(ray._s4[i[3]][i[2]-1][i[1]][i[0]]-
						 ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;      
	      }
	    }
	  }
	}
      }

      else if (dir==3 && rax.ndim==4) {
	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s4[s0[3]][i[2]][i[1]][i[0]] = -ray._s4[s0[3]][i[2]][i[1]][i[0]]*fac;
	      rax._s4[e0[3]][i[2]][i[1]][i[0]] = ray._s4[e0[3]-1][i[2]][i[1]][i[0]]*fac;
#pragma ivdep
	      for (i[3]=s0[3]+1;i[3]<e0[3];i[3]++) {
		rax._s4[i[3]][i[2]][i[1]][i[0]]=(ray._s4[i[3]-1][i[2]][i[1]][i[0]]-
						 ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;      
	      }
	    }
	  }
	}
      }

#if RARR_MAX_NDIM > 4

      else if (dir==0 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
		rax._s5[i[4]][i[3]][i[2]][i[1]][s0[0]] = -ray._s5[i[4]][i[3]][i[2]][i[1]][s0[0]]*fac;
		rax._s5[i[4]][i[3]][i[2]][i[1]][e0[0]] = ray._s5[i[4]][i[3]][i[2]][i[1]][e0[0]-1]*fac;
#pragma ivdep
		for (i[0]=s0[0]+1;i[0]<e0[0];i[0]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]]=(ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]-1]-
							 ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;      
		}
	      }
	    }
	  }
	}
      }
      else if (dir==1 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s5[i[4]][i[3]][i[2]][s0[1]][i[0]] = -ray._s5[i[4]][i[3]][i[2]][s0[1]][i[0]]*fac;
		rax._s5[i[4]][i[3]][i[2]][e0[1]][i[0]] = ray._s5[i[4]][i[3]][i[2]][e0[1]-1][i[0]]*fac;
#pragma ivdep
		for (i[1]=s0[1]+1;i[1]<e0[1];i[1]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]]=(ray._s5[i[4]][i[3]][i[2]][i[1]-1][i[0]]-
							 ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;      
		}
	      }
	    }
	  }
	}
      }
      else if (dir==2 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s5[i[4]][i[3]][s0[2]][i[1]][i[0]] = -ray._s5[i[4]][i[3]][s0[2]][i[1]][i[0]]*fac;
		rax._s5[i[4]][i[3]][e0[2]][i[1]][i[0]] = ray._s5[i[4]][i[3]][e0[2]-1][i[1]][i[0]]*fac;
#pragma ivdep
		for (i[2]=s0[2]+1;i[2]<e0[2];i[2]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]]=(ray._s5[i[4]][i[3]][i[2]-1][i[1]][i[0]]-
							 ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;      
		}
	      }
	    }
	  }
	}
      }
      else if (dir==3 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s5[i[4]][s0[3]][i[2]][i[1]][i[0]] = -ray._s5[i[4]][s0[3]][i[2]][i[1]][i[0]]*fac;
		rax._s5[i[4]][e0[3]][i[2]][i[1]][i[0]] = ray._s5[i[4]][e0[3]-1][i[2]][i[1]][i[0]]*fac;
#pragma ivdep
		for (i[3]=s0[3]+1;i[3]<e0[3];i[3]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]]=(ray._s5[i[4]][i[3]-1][i[2]][i[1]][i[0]]-
							 ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;      
		}
	      }
	    }
	  }
	}
      }

      else if (dir==4 && rax.ndim==5) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s5[s0[4]][i[3]][i[2]][i[1]][i[0]] = -ray._s5[s0[4]][i[3]][i[2]][i[1]][i[0]]*fac;
		rax._s5[e0[4]][i[3]][i[2]][i[1]][i[0]] = ray._s5[e0[4]-1][i[3]][i[2]][i[1]][i[0]]*fac;
#pragma ivdep
		for (i[4]=s0[4]+1;i[4]<e0[4];i[4]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]]=(ray._s5[i[4]-1][i[3]][i[2]][i[1]][i[0]]-
							 ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;      
		}
	      }
	    }
	  }
	}	
      }

#endif // RARR_MAX_NDIM > 4
#endif // > 3
#endif // > 2
#endif // > 1

      else {
	RVLException e;
	e<<"Error: GridAdjDerivFO::operator()\n";
	e<<"  attempt to apply divided diff in direction "<<dir<<" on array of dim "<<rax.ndim<<"\n";
	throw e;
      }
#endif // > 0
      // cerr<<"GridAdjDerivFO::operator() end\n";
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridAdjDerivFO::operator()\n";
      e<<"input type error - not CP<ireal,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridAdjDerivFO::operator()\n";
      throw e;
    }
  }

  GridDerivOp::GridDerivOp(Space<ireal> const & _dom,
			   int _dir, ireal scale)
    : dir(_dir), fac(0), dom(_dom) {
    try {
      ProductSpace<ireal> const * pdom = NULL;
      pdom = dynamic_cast<ProductSpace<ireal> const *>(&dom);
      int n_fac=1;
      if (pdom) n_fac=pdom->getSize();
      Space<ireal> const * sp = NULL;
      for (int j=0; j<n_fac; j++) {
	if (pdom) sp = &((*pdom)[j]);
	else sp = &dom;
	myGridSpace const * gdom = dynamic_cast<myGridSpace const *>(sp);
	if (!gdom) {
	  RVLException e;
	  e<<"Error: GridDerivOp constructor\n";
	  e<<"  factor "<<j<<" of input space is not a GridSpace\n";
	  e<<"  description:\n";
	  sp->write(e);
	  throw e;	  
	}
	// pure out of core: real factors only on rk=0
	if (retrieveGlobalRank()==0) {
	  if (dir < 0 || dir > gdom->getGrid().gdim-1) {
	    RVLException e;
	    e<<"Error: GridDerivOp constructor\n";
	    e<<"  direction index "<<dir<<" out of dimension range [0,"<<gdom->getGrid().gdim-1<<"\n";
	    throw e;
	  }
	  RPNT d;
	  get_d(d,gdom->getGrid());
	  fac.push_back(scale/d[dir]);
	}
	else {
	  fac.push_back(REAL_ZERO);
	}
      }
    }
    catch (RVLException e) {
      e<<"\ncalled from GridDerivOp constructor\n";
      throw e;
    }
  }

  GridDerivOp::GridDerivOp(GridDerivOp const & op)
    : dir(op.dir), fac(op.fac), dom(op.dom) {}

  void GridDerivOp::apply(Vector<ireal> const & x,
			  Vector<ireal> & y) const {
    try {
      Components<ireal> cx(x);
      Components<ireal> cy(y);
      for (int j=0;j<(int)cx.getSize();j++) {
	GridFwdDerivFO f(dir,fac[j]);
	MPISerialFunctionObject<ireal> mpif(f);
    	cy[j].eval(mpif,cx[j]);
      }
#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridDerivOp::apply\n";
      throw e;
    }
  }

  void GridDerivOp::applyAdj(Vector<ireal> const & x,
			     Vector<ireal> & y) const {
    try {
      Components<ireal> cx(x);
      Components<ireal> cy(y);
      for (int j=0;j<(int)cx.getSize();j++) {
	GridAdjDerivFO f(dir,fac[j]);
	MPISerialFunctionObject<ireal> mpif(f);
	cy[j].eval(mpif,cx[j]);
      }
#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridDerivOp::applyAdj\n";
      throw e;
    }
  }

  ostream & GridDerivOp::write(ostream & str) const {
    str<<"GridDerivOp: directional derivative, axis = "<<dir<<"\n";
    str<<"Domain:\n";
    dom.write(str);
    return str;
  }

  void GridFwdExtendFO::operator()(LocalDataContainer<ireal> & x,
				   LocalDataContainer<ireal> const & y) {
    try {
      // cerr<<"GridFwdExtendFO::operator() begin\n";
      ContentPackage< ireal, RARR > const & gy = 
	dynamic_cast<ContentPackage< ireal, RARR > const &>(y);

      ContentPackage< ireal, RARR > & gx = 
	dynamic_cast<ContentPackage< ireal, RARR > &>(x);

      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();	  

      // presumption: these are called only from GridExtendOp, which
      // does all sanity testing - could make them private data
      // exception: number of extended axes can be checked!!
      if (rax.ndim-ray.ndim != n_ext) {
	RVLException e;
	e<<"Error: GridExtendFwdFO::operator()\n"; 
	e<<"  dimension difference of input, output does not match\n";
	e<<"  number of extended axes\n";
	throw e;
      }

      IPNT s0, e0, i;
      ra_a_gse(&rax,s0,e0);

      // zero if internal
      if (!ext) {
	// sanity check - must have zero section in all internal extended axes
	for (int ii=ray.ndim;ii<rax.ndim;ii++) {
	  if (s0[ii] > 0 || e0[ii] < 0) {
	    RVLException e;
	    e<<"Error: GridExtendFwdFO::operator()\n"; 
	    e<<"  index range on axis "<<ii<<" of output extended array\n";
	    e<<"  = ["<<s0[ii]<<","<<e0[ii]<<"], does not contain zero as is\n";
	    e<<"  required for internal extended axes\n";
	    throw e;
	  }
	}
	ra_a_zero(&rax);
      }

      if (ray.ndim==1) {

#if RARR_MAX_NDIM > 1

	if (rax.ndim==2) {
	  if (ext) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s2[i[1]][i[0]] = ray._s0[i[0]];
	      }
	    }
	  }
	  else {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
#pragma ivdep
	      rax._s2[0][i[0]] = fac*ray._s0[i[0]];
	    }
	  }
	}

#if RARR_MAX_NDIM > 2

	if (rax.ndim==3) {
	  if (ext) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  rax._s3[i[2]][i[1]][i[0]] = ray._s0[i[0]];
		}
	      }
	    }
	  }
	  else {
#pragma ivdep
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s3[0][0][i[0]] = fac*ray._s0[i[0]];
	    }
	  }
	}

#if RARR_MAX_NDIM > 3

	if (rax.ndim==4) {
	  if (ext) {
	    for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	      for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		    rax._s4[i[3]][i[2]][i[1]][i[0]] = ray._s0[i[0]];
		  }
		}
	      }
	    }
	  }
	  else {
#pragma ivdep
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s4[0][0][0][i[0]] = fac*ray._s0[i[0]];
	    }
	  }
	}

#if RARR_MAX_NDIM > 4

	if (rax.ndim==5) {
	  if (ext) {
	    for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	      for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
		for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		      rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = ray._s0[i[0]];
		    }
		  }
		}
	      }
	    }
	  }
	  else {
#pragma ivdep
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s5[0][0][0][0][i[0]] = fac*ray._s0[i[0]];
	    }
	  }
	}

#endif
#endif
#endif
#endif
      }

      else if (ray.ndim==2) {

#if RARR_MAX_NDIM > 2

	if (rax.ndim==3) {
	  if (ext) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  rax._s3[i[2]][i[1]][i[0]] = ray._s2[i[1]][i[0]];
		}
	      }
	    }
	  }
	  else {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s3[0][i[1]][i[0]] = fac*ray._s2[i[1]][i[0]];
	      }
	    }
	  }
	}

#if RARR_MAX_NDIM > 3

	if (rax.ndim==4) {
	  if (ext) {
	    for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	      for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		    rax._s4[i[3]][i[2]][i[1]][i[0]] = ray._s2[i[1]][i[0]];
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s4[0][0][i[1]][i[0]] = fac*ray._s2[i[1]][i[0]];
	      }
	    }
	  }
	}

#if RARR_MAX_NDIM > 4

	if (rax.ndim==5) {
	  if (ext) {
	    for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	      for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
		for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		      rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = ray._s2[i[1]][i[0]];
		    }
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s5[0][0][0][i[1]][i[0]] = fac*ray._s2[i[1]][i[0]];
	      }
	    }
	  }
	}

#endif
#endif
#endif

      }

      else if (ray.ndim==3) {

#if RARR_MAX_NDIM > 3

	if (rax.ndim==4) {
	  if (ext) {
	    for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	      for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		    rax._s4[i[3]][i[2]][i[1]][i[0]] = ray._s3[i[2]][i[1]][i[0]];
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  rax._s4[0][i[2]][i[1]][i[0]] = fac*ray._s3[i[2]][i[1]][i[0]];
		}
	      }
	    }
	  }
	}

#if RARR_MAX_NDIM > 4

	if (rax.ndim==5) {
	  if (ext) {
	    for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	      for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
		for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		      rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = ray._s3[i[2]][i[1]][i[0]];
		    }
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  rax._s5[0][0][i[2]][i[1]][i[0]] = fac*ray._s3[i[2]][i[1]][i[0]];
		}
	      }
	    }
	  }
	}

#endif
#endif

      }

      else {
	RVLException e;
	e<<"Error: GridFwdExtendFO::operator()\n";
	e<<"  input rarr dimension not 1, 2, or 3 - only permitted spatial dims\n";
	throw e;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridFwdExtendFO::operator()\n";
      e<<"  input type error - not CP<ireal,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridFwdExtendFO::operator()\n";
      throw e;
    }
  }

  void GridAdjExtendFO::operator()(LocalDataContainer<ireal> & y,
				   LocalDataContainer<ireal> const & x) {
    try {

      ContentPackage< ireal, RARR > const & gx = 
	dynamic_cast<ContentPackage< ireal, RARR > const &>(x);

      ContentPackage< ireal, RARR > & gy = 
	dynamic_cast<ContentPackage< ireal, RARR > &>(y);

      RARR & ray = gy.getMetadata();
      RARR const & rax = gx.getMetadata();	  

      // presumption: these are called only from GridExtendOp, which
      // does all sanity testing - could make them private data
      // exception: number of extended axes can be checked!!
      if (rax.ndim-ray.ndim != n_ext) {
	RVLException e;
	e<<"Error: GridExtendAdjFO::operator()\n"; 
	e<<"  dimension difference of input, output does not match\n";
	e<<"  number of extended axes\n";
	throw e;
      }

      IPNT s0, e0, i;
      ra_a_gse(&rax,s0,e0);

      if (!ext) {
	// sanity check - must have zero section in all internal extended axes
	for (int ii=ray.ndim;ii<rax.ndim;ii++) {
	  if (s0[ii] > 0 || e0[ii] < 0) {
	    RVLException e;
	    e<<"Error: GridExtendAdjFO::operator()\n"; 
	    e<<"  index range on axis "<<ii<<" of output extended array\n";
	    e<<"  = ["<<s0[ii]<<","<<e0[ii]<<"], does not contain zero as is\n";
	    e<<"  required for internal extended axes\n";
	    throw e;
	  }
	}
      }

      // zero output anyway
      ra_a_zero(&ray);

      if (ray.ndim==1) {

#if RARR_MAX_NDIM > 1

	if (rax.ndim==2) {
	  if (ext) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		ray._s0[i[0]] += fac*rax._s2[i[1]][i[0]];
	      }
	    }
	  }
	  else {
#pragma ivdep
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      ray._s0[i[0]] = rax._s2[0][i[0]];
	    }
	  }
	}

#if RARR_MAX_NDIM > 2

	if (rax.ndim==3) {
	  if (ext) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  ray._s0[i[0]]+= fac*rax._s3[i[2]][i[1]][i[0]];
		}
	      }
	    }
	  }
	  else {
#pragma ivdep
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      ray._s0[i[0]]=rax._s3[0][0][i[0]];
	    }
	  }
	}

#if RARR_MAX_NDIM > 3

	if (rax.ndim==4) {
	  if (ext) {
	    for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	      for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		    ray._s0[i[0]] += fac* rax._s4[i[3]][i[2]][i[1]][i[0]] ;
		  }
		}
	      }
	    }
	  }
	  else {
#pragma ivdep
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      ray._s0[i[0]]=rax._s4[0][0][0][i[0]];
	    }
	  }
	}

#if RARR_MAX_NDIM > 4

	if (rax.ndim==5) {
	  if (ext) {
	    for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	      for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
		for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		      ray._s0[i[0]] += fac* rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]];
		    }
		  }
		}
	      }
	    }
	  }
	  else {
#pragma ivdep
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      ray._s0[i[0]]=rax._s5[0][0][0][0][i[0]];
	    }
	  }
	}

#endif
#endif
#endif
#endif
      }

      else if (ray.ndim==2) {

#if RARR_MAX_NDIM > 2

	if (rax.ndim==3) {
	  if (ext) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  ray._s2[i[1]][i[0]] += fac* rax._s3[i[2]][i[1]][i[0]];
		}
	      }
	    }
	  }
	  else {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		ray._s2[i[1]][i[0]] = rax._s3[0][i[1]][i[0]];
	      }
	    }
	  }
	}

#if RARR_MAX_NDIM > 3

	if (rax.ndim==4) {
	  if (ext) {
	    for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	      for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		    ray._s2[i[1]][i[0]] += fac* rax._s4[i[3]][i[2]][i[1]][i[0]];
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		ray._s2[i[1]][i[0]] = rax._s4[0][0][i[1]][i[0]];
	      }
	    }
	  }
	}

#if RARR_MAX_NDIM > 4

	if (rax.ndim==5) {
	  if (ext) {
	    for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	      for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
		for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		      ray._s2[i[1]][i[0]] += fac* rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]];
		    }
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		ray._s2[i[1]][i[0]] = rax._s5[0][0][0][i[1]][i[0]];
	      }
	    }
	  }
	}

#endif
#endif
#endif

      }

      else if (ray.ndim==3) {

#if RARR_MAX_NDIM > 3

	if (rax.ndim==4) {
	  if (ext) {
	    for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	      for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		    ray._s3[i[2]][i[1]][i[0]] += fac* rax._s4[i[3]][i[2]][i[1]][i[0]];
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  ray._s3[i[2]][i[1]][i[0]] = rax._s4[0][i[2]][i[1]][i[0]];
		}
	      }
	    }
	  }
	}

#if RARR_MAX_NDIM > 4

	if (rax.ndim==5) {
	  if (ext) {
	    for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	      for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
		for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		      ray._s3[i[2]][i[1]][i[0]] += fac* rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]];
		    }
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  ray._s3[i[2]][i[1]][i[0]] = rax._s5[0][0][i[2]][i[1]][i[0]];
		}
	      }
	    }
	  }
	}

#endif
#endif

      }

      else {
	RVLException e;
	e<<"Error: GridAdjExtendFO::operator()\n";
	e<<"  input rarr dimension not 1, 2, or 3 - only permitted spatial dims\n";
	throw e;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridAdjExtendFO::operator()\n";
      e<<"  input type error - not CP<ireal,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridAdjExtendFO::operator()\n";
      throw e;
    }
  }

  
  GridExtendOp::GridExtendOp(Space<ireal> const & _dom,
			     Space<ireal> const & _rng)
    : dom(_dom), rng(_rng), n_ext(0), ext(0), fac(0) {
    try {

      // step 1: compatibility of product structures
      ProductSpace<ireal> const * pdom = NULL;
      ProductSpace<ireal> const * prng = NULL;
      pdom = dynamic_cast<ProductSpace<ireal> const *>(&dom);
      prng = dynamic_cast<ProductSpace<ireal> const *>(&rng);
      int n_dom=1;
      int n_rng=1;
      if (pdom) n_dom=pdom->getSize();
      if (prng) n_rng=prng->getSize();
      if (n_dom != n_rng) {
	RVLException e;
	e<<"Error: GridExtendOp constructor\n";
	e<<"  domain, range spaces have differing numbers of factors\n";
	e<<"  domain space:\n";
	dom.write(e);
	e<<"  range space:\n";
	rng.write(e);
	throw e;
      }

      // step 2: every factor of range is really extension of 
      // corresponding factor in domain, and domain factors are
      // all simple spatial grids with no extended axes
      Space<ireal> const * dsp = NULL;
      Space<ireal> const * rsp = NULL;
      for (int j=0; j<n_dom; j++) {
	if (pdom) dsp = &((*pdom)[j]);
	else dsp=&dom;
	if (prng) rsp = &((*prng)[j]); 
	else rsp=&rng;
	myGridSpace const & gdom = dynamic_cast<myGridSpace const &>(*dsp);    
	myGridSpace const & grng = dynamic_cast<myGridSpace const &>(*rsp);    
	if (retrieveGlobalRank()==0) {
	  // temporary limitation: only incore range
	  if (!(grng.isIncore())) {
	    RVLException e;
	    e<<"Error: GridExtendOp constructor\n";
	    e<<"  only incore range space allowed\n";
	    e<<"  this range space:\n";
	    rsp->write(e);
	    throw e;
	  }
	  // no extended axes in domain
	  if (gdom.getGrid().gdim != gdom.getGrid().dim) {
	    RVLException e;
	    e<<"Error: GridExtendOp constructor\n";
	    e<<"  domain factor "<<j<<" has extended axes - not permitted\n";
	    e<<"  this domain factor:\n";
	    dsp->write(e);
	    throw e;
	  }
	  // range is extension of domain
	  if (gdom.getGrid().dim != grng.getGrid().dim) {
	    RVLException e;
	    e<<"Error: GridExtendOp constructor\n";
	    e<<"  spatial dims (keyword dim) of domain, range factors "<<j<<"\n";
	    e<<"  differ - must be same\n";
	    e<<"  domain factor "<<j<<":\n";
	    dsp->write(e);
	    e<<"  range factor "<<j<<":\n";
	    rsp->write(e);
	    throw e;
	  }
	  // find spatial axes, must be in same order, and come before any extended
	  // also must be geometrically same
	  for (int i=0;i<gdom.getGrid().dim;i++) {
	    int idom=-1;
	    int irng=-1;
	    for (int k=0;k<gdom.getGrid().gdim;k++) 
	      if (gdom.getGrid().axes[k].id == i) idom=k;
	    for (int k=0;k<grng.getGrid().gdim;k++) 
	      if (grng.getGrid().axes[k].id == i) irng=k;
	    if (idom<0 || irng<0 || idom > gdom.getGrid().dim-1 || idom != irng ) {
	      RVLException e;
	      e<<"Error: GridExtendOp constructor\n";
	      if (idom<0 || irng<0)
		e<<"  failed to identify axis in domain or range with id="<<i<<"\n";
	      if (idom>gdom.getGrid().dim-1)
		e<<"  spatial grid axes must be ordered first - but axis id="<<i<<" has index "<<idom<<"\n";
	      if (idom != irng)
		e<<"  indices for id="<<i<<" differ: dom="<<idom<<" rng="<<irng<<"\n";
	      throw e;
	    }
	    if ((compare_axis(gdom.getGrid().axes[idom],grng.getGrid().axes[irng]))) {
	      RVLException e;
	      e<<"Error: GridExtendOp constructor\n";
	      e<<"  axes differ for spatial axis id = "<<i<<" dom = rng axis index = "<<irng<<"\n";
	      throw e;
	    }
	  }
	  // determine number of extended axes, whether external or internal, and 
	  // scale factor
	  ireal tmpfac = REAL_ONE;
	  // set tmp external flag by first external axis, if there is one
	  int tmp_n_ext = grng.getGrid().gdim - grng.getGrid().dim;
	  bool tmpext = false;
	  if (tmp_n_ext > 0) {
	    tmpext = (grng.getGrid().axes[grng.getGrid().dim].id < 100);
	  }
	  for (int i=gdom.getGrid().dim;i<grng.getGrid().gdim;i++) {
	    if (grng.getGrid().axes[i].id > 99) {
	      if (tmpext) {
		RVLException e;
		e<<"Error: GridExtendOp constructor\n";
		e<<"  range space mixes internal, external extended axes\n";
		e<<"  not permitted in current design\n";
		//		grng.write(e);
		throw e;
	      }
	      ireal newtmpfac = REAL_ONE;
	      if (ProtectedDivision<float>(tmpfac,grng.getGrid().axes[i].d,newtmpfac)) {
		RVLException e;
		e<<"Error: GridExtendOp constructor\n";
		e<<"  zerodivide by cell vol in axis="<<i<<" id="
		 <<grng.getGrid().axes[i].id<<" of GridSpace:\n";
		//		grng.write(e);
		throw e;
	      }
	      tmpfac*=newtmpfac;
	    }
	    if (grng.getGrid().axes[i].id < 100) {
	      if (!tmpext) {
		RVLException e;
		e<<"Error: GridExtendOp constructor\n";
		e<<"  range space mixes internal, external extended axes\n";
		e<<"  not permitted in current design\n";
		//		grng.write(e);
		throw e;
	      }
	      tmpfac*=grng.getGrid().axes[i].d;
	    }
	  }
	  n_ext.push_back(tmp_n_ext);
	  ext.push_back(tmpext);
	  fac.push_back(tmpfac);
	}
	else {
	  n_ext.push_back(0);
	  ext.push_back(false);
	  fac.push_back(REAL_ZERO);
	}
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridExtendOp constructor\n";
      e<<"  at least one factor in input space is not a GridSpace\n";
      throw e;
    }
  }

  GridExtendOp::GridExtendOp(GridExtendOp const & op)
    : dom(op.dom), rng(op.rng), n_ext(op.n_ext), ext(op.ext), fac(op.fac) {}

  void GridExtendOp::apply(Vector<ireal> const & x,
			   Vector<ireal> & y) const {
    try {
      // note that LinearOp base class takes care of space
      // membership tests, hence sanity tests
      // here x is extended to y
      Components<ireal> cx(x);
      Components<ireal> cy(y);
      for (int j=0;j<(int)cx.getSize();j++) {
	if (n_ext[j] <= 0) {
	  cy[j].copy(cx[j]);
	}
	else {
	  GridFwdExtendFO f(n_ext[j], ext[j], fac[j]);
	  MPISerialFunctionObject<ireal> mpif(f);
	  cy[j].eval(mpif,cx[j]);
	}
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridExtendOp::apply\n";
      throw e;
    }
  }

  void GridExtendOp::applyAdj(Vector<ireal> const & x,
			      Vector<ireal> & y) const {
    try {
      Components<ireal> cx(x);
      Components<ireal> cy(y);
      for (int j=0;j<(int)cx.getSize();j++) {
	if (n_ext[j] <= 0) {
	  cy[j].copy(cx[j]);
	}
	else {
	  GridAdjExtendFO f(n_ext[j], ext[j], fac[j]);
	  MPISerialFunctionObject<ireal> mpif(f);
	  cy[j].eval(mpif,cx[j]);
	}
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridExtendOp::applyAdj\n";
      throw e;
    }
  }

  ostream & GridExtendOp::write(ostream & str) const {
    str<<"GridExtendOp: inject spatial grid into extended grid\n";
    str<<"Domain:\n";
    dom.write(str);
    str<<"Range:\n";
    rng.write(str);
    return str;
  }


  void HelmFO::operator()(LocalDataContainer<ireal> & x,
			  LocalDataContainer<ireal> const & y){
    try{
      float *indata=NULL;
      float *outdata=NULL;
      float *work=NULL;
      integer f2c_n1;
      integer f2c_n2;
      integer lenwork;
      ContentPackage<ireal, RARR>  & gx =
        dynamic_cast<ContentPackage <ireal, RARR>  &> (x);
      ContentPackage<ireal, RARR> const & gy =
        dynamic_cast<ContentPackage <ireal, RARR> const &> (y);
        
      // precondition - metadata are same dimn
      RARR  & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();
      int dimx; int dimy;
      int lendom;
      ra_ndim(&rax,&dimx);
      ra_ndim(&ray,&dimy);
      //cerr << "\n xdim=" << dimx << endl;
      //cerr << "\n ydim=" << dimy << endl;
      if (dimx != dimy) {
	RVLException e;
	e<<"Error: HelmFO::operator()\n";
	e<<"arguments have different dims:\n";
	e<<"dimx="<<dimx<<" dimy="<<dimy<<"\n";
	throw e;
      }
        
      // compute grid params
      IPNT gsx; IPNT gex;
      IPNT gsy; IPNT gey;
      IPNT s; IPNT e;
      ra_a_gse(&rax,gsx,gex);
      ra_a_gse(&ray,gsy,gey);
      //        cerr << "\n===========================\n";
      //        cerr << "\n gsx[0]=" << gsx[0] << endl;
      //        cerr << "\n gex[0]=" << gex[0] << endl;
      //        cerr << "\n gsx[1]=" << gsx[1] << endl;
      //        cerr << "\n gex[1]=" << gex[1] << endl;
      //        cerr << "\n===========================\n";
      //        cerr << "\n gsy[0]=" << gsy[0] << endl;
      //        cerr << "\n gey[0]=" << gey[0] << endl;
      //        cerr << "\n gsy[1]=" << gsy[1] << endl;
      //        cerr << "\n gey[1]=" << gey[1] << endl;
      //        cerr << "\n===========================\n";
      // calculate grid overlap
      for (int ii=0;ii<dimx;ii++)  {
	s[ii]=max(gsy[ii],gsx[ii]);
	e[ii]=min(gey[ii],gex[ii]);
      }
        
      f2c_n1 = n_arr[0];
      f2c_n2 = n_arr[1];
      lendom=f2c_n1*f2c_n2;
      float _scale1=scale1;
      float _scale2=scale2;
      float _power=power;
      float _datum=datum;
      integer iter=0;
        
      // initialize workspace
      lenwork = 6*n_arr[1]*n_arr[0]+3*iwave_max(n_arr[1],2*n_arr[0])+21;
      //        cerr << "\n lenwork=" << lenwork << endl;
      //        cerr << "\n length of data = " << get_datasize_grid(gdom) << endl;
      //        cerr << "\n n_arr[0] = " << n_arr[0] << endl;
      //        cerr << "\n n_arr[1] = " << n_arr[1] << endl;
      //cerr << "\n physical domain size=" << lendom << endl;
      //cerr << "\n retrieveGlobalRank()=" << retrieveGlobalRank() << endl;        
      if (!(work = (float *)malloc(lenwork*sizeof(float)))) {
	RVLException e;
	e<<"Error: HelmOp::apply - failed to allocate " << lenwork << " floats for work buffer\n";
	throw e;
      }
      // allocate data arrays
      if (!(indata = (float *)malloc(lendom*sizeof(float)))) {
	RVLException e;
	e<<"Error: HelmOp::apply - failed to allocate " << lendom << " floats for input data\n";
	throw e;
      }
      if (!(outdata = (float *)malloc(lendom*sizeof(float)))) {
	RVLException e;
	e<<"Error: HelmOp::apply - failed to allocate " << lendom << " floats for output data\n";
	throw e;
      }
      IPNT i;
      integer idx;
#if RARR_MAX_NDIM > 0
      if (dimx==1) {
	for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	  indata[i[0]-s[0]]=ray._s1[i[0]];
	}
        helm_(DirichletSides,&f2c_n1,&f2c_n2,
              &(d_arr[0]),&(d_arr[1]),
              &(_scale1),&(_scale2),
              &_power,&_datum,
              indata,
              outdata,
              work,
              &lenwork,
              &iter);
        fprintf(stderr, "\n indata [100] = %f\n", indata[100]);
        fprintf(stderr, "\n outdata [100] = %f\n", outdata[100]);
	for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	  rax._s1[i[0]]=outdata[i[0]-s[0]];
	}
      }
#endif
#if RARR_MAX_NDIM > 1
      if (dimx==2) {
	for (i[1]=s[1];i[1]<=e[1];i[1]++) {
#pragma ivdep
	  for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	    idx = (i[1]-s[1])*n_arr[0] + i[0]-s[0];
	    indata[idx]=ray._s2[i[1]][i[0]];
	  }
	}
        helm_(DirichletSides,&f2c_n1,&f2c_n2,
              &(d_arr[0]),&(d_arr[1]),
              &(_scale1),&(_scale2),
              &_power,&_datum,
              indata,
              outdata,
              work,
              &lenwork,
              &iter);
        fprintf(stderr, "\n indata [100] = %f\n", indata[100]);
        fprintf(stderr, "\n outdata [100] = %f\n", outdata[100]);
        // copy data back
	for (i[1]=s[1];i[1]<=e[1];i[1]++) {
#pragma ivdep
	  for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	    idx = (i[1]-s[1])*n_arr[0] + i[0]-s[0];
	    rax._s2[i[1]][i[0]]=outdata[idx];
	  }
	}
      }
#endif
#if RARR_MAX_NDIM > 2
      if (dimx==3) {
	//cerr << "\n dim3=" << e[2] << endl;
	for (i[2]=s[2];i[2]<=e[2];i[2]++) {
	  for (i[1]=s[1];i[1]<=e[1];i[1]++) {
#pragma ivdep
	    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	      idx = (i[1]-s[1])*n_arr[0] + i[0]-s[0];
	      indata[idx]=ray._s3[i[2]][i[1]][i[0]];
	    }
	  }
	  helm_(DirichletSides,&f2c_n1,&f2c_n2,
		&(d_arr[0]),&(d_arr[1]),
		&(_scale1),&(_scale2),
		&_power,&_datum,
		indata,
		outdata,
		work,
		&lenwork,
		&iter);
	  // copy data back
	  for (i[1]=s[1];i[1]<=e[1];i[1]++) {
#pragma ivdep
	    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	      idx = (i[1]-s[1])*n_arr[0] + i[0]-s[0];
	      rax._s3[i[2]][i[1]][i[0]]=outdata[idx];
	    }
	  }    
	}
        fprintf(stderr, "\n indata [100] = %f\n", indata[10]);
        fprintf(stderr, "\n outdata [100] = %f\n", outdata[10]);
      }
#endif
      if (dimx<1 || dimx>3) {
	RVLException e;
	e<<"Error: HelmFO::operator()\n";
	e<<"dim = "<<dimx<<" outside of admissible set {1, 2}\n";
	throw e;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"\nError: HelmFO::operator()\n";
      e<<"at least one arg is not ContentPackage<ireal,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from HelmFO::operator()\n";
      throw e;
    }
        
  }
    

  void GridHelmOp::apply(const Vector<float> & x,
			 Vector<float> & y) const {
    try {
      // extract components - fine even if only one!
      Components<float> cx(x);
      Components<float> cy(y);

      // detect product structure
      ProductSpace<ireal> const * pdom = NULL;
      pdom = dynamic_cast<ProductSpace<ireal> const *>(&dom);
      int n_fac=1;
      if (pdom) n_fac=pdom->getSize();
      Space<ireal> const * sp = NULL;
   
      // component loop
      for (int j=0; j<n_fac; j++) {
	if (pdom) sp = &((*pdom)[j]);
	else sp = &dom;

	// late tests
	myGridSpace const * gdom = dynamic_cast<myGridSpace const *>(sp);
	if (!gdom) {
	  RVLException e;
	  e<<"Error: GridHelmOp::apply\n";
	  e<<"  factor "<<j<<" of input space is not a GridSpace\n";
	  e<<"  description:\n";
	  sp->write(e);
	  throw e;	  
	}
        if (retrieveGlobalRank() == 0) {
	if (gdom->getGrid().dim != 2) {
	  RVLException e;
	  e<<"Error: GridHelmOp::apply\n";
	  e<<"  current implementation is 2D only\n";
	  throw e;
	}
        }

	IPNT n_arr;
	RPNT d_arr;
        if (retrieveGlobalRank() == 0) {
	  get_d(d_arr,gdom->getGrid());
	  get_n(n_arr,gdom->getGrid());
	}
	HelmFO fo(n_arr,d_arr,weights[0],weights[1],power,datum,DirichletSides);
	MPISerialFunctionObject<float> mpifo(fo);
	cy[j].eval(mpifo,cx[j]);    
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled in GridHelmOp::apply\n";
      throw e;
    }
            
  }
        
  void GridHelmOp::applyAdj(const Vector<float> & x,
			    Vector<float> & y) const {
    try {
      apply(x,y);
    }
    catch (RVLException & e) {
      e<<"\ncalled in GridHelmOp::applyAdj\n";
      throw e;
    }
  }

}
