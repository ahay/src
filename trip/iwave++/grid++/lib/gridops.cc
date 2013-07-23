#include "gridops.hh"

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

#ifdef IWAVE_USE_MPI
  typedef MPIGridSpace myGridSpace; 
#else
  typedef GridSpace myGridSpace; 
#endif

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
	for (i[0]=s0[0];i[0]<e0[0];i[0]++) 
	  rax._s1[i[0]] = (ray._s1[i[0]+1]-ray._s1[i[0]])*fac;
      }

#if RARR_MAX_NDIM > 1

      else if (dir==0 && rax.ndim>1) {
	for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	  rax._s2[i[1]][e0[0]] = REAL_ZERO;
	  for (i[0]=s0[0];i[0]<e0[0];i[0]++) 
	    rax._s2[i[1]][i[0]] = (ray._s2[i[1]][i[0]+1]-ray._s2[i[1]][i[0]])*fac;
	}
      }

      else if (dir==1 && rax.ndim==2) {
	for (i[1]=s0[1];i[1]<e0[1];i[1]++) {
	  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	    rax._s2[e0[1]][i[0]] = REAL_ZERO;
	    rax._s2[i[1]][i[0]] = (ray._s2[i[1]+1][i[0]]-ray._s2[i[1]][i[0]])*fac;
	  }
	}
      }

#if RARR_MAX_NDIM > 2

      else if (dir==1 && rax.ndim > 2) {
	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	  for (i[1]=s0[1];i[1]<e0[1];i[1]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
	      rax._s3[i[2]][e0[1]][i[0]] = REAL_ZERO;
	      rax._s3[i[2]][i[1]][i[0]] = (ray._s3[i[2]][i[1]+1][i[0]]-ray._s3[i[2]][i[1]][i[0]])*fac;
	    }	
	  }
	}
      }

      else if (dir==2 && rax.ndim == 3) {
	for (i[2]=s0[2];i[2]<e0[2];i[2]++) {
	  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s3[e0[2]][i[1]][i[0]] = REAL_ZERO;
	      rax._s3[i[2]][i[1]][i[0]] = (ray._s3[i[2]+1][i[1]][i[0]]-ray._s3[i[2]][i[1]][i[0]])*fac;
	    }
	  }
	}
      }

#if RARR_MAX_NDIM > 3
	
      else if (dir==2 && rax.ndim>3) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
		rax._s4[i[3]][e0[2]][i[1]][i[0]] = REAL_ZERO;
		rax._s4[i[3]][i[2]][i[1]][i[0]] = (ray._s4[i[3]][i[2]+1][i[1]][i[0]]-
						   ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;
	      }	
	    }
	  }
	}
      }

      else if (dir==3 && rax.ndim==4) {
	for (i[3]=s0[3];i[3]<e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
		rax._s4[e0[3]][i[2]][i[1]][i[0]] = REAL_ZERO;
		rax._s4[i[3]][i[2]][i[1]][i[0]] = (ray._s4[i[3]+1][i[2]][i[1]][i[0]]-
						   ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;
	      }	
	    }
	  }
	}
      }

#if RARR_MAX_NDIM > 4

      else if (dir==3 && rax.ndim>4) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<e0[3];i[3]++) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
		  rax._s4[i[4]][e0[3]][i[2]][i[1]][i[0]] = REAL_ZERO;
		  rax._s4[i[4]][i[3]][i[2]][i[1]][i[0]] = (ray._s4[i[4]][i[3]+1][i[2]][i[1]][i[0]]-
							   ray._s4[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;
		}	
	      }
	    }
	  }
	}
      }

      else if (dir==4 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
		  rax._s4[e0[4]][i[3]][i[2]][i[1]][i[0]] = REAL_ZERO;
		  rax._s4[i[4]][i[3]][i[2]][i[1]][i[0]] = (ray._s4[i[4]+1][i[3]][i[2]][i[1]][i[0]]-
							   ray._s4[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;
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
	for (i[0]=s0[0]+1;i[0]<e0[0];i[0]++) 
	  rax._s1[i[0]]=(ray._s1[i[0]-1]-ray._s1[i[0]])*fac;      }

#if RARR_MAX_NDIM > 1

      else if (dir==0 && rax.ndim>1) {
	for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	  for (i[0]=s0[0]+1;i[0]<e0[0];i[0]++) {
	    rax._s2[i[1]][s0[0]] = -ray._s2[i[1]][s0[0]]*fac;
	    rax._s2[i[1]][e0[0]] = ray._s2[i[1]][e0[0]-1]*fac;
	    rax._s2[i[1]][i[0]]=(ray._s2[i[1]][i[0]-1]-ray._s2[i[1]][i[0]])*fac;      
	  }
	}
      }
      

      else if (dir==1 && rax.ndim==2) {
	for (i[1]=s0[1]+1;i[1]<e0[1];i[1]++) {
	  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	    rax._s2[s0[1]][i[0]] = -ray._s2[s0[1]][i[0]]*fac;
	    rax._s2[e0[1]][i[0]] = ray._s2[e0[1]-1][i[0]]*fac;
	    rax._s2[i[1]][i[0]]=(ray._s2[i[1]-1][i[0]]-ray._s2[i[1]][i[0]])*fac;      
	  }
	}
      }

#if RARR_MAX_NDIM > 2

      else if (dir==1 && rax.ndim > 2) {
	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	  for (i[1]=s0[1]+1;i[1]<e0[1];i[1]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s3[i[2]][s0[1]][i[0]] = -ray._s3[i[2]][s0[1]][i[0]]*fac;
	      rax._s3[i[2]][e0[1]][i[0]] = ray._s3[i[2]][e0[1]-1][i[0]]*fac;
	      rax._s3[i[2]][i[1]][i[0]]=(ray._s3[i[2]][i[1]-1][i[0]]-ray._s3[i[2]][i[1]][i[0]])*fac;      
	    }
	  }
	}
      }

      else if (dir==2 && rax.ndim == 3) {
	for (i[2]=s0[2]+1;i[2]<e0[2];i[2]++) {
	  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s3[s0[2]][i[1]][i[0]] = -ray._s3[s0[2]][i[1]][i[0]]*fac;
	      rax._s3[e0[2]][i[1]][i[0]] = ray._s3[e0[2]-1][i[1]][i[0]]*fac;
	      rax._s3[i[2]][i[1]][i[0]]=(ray._s3[i[2]-1][i[1]][i[0]]-ray._s3[i[2]][i[1]][i[0]])*fac;      
	    }
	  }
	}
      }

#if RARR_MAX_NDIM > 3
	
      else if (dir==2 && rax.ndim>3) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[2]=s0[2]+1;i[2]<e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s4[i[3]][s0[2]][i[1]][i[0]] = -ray._s4[i[3]][s0[2]][i[1]][i[0]]*fac;
		rax._s4[i[3]][e0[2]][i[1]][i[0]] = ray._s4[i[3]][e0[2]-1][i[1]][i[0]]*fac;
		rax._s4[i[3]][i[2]][i[1]][i[0]]=(ray._s4[i[3]][i[2]-1][i[1]][i[0]]-
						 ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;      
	      }
	    }
	  }
	}
      }

      else if (dir==3 && rax.ndim==4) {
	for (i[3]=s0[3]+1;i[3]<e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s4[s0[3]][i[2]][i[1]][i[0]] = -ray._s4[s0[3]][i[2]][i[1]][i[0]]*fac;
		rax._s4[e0[3]][i[2]][i[1]][i[0]] = ray._s4[e0[3]-1][i[2]][i[1]][i[0]]*fac;
		rax._s4[i[3]][i[2]][i[1]][i[0]]=(ray._s4[i[3]-1][i[2]][i[1]][i[0]]-
						 ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;      
	      }
	    }
	  }
	}
      }

#if RARR_MAX_NDIM > 4

      else if (dir==3 && rax.ndim>4) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[3]=s0[3]+1;i[3]<e0[3];i[3]++) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  rax._s5[i[4]][s0[3]][i[2]][i[1]][i[0]] = -ray._s5[i[4]][s0[3]][i[2]][i[1]][i[0]]*fac;
		  rax._s5[i[4]][e0[3]][i[2]][i[1]][i[0]] = ray._s5[i[4]][e0[3]-1][i[2]][i[1]][i[0]]*fac;
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]]=(ray._s5[i[4]][i[3]-1][i[2]][i[1]][i[0]]-
						   ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;      
		}
	      }
	    }
	  }
	}
      }

      else if (dir==4 && rax.ndim==5) {
	for (i[4]=s0[4]+1;i[4]<e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  rax._s5[s0[4]][i[3]][i[2]][i[1]][i[0]] = -ray._s5[s0[4]][i[3]][i[2]][i[1]][i[0]]*fac;
		  rax._s5[e0[4]][i[3]][i[2]][i[1]][i[0]] = ray._s5[e0[4]-1][i[3]][i[2]][i[1]][i[0]]*fac;
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
	if (n_fac > 1) sp = &(pdom[j]);
	else sp=&dom;
	myGridSpace const & gdom = dynamic_cast<myGridSpace const &>(*sp);    
	// pure out of core: real factors only on rk=0
	if (retrieveGlobalRank()==0) {
	  if (dir < 0 || dir > gdom.getGrid().gdim-1) {
	    RVLException e;
	    e<<"Error: GridDerivOp constructor\n";
	    e<<"  direction index "<<dir<<" out of dimension range [0,"<<gdom.getGrid().gdim-1<<"\n";
	    throw e;
	  }
	  RPNT d;
	  get_d(d,gdom.getGrid());
	  fac.push_back(scale/d[dir]);
	}
	else {
	  fac.push_back(REAL_ZERO);
	}
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridDerivOp constructor\n";
      e<<"  at least one factor in input space is not a GridSpace\n";
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
}
