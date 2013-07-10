#include "sigrho.hh"

namespace ASG {
  
  void SigRhoOp::apply(Vector<float> const & x,
		       Vector<float> & y) const {
    try {
      // sanity test - dom & rng have all necessary str
      SpaceTest<float>(dom,x,"SigRhoOp::apply - dom");
      SpaceTest<float>(rng,y,"SigRhoOp::apply - rng");

      // divvy up both input and output into components
      Components<float> cx(x);
      Components<float> cy(y);

      // now use ScalarFOs
      // cy[0]=kappa
      cy[0].eval(kappa_sigrho,cx[0],cx[1]);
      // cy[1]=buoy
      cy[1].eval(buoy_sigrho,cx[0],cx[1]);
      // all done!!
    }
    catch (RVLException &e) {
      e<<"\ncalled from SigRhoOp::apply\n";
      throw e;
    }
  }

  void SigRhoOp::applyDeriv(Vector<float> const & x,
			    Vector<float> const & dx,
			    Vector<float> & dy) const {
    try {
      // sanity test
      SpaceTest<float>(dom,x,"SigRhoOp::applyDeriv - dom");
      SpaceTest<float>(dom,dx,"SigRhoOp::applyDeriv - dom pert");
      SpaceTest<float>(rng,dy,"SigRhoOp::applyDeriv - rng pert");

      // workspace for building up total derivs
      Vector<float> tmp(dy.getSpace());

      // get components
      Components<float> cx(x);
      Components<float> cdx(dx);
      Components<float> cdy(dy);
      Components<float> ctmp(tmp);

      // first compute dkappa/dsig and store in cdy[0]
      cdy[0].eval(dkappa_dsig,cx[0],cx[1],cdx[0]);
      // compute dkappa/drho, store in tmp[0]
      ctmp[0].eval(dkappa_drho,cx[0],cx[1],cdx[1]);
      // overwrite cdy[0] (dkappa) with sum of cdy[0] and tmp[0]
      // note y.linComb(a,x) does y <- ax+y (axpy)
      cdy[0].linComb(1.0f,ctmp[0]);

      // next compute dbuoy/dsig and store in cdy[1]
      cdy[1].eval(dbuoy_dsig,cx[0],cx[1],cdx[0]);
      // compute dkappa/drho, store in tmp[0]
      ctmp[1].eval(dbuoy_drho,cx[0],cx[1],cdx[1]);
      // overwrite cdy[1] (dbuoy) with sum of cdy[1] and tmp[1]
      cdy[1].linComb(1.0f,ctmp[1]);
	
      // all done!
    }
    catch (RVLException &e) {
      e<<"\ncalled from SigRhoOp::applyDeriv\n";
      throw e;
    }
  }
    
  void SigRhoOp::applyAdjDeriv(Vector<float> const & x,
			       Vector<float> const & dy,
			       Vector<float> & dx) const {
    try {
       // sanity test
      SpaceTest<float>(dom,x,"SigRhoOp::applyDeriv - dom");
      SpaceTest<float>(dom,dx,"SigRhoOp::applyDeriv - dom pert");
      SpaceTest<float>(rng,dy,"SigRhoOp::applyDeriv - rng pert");

      // workspace for building up total derivs
      Vector<float> tmp(dx.getSpace());

      // get components
      Components<float> cx(x);
      Components<float> cdx(dx);
      Components<float> cdy(dy);
      Components<float> ctmp(tmp);

      // first compute dkappa/dsig and store in cdx[0]
      cdx[0].eval(dkappa_dsig,cx[0],cx[1],cdy[0]);
      // compute dbuoy/drho, store in tmp[0]
      ctmp[0].eval(dbuoy_dsig,cx[0],cx[1],cdy[1]);
      // overwrite cdx[0] (dkappa) with sum of cdx[0] and tmp[0]
      // note y.linComb(a,x) does y <- ax+y (axpy)
      cdx[0].linComb(1.0f,ctmp[0]);

      // next compute dkappa/drho and store in cdx[1]
      cdx[1].eval(dkappa_drho,cx[0],cx[1],cdy[0]);
      // compute dbuoy/drho, store in tmp[0]
      ctmp[1].eval(dbuoy_drho,cx[0],cx[1],cdy[1]);
      // overwrite cdx[1] (dbuoy) with sum of cdx[1] and tmp[1]
      cdx[1].linComb(1.0f,ctmp[1]);
	
      // all done!
    }
    catch (RVLException &e) {
      e<<"\ncalled from SigRhoOp::applyAdjDeriv\n";
      throw e;
    }
  }

  SigRhoOp::SigRhoOp(Space<float> const & _dom,
		     Space<float> const & _rng)
    : dom(_dom), rng(_rng) {
    try {
      // check that these are both 2-component product spaces by
      // dynamic cast. Don't bother to test other attributes - could
      // test that both components are GridSpaces or MPIGridSpaces,
      // and even that the components have the right physical types,
      // but not for now
      ProductSpace<float> const & td 
	= dynamic_cast<ProductSpace<float> const &>(dom);
      if (td.getSize() != 2) {
	RVLException e;
	e<<"Error: SigRhoOp constructor\n";
	e<<"domain is not 2 component ProductSpace\n";
	throw e;
      }
      ProductSpace<float> const & tr
	= dynamic_cast<ProductSpace<float> const &>(rng);
      if (tr.getSize() != 2) {
	RVLException e;
	e<<"Error: SigRhoOp constructor\n";
	e<<"range is not 2 component ProductSpace\n";
	throw e;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: SigRhoOp constructor\n";
      e<<"either domain or range is not ProductSpace\n";
      throw e;
    }
  }
}
