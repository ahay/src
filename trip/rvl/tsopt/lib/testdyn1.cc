#include "tests.hh"

namespace TSOpt {

  Dyn1::Dyn1(ContentPackage<float,size_t> const & _c, 
	     float _dt, 
	     bool _verbose, 
	     ostream & _fpout)
  :  StdDFTStep<RnState>(), fpout(_fpout), dt(_dt), verbose(_verbose) {
    try {
      size_t n=_c.getMetadata();
      mystate.initialize(n,n,0);
      memcpy(mystate.getrn().u,_c.getData(),n*sizeof(float));
      memcpy(mystate.getrn().c,_c.getData(),n*sizeof(float));	
      
      //    cerr<<"in Dyn1\n";
      //    cerr<<"rn mystate="<<&(mystate.getrn())<<endl;
      //    cerr<<"rn mystate u = "<<mystate.getrn().u<<endl;
      //    cerr<<"rn mystate c  = "<<mystate.getrn().c<<endl;
      
      StdDiscreteTime t0;
      t0=0;
      mystate.setTime(t0);
      

    }
    catch (RVLException & e) {
      e<<"\ncalled from Dyn1 constructor\n";
      throw e;
    }
  }

  Dyn1::Dyn1(float _dt, 
	     bool _verbose, 
	     ostream & _fpout)
  : StdDFTStep<RnState>(), fpout(_fpout), dt(_dt), verbose(_verbose) {
    
    //    cerr << "calling Dyn1 alt constr" << endl;
  }

  Dyn1::Dyn1(Dyn1 const & t) 
    : StdDFTStep<RnState>(t), fpout(t.fpout), dt(t.dt), verbose(t.verbose) {

  }


  void Dyn1::dyn1fwd(rn const & x, float dt) {
    for (size_t i=0;i<x.nu;i++)
      x.u[i] = x.u[i] + dt*(1.0 - x.u[i]*x.u[i]);    
    
    ++(x.it);
    StdDiscreteTime t;
    t = x.it;
    mystate.setTime(t);
  }

  void Dyn1::run() {
    try {

      if (verbose) {
	fpout<<"------------------------------\n";
	fpout<<"Dyn1 State - before \n";
	mystate.write(fpout);
	fpout<<"------------------------------\n";
      }
      dyn1fwd(mystate.getrn(),dt);
      if (verbose) {
	fpout<<"------------------------------\n";
	fpout<<"Dyn1 State - after \n";
	mystate.write(fpout);
	fpout<<"------------------------------\n";
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from Dyn1::run\n";
      throw e;
    }
  }

  ostream & Dyn1::write(ostream & str) const {
    str<<"Dyn1\n";
    return str;
  }

  /*    
  void D_Dyn1InitFO::operator () (LocalDataContainer<float> const & x) {
    try {
      // initialize dimension - it=0
      s.initialize(x.getSize(),x.getSize());
      memcpy(s.getrn().c,x.getData(),x.getSize()*sizeof(float));
      memcpy(s.getrn().u,x.getData(),x.getSize()*sizeof(float));
      DiscreteTime t;
      t=0;
      s.setTime(t);
    }
    catch (RVLException & e) {
      e<<"\ncalled from Dyn1InitFO::operator()\n";
      throw e;
    } 
  }
  */

  D_Dyn1::D_Dyn1(ContentPackage<float,size_t> const & _dc, 
		 RnState & _ref, 
		 float _dt, 
		 bool _verbose, 
		 ostream & _fpout)
  : StdDFTStep<RnState>(),
	 fpout(_fpout),  
	 dt(_dt), 
	 verbose(_verbose) ,
	 ref(_ref)
  {

    size_t n = ref.getrn().nu;
    if (n != _dc.getMetadata()) {
      RVLException e; 
      e<<"Error: D_Dyn1 main constructor: incompatible control pert, ref\n";
      throw e;
    }
    mystate.initialize(n,n,0);
    memcpy(mystate.getrn().u,_dc.getData(),n*sizeof(float));
    memcpy(mystate.getrn().c,_dc.getData(),n*sizeof(float));	

    //    cerr<<"in D_Dyn1\n";
    //    cerr<<"rn mystate="<<&(mystate.getrn())<<endl;
    //    cerr<<"rn mystate u = "<<mystate.getrn().u<<endl;
    //    cerr<<"rn mystate c  = "<<mystate.getrn().c<<endl;

    StdDiscreteTime t0;
    t0 = 0;
    mystate.setTime(t0);


  
  }

  D_Dyn1::D_Dyn1(RnState & _ref, 
		 float _dt, 
		 bool _verbose, 
		 ostream & _fpout)
  : StdDFTStep<RnState>(), 
    fpout(_fpout),
    dt(_dt), 
    verbose(_verbose), 
    ref(_ref)
  {

  }

  D_Dyn1::D_Dyn1(D_Dyn1 const & t) 
    : StdDFTStep<RnState>(t), 
    fpout(t.fpout),
      dt(t.dt), 
      verbose(t.verbose), 
      ref(t.ref) 
  {
  
  }

  void D_Dyn1::dyn1lin(rn const & du, rn const & u, 
		       float dt) {
    try {
      if (du.nu != u.nu) {
	RVLException e;
	e<<"Error Dyn::dyn1lin - barf!\n";
	throw e;
      }

      for (size_t i=0;i<u.nu;i++) 
	du.u[i]*=(1.0-2.0*dt*u.u[i]);
      ++(du.it);
      
      StdDiscreteTime t;
      t = du.it;
      mystate.setTime(t);

    }

    catch (bad_cast) {
      RVLException e;
      e<<"Error Dyn1::linact - args not RnState instances\n";
      throw e;
    }      
  }
  
  void D_Dyn1::run() {
    try {


      if (verbose) {
	fpout<<"------------------------------\n";
	fpout<<"D_Dyn1 before State = \n";
	mystate.write(fpout);
	fpout<<"------------------------------\n";
      }

      dyn1lin(mystate.getrn(), ref.getrn(),dt);

      if (verbose) {
	fpout<<"------------------------------\n";
	fpout<<"D_Dyn1 after State = \n";
	mystate.write(fpout);
	fpout<<"------------------------------\n";
      }

    }

    catch (RVLException & e) {
      e<<"\ncalled from D_Dyn1::run\n";
      throw e;
    }

  }

  ostream & D_Dyn1::write(ostream & str) const {
    str<<"D_Dyn1\n";
    return str;
  }


  /** Implementations for the adjoint evolution */
  A_Dyn1::A_Dyn1(ContentPackage<float,size_t> const & _ac, 
		 RnState & _ref, 
		 float _dt, 
		 bool _verbose, 
		 ostream & _fpout)
  : StdDBTStep<RnState>(), 
	 fpout(_fpout),
	 dt(_dt), 
	 verbose(_verbose), 
	 ref(_ref)
  {
    //    cerr<<"======================================="<<endl;
    //    cerr<<"A_Dyn1 constructor = "<<this<<"\n";

    size_t n = ref.getrn().nu;
    if (n != _ac.getMetadata()) {
      RVLException e; 
      e<<"Error: A_Dyn1 main constructor: incompatible control pert, ref\n";
      throw e;
    }

    mystate.initialize(n,n,0);
    for (size_t i=0;i<n;i++) {
      mystate.getrn().u[i]=_ac.getData()[i];
      mystate.getrn().c[i]=_ac.getData()[i];
    }
    //    memcpy(mystate.getrn().u,_ac.getData(),n*sizeof(float));
    //    memcpy(mystate.getrn().c,_ac.getData(),n*sizeof(float));	

    //    cerr<<"rn mystate="<<&(mystate.getrn())<<endl;
    //    cerr<<"rn mystate u = "<<mystate.getrn().u<<endl;
    //    cerr<<"rn mystate c  = "<<mystate.getrn().c<<endl;

    StdDiscreteTime t0;
    t0 = 9;  // <ME?> set to final sim time? 
             //(how to pass this info without setting it manually?)
    mystate.setTime(t0);
     
    //    cerr<<"======================================="<<endl;        

  }

  A_Dyn1::A_Dyn1(RnState & _ref, 
		 float _dt, 
		 bool _verbose, 
		 ostream & _fpout)
    : StdDBTStep<RnState>(), 
	 fpout(_fpout),
	 dt(_dt),
	 verbose(_verbose), 
	 ref(_ref) 
  {
    //    cerr<<"in A_Dyn1 2nd constr\n";
    //    cerr<<"rn mystate="<<&(mystate.getrn())<<endl;
    //    cerr<<"rn mystate u = "<<mystate.getrn().u<<endl;
    //    cerr<<"rn mystate c  = "<<mystate.getrn().c<<endl;   
  }

  A_Dyn1::A_Dyn1(A_Dyn1 const & t) 
    : StdDBTStep<RnState>(t), 
	 fpout(t.fpout),
	 dt(t.dt), 
	 verbose(t.verbose), 
	 ref(t.ref) 
  {
    //    cerr<<"A_Dyn1 alt constr\n";
  }

  void A_Dyn1::dyn1adj(rn const & au, rn const & u, 
		       float dt) {

     try {
      if (au.nu != u.nu) {
	RVLException e;
	e<<"Error Dyn::dyn1lin - barf!\n";
	throw e;
      }

      // for each element in the state field, update
      for (size_t i=0;i<u.nu;i++) 
	au.u[i]*=(1.0-2.0*dt*u.u[i]);

      //decrement time index
      --(au.it);
  
      StdDiscreteTime t;
      t = au.it;
      mystate.setTime(t);

    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error Dyn1::linact - args not RnState instances\n";
      throw e;
    }      
  }
  
  void A_Dyn1::run() {

    //cout << endl << "*** in A_DYN1::run() ***" << endl;

    try {

      if (verbose) {
	fpout<<"------------------------------\n";
	fpout<<"A_Dyn1 before State = "<<this<<"\n";
	mystate.write(fpout);
	
	fpout<<"\nRef before State = \n";
	ref.write(fpout);
	
	fpout<<"------------------------------\n";
      }
      dyn1adj(mystate.getrn(), ref.getrn(),dt);
      if (verbose) {
	fpout<<"------------------------------\n";
	fpout<<"A_Dyn1 after State = \n";
	mystate.write(fpout);

	//    cerr<<"rn mystate="<<&(mystate.getrn())<<endl;
	//    cerr<<"rn mystate u = "<<mystate.getrn().u<<endl;
	//    cerr<<"rn mystate c  = "<<mystate.getrn().c<<endl;   
	
	fpout<<"\nRef after State = \n";
	ref.write(fpout);

	fpout<<"------------------------------\n";
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from A_Dyn1::run\n";
      throw e;
    }

  }

  ostream & A_Dyn1::write(ostream & str) const {
    str<<"======================================="<<endl;
    str<<"A_Dyn1 object = "<<this<<"\n";
    mystate.write(str);
    str<<"======================================="<<endl;
    return str;
  }





}




/*			 
void testdyn1::derStep(LocalDataContainer<double> & u,
			 LocalDataContainer<double> & du,
			 LocalDataContainer<double> const & c,
			 LocalDataContainer<double> const & dc,
			 Clock<double> const & clk) {
	double a = clk.getTimeStep();
	if (verbose) {
		cerr<<endl<<"//////////////////////////////////////"<<endl;
		cerr<<"drhs: a = "<<a<<" t = "<<clk.getTime()<<endl;
	cerr<<"before:"<<endl;
		cerr<<"c  = "<<c.getData()[0]<<endl;
		cerr<<"u  = "<<u.getData()[0]<<endl;
		cerr<<"dc = "<<dc.getData()[0]<<endl;
		cerr<<"du = "<<du.getData()[0]<<endl;
	}
	int n = u.getSize();
	double * uc = u.getData();
	double * duc = du.getData();
	for (int i=0;i<n;i++) {
		duc[i] = duc[i]*(1.0 - 2.0*a*uc[i]);
	//	uc[i] = uc[i] + a*(1.0 - uc[i]*uc[i]);
	}
	if (verbose) {
		cerr<<"after:"<<endl;
		cerr<<"c  = "<<c.getData()[0]<<endl;
		cerr<<"u  = "<<u.getData()[0]<<endl;
		cerr<<"dc = "<<dc.getData()[0]<<endl;
		cerr<<"du = "<<du.getData()[0]<<endl;
		cerr<<"//////////////////////////////////////"<<endl;
	}
	
}
void testdyn1::adjStep(LocalDataContainer<double> & u,
			 LocalDataContainer<double> & du,
			 LocalDataContainer<double> const & c,
			 LocalDataContainer<double> & dc,
			 Clock<double> const & clk) {
	double a = clk.getTimeStep();
	if (verbose) {
		cerr<<endl<<"//////////////////////////////////////"<<endl;
		cerr<<"arhs: a = "<<a<<" t = "<<clk.getTime()<<endl;
		cerr<<"before:"<<endl;
		cerr<<"c  = "<<c.getData()[0]<<endl;
		cerr<<"u  = "<<u.getData()[0]<<endl;
		cerr<<"dc = "<<dc.getData()[0]<<endl;
		cerr<<"du = "<<du.getData()[0]<<endl;
	}
	int n = u.getSize();
	for (int i=0;i<n;i++) {
		du.getData()[i] = du.getData()[i]*(1.0 - 2.0*a*u.getData()[i]);
	}
	if (verbose) {
		cerr<<"after:"<<endl;
		cerr<<"c  = "<<c.getData()[0]<<endl;
		cerr<<"u  = "<<u.getData()[0]<<endl;
		cerr<<"dc = "<<dc.getData()[0]<<endl;
		cerr<<"du = "<<du.getData()[0]<<endl;
		cerr<<"//////////////////////////////////////"<<endl;
	}
	
}
   //void query(LocalDataContainer<double> const & x) {}
 
}











*/
