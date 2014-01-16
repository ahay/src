#include "rns.hh"

namespace TSOpt {

  int create_rn(rn & s) {
    s.it=0;
    s.nu=0;
    s.nc=0;
    s.u=NULL;
    s.c=NULL;
    return 0;
  }

  int init_rn(rn & s, int _it, size_t _nu, size_t _nc) {
    if (_nu<1 || _nc<1 || s.u || s.c) { 
      cout << "_nu = " << _nu << endl;
      cout << "_nc = " << _nc << endl;
      cout << "s.u = " << s.u[0] << endl;
      cout << "s.c =  " << s.c[0] << endl;
      return 1;
    }
    s.it=_it; s.nu=_nu; s.nc=_nc; 
    //    s.u = new (nothrow) float[s.nu];
    //    s.c = new (nothrow) float[s.nc];
    s.u = new float[s.nu];
    s.c = new float[s.nc];
    //    cerr<<"rn init this = "<<&s<<endl;
    //    cerr<<"rn init u="<<s.u<<endl;
    //    cerr<<"rn init c="<<s.c<<endl;
    
    return 0;
  }

  int copy_rn(rn & s, rn const & t) {
    
    
    if (t.nu<1||t.nc<1) {
      cout << "t.nu = " << t.nu << endl;
      cout << "t.nc = " << t.nc << endl;
      return 1;
    }
    
    //printf("t.it = %d, t.nu = %d,  t.nc = %d", t.it, t.nu, t.nc); 
    
    s.it=t.it; 
    if (s.nu==t.nu && s.nc==t.nc && s.u && s.c && t.u && t.c) {
    /*
    if (s.u) {delete [] s.u; s.u = NULL;}
    if (s.c) {delete [] s.c; s.c = NULL;}

    s.u = new float[s.nu];
    s.c = new float[s.nc];
    */

      for (size_t i=0;i<s.nu;i++) s.u[i]=t.u[i];
      for (size_t i=0;i<s.nc;i++) s.c[i]=t.c[i];
      
      //      cerr<<"rn copy this = "<<&s<<endl;
      //      cerr<<"rn copy u="<<s.u<<endl;
      //      cerr<<"rn copy c="<<s.c<<endl;
      return 0;
    }
    else {
      return 1;
    }
  }
  
  void delete_rn(rn & s) { 
    //    cerr<<"rn delete this="<<&s<<endl;
    if (s.u) { 
      //cerr<<"rn delete u = "<<s.u<<endl; 
      delete [] s.u;  s.u = NULL; }
    if (s.c) { 
      //cerr<<"rn delete c = "<<s.c<<endl; 
      delete [] s.c;  s.c = NULL; }
    create_rn(s);
  }

  void fprint_rn(rn const & s, ostream & fp) {
    fp<<"rn struct: nu="<<s.nu<<" nc="<<s.nc<<" it="<<s.it<<"\n";
  
    if (s.u) {
      for (size_t i=0;i<s.nu;i++) fp<<"i="<<i<<" u="<<s.u[i]<<"\n";
    }
    else {
      fp<<"data uninitialized\n";
    }
    
    if (s.c) {
      for (size_t i=0;i<s.nc;i++) fp<<"i="<<i<<" c="<<s.c[i]<<"\n";
    }
    else {
      fp<<"control uninitialized\n";
    }
    
    
  }

  RnState::RnState() {
    create_rn(s);
    dt=0;
  }

  RnState::RnState(RnState const & t) { 
    create_rn(s);
    init_rn(s,t.getrn().it,t.getrn().nu,t.getrn().nc);
    if (copy_rn(s,t.s)) {
      RVLException e;
      e<<"Error: RnState::copy constr from copy_rn\n";
      throw e;
    } 
    dt=t.dt; 
    //    cerr<<"-----------------------"<<endl;
    //    cerr<<"RnState copy constr\n";
    //    this->write(cerr);
    //    cerr<<"-----------------------"<<endl;
  }

  RnState::~RnState()  { 
    //    cerr<<"-----------------------"<<endl;
    //    cerr<<"RnState destructor\n";
    //    this->write(cerr);
    //    cerr<<"-----------------------"<<endl;
    delete_rn(s); 
  }

  void RnState::initialize(size_t nu, size_t nc, int it) {
    if (init_rn(s,it,nu,nc)) {
      RVLException e;
      e<<"Error: RnState::initialize from init_rn\n";
      throw e;
    }
    dt=it;
    //    cerr<<"-----------------------"<<endl;
    //    cerr<<"RnState initialize\n";
    //    this->write(cerr);
    //    cerr<<"-----------------------"<<endl;
  }
  
  void RnState::setTime(Time const & t) {
    try {
      dt=t;
      s.it=dt.getint();
    }
    catch (RVLException & e) {
      e<<"\ncalled from RnState::setTime\n";
      throw e;
    }
  }

  Time const & RnState::getTime() const { 
    dt = s.it;
    return dt;
  }

  void RnState::copy(RnState const & x) { 
    
    
    //cout << "in RnState::copy" << endl;
    //fprint_rn(s, cout);
    //fprint_rn(x.s,cout);

    try {
      if (copy_rn(s,x.s)) {
	RVLException e;
	e<<"Error: RnState::copy from copy_rn\n";
	throw e;
      } 
      dt=x.dt;
      //      cerr<<"-----------------------"<<endl;
      //      cerr<<"RnState copy \n";
      //      this->write(cerr);
      //      cerr<<"-----------------------"<<endl;
    }
    catch (RVLException & e) {
      e<<"\ncalled from RnState::copy\n";
      throw e;
    }
  }





}
