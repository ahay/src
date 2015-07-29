#ifndef __RVL_FCN_FO__
#define __RVL_FCN_FO__

#include "local.hh"

/** Standard constructions of Function Objects from C-style
    functions. Each template accepts a type, and a function type for
    functions returning an object of that type from one to five
    arguments of the same type. Each number of arguments requires a
    different template, unfortunately. 
 */

namespace RVL {

  /** unary test function */
  float test1(float x) { return x*x; }

  /** binary test function */
  double test2(double x, double y) { return x*y; }

  /** common sanity test - if successful, returns length of output, else throws exception */
  template<typename T>
  int ScalarFOSanity(int na,
		      LocalDataContainer<T> & target,
		      vector<LocalDataContainer<T> const *> & sources) {
    try {
      int n=sources.size();
      if (n!=na) {
	RVLException e;
	e<<"Error: ScalarFO::operator\n";
	e<<"size of sources = "<<n<<" != na\n";
	throw e;
      }
      int nt=target.getSize();
      for (int i=0;i<n;i++) {
	if ((sources.at(i))->getSize() != nt) {
	  RVLException e;
	  e<<"Error: ScalarFO::operator\n";
	  e<<"input data len of sources "<<i<<" = "<<sources.at(i)->getSize()<<" not = target data len= "<<nt<<"\n";
	  throw e;
	}
      }	  
      return nt;
    }
    catch (out_of_range) {
      RVLException e;
      e<<"Error: ScalarFOSanity\n";
      e<<"dimensional mismatch - most likely in number of args\n";
      throw e;
    }
  }

  /** Unary template */
  template<typename T, T f(T)>
  class ScalarFO1: public LocalFunctionObject<T> {
    
  public:
    
    virtual void operator()(LocalDataContainer<T> & target,
			    vector<LocalDataContainer<T> const *> & sources) {
      try {
	int nt=ScalarFOSanity(1,target,sources);
	for (int it=0;it<nt;it++) 
	  target.getData()[it] = f(sources.at(0)->getData()[it]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ScalarFO1\n";
      }
    }
    
    string getName() const { string tmp="ScalarFO1"; return tmp; }
  };

  /** Binary template */
  template<typename T, T f(T,T)>
  class ScalarFO2: public LocalFunctionObject<T> {
    
  public:
    
    virtual void operator()(LocalDataContainer<T> & target,
			    vector<LocalDataContainer<T> const *> & sources) {
      try {
	int nt=ScalarFOSanity(2,target,sources);	
	for (int it=0;it<nt;it++) 
	  target.getData()[it] = f(sources.at(0)->getData()[it],
				   sources.at(1)->getData()[it]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ScalarFO3\n";
      }
    }
    
    string getName() const { string tmp="ScalarFO"; return tmp; }
  };
  
  /** Ternary template */
  template<typename T, T f(T,T,T)>
  class ScalarFO3: public LocalFunctionObject<T> {
    
  public:
    
    virtual void operator()(LocalDataContainer<T> & target,
			    vector<LocalDataContainer<T> const *> & sources) {
      try {
	
	int n=sources.size();
	if (n!=3) {
	  RVLException e;
	  e<<"Error: ScalarFO::operator\n";
	  e<<"size of sources = "<<n<<" != 3\n";
	  throw e;
	}
	int nt=target.getSize();
	for (int i=0;i<n;i++) {
	  if ((sources.at(i))->getSize() != nt) {
	    RVLException e;
	    e<<"Error: ScalarFO::operator\n";
	    e<<"input data len of sources "<<i<<" = "<<sources.at(i)->getSize()<<" not = target data len= "<<nt<<"\n";
	    throw e;
	  }
	}	  
	for (int it=0;it<nt;it++) 	
	  target.getData()[it] = f(sources.at(0)->getData()[it],
				   sources.at(1)->getData()[it],
				   sources.at(2)->getData()[it]);
      }
      
      catch (out_of_range) {
	RVLException e;
	e<<"Error: ScalarFO::operator()\n";
	e<<"dimensional mismatch - most likely in number of args\n";
	throw e;
      }
    }
    
    string getName() const { string tmp="ScalarFO"; return tmp; }
  };

  template<typename T, T f(T,T,T,T)>
  class ScalarFO4: public LocalFunctionObject<T> {
    
  public:
    
    virtual void operator()(LocalDataContainer<T> & target,
			    vector<LocalDataContainer<T> const *> & sources) {
      try {
	
	int n=sources.size();
	if (n!=4) {
	  RVLException e;
	  e<<"Error: ScalarFO::operator\n";
	  e<<"size of sources = "<<n<<" != 4\n";
	  throw e;
	}
	int nt=target.getSize();
	for (int i=0;i<n;i++) {
	  if ((sources.at(i))->getSize() != nt) {
	    RVLException e;
	    e<<"Error: ScalarFO::operator\n";
	    e<<"input data len of sources "<<i<<" = "<<sources.at(i)->getSize()<<" not = target data len= "<<nt<<"\n";
	    throw e;
	  }
	}	  
	for (int it=0;it<nt;it++) 	
	  target.getData()[it] = f(sources.at(0)->getData()[it],
				   sources.at(1)->getData()[it],
				   sources.at(2)->getData()[it],
				   sources.at(3)->getData()[it]);
      }
      
      catch (out_of_range) {
	RVLException e;
	e<<"Error: ScalarFO::operator()\n";
	e<<"dimensional mismatch - most likely in number of args\n";
	throw e;
      }
    }
    
    string getName() const { string tmp="ScalarFO"; return tmp; }
  };

  template<typename T, T f(T,T,T,T,T)>
  class ScalarFO5: public LocalFunctionObject<T> {
    
  public:
    
    virtual void operator()(LocalDataContainer<T> & target,
			    vector<LocalDataContainer<T> const *> & sources) {
      try {
	
	int n=sources.size();
	if (n!=5) {
	  RVLException e;
	  e<<"Error: ScalarFO::operator\n";
	  e<<"size of sources = "<<n<<" != 5\n";
	  throw e;
	}
	int nt=target.getSize();
	for (int i=0;i<n;i++) {
	  if ((sources.at(i))->getSize() != nt) {
	    RVLException e;
	    e<<"Error: ScalarFO::operator\n";
	    e<<"input data len of sources "<<i<<" = "<<sources.at(i)->getSize()<<" not = target data len= "<<nt<<"\n";
	    throw e;
	  }
	}	  
	for (int it=0;it<nt;it++) 	
	  target.getData()[it] = f(sources.at(0)->getData()[it],
				   sources.at(1)->getData()[it],
				   sources.at(2)->getData()[it],
				   sources.at(3)->getData()[it],
				   sources.at(4)->getData()[it]);	
      }
      
      catch (out_of_range) {
	RVLException e;
	e<<"Error: ScalarFO::operator()\n";
	e<<"dimensional mismatch - most likely in number of args\n";
	throw e;
      }
    }
    
    string getName() const { string tmp="ScalarFO"; return tmp; }
  };

  template<typename T, T f(T,T,T,T,T,T)>
  class ScalarFO6: public LocalFunctionObject<T> {
    
  public:
    
    virtual void operator()(LocalDataContainer<T> & target,
			    vector<LocalDataContainer<T> const *> & sources) {
      try {
	
	int n=sources.size();
	if (n!=6) {
	  RVLException e;
	  e<<"Error: ScalarFO::operator\n";
	  e<<"size of sources = "<<n<<" != 6\n";
	  throw e;
	}
	int nt=target.getSize();
	for (int i=0;i<n;i++) {
	  if ((sources.at(i))->getSize() != nt) {
	    RVLException e;
	    e<<"Error: ScalarFO::operator\n";
	    e<<"input data len of sources "<<i<<" = "<<sources.at(i)->getSize()<<" not = target data len= "<<nt<<"\n";
	    throw e;
	  }
	}	  
	for (int it=0;it<nt;it++) 	
	  target.getData()[it] = f(sources.at(0)->getData()[it],
				   sources.at(1)->getData()[it],
				   sources.at(2)->getData()[it],
				   sources.at(3)->getData()[it],
				   sources.at(4)->getData()[it],
				   sources.at(5)->getData()[it]);      
      }
      catch (out_of_range) {
	RVLException e;
	e<<"Error: ScalarFO::operator()\n";
	e<<"dimensional mismatch - most likely in number of args\n";
	throw e;
      }
    }
    
    string getName() const { string tmp="ScalarFO"; return tmp; }
  };

}

#endif
