#include "seqdc.hh"

using namespace RVL;

void SeqDC::eval(FunctionObject & f,
		 vector<DataContainer const *> & x) {
  try {
    // make no mix-and-match pretense - all DCs must be SeqDCs
    vector<SeqDC const *> sx(x.size(),NULL);
    vector<DataContainer const *>::const_iterator dc = x.begin();
    vector<SeqDC const *>::iterator sdc = sx.begin();
    while (dc != x.end()) {
      *sdc=dynamic_cast<SeqDC const *>(*dc);
      ++sdc; ++dc;
    }

    // first the LocalFunctionObject case - turn everything into
    // CPs, evaluate
    LocalFunctionObject<double> * lf = NULL;
    if ((lf = dynamic_cast<LocalFunctionObject<double> *>(&f))) {

      // determine max size
      size_t n = this->datalist.size();
      sdc = sx.begin();
      while (sdc != sx.end())
	n = max(n,(*sdc)->datalist.size());

      // turn this into an LDC
      ContentPackage<double,size_t> thisldc;
      thisldc.initialize(n);
      size_t i=0;
      list<double>::const_iterator cp0 = this->datalist.begin();
      while (cp0 != this->datalist.end()) {
	thisldc.getData()[i]=*cp0;
	i++; cp0++;
      }
      while (i<n) {
	thisldc.getData()[i]=0.0;
	i++;
      }

      // turn each arg into an LDC
      vector<DataContainer const *> ldcx(x.size());
      vector<DataContainer const *>::iterator ldcp = 
	ldcx.begin();
      vector<ContentPackage<double,size_t> *> ldcxtmp(x.size());
      vector<ContentPackage<double,size_t> *>::iterator ldcptmp = 
	ldcxtmp.begin();
      sdc = sx.begin();
      while (ldcp != ldcx.end()) {
	*ldcptmp = new ContentPackage<double,size_t>;
	(*ldcptmp)->initialize(n);
	i=0;
	cp0 = ((*sdc)->datalist).begin();
	while (cp0 != (*sdc)->datalist.end()) {
	  (*ldcptmp)->getData()[i]=*cp0;
	  i++; cp0++;
	}
	while (i<n) {
	  (*ldcptmp)->getData()[i]=0.0;
	  i++;
	}
	// assign the const ptr
	*ldcp = *ldcptmp;
	++ldcp;
	++ldcptmp;
	++sdc;
      }
      
      // evaluate
      thisldc.eval(*lf,ldcx);

      // transfer data from ldc back to list
      i=0;
      list<double>::iterator p0 = this->datalist.begin();
      p0 = this->datalist.begin();
      while (p0 != this->datalist.end()) {
	*p0 = thisldc.getData()[i];
	i++; p0++;
      }

      // clean up
      ldcp=ldcx.begin();
      ldcptmp=ldcxtmp.begin();
      while (ldcp != ldcx.end()) {
	delete *ldcp;
	delete *ldcptmp;
	ldcp++;
	ldcptmp++;
      }

      // that's it folks
      return;
    }

    // unary operator on sequences
    UnarySeqFO * uf = NULL;
    if ((uf = dynamic_cast<UnarySeqFO *>(&f))) {
      (*uf)(*this);
      return;
    }

    // binary operator on sequences
    BinarySeqFO * bf = NULL;
    if ((bf = dynamic_cast<BinarySeqFO *>(&f))) {
      if (x.size()<1) {
	RVLException e;
	e<<"Error: SeqDC::eval(FO) - binary sequence op\n";
	e<<"second arg not supplied\n";
	throw e;
      }
      (*bf)(*this,*(sx[0]));
      return;
    }

    // huh?
    RVLException e;
    e<<"Error: SeqDC::eval(FO)\n";
    e<<"type of FO unknown to SeqDC\n";
    throw e;
  
  }
  catch (bad_cast const&) {
    RVLException e;
    e<<"Error: SeqDC::eval(FO)\n";
    e<<"at least one input argument not SeqDC\n";
    throw e;
  }
  catch (RVLException & e) {
    e<<"\ncalled from SeqDC::eval(FO)\n";
    throw e;
  }
}

void SeqDC::eval(FunctionObjectConstEval & f,
		 vector<DataContainer const *> & x) const {
  try {
    // make no mix-and-match pretense - all DCs must be SeqDCs
    vector<SeqDC const *> sx(x.size(),NULL);
    vector<DataContainer const *>::const_iterator dc = x.begin();
    vector<SeqDC const *>::iterator sdc = sx.begin();
    while (dc != x.end()) {
      *sdc=dynamic_cast<SeqDC const *>(*dc);
      ++sdc; ++dc;
    }

    // first the LocalFunctionObject case - turn everything into
    // CPs, evaluate
    LocalConstEval<double> * lf = NULL;
    if ((lf = dynamic_cast<LocalConstEval<double> *>(&f))) {

      // determine max size
      size_t n = this->datalist.size();
      sdc = sx.begin();
      while (sdc != sx.end())
	n = max(n,(*sdc)->datalist.size());

      // turn this into an LDC
      ContentPackage<double,size_t> thisldc;
      thisldc.initialize(n);
      size_t i=0;
      list<double>::const_iterator cp0 = this->datalist.begin();
      while (cp0 != this->datalist.end()) {
	thisldc.getData()[i]=*cp0;
	i++; cp0++;
      }
      while (i<n) {
	thisldc.getData()[i]=0.0;
	i++;
      }

      // turn each arg into an LDC
      vector<DataContainer const *> ldcx(x.size());
      vector<DataContainer const *>::iterator ldcp = 
	ldcx.begin();
      vector<ContentPackage<double,size_t> *> ldcxtmp(x.size());
      vector<ContentPackage<double,size_t> *>::iterator ldcptmp = 
	ldcxtmp.begin();
      sdc = sx.begin();
      while (ldcp != ldcx.end()) {
	*ldcptmp = new ContentPackage<double,size_t>;
	(*ldcptmp)->initialize(n);
	i=0;
	cp0 = ((*sdc)->datalist).begin();
	while (cp0 != (*sdc)->datalist.end()) {
	  (*ldcptmp)->getData()[i]=*cp0;
	  i++; cp0++;
	}
	while (i<n) {
	  (*ldcptmp)->getData()[i]=0.0;
	  i++;
	}
	// assign the const ptr
	*ldcp = *ldcptmp;
	++ldcp;
	++ldcptmp;
	++sdc;
      }
      
      // evaluate
      thisldc.eval(f,ldcx);

      // clean up
      ldcp=ldcx.begin();
      ldcptmp=ldcxtmp.begin();
      while (ldcp != ldcx.end()) {
	delete *ldcp;
	delete *ldcptmp;
	ldcp++;
	ldcptmp++;
      }

      // that's it folks
      return;
    }

    // unary operator on sequences
    UnarySeqFOR * uf = NULL;
    if ((uf = dynamic_cast<UnarySeqFOR *>(&f))) {
      (*uf)(*this);
      return;
    }

    // binary operator on sequences
    BinarySeqFOR * bf = NULL;
    if ((bf = dynamic_cast<BinarySeqFOR *>(&f))) {
      if (x.size()<1) {
	RVLException e;
	e<<"Error: SeqDC::eval(FOR) - binary sequence op\n";
	e<<"second arg not supplied\n";
	throw e;
      }
      (*bf)(*this,*(sx[0]));
      return;
    }

    // huh?
    RVLException e;
    e<<"Error: SeqDC::eval(FOR)\n";
    e<<"type of FOR unknown to SeqDC\n";
    throw e;
  
  }
  catch (bad_cast const&) {
    RVLException e;
    e<<"Error: SeqDC::eval(FOR)\n";
    e<<"at least one input argument not SeqDC\n";
    throw e;
  }
  catch (RVLException & e) {
    e<<"\ncalled from SeqDC::eval(FOR)\n";
    throw e;
  }
}

ostream & SeqDC::write(ostream & str) const {
  str<<"RVL Sequence Data Container of "<<int(this->datalist.size())<<" doubles\n";
  str<<" *** begin list:\n";
  int i = 0;
  std::list<double>::const_iterator p = this->datalist.begin();
  while (p != this->datalist.end()) {
    str<<"  item "<<i<<" = "<<*p<<"\n";
    i++; p++;
  }
  str<<" *** end list\n";
  return str;
}

