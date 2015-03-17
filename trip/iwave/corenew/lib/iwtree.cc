#include "iwtree.hh"

namespace TSOpt {

  IWaveTree::IWaveTree(std::vector<IWAVE *> sv, IWaveInfo const & _ic)
    : own(false), sa(iwave_max(sv.size()/2,0)),
      rd(iwave_max(sv.size()/2,0)), ic(_ic) {
    // cerr<<"iwavetree private constructor\n";
    if (sv.size()>1) {
      for (size_t i=0;i<sv.size()/2;i++) {
	sa[i]=sv[i];
	rd[i]=&((sv[i]->model).ld_c);
      }
    }
    else {
      RVLException e;
      e<<"Error: IWaveTree private constructor\n";
      e<<"  input vector of size "<<sv.size()<<"\n";
      throw e;
    }
    ref=NULL;
    // cerr<<"exit iwavetree private constructor\n";
  }

  IWaveTree::IWaveTree(PARARRAY & _pars, FILE * _stream, IWaveInfo const & _ic,
		       int order)
    : own(true), sa(pow2(order)), rd(pow2(order)), ic(_ic) {
    try {
      // cerr<<"iwavetree constructor\n";
      for (size_t i=0; i<sa.size(); i++) {
	sa[i]=new IWAVE;
	int err=0;
	// cerr<<"i="<<i<<" call iwave_construct\n";
	if ((err=iwave_construct(sa[i],&(_pars),_stream,ic))) {
	  RVLException e;
	  e<<"Error: IWaveTree main constructor, IWAVE["<<i<<"]\n";
	  e<<"  returning error "<<err<<" from iwave_construct\n";
	  throw e;
	}
	// clean it up
	// cerr<<"-> rd_a_zero\n";
	if ((err=rd_a_zero(&((sa[i]->model).ld_a)))) {
	  RVLException e;
	  e<<"Error: IWaveTree main constructor, IWAVE["<<i<<"]\n";
	  e<<"  returning error "<<err<<" from rd_a_zero\n";
	  throw e;
	}
	rd[i]=&((sa[i]->model).ld_c);
      }
      // create reference state out of first half of
      // array, without allocating any new memory,
      // using private constructor
      // cerr<<"priv constr\n";
      if (order>0) ref = new IWaveTree(sa,ic);
      else ref=NULL;
      // extract time step, set up step index 
      // cerr<<"exit IWaveTree constructor\n";
    }
    catch (RVLException & e) {
      e<<"\ncalled from IWaveTree constructor\n";
      throw e;
    }
  }

  IWaveTree::~IWaveTree() {
    if (ref) delete ref;
    if (own) 
      for (size_t i=0; i<sa.size(); i++) {
	iwave_destroy(sa[i],ic.get_mdest());
	delete sa[i];
      }    
  }
   
}
