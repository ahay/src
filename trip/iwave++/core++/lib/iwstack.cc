#include "iwstack.hh"

namespace TSOpt{
  
  using RVL::RVLException;

  /** IWaveTimeStack */
  
  /** check */
  bool IWaveTimeStack::iscurrentState(IWaveState const & _iw){
    int it = _iw.IWaveState::getTSIndex().getCstruct().it;
    int iv = _iw.IWaveState::getTSIndex().getCstruct().iv;
    int cit = iw.IWaveState::getTSIndex().getCstruct().it;
    int civ = iw.IWaveState::getTSIndex().getCstruct().iv;
    if((it==cit) && (iv==civ)) return true;
    else return false;
  }
    
  /** place current state's time at the back of the stack */
  void IWaveTimeStack::push_back(){
    int err = 0;
    int it = iw.IWaveState::getTSIndex().getCstruct().it;
    int iv = iw.IWaveState::getTSIndex().getCstruct().iv;
    // cout<<"IWaveTimeStack::push_back, it = "<<it<<endl;
    /* push_back in two cases:
       1. max_ic_size <= 0: store everything outcore; 
       2. max_ic_size > 0: 
       2(a). ic_size < max_ic_size: store state incore;
       2(b). ic_size == max_ic_size: move incore stack forwards, i.e., 
       move ic_front to outcore, 
       update incore/outcore stack,
       push back cur-state into incore stack           
    */
#ifdef _IWSTACK_DEBUG
    write_IWStack();
#endif
    if(iv==0) {
      IWAVE & iwstate = iw.IWaveState::getIWAVE();
      
      if(max_ic_size <= 0) {
	/*case 1: store dynamic fields into outcore checkfiles */
	err = iwave_dynamic_takeshot(&iwstate,it-itorigin); 
	if (err) {
	  RVLException e;
	  e<<"Error: IWaveTimeStack::push_back()\n";
	  e<<"iwave_dynamic_takeshot return error \n";
	  throw e;
	}
	oc_itstack.push_back(it);
      } 
      else { /* case 2 starts: */
	  if((int) ic_itstack.size() < max_ic_size) {
	  /* case 2a: store dynamic fields into incore stack */
	  
	  /* allocate storage */
	  RDOM * tmp = new RDOM;
	  rd_a_setnull(tmp);
	  tmp->narr = iwstate.model.ld_a.narr;	  	  	  
	  /* copy each dynamic array into tmp */
	  for (int iarr=0; iarr < tmp->narr; iarr++) {
	    if (fd_isdyn(&(iwstate.model),iarr)) {
	      /* copy dynamic arrays of ld_a to tmp */
	      int _ndim = iwstate.model.ld_a._s[iarr].ndim;
	      IPNT gs, ge;
	      for (int i=0; i< _ndim; ++i){
		gs[i] = iwstate.model.ld_a._s[iarr]._dims[i].gs0;
		ge[i] = iwstate.model.ld_a._s[iarr]._dims[i].ge0;
	      }
	      err = ra_create(tmp->_s + iarr,_ndim,gs,ge);
	      if (err) {
		RVLException e;
		e<<"Error: IWaveTimeStack::push_back(), copy dynamic array "<<iarr<<" into stack buffer\n";
		e<<" ra_create() return err "<<err<<" \n" ;
		throw e;
	      }
	    
	      ra_copy(tmp->_s + iarr, iwstate.model.ld_a._s + iarr);		
	    }
	  }
	  rdomvec.push_back(tmp);	
	  ic_itstack.push_back(it);
	}
	else {
	  /* case 2b: move incore stack forwards */ 
	  
	  /* store dyn-fields pointed by rdomvec.front() into outcore files, and reuse this space to store current dyn fields */    	      
	  int it0 = ic_itstack.front();
	  RDOM *tmp = rdomvec.front();
	    
	  for (int iarr=0; iarr < tmp->narr; iarr++) {
	     if (fd_isdyn(&(iwstate.model),iarr)) {
	       char * filename = iwave_get_a_checkfilename(iarr,it0-itorigin);
	      if (!filename) { 
		RVLException e;
		e<<"Error: IWaveTimeStack::push_back, filename is NULL for array" <<iarr<<" \n";
		throw e;		  
	      }
	      err = ra_fwrite(tmp->_s + iarr,filename);  
	      if (err) {
		RVLException e;
		e<<"Error: IWaveTimeStack::push_back()\n";
		e<<"ra_fwrite return err"<<err<<" for array "<<iarr<<" \n";
		throw e;  
	      }
	      if(filename) free(filename);
	      
	      /* store current dyn fields into the storage pointed by tmp */    	      
	      ra_copy(tmp->_s + iarr, iwstate.model.ld_a._s + iarr);
	      
	    }	     
	  }
	  
	  /* update incore/outcore stacks */
	  oc_itstack.push_back(it0);
	  ic_itstack.pop_front();
	  rdomvec.pop_front();
	  rdomvec.push_back(tmp);
	  ic_itstack.push_back(it);
	} /* case 2(b) ended */
	
      } /*case 2 ended */
      
    } // pair with if (iv==0)

#ifdef _IWSTACK_DEBUG
    write_IWStack();
#endif

  }

  /** place some state's time at the back of the stack */
  void IWaveTimeStack::push_back(IWaveState const & _iw){
    if (iscurrentState(_iw)) push_back();
    else {
      RVLException e;
      e<<"Error: IWaveTimeStack::push_back(_iw)\n";
      throw e;
    }
  }
  
  /** Pop the state element at the top of the stack */
  void IWaveTimeStack::pop_back(){
    /* pop_back in two cases:
       1. ic_size == 0: delete check files and update outcore stack;
       2. ic_size > 0: 
       2(a). oc_size == 0: incore stack pop_back and free corresponding memory; 
       2(b). oc_size > 0: move incore stack backwards, i.e., 
       move oc_back to incore,
       delete ic_back,
       update incore/outcore stack		   
    */
#ifdef _IWSTACK_DEBUG
    write_IWStack();
#endif

    int err = 0;
    IWAVE & iwstate = iw.IWaveState::getIWAVE();
    if(ic_itstack.size() == 0) { /* case 1 starts: */
      if(oc_itstack.size() > 0) {
	err = iwave_remove_checkfile(&iwstate, oc_itstack.back()-itorigin);
	if (err){ 
	  RVLException e;
	  e<< "Error: IWaveTimeStack::pop_back, iwave_remove_checkfile return err = "<< err<<"\n";
	  throw e;
	}
	oc_itstack.pop_back();
      }	
    } /* case 1 ended */
    else { /* case 2 starts */
      if (oc_itstack.size() == 0) { /* case 2(a): ic_pop_back and free memory */
	ic_itstack.pop_back();
	RDOM * tmp = rdomvec.back();
	rdomvec.pop_back();
	rd_a_destroy(tmp);
	delete tmp;
      }
      else { /* case 2(b): move incore stack backwards */
	RDOM * tmp = rdomvec.back();
	/* read corresponding check files and into incore buffer*/ 
	int ite = oc_itstack.back();	  
	for (int iarr=0; iarr< tmp->narr; iarr++) {
	  if (fd_isdyn(&(iwstate.model),iarr)) {
	    char* filename = iwave_get_a_checkfilename(iarr, ite - itorigin);
	    if (!filename) { 
	      RVLException e;
	      e<<"Error: IWaveTimeStack::pop_back, filename is NULL for array" <<iarr<<" \n";
	      throw e;		  
	    }
	    err = ra_fread(&(tmp->_s[iarr]), filename);
	    if (err){ 
	      RVLException e;
	      e<< "Error: IWaveTimeStack::pop_back, ra_fread for iarr = "<<iarr<<", return err ="<<err<<", for file "<<filename<<"\n";
	      throw e;
	    }
	    remove(filename);
	    if(filename) free(filename);
	  }
	}    
	
	/*update incore/outcore stacks*/
	oc_itstack.pop_back();
	ic_itstack.pop_back();
	rdomvec.pop_back();
	rdomvec.push_front(tmp);
	ic_itstack.push_front(ite);
      } /* case 2(b) ended */
      
    } /* case 2 ended */

#ifdef _IWSTACK_DEBUG    
    write_IWStack();
#endif

  }
  
  /** Allows access to a specific state at time associated with index idx */
  void IWaveTimeStack::getAt(IWaveState & _iw, int idx){
#ifdef _IWSTACK_DEBUG
    write_IWStack();
#endif
    if ( (idx >= this->size()) || (idx<0)) {
      RVLException e;
      e<<"Error: IWaveTimeStack::getAt(IWaveState &, int)\n";
      e<<"idx = "<<idx<<" out of range [0,"<<this->size()-1<<"]\n";
      throw e;
    }
    
    IWAVE & iwstate = _iw.IWaveState::getIWAVE();      
    int it;

    int ic_idx = idx - oc_itstack.size();
    if(ic_idx >= 0) {
      /* get data from incore buffer */
      it = ic_itstack.at(ic_idx);
      RDOM * tmp = rdomvec[ic_idx];
      for (int iarr=0; iarr < tmp->narr; iarr++) {
	if (fd_isdyn(&(iwstate.model),iarr)) {
	  /* copy dynamic arrays of tmp to ld_a*/
	  int _ndim = iwstate.model.ld_a._s[iarr].ndim;
	  long size = 1L;
	  for (int d =0; d<_ndim; ++d) {
	    size *= (long)(tmp->_s->_dims[d].n0);
	  }
	  if (size > 0L) {
	    if (tmp->_s[iarr]._s != NULL){
	      ra_copy(iwstate.model.ld_a._s + iarr, tmp->_s + iarr); 
	    }
	  }
	}
      }
    }
    else {
      /* get data from outcore files */
      it = oc_itstack.at(idx);
      int err = 0;
      err = giwave_dynamic_init(&iwstate,it,it-itorigin);
      if (err) {
	RVLException e;
	e<<"Error: IWaveTimeStack::getAt, giwave_dynamic_init return error \n";
	throw e;
      }
    }
    /* transfer check time to state */
    iwstate.model.tsind.it=it;
    iwstate.model.tsind.iv=0;    

#ifdef _IWSTACK_DEBUG
    write_IWStack();
#endif
  }
  
   
  IWaveState & IWaveTimeStack::at(int idx){
    getAt(iw,idx); 
    return iw;
  }
  
  /** Returns reference at the head of the stack */
  IWaveState & IWaveTimeStack::front(){
    if(oc_itstack.size() > 0) 
      return at(oc_itstack.front());
    else 
      return at(ic_itstack.front());
  }
    
  /** Returns reference at the tail of the stack */
  IWaveState & IWaveTimeStack::back(){
    if(ic_itstack.size() > 0) 
      return at(ic_itstack.back());
    else 
      return at(oc_itstack.back());
  }
  
}




