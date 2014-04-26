#ifndef __IWAVEPP_STACK
#define __IWAVEPP_STACK

#include "state.hh"
#include "stackBase.hh"

//#define _IWSTACK_DEBUG

namespace TSOpt{
  
  using RVL::RVLException;

  /** IWaveTimeStack */
  class IWaveTimeStack : public stackBase<IWaveState>{

  private:

    IWaveState & iw;          /* refer to the state object in sim */
    int itorigin;             /* starting time index */
    const int max_ic_size;    /* maximum incore buffer size */
    
    /* it_stack consists of (outcore_itstack + incore_itstack) */
    std::deque<int> ic_itstack;  /* incore it stack */
    std::vector<int> oc_itstack; /* outcore it stack */
    std::deque<RDOM *> rdomvec;  /* incore RDOM buffer stack */
    
    IWaveTimeStack();
    IWaveTimeStack(IWaveTimeStack const &);
  
  public: 

    IWaveTimeStack(IWaveState & _iw, int _startit = 0, int _max_incore_size = 0): 
      iw(_iw), itorigin(_startit), max_ic_size(_max_incore_size),  ic_itstack(), oc_itstack(){
#ifdef _IWSTACK_DEBUG
      write_IWStack();
#endif
    }
    
    ~IWaveTimeStack(){
      while(this->size()) pop_back();
    }
    
    /** set starting time index */
    void setOrigin(int _itorigin) { itorigin = _itorigin;}

    /** get starting time index*/
    int getOrigin() { return itorigin;}

    /** check */
    bool iscurrentState(IWaveState const & _iw); 
    
    /** place current state's time at the back of the stack */
    /* push_back in two cases:
       1. max_ic_size <= 0: store everything outcore; 
       2. max_ic_size > 0: 
         2(a). ic_size < max_ic_size: store state incore;
         2(b). ic_size == max_ic_size: move incore stack forwards, i.e., 
             move ic_front to outcore, 
             update incore/outcore stack,
             push back cur-state into incore stack           
    */
    void push_back(); 


    /** place some state's time at the back of the stack */
    void push_back(IWaveState const & _iw);
      
    /** Pop the state element at the top of the stack */
    /* pop_back in two cases:
       1. ic_size == 0: delete check files and update outcore stack;
       2. ic_size > 0: 
          2(a). oc_size == 0: incore stack pop_back and free corresponding memory; 
          2(b). oc_size > 0: move incore stack backwards, i.e., 
                  move oc_back to incore,
                  delete ic_back,
                  update incore/outcore stack		   
    */
    void pop_back();

    /** Returns stack size */
    int size(){
      return ((int)ic_itstack.size() + (int)oc_itstack.size());
    } 
    
    /** Allows access to a specific state at time associated with index idx */
    void getAt(IWaveState & _iw, int idx);
   
    IWaveState & at(int idx);
    
    /** Returns reference at the head of the stack */
    IWaveState & front();
    
    /** Returns reference at the tail of the stack */
    IWaveState & back();
    
    /** Debugging func to tell if TimeStack behaves correctly */
#ifdef _IWSTACK_DEBUG    
    void write_IWStack() {
      bool is_right;
      
      if (oc_itstack.size()>0)
	is_right = (ic_itstack.size() == max_ic_size);
      else 
	is_right = (ic_itstack.size() <= max_ic_size);
      
      if(!is_right) {
	cerr<<"--------- IWaveTimeStack::write_IWStack -------------\n";
	cerr<<"---- Unexpected IWStack : ic_itstack.size = "<<ic_itstack.size()<<", oc_itstack.size = "<<oc_itstack.size()<<", stack_size = "<<size()<<", max_ic_size = "<<max_ic_size<<" ---- \n";
      }	
    }
#endif

  };
}

#endif
