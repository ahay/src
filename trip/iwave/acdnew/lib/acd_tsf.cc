#include "acd.hh"

// defined in fd_acd.cc
int acd_step(RDOM *, int, void*);

// interface to IWaveSim::run() 
void acd_timestep(std::vector<RDOM *> iw, bool fwd, int iv, void *fdpars){
  int n  = iw.size();

  if ( iv != 0 ){
    RVLException e;
    e<<"Error: acd_tsf(). ACD is a single step method. Input iv = " << n << " should be 0.\n";
    throw e;
  }
  switch (n){
  case 1:{
    if(fwd == true) {
      acd_step(iw[0], iv, fdpars);
      return;
    }
    if(fwd == false){
      RVLException e;
      e<<"Error: acd_timestep().  iw.size() = "<< n << " has no adjoint!";
      throw e;
    } 
  }
  case 2:{
    if(fwd == true) {
      acd_tsfm(iw[1], iw[0], iv, fdpars);
      return;
    }
    else {
      acd_tsam(iw[1], iw[0], iv, fdpars);
      return;
    }
  }
  case 4:{
    if(fwd == true) {
      acd_tsfm2(iw[3], iw[2], iw[1], iw[0],
		iv, fdpars);
      return;
    }
    else {
      acd_tsam2(iw[3], iw[2], iw[1], iw[0],
		iv, fdpars);
      return;
    }
            
  }
  default:
    RVLException e;
    e<<"Error: acd_tsf(). Only zero, 1st & 2nd derivatives are implemented.\n";
    throw e;
  }
}

