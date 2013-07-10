#include "srcwdippw.hh"
#define VERBOSE_PSrc
namespace TSOpt {

  SrcWDipPwOp::SrcWDipPwOp(string wavelet,
			   SrcData const & _sd,
			   float _dip) 
    : sd(_sd),
      rng(sd.getCRG()),
      dom(sd.get_nsrc()),
      dip(_dip) {

    if (retrieveRank()==0) {

      FILE * fpw; /* file for wavelet */

      /* open wavelet file*/    
      char const * srcname= NULL;
      if (sd.getSrc().size()) srcname=sd.getSrc().c_str();
      if (!(fpw=iwave_const_fopen(wavelet.c_str(),"r",srcname,stderr))) {
	RVLException e;
	e<<"Error: SrcWPW constructor\n";
	e<<"failed to open wavelet file "<<wavelet<<"\n";
	throw e;
      }

      /* read it */
      if (!fgettr(fpw,&w)) {
	RVLException e;
	e<<"Error: SrcWPW constructor\n";
	e<<"failed to read wavelet trace on file "<<wavelet<<"\n";
	throw e;
      }

      /* sanity check */
#ifdef VERBOSE_PSrc
      cerr<<"SrcWDipPwOp constructor\n";
      cerr<<"wavelet nt = "<<w.ns<<" sd nt = "<<sd.get_ntsrc()<<"\n";
      cerr<<"wavelet dt = "<<w.dt<<" sd dt = "<<sd.get_dtsrc()<<"\n";
      cerr<<"wavelet t0 = "<<w.delrt<<" sd t0 = "<<sd.get_t0src()<<"\n";
#endif
      /* Note: wavelet_nt may not be the same as sd_nt;
	 to accommodate plane-wave with any angle, sd_nt depends on the dipping angle, i.e.,
	 sd_nt_dip = sd_nt_0 + ceil(fabs((sd.get_ntr()-1)*dip)/sd_dt)
       */
      if ((w.dt    != sd.get_dtsrc()) ||
	  (w.delrt != sd.get_t0src())) {
	RVLException e;
	e<<"Error: SrcWDipPwOp constructor\n";
	e<<"wavelet "<<wavelet<<" time data incompatible with SrcData input\n";
	e<<"wavelet dt = "<<w.dt<<" sd dt = "<<sd.get_dtsrc()<<"\n";
	e<<"wavelet t0 = "<<w.delrt<<" sd t0 = "<<sd.get_t0src()<<"\n";
	throw e;
      }
      
      iwave_fclose(fpw);
    }
    
  }
  
  void SrcWDipPwOp::apply(Vector<float> const & x,
		          Vector<float> & y) const {
    try {
      
      if (retrieveRank()==0) {
	
	/* check domain and range */
	if (x.getSpace() != dom) {
	  RVLException e;
	  e<<"Error: SrcWDipPwOp::applyOp\n";
	  e<<"input vector not member of domain\n";
	  throw e;
	}
	if (y.getSpace() != rng) {
	  RVLException e;
	  e<<"Error: SrcWDipPwOp::applyOp\n";
	  e<<"output vector not member of range\n";
	  throw e;
	}

	LocalVector<float> lx(x);

	// the next block of code extracts the filename of the 
	// target SEGYSpace vector

	PARARRAY par;
	ps_setnull(&par);

	AssignParams apg(par,"CRG");
	y.eval(apg);

	//      ps_printall(par,stderr);

	char * gstr;
	if (ps_ffcstring(par,"CRG",&gstr)) {
	  RVLException e;
	  e<<"Error: SrcWDipPwOp::apply\n";
	  e<<"failed to extract gather (output) filename\n";
	  throw e;
	}

	ps_destroy(&par);

	// open header file from SEGYSpace
	FILE * fph = NULL;
	SEGYDCF const & dcf = dynamic_cast< SEGYDCF const &>(rng.getDCF());
	if (!(fph=iwave_const_fopen((dcf.getFilename()).c_str(),"r",NULL,stderr))) {
	  RVLException e;
	  e<<"Error: SrcWDipPwOp::apply\n";
	  e<<"failed to open header filename of space = "<<dcf.getFilename()<<"\n";
	  throw e;
	}

	// open target file 
	FILE * fpg = NULL;
	if (!(fpg=iwave_const_fopen(gstr,"w",sd.getCRG().c_str(),stderr))) {
	  RVLException e;
	  e<<"Error: SrcWPWOp::apply\n";
	  e<<"failed to open output source filename = "<<gstr<<"\n";
	  throw e;
	}	
	if (fseeko(fpg,0L,SEEK_SET)) {
	  RVLException e;
	  e<<"Error: SrcWPWOp::apply\n";
	  e<<"failed to reset output source filename = "<<gstr<<"\n";
	  throw e;
	}	

	segy tr;
	int itr = 0;
	int etr = 0;
	if(dip < 0.0)  etr = sd.get_ntr()-1;
	
	// cerr<<"SrcWDipPwOp::applyOp, etr ="<<etr<<endl;
	
 	float msdt = 0.001*sd.get_dtsrc();
        // cerr<<"msdt ="<<msdt<<endl; 
	while (fgettr(fph,&tr)) {
	  float fit0 = (itr - etr)*dip/msdt;
	  int itl0 = static_cast<int>(fit0);
	  float alpha = fit0 - itl0;
	  /* initialize the trace */
	  for (int i=0; i < tr.ns; i++) tr.data[i]=0.0;
	  /* generate trace from wavelet w and weight vector lx*/
	  for (int it=0; it < tr.ns; it++) {
	    int itless = it - itl0;
	    if (itless < 0) tr.data[it]=0.0;
	    else if (itless < w.ns){      
	      /* linear interpolation */
	      tr.data[it] += lx.getData()[itr]*w.data[itless]*(1.0 - alpha);	 
	      if (it < tr.ns - 1 )
		tr.data[it+1] += lx.getData()[itr]*w.data[itless]*alpha;		 
	    }
	    else tr.data[it] = 0.0;
	    //cerr<<"itr ="<<itr<<" tr["<<it<<"] ="<<tr.data[it]<<"  "; 
	  } 
	  //cerr<<"itr ="<<itr<<", itl0 ="<<itl0<<", alpha ="<<alpha<<", tr[itl0] ="<<tr.data[itl0]<<endl;
	  
	  fputtr(fpg,&tr);
	  itr++;
	}
	fflush(fpg);
	
	iwave_fclose(fph);
	iwave_fclose(fpg);
      }
      
    }
    catch (RVLException & e) {
      e<<"\ncalled from SrcWDipPwOp::applyOp\n";
      throw e;
    }
  }
	
  void SrcWDipPwOp::applyAdj(Vector<float> const & y,
			     Vector<float> & x) const {
  
    try {
      
      if (retrieveRank()==0) {
	
	/* check domain and range */
	if (x.getSpace() != dom) {
	  RVLException e;
	  e<<"Error: SrcWPWOp::applyAdjOp\n";
	  e<<"output vector not member of domain\n";
	  throw e;
	}
	if (y.getSpace() != rng) {
	  RVLException e;
	  e<<"Error: SrcWPWOp::applyAdjOp\n";
	  e<<"input vector not member of range\n";
	  throw e;
	}

	LocalVector<float> lx(x);
	lx.zero();

	// the next block of code extracts the filename of the 
	// target SEGYSpace vector

	PARARRAY par;
	ps_setnull(&par);

	AssignParams apg(par,"CRG");
	y.eval(apg);

	//      ps_printall(par,stderr);

	char * gstr;
	if (ps_ffcstring(par,"CRG",&gstr)) {
	  RVLException e;
	  e<<"Error: SrcWPWOp::apply\n";
	  e<<"failed to extract gather (output) filename\n";
	  throw e;
	}

	ps_destroy(&par);

	// open target file 
	FILE * fpg = NULL;
	if (!(fpg=iwave_const_fopen(gstr,"r",sd.getCRG().c_str(),stderr))) {
	  RVLException e;
	  e<<"Error: SrcWPWOp::applyAdj\n";
	  e<<"failed to open input source filename = "<<gstr<<"\n";
	  throw e;
	}	

	// reset file
	if (fseeko(fpg,0L,SEEK_SET)) {
	  RVLException e;
	  e<<"Error: SrcWPWOp::apply\n";
	  e<<"failed to reset input source filename = "<<gstr<<"\n";
	  throw e;
	}	

	segy tr;
	int itr=0;

	int etr = 0;
	if(dip < 0.0)  etr = sd.get_ntr()-1;
 	float msdt = 0.001*sd.get_dtsrc();

	while (fgettr(fpg,&tr)) {
	  float fit0 = (itr - etr)*dip/msdt;
	  int itl0 = static_cast<int>(fit0);
	  float alpha = fit0 - itl0;
	  lx.getData()[itr] = 0.0;
	  for (int it=0;it<tr.ns;it++){
	    int itless = it - itl0;
	    if (itless >= 0 && itless < w.ns){
	      lx.getData()[itr] += tr.data[it] * w.data[itless] * (1.0 - alpha);
	      if (it < tr.ns - 1)
		lx.getData()[itr] += tr.data[it + 1] * w.data[itless] * alpha;
	    }
	  }
	  lx.getData()[itr]*=rng.getDt();
	  //	cerr<<" lx["<<itr<<"]="<<lx.getData()[itr]<<endl;
	  itr++;
	}

	/* copy and scale output */
	x.copy(lx);
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from SrcWDipPwOp::applyAdjOp\n";
      throw e;
    }

  }

  ostream & SrcWDipPwOp::write(ostream & str) const {
    str<<"Weighted Dipping plane wave source op\n";
    return str;
  }

}

