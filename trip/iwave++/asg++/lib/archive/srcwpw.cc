#include "srcwpw.hh"

namespace TSOpt {

  SrcWPWOp::SrcWPWOp(string wavelet,
		     SrcData const & _sd) 
    : sd(_sd),
      rng(sd.getCRG()),
      dom(sd.get_nsrc()) {

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
      if ((w.ns    != sd.get_ntsrc()) ||
	  (w.dt    != sd.get_dtsrc()) ||
	  (w.delrt != sd.get_t0src())) {
	RVLException e;
	e<<"Error: SrcWPW constructor\n";
	e<<"wavelet "<<wavelet<<" time data incompatible with SrcData input\n";
	e<<"wavelet nt = "<<w.ns<<" sd nt = "<<sd.get_ntsrc()<<"\n";
	e<<"wavelet dt = "<<w.dt<<" sd dt = "<<sd.get_dtsrc()<<"\n";
	e<<"wavelet t0 = "<<w.delrt<<" sd t0 = "<<sd.get_t0src()<<"\n";
	throw e;
      }
      iwave_fclose(fpw);
    }
    
  }
    
  void SrcWPWOp::apply(Vector<float> const & x,
		       Vector<float> & y) const {
    try {

      if (retrieveRank()==0) {

	/* check domain and range */
	if (x.getSpace() != dom) {
	  RVLException e;
	  e<<"Error: SrcWPWOp::applyOp\n";
	  e<<"input vector not member of domain\n";
	  throw e;
	}
	if (y.getSpace() != rng) {
	  RVLException e;
	  e<<"Error: SrcWPWOp::applyOp\n";
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
	  e<<"Error: SrcWPWOp::apply\n";
	  e<<"failed to extract gather (output) filename\n";
	  throw e;
	}

	ps_destroy(&par);

	// open header file from SEGYSpace
	FILE * fph = NULL;
	SEGYDCF const & dcf = dynamic_cast< SEGYDCF const &>(rng.getDCF());
	if (!(fph=iwave_const_fopen((dcf.getFilename()).c_str(),"r",NULL,stderr))) {
	  RVLException e;
	  e<<"Error: SrcWPWOp::apply\n";
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
	int itr=0;
      
	while (fgettr(fph,&tr)) {
	  for (int it=0;it<tr.ns;it++) 
	    tr.data[it]=lx.getData()[itr]*w.data[it];
	  fputtr(fpg,&tr);
	  itr++;
	}
	fflush(fpg);

	iwave_fclose(fph);
	iwave_fclose(fpg);
      }

    }
    catch (RVLException & e) {
      e<<"\ncalled from SrcWPW::applyOp\n";
      throw e;
    }
  }
	
  void SrcWPWOp::applyAdj(Vector<float> const & y,
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

	while (fgettr(fpg,&tr)) {
	  //	cerr<<"tr.ns = "<<tr.ns<<" w.ns = "<<w.ns;
	  for (int it=0;it<tr.ns;it++)
	    lx.getData()[itr]+=tr.data[it]*w.data[it];
	  lx.getData()[itr]*=rng.getDt();
	  //	cerr<<" lx["<<itr<<"]="<<lx.getData()[itr]<<endl;
	  itr++;
	}

	/* copy and scale output */
	x.copy(lx);
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from SrcWPW::applyOp\n";
      throw e;
    }

  }

  ostream & SrcWPWOp::write(ostream & str) const {
    str<<"Weighted plane wave source op\n";
    return str;
  }

}

