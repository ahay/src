#include "srcstk.hh"

namespace TSOpt {

  SrcData::SrcData(string _line, 
		   string _csg, 
		   string _crg, 
		   string _time, 
		   float tmute)
    : nt(0), ntr(0), nsrc(0), kmute(0), isfix(true), 
      line(_line), csg(_csg), crg(_crg), time(_time),
      fpl(NULL), fpg(NULL), fpr(NULL), fpt(NULL), 
      nsmod(0), dtmod(0), t0mod(0) {

    segy tr;
    int sz,sx,sy;
    int itr=0;
    int ir=0;
    
    int gx[MAX_TRACES];
    int gy[MAX_TRACES];
    int gz[MAX_TRACES];
      
    ntr=0;
    nsrc=0;

    if (time.size()) {
      // no prototype
      //if (!(fpt=fopen(time.c_str(),"r"))) {
      if (!(fpt=iwave_const_fopen(time.c_str(),"r",NULL,stderr))) {
	RVLException e;
	e<<"Error: SrcData constructor\n";
	e<<"failed to open time data file "<<time<<"\n";
	throw e;
      }
      if (!(fgettr(fpt,&tr))) {
	RVLException e;
	e<<"Error: SrcData constructor\n";
	e<<"failed to read trace from time data file "<<time<<"\n";
	throw e;
      }      
      nsmod=tr.ns;
      dtmod=tr.dt;
      t0mod=tr.delrt;
    }
      
    if (!(fpl=iwave_const_fopen(line.c_str(),"r",NULL,stderr))) {
      RVLException e;
      e<<"Error: SrcData constructor\n";
      e<<"failed to open source file "<<line<<"\n";
      throw e;
    }
    if (csg.size()) {
      if (!(fpg=iwave_const_fopen(csg.c_str(),"w",NULL,stderr))) {
	RVLException e;
	e<<"Error: SrcData constructor\n";
	e<<"failed to open target CSG file "<<csg<<"\n";
	throw e;
      }
    }
    if (crg.size()) {
      if (!(fpr=iwave_const_fopen(crg.c_str(),"w",NULL,stderr))) {
	RVLException e;
	e<<"Error: SrcData constructor\n";
	e<<"failed to open target CRG file "<<crg<<"\n";
	throw e;
      }
    }
    /* find first source coord vector on first trace */
    if (fgettr(fpl,&tr)) {
      if (fpg) fputtr(fpg,&tr);
      itr=1;
      nt=tr.ns;
      sz=tr.selev;
      sx=tr.sx;
      sy=tr.sy;
      gz[itr-1]=tr.gelev;
      gx[itr-1]=tr.gx;
      gy[itr-1]=tr.gy;
      // calculate kmute
      kmute=(int)(1000*tmute/((float)tr.dt));

      /* 010510: require that sampling of source and seismic data
	 is same - in practice, this means resampling source files
	 offline */

      if ((time.size()) && 
	  ((nsmod != nt) ||
	   (dtmod != tr.dt) ||
	   (t0mod != tr.delrt))) {
	RVLException e;
	e<<"Error: SrcData constructor\n";
	e<<"source sampling (file = "<<time<<") differs from data sampling (file = "<<line<<")\n";
	throw e;
      }


      ntr++;
      nsrc++;

      // now overwrite receiver data

      if (fpr) {
	tr.sx=0;
	tr.sy=0;
	tr.selev=0;
	tr.tracl=nsrc;
	tr.tracr=nsrc;
	tr.tracf=nsrc;
	tr.gelev=sz;
	tr.gx=sx;
	tr.gy=sy;
	if (fpt) {
	  tr.ns=nsmod;
	  tr.dt=dtmod;
	  tr.delrt=t0mod;
	}
	fputtr(fpr,&tr);
      }

      //      cerr<<"itr="<<itr<<" ntr="<<ntr<<" nsrc="<<nsrc<<endl;
    }
    /* read through the rest of first csg, determine ntr */
    while ((ir=fgettr(fpl,&tr)) && 
	   (sz==tr.selev) && 
	   (sx==tr.sx) && 
	   (sy==tr.sy)) {
      if (fpg) fputtr(fpg,&tr);
      itr++;
      ntr++;
      gz[itr-1]=tr.gelev;
      gx[itr-1]=tr.gx;
      gy[itr-1]=tr.gy;
      //      cerr<<"itr="<<itr<<" ntr="<<ntr<<" nsrc="<<nsrc<<endl;
    }      
    
    /* the first trace of the second csg has been read */
    if (ir) {
      itr=1;
      sz=tr.selev;
      sx=tr.sx;
      sy=tr.sy;
      if ((gz[itr-1]!=tr.gelev) ||
	  (gx[itr-1]!=tr.gx) ||
	  (gy[itr-1]!=tr.gy)) isfix=false;
      nsrc++;
      // now overwrite receiver data
      if (fpr) {
	tr.sx=0;
	tr.sy=0;
	tr.selev=0;
	tr.tracl=nsrc;
	tr.tracr=nsrc;
	tr.tracf=nsrc;
	tr.gelev=sz;
	tr.gx=sx;
	tr.gy=sy;
	if (fpt) {
	  tr.ns=nsmod;
	  tr.dt=dtmod;
	  tr.delrt=t0mod;
	}
	fputtr(fpr,&tr);
      }
      //      cerr<<"itr="<<itr<<" ntr="<<ntr<<" nsrc="<<nsrc<<endl;
      
      while (ir) {
	ir=fgettr(fpl,&tr);
	if (ir) {
	  if ((sz==tr.selev) && 
	      (sx==tr.sx) && 
	      (sy==tr.sy)) {
	    itr++;
	    if ((gz[itr-1]!=tr.gelev) ||
		(gx[itr-1]!=tr.gx) ||
		(gy[itr-1]!=tr.gy)) isfix=false;
	    //	    cerr<<"itr="<<itr<<" ntr="<<ntr<<" nsrc="<<nsrc<<endl;    
	  }
	  else {
	    if (ntr!=itr) {
	      RVLException e;
	      e<<"Error: SrcData constructor\n";
	      e<<"csg "<<nsrc<<" has "<<itr<<" traces, which"<<"\n";
	      e<<"differs from the "<<ntr<<" traces in csg 1\n";
	      throw e;
	    }
	    nsrc++;
	    itr=1;
	    sz=tr.selev;
	    sx=tr.sx;
	    sy=tr.sy;
	    if ((gz[itr-1]!=tr.gelev) ||
		(gx[itr-1]!=tr.gx) ||
		(gy[itr-1]!=tr.gy)) isfix=false;
	    //	    cerr<<"itr="<<itr<<" ntr="<<ntr<<" nsrc="<<nsrc<<endl;
	    // now overwrite receiver data
	    if (fpr) {
	      tr.sx=0;
	      tr.sy=0;
	      tr.selev=0;
	      tr.tracl=nsrc;
	      tr.tracr=nsrc;
	      tr.tracf=nsrc;
	      tr.gelev=sz;
	      tr.gx=sx;
	      tr.gy=sy;
	      if (fpt) {
		tr.ns=nsmod;
		tr.dt=dtmod;
		tr.delrt=t0mod;
	      }
	      fputtr(fpr,&tr);
	    }
	  }
	}
      }
      
    }

    if (fpg) fflush(fpg);
    if (fpr) fflush(fpr);

    iwave_fclose(fpl);
    if (fpr) iwave_fclose(fpr);
    if (fpg) iwave_fclose(fpg);
    if (fpt) iwave_fclose(fpt);

  }

  SrcData::~SrcData() {  }

  ostream & SrcData::write(ostream & str) const {
    str<<"SrcData object: src, rec gathers\n";
    str<<"nt     (csg) = "<<nt<<"\n";
    str<<"ntr    (csg) = "<<ntr<<"\n";
    str<<"nsrc         = "<<nsrc<<"\n";
    str<<"kmute        = "<<kmute<<"\n";
    str<<"fixed spread = "<<isfix<<"\n";
    str<<"line source  = "<<line<<"\n";
    str<<"CSG headers  = "<<csg<<"\n";
    str<<"CRG headers  = "<<crg<<"\n";
    str<<"nt     (crg) = "<<nsmod<<"\n";
    str<<"dt (ms)(crg) = "<<0.001*(float)(dtmod)<<"\n";
    str<<"t0     (crg) = "<<t0mod<<"\n";
    return str;
  }

  void SrcOp::apply(Vector<float> const & x,
		    Vector<float> & y) const {
    try {

      FILE * fpg;

      /* check domain and range */
      if (x.getSpace() != dom) {
	RVLException e;
	e<<"Error: SrcOp::applyOp\n";
	e<<"input vector not member of domain\n";
	throw e;
      }
      if (y.getSpace() != rng) {
	RVLException e;
	e<<"Error: SrcOp::applyOp\n";
	e<<"output vector not member of range\n";
	throw e;
      }

      PARARRAY par;
      ps_setnull(&par);

      AssignParams apg(par,"gather");
      y.eval(apg);

      //      ps_printall(par,stderr);

      LocalVector<float> lx(x);
      lx.copy(x);

      //      cerr<<"\n****** applyOp: x="<<endl;
      //      lx.write(cerr);

      char * gstr;
      if (ps_ffcstring(par,"gather",&gstr)) {
	RVLException e;
	e<<"Error: SrcOp::apply\n";
	e<<"failed to extract gather (output) filename\n";
	throw e;
      }

      ps_destroy(&par);

      segy trl;
      
      segy trg; /* workspace for line input */

      float * bufg = new float[sd.get_ntr()*sd.get_nt()];     
      for (int i=0;i<sd.get_ntr()*sd.get_nt();i++) bufg[i]=0.0;

      if (fseeko(fpl,0L,SEEK_SET)) {
	RVLException e;
	e<<"Error: SrcOp::apply\n";
	e<<"failed to seek to start of line file\n";
	throw e;
      }
    
      //      cerr<<"\n***** applyOp: before loop, nsrc="<<sd.get_nsrc()<<" ntr="<<sd.get_ntr()<<" nt="<<sd.get_nt()<<endl;

      /* accumulation loop */
      for (int j=0;j<sd.get_nsrc();j++) {

	//	cerr<<"***** applyOp: src="<<j<<" lx="<<lx.getData()[j]<<endl;
	for (int i=0;i<sd.get_ntr();i++) {

	  //	  cerr<<"trc="<<i<<endl;
	  //	  cerr<<"read line\n";
	  if (!(fgettr(fpl,&trl))) {
	    RVLException e;
	    e<<"Error: SrcOp::apply\n";
	    e<<"failed to read trace "<<i+j*sd.get_ntr()<<" of output file\n";
	    throw e;
	  }      
	  //	  cerr<<"addemup\n";
	  //	  float jnk=0.0f;
	  for (int k=0;k<sd.get_nt();k++) {
	    bufg[k+i*sd.get_nt()]+=lx.getData()[j]*trl.data[k];
	    //	    jnk+=bufg[k+i*sd.get_nt()]*bufg[k+i*sd.get_nt()];
	  }
	  //	  cerr<<"------ trc="<<i<<"ms accum="<<jnk<<"\n";

	}
      }  
      
      //      fprintf(stderr,"open gather file %s\n",gstr);
      if (!(fpg=iwave_const_fopen(gstr,"w",sd.getCSG().c_str(),stderr))) {
	RVLException e;
	e<<"Error: SrcOp::apply\n";
	e<<"failed to open gather (output) file\n";
	throw e;
      }    

      fseeko(fpl,0L,SEEK_SET);
      fseeko(fpg,0L,SEEK_SET);

      //      cerr<<"writeback loop\n";
      for (int i=0;i<sd.get_ntr();i++) {
	//	cerr<<"read line trace "<<i<<endl;
	if (!(fgettr(fpl,&trg))) {
	  RVLException e;
	  e<<"Error: SrcOp::apply\n";
	  e<<"failed to read trace "<<i<<" of output file\n";
	  throw e;
	}      	
	//	cerr<<"addemup\n";
	for (int k=0;k<sd.get_kmute();k++) 
	  trg.data[k]=0.0;
	for (int k=sd.get_kmute();k<sd.get_nt();k++) 
	  trg.data[k]=bufg[k+i*sd.get_nt()];
	//	cerr<<"write to gather trace "<<i<<endl;
	fputtr(fpg,&trg);
      }


      fflush(fpg);
      iwave_fclose(fpg);

    }
    catch (RVLException e) {
      e<<"\ncalled from SrcOp::applyOp\n";
      throw e;
    }
  }

  void SrcOp::applyAdj(Vector<float> const & x,
		       Vector<float> & y) const {
    try {

      FILE * fp;

      /* check domain and range */
      if (y.getSpace() != dom) {
	RVLException e;
	e<<"Error: SrcOp::applyOp\n";
	e<<"output vector not member of domain\n";
	throw e;
      }
      if (x.getSpace() != rng) {
	RVLException e;
	e<<"Error: SrcOp::applyOp\n";
	e<<"input vector not member of range\n";
	throw e;
      }

      segy tr; /* workspace for line input */

      //      float * bufl = new float[sd.get_ntr()*sd.get_nt()];
      float * bufg = new float[sd.get_ntr()*sd.get_nt()];
      for (int i=0;i<sd.get_ntr()*sd.get_nt();i++) {
	//	bufl[i]=0.0;
	bufg[i]=0.0;
      }

      PARARRAY par;
      ps_setnull(&par);

      AssignParams apg(par,"gather");
      x.eval(apg);

      LocalVector<float> ly(y);
      ly.zero();

      char * gstr;
      if (ps_ffcstring(par,"gather",&gstr)) {
	RVLException e;
	e<<"Error: SrcOp::applyAdj\n";
	e<<"failed to extract gather (input) filename\n";
	throw e;
      }

      ps_destroy(&par);

      string xxx=gstr;

      if (!(fp=iwave_const_fopen(gstr,"r",sd.getCSG().c_str(),stderr))) {
	RVLException e;
	e<<"Error: SrcOp::applyAdj\n";
	e<<"failed to open gather (input) file "<<xxx<<"\n";
	throw e;
      }    

      fseeko(fp,0L,SEEK_SET);
    
      for (int i=0;i<sd.get_ntr();i++) {
	if (!(fgettr(fp,&tr))) {
	  RVLException e;
	  e<<"Error: SrcOp::applyAdj\n";
	  e<<"failed to read trace "<<i<<" of input file "<<xxx<<"\n";
	  throw e;
	}	
	for (int k=sd.get_kmute();k<sd.get_nt();k++)
	  bufg[k+i*sd.get_nt()]=tr.data[k];
      }
      iwave_fclose(fp);

      SEGYSpace const & xsp = dynamic_cast<SEGYSpace const &>(x.getSpace());
      SEGYDCF const & f = dynamic_cast<SEGYDCF const &>(xsp.getDCF());
      float dt=f.getDt();

      fseeko(fpl,0L,SEEK_SET);

      for (int j=0;j<sd.get_nsrc();j++) {
	for (int i=0;i<sd.get_ntr();i++) {
	  if (!(fgettr(fpl,&tr))) {
	    RVLException e;
	    e<<"Error: SrcOp::applyAdjOp\n";
	    e<<"failed to read trace "<<i+j*sd.get_ntr()<<" of line data - file = "<<sd.getLine()<<"\n";
	    throw e;
	  }		
	  for (int k=0;k<sd.get_nt();k++)
	    //	    bufl[k+i*sd.get_nt()]=tr.data[k];
	    //	}
	    //	for (int i=0;i<sd.get_ntr();i++) {
	    //	  for (int k=0;k<sd.get_nt();k++) {
	    //	    ly.getData()[j]+=bufg[k+i*sd.get_nt()]*bufl[k+i*sd.get_nt()];
	    ly.getData()[j]+=bufg[k+i*sd.get_nt()]*tr.data[k];
	
	  //	  cerr<<"***** applyAdjOp: src="<<j<<" trc="<<i<<" dot="<<ly.getData()[j]<<endl;
	}
	ly.getData()[j]*=dt;
	//	cerr<<"***** applyAdjOp: src="<<j<<"lx="<<ly.getData()[j]<<endl;
      }

      iwave_fclose(fp);

      /* copy and scale by dt to match IP in SEGYSpace */
      y.copy(ly);
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: SrcOp::applyAdjOp - bad cast\n";
      e<<"somehow input space is not SEGYSpace\n";
      throw e;
    }
    catch (RVLException e) {
      e<<"\ncalled from SrcOp::applyAdjOp\n";
      throw e;
    }
  }

  ostream & SrcOp::write(ostream & str) const {
    str<<"Data Synthesis Operator\n";
    str<<"number of time samples = "<<sd.get_nt()<<"\n";
    str<<"number of traces       = "<<sd.get_ntr()<<"\n";
    str<<"number of gathers      = "<<sd.get_nsrc()<<"\n";
    str<<"line filename          = "<<sd.getLine()<<"\n";
    str<<"CSG filename           = "<<sd.getCSG()<<"\n";
    str<<"CRG filename           = "<<sd.getCRG()<<"\n";
    return str;
  }
}
