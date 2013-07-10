#ifndef __IWAVE_INV_OP
#define __IWAVE_INV_OP

//#include "umin.hh"
//#include "linearsolver.hh"
#include "LBFGSBT.hh"
#include "cgalg.hh"

#include "iwop.hh"
#include "ls.hh"
#include "blockop.hh"

#include "gridops.hh"
#include "griddiffops.hh"

#include "segyops.hh"
#include "localspace.hh"

#include "fftwsegyops.hh"

/** Operator for constructing the mapping from data points to models via solving least-squares optimization (inversion) driven by simulation (wrapped in IWaveOp).
 *  In fact, this mapping is defined by the implicit equation derived from the first order necessary condition of LS optimization 
*/

namespace TSOpt{
  using namespace RVLUmin;
  using namespace RVLAlg;
  using namespace TSOpt;
  using namespace RVL;
   
  template
  <
    class FwdSamplerPolicy,
    class LinSamplerPolicy,
    class AdjSamplerPolicy,
    class FwdSimPolicy,
    class LinFwdSimPolicy,
    class LinSimPolicy,
    class AdjFwdSimPolicy,
    class AdjSimPolicy,
    class Scalar = ireal
  >
  class IWaveInvOp: public Operator<Scalar>,
		    public FwdSamplerPolicy,
		    public LinSamplerPolicy,
		    public AdjSamplerPolicy,
		    public FwdSimPolicy,
		    public LinFwdSimPolicy,
		    public LinSimPolicy,
		    public AdjFwdSimPolicy,
		    public AdjSimPolicy
  {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  private:

    Space<Scalar> const & dom;
    Space<Scalar> const & rng; 
    
    mutable FILE * stream;              /* output stream            */
    mutable PARARRAY lpars;             /* parameter array worspace */
    mutable PARARRAY pars;              /* parameter array ref copy */

    /** reference of initial model (no-window case) or background model (window case)*/
    Vector<Scalar> const & xinit; 
    
    /** reference to bound vectors (for iwop) */ 
    Vector<Scalar> const & lb;
    Vector<Scalar> const & ub;
    
    /** pxfnl stores the pointer to the resulting model from apply() */
    mutable Vector<Scalar> * pxfnl;
    
    /** operator for FD simulation */
    IWaveOp<
      FwdSamplerPolicy,
      LinSamplerPolicy,
      AdjSamplerPolicy,
      FwdSimPolicy,
      LinFwdSimPolicy,
      LinSimPolicy,
      AdjFwdSimPolicy,
      AdjSimPolicy
      >  iwop;
    
    /** operator for processing data output from iwop (e.g., muting op) */
    //    mutable CloneHandle<Operator<Scalar> > postiwop;

    /** simulation operator defined via the composition of postiwop and iwop */
    mutable OpComp<Scalar> simop;    // mutable CloneHandle<Operator<Scalar> > simop;  

    /** operator for manipulating model input for simop (e.g,, windowing, bulkonly projecting) */
    mutable CloneHandle<Operator<Scalar> > presimop;
   
    /** functional for implementing different types of regularization strategies */
    mutable CloneHandle<Functional<Scalar> > regfcnl;

    /** the name vector of model parameters for iwop (for storing x, i.e., m_inv)*/
    std::vector<string> m_names;
    bool is_set;

    /** verbosity control */
    int dump_pars;
    /* inversion history output related flags */
    int dump_model_history;
    int dump_data_history;
    int dump_grad_history;
    
    IWaveInvOp();

    ///////////
    void dump_subinv_model(Vector<Scalar> const & vec, int Isub) const{
      try{
	// SpaceTest(rng,vec,"TSOpt::IWaveInvOp::dump_subinv_model (vec not in rng(IWaveInvOp))");
	SpaceTest(presimop.get().getRange(),vec,"TSOpt::IWaveInvOp::dump_subinv_model (vec not in rng(presimop))");

	/* index string (used to generate history filenames) */
	std::ostringstream sind;
	sind << Isub;
	std::string indstr = sind.str();
	try{      
	  // ProductSpace<Scalar> const & prng = dynamic_cast<ProductSpace<Scalar> const &>(rng);
	  ProductSpace<Scalar> const & prng = dynamic_cast<ProductSpace<Scalar> const &>(presimop.get().getRange());
	  int size = prng.getSize();
	  
	  // Vector<Scalar> xinvcur(rng);
	  Vector<Scalar> xinvcur(presimop.get().getRange());
	  Components<Scalar> cxinvcur(xinvcur);
	  for(int i=1; i<=size; i++){
	    std::ostringstream smind;
	    smind << i;
	    std::string smindstr = smind.str();
	    AssignFilename mfn("m"+smindstr+"_fnl"+indstr+".rsf");
	    cxinvcur[i-1].eval(mfn);
	  }
	  xinvcur.copy(vec);	
	}
	catch (std::bad_cast &bd){
	  RVLException e;
	  e<<bd.what()<<"\n";
	  e<< "IWaveInvOp::dump_subinv_model at subinv "<<Isub<<", rng is not compatiable with ProductSpace\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from IWaveInvOp::dump_subinv_model() \n";
	throw e;
      }
    }

    void dump_subinv_data(Vector<Scalar> const & data, Vector<Scalar> const & fnldata, int Isub) const{
      try{
	SpaceTest(dom,data,"IWaveInvOp::dump_subinv_data (data not in dom(IWaveInvOp))");
	SpaceTest(dom,fnldata,"IWaveInvOp::dump_subinv_data (fnldata not in dom(IWaveInvOp))");
	/* index string (used to generate history filenames) */
	std::ostringstream sind;
	sind << Isub;
	std::string indstr = sind.str();
	
	Vector<Scalar> yc(dom);
	Vector<Scalar> yfc(dom);
	AssignFilename tfn_cur("data_cur"+indstr+".su");
	AssignFilename tfn_fnlcur("data_fnl"+indstr+".su");
	yc.eval(tfn_cur);
	yfc.eval(tfn_fnlcur);
	
	yc.copy(data);     
	yfc.copy(fnldata);	

      }
      catch (RVLException & e) {
	e<<"\ncalled from IWaveInvOp::dump_subinv_data() \n";
	throw e;
      }
    }
    //////////////
  protected:
    
    void apply(const Vector<Scalar> & y,
	       Vector<Scalar> & x) const {
      try{
	
	int err; err = 0;
	cerr<<"IWaveInvOp::APPLY \n";
	
	SpaceTest(this->getDomain(),y,"TSOpt::IWaveInvOp::apply (y not in dom(IWaveInvOp) )");
	SpaceTest(this->getRange(),x,"TSOpt::IWaveInvOp::apply (x not in rng(IWaveInvOp) )");
	
	if (pxfnl != NULL && pxfnl != &x) {
	  RVLException e;
	  e<<"Error: IWaveInvOp::APPLY \n";
	  e<<"Error: pxfnl should equal to &x \n";
	  e<<"Reason: pxfnl is set to NULL in all constructors, and only to non-NULL in this->apply(),\n";
	  e<<"        which can only be called through OperatorEvaluation::getValue() or OpComp::apply/applyDeriv/applyAdjDeriv \n";
	  throw e;	  
	}
	else if (pxfnl == NULL) {
	  cerr<<"IWaveInvOp::APPLY, pxfnl == NULL, inversion procedure starts \n";
	  /*---  data input for current sub-inversion ---*/
	  Vector<Scalar> ycur(this->getDomain()); 
	  ycur.copy(y);
	  
	  /////////////////////////////////////////////////////
	  //         Set up inversion environment            //
	  /////////////////////////////////////////////////////
	  
	  /*------ assign info for inversion with freq continuation ------*/
	  /* Number of frequency continuation sub inversions (default value 1)*/
	  int cfinvnum = 1;
	  if(ps_ffint(pars,"num_subinv",&cfinvnum)){
	    cerr<<"----NOTE: IWaveInvOp::APPLY, failed to extract value for key = num_subinv from param table\n";
	    cerr<<"---------- using default value: 1\n";
	  }
	  if (cfinvnum < 1) {
	    RVLException e;
	    e<<"Error: IWaveInvOp::APPLY \n";
	    e<<"Error: num_subinv must be at least 1 (current value = "<<cfinvnum<<")\n";
	    throw e;
	  }
	  /* power for defining trapezoid filters (default value 2.0)*/
	  float filterpwr = 2.;
	  if(ps_fffloat(pars,"filter_power",&filterpwr)){
	    cerr<<"----NOTE: IWaveInvOp::APPLY \n";
	    cerr<<"---------- failed to extract value for key = filter_power from param table\n";
	    cerr<<"---------- using default value: 2.0\n";
	  }
	  /* frequency bands for defining trapezoid filters*/
	  float *freqb3 = new float[cfinvnum];
	  float *freqb4 = new float[cfinvnum];
	  freqb3[0] = 3.;
	  freqb4[0] = 5.;
	  for (int i = 1; i<= cfinvnum; i++){
	    std::ostringstream sind;
	    sind << i;
	    string key = "freqbd3_inv"+ sind.str();
	    if( ps_fffloat(pars, key.c_str(), freqb3 + (i-1)) && (cfinvnum - 1)){
	      RVLException e;
	      e<<"Error: IWaveInvOp::APPLY\n";
	      e<<"Error: failed to extract value for key = "<<key<<" from param table\n";
	      throw e;
	    }
	    //	cerr<<"freqb3["<<(i-1)<<"]="<<freqb3[i-1]<<endl;
	    key = "freqbd4_inv"+sind.str();
	    if((cfinvnum - 1) && ps_fffloat(pars, key.c_str(), freqb4 + (i-1))){
	      RVLException e;
	      e<<"Error: IWaveInvOp::APPLY \n";
	      e<<"Error: failed to extract value for key = "<<key<<" from param table\n";
	      throw e;
	    }
	    //	cerr<<"freqb4["<<(i-1)<<"]="<<freqb4[i-1]<<endl;
	  }
	  
	  /* build ULBoundsTest */
	  RVLMin<Scalar> mn;
#ifdef IWAVE_USE_MPI
	  MPISerialFunctionObjectRedn<float,float> mpimn(mn);
	  ULBoundsTest<float> ultest(lb,ub,mpimn);
#else
	  ULBoundsTest<float> ultest(lb,ub,mn);
#endif
	
	  /* employ filter (to do freq-continuous inversion) */	
	  //cerr<<" IWaveInvOp::APPLY, employ filter \n";
	  SEGYBandFilter filt;    
	  LinearOpFO<Scalar> filtop(dom,dom,filt,filt);
	  LNLOperator<Scalar> nfiltop(filtop);
	  OpComp<Scalar> fwdop(simop,nfiltop);
	  
	  /*------ build lsobjfcnl objective ------*/
	  //cerr<<" IWaveInvOp::APPLY, build lsobjfcnl \n";
	  CloneHandle<Functional<Scalar> > lsobjfcnl;
	  /* build fls (standard least-squares fcnl) */	
	  StdLeastSquaresFcnlGN<Scalar> fls(fwdop,ycur);
	  /* employ regularization */
	  //if(regflag){
	  LinCombFunctional<Scalar> regfls(ScalarFieldTraits<Scalar>::One(),fls,ScalarFieldTraits<Scalar>::One(),regfcnl.get());
	  lsobjfcnl.set(regfls);
	  //}
	  //else{
	  //  lsobjfcnl.set(fls);
	  //}
	
	  /*------ build ls-objective with bound constrains ------*/
	  //cerr<<" IWaveInvOp::APPLY, build fobjbd \n";
	  FunctionalBd<Scalar, ULBoundsTest<Scalar> > fobjbd(lsobjfcnl.get(),ultest);
	 
	  /*------ build final ls-objective ------*/	
	  cerr<<" IWaveInvOp::APPLY, build f \n";
	  FcnlOpComp<Scalar> f(fobjbd,presimop.get());
	
	  /*------ inversion with frequency continuation ------*/
	  
	  /*--- initial model (iterator) for sub-inversions ---*/
	  //Vector<Scalar> x0(presimop->getDomain());
	  //x0.copy(xinit);
	  
	  x.copy(xinit);
	  
	  /** NOTE: replace xfnl with xinit, which is defined in driver, 
	      must ensure the xinit in driver lives longer than any objects of this operator */
	  /** D.S. 03/28/11: build umin inside sub-inv loop to configure sub-inversions seperately */
	  /*------ minimization pars ------*/
	  // string uminparname;
	  // char * cbuf = NULL;
	  // ps_ffcstring(pars,"uminpar",&cbuf);
	  // if (!cbuf) {
	  //    RVLException e;
	  //    e<<"Error: IWaveInvOp::Apply \n";
	  //    e<<"failed to extract value for key = uminpar from param table\n";
	  //    throw e;
	  // }
	  // else {
	  //    uminparname = cbuf;
	  //    free(cbuf);
	  //    cbuf=NULL;
	  // }
	  // UMinTable<float> uminpar(uminparname);
	  // cerr<<" IWaveInvOp::APPLY, build umin \n";
	  // UMinMethod<Scalar> umin(f,x,uminpar,cerr);
	  //---------------------------------------------------------- 
		
	  // cerr<<" IWaveInvOp::APPLY, build presimeval \n";
	  OperatorEvaluation<Scalar> presimeval(presimop.get(),x); //only used for dump history
	  
	  // cerr<<" IWaveInvOp::APPLY, start inversion loop \n";
	  for (int i = 1; i <= cfinvnum; i++){
	    std::ostringstream sind_inv;
	    sind_inv << i;
	    string parkey = "uminpar_inv"+ sind_inv.str();
	    string uminparname;
	    char * cbuf = NULL;
	    ps_ffcstring(pars,parkey.c_str(),&cbuf);
	    if (!cbuf){
	      ps_ffcstring(pars,"uminpar",&cbuf);
	      if (!cbuf) {
		RVLException e;
		e<<"Error: IWaveInvOp::Apply \n";
		e<<"failed to extract value for key = "<<parkey<<" and uminpar from param table\n";
		throw e;
	      }
	      else {
		uminparname = cbuf;
		free(cbuf);
		cbuf=NULL;
		// cerr<<" IWaveInvOp::APPLY, specific uminpar-file "<<parkey<<" is not provided; use default file "<<uminparname<<" for all subinversions"<<endl;
	      }
	    }
	    else {
	      uminparname = cbuf;
	      free(cbuf);
	      cbuf=NULL;
	      cerr<<"----- IWaveInvOp::APPLY, specific uminpar-file "<<parkey<<" is provided as "<<uminparname<<" for sub-inversion "<<i<<endl;
	    }
	    
	    //UMinTable<float> uminpar(uminparname);
	    Table uminpar(uminparname);

	    // cerr<<" IWaveInvOp::APPLY, subinv #"<<i<<", build umin"<<endl;
	    //UMinMethod<Scalar> umin(f,x,uminpar,cerr);
	    LBFGSBT<float> umin(f,x,uminpar,cerr);

	    /* set-up band-filter */
	    if (i == cfinvnum) {
	      filt.set(0.,0.,freqb3[i-1],freqb4[i-1],1.,1.,1.,1.,filterpwr);   
	    }
	    else{
	      filt.set(0.,0.,freqb3[i-1],freqb4[i-1],1.,1.,1.,0.,filterpwr);   
	    }
	    
	    /* get data input for current sub-inversion */
	    ycur.eval(filt,y); 
	    
	    //cerr<<" IWaveInvOp::APPLY, subinverion "<<i<<", umin.run starts \n"; 
	    umin.run();
	    // cerr<<" IWaveInvOp::APPLY, subinverion "<<i<<", umin.run ended \n";  
	    
	    /** Note: now xfnl stores the final model for current sub-inversion  
	     *   which'll be  the initial model for the next sub-inversion 
	     */
	    
	    
	    /*------ store inversion history into files 
	      (finalmodel and finaldata for sub-inversion)------*/
	    /* output model history, when "dump_model_history != 0" */ 	 
	    // set output filenames to "m{1,...,NPAR}_fnl{$Isub}.rsf" 
	    // data output only allowed when dumping model history also requested 
	    // in addition, output data history, when "dump_data_history != 0" 
	    // set output filenames to "data_cur${Isub}.su" and "data_fnl${Isub}.su" 
	    
	    // NOTE: store MODEL, not MODEL INCREMENT(or difference with background)!!!
	    
	    if(dump_model_history){
	      dump_subinv_model(presimeval.getValue(),i);
	      if(dump_data_history){
		OperatorEvaluation<Scalar> final_opeval(fwdop,presimeval.getValue());
		dump_subinv_data(ycur,final_opeval.getValue(),i);
	      }
	    }
	    
	  }//end of inversion loop
	
	  /* resulting model from inversion*/
	  // x.copy(presimeval.getValue());
	  pxfnl = &x;
	
	  /* store inversion results into disk files */
	  if(dump_model_history){
	    Vector<Scalar> x_s(rng);
	    Components<Scalar> cx_s(x_s);
	    int c_size = cx_s.getSize();
	    for (int i = 1; i<= c_size; i++){
	      std::ostringstream smind;
	      smind << i;
	      std::string smindstr = smind.str();
	      AssignFilename mfn("m"+smindstr+"_invfnl.rsf");
	      cx_s[i-1].eval(mfn);	
	    }
	    x_s.copy(x);
	  }
	  /* set applied indicator to be true */
	  // applied = true;
	  
	  /* clean up */
	  delete [] freqb3;
	  delete [] freqb4;
	} // if(pxfnl == NULL)
	cerr<<"IWaveInvOp::apply, EXIT \n";
	
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSOpt::IWaveInvOp::apply\n";
	throw e;
      } 
      
    }
    
    void applyDeriv(const Vector<Scalar> & y,
		    const Vector<Scalar> & dy,
		    Vector<Scalar> &dx) const{
      try {
	int err; err = 0;
	cerr<<" IWaveInvOp::applyDeriv \n";
	////////// In Progress ////////////
	
	SpaceTest(this->getDomain(),y,"TSOpt::IWaveInvOp::applyDeriv (y ?in dom)");
	SpaceTest(this->getDomain(),dy,"TSOpt::IWaveInvOp::applyDeriv (dy ?in dom)");
	SpaceTest(this->getRange(),dx,"TSOpt::IWaveInvOp::applyDeriv (dx ?in rng)");

	/* set up temp copy x */
	Vector<Scalar> x(presimop.get().getDomain());
	
	/*------ get current xfnl (without duplicated inversion) ------*/
	if (pxfnl == NULL) {
	  if (is_set) {
	    int size = m_names.size();
	    if (size == 1) {
	      AssignFilename tmpfn(m_names[0]);
	      x.eval(tmpfn);	      
	    }
	    else if (size >1) {
	      Components<Scalar> cx(x);
	      for (int j =0; j < size; j++) {
		AssignFilename tmpfn(m_names[j]);
		cx[j].eval(tmpfn);
	      }
	    }
	    else {
	      RVLException e;
	      e<<"Error: IWaveInvOp::APPLYDeriv \n";
	      e<<"Error: pxfnl == NULL and is_set == true, but m_names provides inconsistent information \n";
	      throw e;
	    }
	  }
	  else{
	    RVLException e;
	    e<<"Error: IWaveInvOp::APPLYDeriv \n";
	    e<<"Error: pxfnl == NULL and is_set == false, applyDeriv must be called after OperatorEvaluation::getValue() \n";
	    throw e;	  
	  }
	}
	else { 
	  cerr<<"IWaveInvOp::applyDeriv: xfnl is the latest version; donot redo inversion in this->apply()\n";
	  x.copy(*pxfnl);
	}
	
	/** Approach II: directly solving Normal Equation
	          
	          \f$ H \delta m = DF[m]^T \delta d \f$,
	    
	    where 
	    \f$ H = DF[m]^T DF[m] + D^2 R \f$,
	    \f$ m = x = opeval(IWaveInvOp,y), dy = \delta d, dx = \delta m \f$	  
	*/
	
	/* initial iterator */
	dx.zero();
	
	// donot need the following copy procedure, if change constructors in ls.hh
	// Vector<Scalar> dycur(dom);
	// dycur.copy(dy);

	/** The following approach is only for testing purpose. 
	 *  Later, it will be replaced by a more efficient optimization process,
	 *  which makes use of a more sophisticated linear solver.
	 */

	/*------ read parameters  ------*/
	int _cgmaxiter = 10;
	atype _cgtol = 1e-8;
	atype _cgmaxstep = numeric_limits<atype>::max();
	if(!parse<int>(pars,"cg_maxiter",_cgmaxiter)) {
	  cerr<<"----NOTE: IWaveInvOp::APPLYDERIV, failed to extract value for key = cg_maxiter from param table\n";
	  cerr<<"---------- using default value: "<<_cgmaxiter<<"\n";
	}
	if(!parse<atype>(pars,"cg_tol",_cgtol)) {
	  cerr<<"----NOTE: IWaveInvOp::APPLYDERIV, failed to extract value for key = cg_tol from param table\n";
	  cerr<<"---------- using default value: "<<_cgtol<<"\n";
	}
	if(!parse<atype>(pars,"cg_maxstep",_cgmaxstep)) {
	  cerr<<"----NOTE: IWaveInvOp::APPLYDERIV, failed to extract value for key = cg_maxstep from param table\n";
	  cerr<<"---------- using default value: "<<_cgmaxstep<<"\n";
	}

	/*------ build linear operator and right hand side ------*/
	OpComp<Scalar> Fbar(presimop.get(),simop);

	cerr<<" IWaveInvOp::applyDeriv, build opeval \n";
	OperatorEvaluation<Scalar> opeval(Fbar,x);

	//-- Hessian of regularization term ----
	FcnlOpComp<Scalar> regfcnlbar(regfcnl.get(),presimop.get());	
	FunctionalEvaluation<Scalar> regeval(regfcnlbar,x);
	HessianEvaluation<Scalar> regH(regeval);

	cerr<<" IWaveInvOp::applyDeriv, build (DFbar^T DFbar) \n";
	NormalLinearOp<Scalar> NormOpFbar(opeval.getDeriv());
	
	LinCombLinearOp<Scalar> NormOp(1.,NormOpFbar,1.,regH);
	
	cerr<<" IWaveInvOp::applyDeriv, build (DFbar^T dy) \n";
	Vector<Scalar> dx_rt(presimop.get().getDomain());
	// output intermediate results
	Components<Scalar> cdx_rt(dx_rt);
	for (int j =0; j < cdx_rt.getSize(); j++) {
	  std::ostringstream cdxind;
	  cdxind << j;
	  std::string cdxstr = cdxind.str();
	  AssignFilename tmpfn("m"+cdxstr+"_iwinvop_dev_dump.rsf");
	  cdx_rt[j].eval(tmpfn);
	}
	// -------------------
	opeval.getDeriv().applyAdjOp(dy,dx_rt);
	
	
	// /*------ build CG-linsolver ------*/
	//CGLinearSolver<Scalar> sglinsolv(_cgmaxiter,_cgtol,_cgmaxstep,cerr);
	//sglinsolv.setSystem(NormOp, dx, dx_rt);
	//cerr<<" IWaveInvOp::applyDeriv, bulkonly-case, sglinsolv.run \n";
	//sglinsolv.run();

	/*------ build CG-Alg ------*/
	atype rnormsq = dx_rt.normsq();
	CGAlg<Scalar> sglinsolv(dx, NormOp, dx_rt, rnormsq, _cgtol,_cgmaxiter, _cgmaxstep, cerr);
	cerr<<" IWaveInvOp::applyDeriv, bulkonly-case, sglinsolv.run \n";
	sglinsolv.run();
	
	//cerr<<" IWaveInvOp::applyDeriv, bulkonly-case, build presimeval \n";
	// OperatorEvaluation<Scalar> presimeval(presimop.get(),xfnl);
	// presimeval.getDeriv().applyOp(dx0,dx);

	cerr<<" IWaveInvOp::applyDeriv, EXIT \n";	
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSOpt::IWaveInvOp::applyDeriv\n";
	throw e;
      } 
      
    }
    
    void applyAdjDeriv(const Vector<Scalar> & y,
		       const Vector<Scalar> & dx,
		       Vector<Scalar> & dy) const {
      try{
	int err; err =0;
	cerr<<" IWaveInvOp::applyAdjDeriv \n";
	/////////// In Progress ///////////////

	SpaceTest(this->getDomain(),y,"TSOpt::IWaveInvOp::applyDeriv (y ?in dom)");
	SpaceTest(this->getDomain(),dy,"TSOpt::IWaveInvOp::applyDeriv (dy ?in dom)");
	SpaceTest(this->getRange(),dx,"TSOpt::IWaveInvOp::applyDeriv (dx ?in rng)");

	/* set up temp copy x */
	Vector<Scalar> x(presimop.get().getDomain());
	
	/*------ get current xfnl (without duplicated inversion) ------*/
	if (pxfnl == NULL) {
	  if (is_set) {
	    int size = m_names.size();
	    if (size == 1) {
	      AssignFilename tmpfn(m_names[0]);
	      x.eval(tmpfn);	      
	    }
	    else if (size >1) {
	      Components<Scalar> cx(x);
	      for (int j =0; j < size; j++) {
		AssignFilename tmpfn(m_names[j]);
		cx[j].eval(tmpfn);
	      }
	    }
	    else {
	      RVLException e;
	      e<<"Error: IWaveInvOp::ApplyAdjDeriv \n";
	      e<<"Error: pxfnl == NULL and is_set == true, but m_names provides inconsistent information \n";
	      throw e;
	    }
	  }
	  else{
	    RVLException e;
	    e<<"Error: IWaveInvOp::ApplyAdjDeriv \n";
	    e<<"Error: pxfnl == NULL and is_set == false, applyAdjDeriv must be called after OperatorEvaluation::getValue() \n";
	    throw e;	  
	  }
	}
	else { 
	  cerr<<"IWaveInvOp::applyAdjDeriv: xfnl is the latest version; donot redo inversion in this->apply()\n";
	  x.copy(*pxfnl);
	}

	/** Approach II: directly solving Normal Equation
	          
	          \f$ H \delta m_1 = \delta m \f$,
	     then compute
	          \f$ \delta d = DF[m] \delta m_1 \f$,
	    where 
	    \f$ H = DF[m]^T DF[m] + D^2 R \f$,
	    \f$ m = x = opeval(IWaveInvOp,y), dy = \delta d, dx = \delta m \f$	  
	*/

	/* initial iterator */
	dy.zero();

	// donot need the following copy procedure, if change constructors in ls.hh
	//Vector<Scalar> adxc(rng);
	//adxc.copy(dx);

	/** The following approach is only for testing purpose. 
	 *  Later, it will be replaced by a more efficient optimization process,
	 *  which makes use of a more sophisticated linear solver.
	 */

	/*------ read parameters  ------*/
	int _cgmaxiter = 10;
	atype _cgtol = 1e-8;
	atype _cgmaxstep = numeric_limits<atype>::max();
	if(!parse<int>(pars,"cg_maxiter",_cgmaxiter)) {
	  cerr<<"----NOTE: IWaveInvOp::APPLYADJDERIV, failed to extract value for key = cg_maxiter from param table\n";
	  cerr<<"---------- using default value: "<<_cgmaxiter<<"\n";
	}
	if(!parse<atype>(pars,"cg_tol",_cgtol)) {
	  cerr<<"----NOTE: IWaveInvOp::APPLYADJDERIV, failed to extract value for key = cg_tol from param table\n";
	  cerr<<"---------- using default value: "<<_cgtol<<"\n";
	}
	if(!parse<atype>(pars,"cg_maxstep",_cgmaxstep)) {
	  cerr<<"----NOTE: IWaveInvOp::APPLYADJDERIV, failed to extract value for key = cg_maxstep from param table\n";
	  cerr<<"---------- using default value: "<<_cgmaxstep<<"\n";
	}
	
	/*------ build linear operator and right hand side ------*/
	OpComp<Scalar> Fbar(presimop.get(),simop);

	cerr<<" IWaveInvOp::applyAdjDeriv, build opeval \n";
	OperatorEvaluation<Scalar> opeval(Fbar,x);

	//-- Hessian of regularization term ----
	FcnlOpComp<Scalar> regfcnlbar(regfcnl.get(),presimop.get());	
	FunctionalEvaluation<Scalar> regeval(regfcnlbar,x);
	HessianEvaluation<Scalar> regH(regeval);

	cerr<<" IWaveInvOp::applyAdjDeriv, build (DFbar^T DFbar) \n";
	NormalLinearOp<Scalar> NormOpFbar(opeval.getDeriv());

	LinCombLinearOp<Scalar> NormOp(1.,NormOpFbar,1.,regH);	

	cerr<<" IWaveInvOp::applyAdjDeriv, build (Hinv dx) \n";
	Vector<Scalar> dxc(rng);
	dxc.copy(dx);
	Vector<Scalar> dx_tmp(rng);
	// output intermediate results
	Components<Scalar> cdx_tmp(dx_tmp);
	for (int j =0; j < cdx_tmp.getSize(); j++) {
	  std::ostringstream cdxind;
	  cdxind << j;
	  std::string cdxstr = cdxind.str();
	  AssignFilename tmpfn("m"+cdxstr+"_iwinvop_adjdev_dump.rsf");
	  cdx_tmp[j].eval(tmpfn);
	}
	// -------------------
	
	dx_tmp.zero();
	
	///*------ build CG-linsolver ------*/
	//CGLinearSolver<Scalar> sglinsolv(_cgmaxiter,_cgtol,_cgmaxstep,cerr);
	//sglinsolv.setSystem(NormOp, dx_tmp, dxc);

	//cerr<<" IWaveInvOp::applyAdjDeriv, bulkonly-case, sglinsolv.run \n";
	
	//sglinsolv.run();

	/*------ build CG-Alg ------*/
	atype rnormsq = dxc.normsq();
	CGAlg<Scalar> sglinsolv(dx_tmp, NormOp, dxc, rnormsq, _cgtol,_cgmaxiter, _cgmaxstep, cerr);
	cerr<<" IWaveInvOp::applyDeriv, bulkonly-case, sglinsolv.run \n";
	sglinsolv.run();
	
	cerr<<" IWaveInvOp::applyAdjDeriv, compute DFbar dx_tmp \n";
	opeval.getDeriv().applyOp(dx_tmp,dy);

	cerr<<" IWaveInvOp::applyAdjDeriv, EXIT \n";
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSOpt::IWaveInvOp::applyAdjDeriv\n";
	throw e;
      } 
    }
    
    Operator<Scalar> * clone() const { 
      // cerr<<"IWaveInvOp::clone \n";
      return new IWaveInvOp< 
                     FwdSamplerPolicy,
		     LinSamplerPolicy,
		     AdjSamplerPolicy,
		     FwdSimPolicy,
		     LinFwdSimPolicy,
		     LinSimPolicy,
		     AdjFwdSimPolicy,
		     AdjSimPolicy,
		     Scalar
		     >(*this);
      // cerr<<"exit IWaveInvOp::clone \n";
    }
   
    
  public:
    
    IWaveInvOp(Vector<Scalar> const & _xinit,
	       Vector<Scalar> const & _lb,
	       Vector<Scalar> const & _ub,
	       PARARRAY _pars,
	       FILE * _stream,
	       int (*_minit)(struct GFD_MODEL *gfdm),
	       Operator<Scalar> const & _presimop,
	       Operator<Scalar> const & _postiwop, 
	       Functional<Scalar> const & _regfcnl
	       )
      :dom(_postiwop.getRange()), rng(_presimop.getDomain()), xinit(_xinit), 
       lb(_lb), ub(_ub), pxfnl(NULL), stream(_stream),
       iwop(_presimop.getRange(),_postiwop.getDomain(),_pars,stream,_minit), 
       simop(iwop,_postiwop), presimop(_presimop), regfcnl(_regfcnl)
    {
      try{
	int err=0;
	/* Space Test */
	SpaceTest(presimop.get().getDomain(),xinit,"IWaveInvOp Constructor (xinit not in dom(presimop))");
	SpaceTest(iwop.getDomain(),lb,"IWaveInvOp Constructor (lb not in dom(iwop))");
	SpaceTest(iwop.getDomain(),ub,"IWaveInvOp Constructor (ub not in dom(iwop))");
	
	/* copy input par array to data member */
	ps_setnull(&pars);
	if (err=ps_copy(_pars,&pars)) {
	  RVLException e;
	  e<<"Error: IWaveInvOp constructor from ps_copy, err="<<err<<"\n";
	  throw e;
	}
	
	is_set=false;
	
	/* set dump controls */
	dump_pars = dump_model_history = dump_data_history = dump_grad_history = 0;
	ps_ffint(pars,"dump_pars",&dump_pars);
	ps_ffint(pars,"dump_model_history",&dump_model_history);
	ps_ffint(pars,"dump_data_history",&dump_data_history);
	ps_ffint(pars,"dump_grad_history",&dump_grad_history);
	
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSOpt::IWaveInvOp Constructor\n";
	throw e;
      } 
    }

    IWaveInvOp(IWaveInvOp
	       <
	       FwdSamplerPolicy,
	       LinSamplerPolicy,
	       AdjSamplerPolicy,
	       FwdSimPolicy,
	       LinFwdSimPolicy,
	       LinSimPolicy,
	       AdjFwdSimPolicy,
	       AdjSimPolicy,
	       Scalar
	       > const &_iwinv)
      :dom(_iwinv.dom), rng(_iwinv.rng), xinit(_iwinv.xinit), lb(_iwinv.lb), ub(_iwinv.ub),
       pxfnl(_iwinv.pxfnl), stream(_iwinv.stream), iwop(_iwinv.iwop), 
       simop(_iwinv.simop), presimop(_iwinv.presimop), regfcnl(_iwinv.regfcnl),
       m_names(_iwinv.m_names), is_set(_iwinv.is_set),
       dump_pars(_iwinv.dump_pars), dump_model_history(_iwinv.dump_model_history),
       dump_data_history(_iwinv.dump_data_history), dump_grad_history(_iwinv.dump_grad_history)
    {
      int err=0;
      //cerr<<"IWaveInvOp copy constructor \n";
      
      /*------ set pars ------*/
      ps_setnull(&pars);
      if (err=ps_copy(_iwinv.pars,&pars)) {
	RVLException e;
	e<<"Error: IWInvOP copy constructor from ps_copy, err="<<err<<"\n";
	throw e;
      }
      
    }
  
  
    ~IWaveInvOp(){ }  //cerr<<"exit IWaveInvOp distructor \n";
    
    void set_m_names(int num, string * _mnames){
      for (int i=0; i<num; i++) 
	m_names.push_back(_mnames[i]);
      is_set = true;
    }

    const Space<Scalar> & getDomain() const { return dom; }
    const Space<Scalar> & getRange() const { return rng; }
    
    PARARRAY & getPar() { return pars; }
    PARARRAY const & getPar() const { return pars; }
    
    ostream & write(ostream & str) const { return str; }
    
  };
  
}
#endif
