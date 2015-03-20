#include "acd_defn.hh"
#include "grid.h"
#include "gridpp.hh"
#include "gridops.hh"
#include "segyops.hh"
#include "iwop.hh"
#include "functions.hh"
#include "op.hh"
#include "ls.hh"
#include "blockop.hh"
#include "cgnealg.hh"
#include "TRGNAlg.hh"
#include "LBFGSBT.hh"
#include "acdfwi_selfdoc.h"
#include <par.h>

//#define GTEST_VERBOSE
IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq",    0, true,  true },
  {"data",   1, false, true },
  {"source", 1, true,  false},
  {"",       0, false, false}
};

using RVL::valparse;
using RVL::RVLException;
using RVL::Vector;
using RVL::Components;
using RVL::Operator;
using RVL::OperatorEvaluation;
using RVL::LinearOp;
using RVL::LinearOpFO;
using RVL::OpComp;
using RVL::StdLeastSquaresFcnlGN;
using RVL::SymmetricBilinearOp;
using RVL::AssignFilename;
using RVL::AssignParams;
using RVL::RVLRandomize;
using RVL::ScaleOpFwd;
using RVL::ResidualOperator;
using RVL::TensorOp;
using RVLUmin::CGNEPolicy;
using RVLUmin::TRGNAlg;
using RVLUmin::LBFGSBT;
using TSOpt::IWaveEnvironment;
using TSOpt::IWaveTree;
using TSOpt::IWaveSampler;
using TSOpt::IWaveSim;
using TSOpt::TASK_RELN;
using TSOpt::IOTask;
using TSOpt::IWaveOp;
using TSOpt::SEGYTaperMute;
using TSOpt::GridMaskOp;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
using TSOpt::MPISEGYSpace;
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
#endif
using TSOpt::GridExtendOp;
using TSOpt::GridDerivOp;
//using TSOpt::GridHelmOp;

int xargc;
char **xargv;

/** this version requires that two model files be present:
    csqext - extended version of csq, with [d=spatial dimn]
    gdim=d+1
    dim=d
    n[d+1]=number of shot records, 
    d[d+1]=1
    o[d+1]=0
    id[d+1]=d+1
    csq - physical space version of csq, with dim=gdim=d
    Uses these keywords for geometry only - data content irrelevant

    [intended near-term update: require only physical space material
    params, generate extended version internally, so automatic 
    compatibility with data]
*/
int main(int argc, char ** argv) {

  try {

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif

    PARARRAY * pars = NULL;
    FILE * stream = NULL;
    IWaveEnvironment(argc, argv, 0, &pars, &stream);

    if (retrieveGlobalRank()==0 && argc<2) {
      pagedoc();
      exit(0);
    }

#ifdef IWAVE_USE_MPI
    if (retrieveGroupID() == MPI_UNDEFINED) {
      fprintf(stream,"NOTE: finalize MPI, cleanup, exit\n");
      fflush(stream);
    }
    else {
#endif

      /* the Op */
      IWaveOp iwop(*pars,stream);

      /* generate physical model space */
#ifdef IWAVE_USE_MPI
      MPIGridSpace sp(valparse<std::string>(*pars,"dom"),"notype",true);
#else
      GridSpace sp(valparse<std::string>(*pars,"dom"),"notype",true);
#endif
      // make it a product, so it's compatible with domain of op
      StdProductSpace<ireal> dom(sp); 

      //      SEGYTaperMute(float _s=0.0f, float _tm=0.0f, float _w=0.0f, int _type = 0, float _taper_min=0.0f, float _taper_max=0.0f, float _width=0.0f, int _tapertype=0, float _tw=0.0f)
        
      SEGYTaperMute tnm(valparse<float>(*pars,"mute_slope",0.0f),
                        valparse<float>(*pars,"mute_zotime",0.0f),
                        valparse<float>(*pars,"mute_width",0.0f),
			valparse<int>(*pars,"mute_type",0),
                        valparse<float>(*pars,"taper_min",-numeric_limits<float>::max()),
                        valparse<float>(*pars,"taper_max",numeric_limits<float>::max()),
                        valparse<float>(*pars,"taper_width",0.0f),
                        valparse<int>(*pars,"taper_type",0),
                        valparse<float>(*pars,"time_width",0.0f));

      // time width = width of final time taper
      LinearOpFO<float> tnmop(iwop.getRange(),iwop.getRange(),tnm,tnm);
      OpComp<float> mop(iwop,tnmop);
    
      // vel-squared reference model
      Vector<ireal> m0(iwop.getDomain());
      AssignFilename mf0n(valparse<std::string>(*pars,"csq_bg"));
      Components<ireal> cm0(m0);
      cm0[0].eval(mf0n);

      // spatial window taper vec
      RPNT wind;
      RASN(wind,RPNT_0);
      wind[0]=valparse<float>(*pars,"wind0",0.0f);
      wind[1]=valparse<float>(*pars,"wind1",0.0f);
      wind[2]=valparse<float>(*pars,"wind2",0.0f);
      
      // ceate window op, sim
      GridWindowOp wop(dom,m0,wind);
      OpComp<float> op(wop,mop);

      // velocity update vector
      Vector<ireal> m(op.getDomain());
      // make archival
      AssignFilename mfn(valparse<std::string>(*pars,"csq_update"));
      Components<ireal> cm(m);
      cm[0].eval(mfn);
      m.zero();
    
      // muted data
      Vector<ireal> mdd(op.getRange());
      std::string mddnm = valparse<std::string>(*pars,"datamut","");
      if (mddnm.size()>0) {
        AssignFilename mddfn(mddnm);
        mdd.eval(mddfn);
      }
      { // read data - only needed to define muted data
	Vector<ireal> dd(op.getRange());
	AssignFilename ddfn(valparse<std::string>(*pars,"data"));
	Components<ireal> cdd(dd);
	cdd[0].eval(ddfn);
	tnmop.applyOp(dd,mdd);
      }
      // residual operator - used in TRCG, zero cost so just build it
      ResidualOperator<float> rop(op,mdd);

      // choice of GridDerivOp for semblance op is placeholder - extended axis is dim-th, with id=dim+1
      // note default weight of zero!!!
      // Note added 10.03.14: getGrid is not usably implemented for MPIGridSpace at this time, so 
      // must fish the derivative index out by hand and bcast
      // THIS SHOULD ALL BE FIXED! (1) getGrid properly implemented in parallel, (2) proper
      // defn of DSOp to include internal extd axis case (space/time shift)
      int dsdir = INT_MAX;
      if (retrieveGlobalRank()==0) dsdir=csqsp.getGrid().dim;
#ifdef IWAVE_USE_MPI
      if (MPI_Bcast(&dsdir,1,MPI_INT,0,retrieveGlobalComm())) {
	RVLException e;
	e<<"Error: acdiva, rank="<<retrieveGlobalRank()<<"\n";
	e<<"  failed to bcast dsdir\n";
	throw e;
      }
#endif
      // does not appear to work properly
      // assign window widths - default = 0;
      
      //      StdLeastSquaresFcnlGN<float> f(cop,mdd);
      StdLeastSquaresFcnlGN<float> f(op,mdd);

      // choice of preop is placeholder
      // ScaleOpFwd<float> preop(op.getDomain(),1.0f);

      // lower, upper bds for csq
      Vector<float> lb(f.getDomain());
      Vector<float> ub(f.getDomain());
      float lbcsq=valparse<float>(*pars,"cmin");
      lbcsq=lbcsq*lbcsq;
      RVLAssignConst<float> asl(lbcsq);
      lb.eval(asl);
      float ubcsq=valparse<float>(*pars,"cmax");
      ubcsq=ubcsq*ubcsq;
      RVLAssignConst<float> asu(ubcsq);
      ub.eval(asu);
      RVLMin<float> mn;
#ifdef IWAVE_USE_MPI
      MPISerialFunctionObjectRedn<float,float> mpimn(mn);
      ULBoundsTest<float> ultest(lb,ub,mpimn);
#else
      ULBoundsTest<float> ultest(lb,ub,mn);
#endif

      FunctionalBd<float> fbd(f,ultest);

      /* output stream */

      std::stringstream outstr;
      outstr<<scientific;
      std::ostream * optr = &outstr;
      std::string outfile = valparse<std::string>(*pars,"outfile","");
      if ((outfile.size()==0) && (retrieveRank()==0)) optr = &cout;

      // switch on method
      Algorithm * alg = NULL;
      std::string optmethod = valparse<std::string>(*pars,"OptMethod","lbfgs");
      if (optmethod == "lbfgs") {
	alg = new LBFGSBT<float>
	  (fbd, m,
	   valparse<float>(*pars,"InvHessianScale",1.0f),
	   valparse<int>(*pars,"MaxInvHessianUpdates",5),
	   valparse<int>(*pars,"MaxSubSteps",10),
	   valparse<bool>(*pars,"VerboseDisplay",true), 
	   valparse<float>(*pars,"InitStepBound",1.0f), 
	   valparse<float>(*pars,"MinDecrease",0.1f),
	   valparse<float>(*pars,"GoodDecrease",0.9f), 
	   valparse<float>(*pars,"StepDecrFactor",0.5f), 
	   valparse<float>(*pars,"StepIncrFactor",1.8f),
	   valparse<float>(*pars,"MaxFracDistToBdry",1.0), 
	   valparse<float>(*pars,"LSMinStepFrac",1.e-06),
	   valparse<int>(*pars,"MaxSteps",10), 
	   valparse<float>(*pars,"AbsGradThresh",0.0), 
	   valparse<float>(*pars,"RelGradThresh",1.e-2), 
	   *optr);
      }
      else if (optmethod == "trcg") {
	TRGNAlg<float, CGNEPolicy<float> > * tralg = new TRGNAlg<float, CGNEPolicy<float> >
	  (rop, m,
	   valparse<int>(*pars,"MaxSteps",10),             // _maxcount,
	   valparse<float>(*pars,"ResidualTol",0.0f),       // _jtol,
	   valparse<float>(*pars,"AbsGradThresh",0.0f),    // _agtol,
	   valparse<float>(*pars,"RelGradThresh",1.0e-2),  // _rgtol,
	   valparse<float>(*pars,"MinDecrease",0.1f),      // _eta1
	   valparse<float>(*pars,"GoodDecrease",0.9f),     // _eta2
	   valparse<float>(*pars,"StepDecrFactor",0.5f),   // _gamma1
	   valparse<float>(*pars,"StepIncrFactor",1.8f),   // _gamma2
	   *optr);

	// assign CG params
	tralg->assign
	  (valparse<float>(*pars,"CGNE_ResTol",0.0f),      // rtol
	   valparse<float>(*pars,"CGNE_GradTol",0.001f),   // nrtol,
	   valparse<float>(*pars,"InitStepBound",1.0f),    // Delta
	   valparse<int>(*pars,"MaxSubSteps",10),          // maxcount
	   valparse<bool>(*pars,"VerboseDisplay",true));   // verbose
	alg=tralg;
      }
      else {
	RVLException e;
	e<<"Error: acdfwi\n";
	e<<"  failed to specify legit opt method choice\n";
	e<<"  currently available: lbfgs or trcg\n";
	throw(e);
      }

      if (valparse<int>(*pars,"MaxLBFGSIter",3) <= 0) {
	FunctionalEvaluation<float> feval(f,m);
	//	Vector<float> grad(alg.getFunctionalEvaluation().getDomain());
	Vector<float> grad(feval.getDomain());
	AssignFilename gfn("grad.rsf");
	Components<float> cgrad(grad);
	cgrad[0].eval(gfn);
	cerr<<"fcnl gradient"<<endl;
	grad.copy(feval.getGradient());
	cerr<<"grad norm from fcnl = "<<grad.norm()<<endl;
	cerr<<"fcnl value"<<endl;
	cerr<<"val = "<<feval.getValue()<<endl;
	Vector<float> din(op.getRange());
	OperatorEvaluation<float> opeval1(op,m);
	AssignFilename dinfn("din.su");
	din.eval(dinfn);
	din.copy(mdd);
	cerr<<"norm of simulated data = "<<opeval1.getValue().norm()<<endl;
	din.linComb(-1.0,opeval1.getValue());
	cerr<<"norm of residual = "<<din.norm()<<endl;
	Vector<float> grad1(op.getDomain());
	AssignFilename gfn1("grad1.rsf");
	Components<float> cgrad1(grad1);
	cgrad1[0].eval(gfn1);
	grad1.zero();
	opeval1.getDeriv().applyAdjOp(din,grad1);
	cerr<<"grad norm from basic = "<<grad1.norm()<<endl;
      }
      else {
        alg->run();
      }      
      
      std::string dataest = valparse<std::string>(*pars,"dataest","");
      std::string datares = valparse<std::string>(*pars,"datares","");
      if (dataest.size()>0) {
	OperatorEvaluation<float> opeval(op,m);
	AssignFilename mdlfn("model.rsf");
	Vector<float> mdl(op.getDomain());
	Components<float> cmdl(mdl);
	cmdl[0].eval(mdlfn);
	mdl.copy(opeval.getPoint());

	Vector<float> est(op.getRange());
	AssignFilename estfn(dataest);
	est.eval(estfn);
	est.copy(opeval.getValue());
	if (datares.size()>0) {
	  Vector<float> res(op.getRange());
	  AssignFilename resfn(datares);
	  res.eval(resfn);
	  res.copy(est);
	  res.linComb(-1.0f,mdd);
	}
      }

      if ((outfile.size()>0) && (retrieveRank()==0)) {
	ofstream outf(outfile.c_str());
	outf<<outstr.str();
	outf.close();
      }

      if (alg) delete(alg);
    
#ifdef IWAVE_USE_MPI
    }
    MPI_Finalize();
#endif
  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
  
}
