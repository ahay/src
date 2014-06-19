#include "acd_defn.hh"
#include "grid.h"
#include "gridpp.hh"
#include "gridops.hh"
#include "istate.hh"
#include "iwop.hh"
#include "functions.hh"
#include "op.hh"
#include "blockop.hh"
#include "cgnealg.hh"
#include "LBFGSBT.hh"
#include "LinFitLS.hh"
#include "acdds_grad_selfdoc.hh"
#include <par.h>
#include "gradtest.hh"
#include "scantest.hh"
#include "segyops.hh"

//#define GTEST_VERBOSE
IOKEY IWaveInfo::iwave_iokeys[]
= {
    {"csqext",    0, true,  true },
    {"data",   1, false, true },
    {"source", 1, true,  false},
    {"",       0, false, false}
};

using RVL::parse;
using RVL::RVLException;
using RVL::Vector;
using RVL::Components;
using RVL::Operator;
using RVL::OperatorEvaluation;
using RVL::LinearOp;
using RVL::LinearOpFO;
using RVL::OpComp;
using RVL::SymmetricBilinearOp;
using RVL::AssignFilename;
using RVL::AssignParams;
using RVL::RVLRandomize;
using RVL::ScaleOpFwd;
using RVL::TensorOp;
using RVL::GradientTest;
using RVL::Scan;
using RVLUmin::CGNEPolicy;
using RVLUmin::CGNEPolicyData;
using RVLUmin::LBFGSBT;
using RVLUmin::LinFitLS;
using RVLUmin::CGNEAlg;
using TSOpt::IWaveEnvironment;
using TSOpt::IWaveTree;
using TSOpt::IWaveSampler;
using TSOpt::IWaveSim;
using TSOpt::TASK_RELN;
using TSOpt::IOTask;
using TSOpt::IWaveOp;
using TSOpt::SEGYTaperMute;
using TSOpt::GridWindowOp;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
using TSOpt::MPISEGYSpace;
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
#endif
using TSOpt::GridExtendOp;
using TSOpt::GridDerivOp;

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
            
            // the Op - note that it comes equpped with domain, range spaces!
            IWaveOp iwop(*pars,stream);
            
            SEGYTaperMute tnm(valparse<float>(*pars,"mute_slope",0.0f),
                              valparse<float>(*pars,"mute_zotime",0.0f),
                              valparse<float>(*pars,"mute_width",0.0f),0,
                              valparse<float>(*pars,"min_gx",0.0f),
                              valparse<float>(*pars,"max_gx",numeric_limits<float>::max()),
                              valparse<float>(*pars,"taper_width",0.0f),0);
            
            LinearOpFO<float> tnmop(iwop.getRange(),iwop.getRange(),tnm,tnm);
        
            OpComp<float> op(iwop,tnmop);
            
            /* generate physical model space - a priori distinct from extended space
             without even the same grid */
#ifdef IWAVE_USE_MPI
            MPIGridSpace csqsp(valparse<std::string>(*pars,"csq"),"notype",true);
#else
            GridSpace csqsp(valparse<std::string>(*pars,"csq"),"notype",true);
#endif
            // make it a product, so it's compatible with domain of op
            StdProductSpace<ireal> dom(csqsp);
            // build grid extension op - checks compatibility between physical,
            // extended model-space specs
            
            GridExtendOp gext(dom,op.getDomain());
            
            // reference vel-squared model - assumed to be data used to spec
            // physical model space
            Vector<ireal> m(dom);
            AssignFilename mfn(valparse<std::string>(*pars,"csq"));
            Components<ireal> cm(m);
            cm[0].eval(mfn);
            
            // muted data - optionally archived
            Vector<ireal> mdd(op.getRange());
            std::string mddnm = valparse<std::string>(*pars,"datamut","");
            if (mddnm.size()>0) {
                AssignFilename mddfn(mddnm);
                mdd.eval(mddfn);
            }
            {
                // read data - only needed to define muted data
                Vector<ireal> dd(op.getRange());
                AssignFilename ddfn(valparse<std::string>(*pars,"data"));
                Components<ireal> cdd(dd);
                cdd[0].eval(ddfn);
                tnmop.applyOp(dd,mdd);
            }
            
            /* output stream */
            std::stringstream res,strgrad;
            res<<scientific;
            //strgrad<<scientific;
            
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
            
            // assign window widths - default = 0;
            RPNT wind;
            RASN(wind,RPNT_0);
            wind[0]=valparse<float>(*pars,"windw1",0.0f);
            wind[1]=valparse<float>(*pars,"windw2",0.0f);
            wind[2]=valparse<float>(*pars,"windw3",0.0f);
            
            GridDerivOp dsop0(op.getDomain(),dsdir,valparse<float>(*pars,"DSWeight",0.0f));
            
            // need to read in model space for bg input to GridWindowOp
            Vector<ireal> m_in(op.getDomain());
            AssignFilename minfn(valparse<std::string>(*pars,"csqext"));
            Components<ireal> cmin(m_in);
            cmin[0].eval(minfn);
            GridWindowOp wop(op.getDomain(),m_in,wind);
            
            OpComp<float> dsop(wop,dsop0);
            
            TensorOp<float> top(op,dsop);
            // create RHS of block system
            Vector<float> td(top.getRange());
            Components<float> ctd(td);
            ctd[0].copy(mdd);
            ctd[1].zero();
            
            // choice of preop is placeholder
            ScaleOpFwd<float> preop(top.getDomain(),1.0f);
            CGNEPolicyData<float> pd(valparse<float>(*pars,"ResidualTol",100.0*numeric_limits<float>::epsilon()),
                                     valparse<float>(*pars,"GradientTol",100.0*numeric_limits<float>::epsilon()),
                                     valparse<float>(*pars,"MaxStep",numeric_limits<float>::max()),
                                     valparse<int>(*pars,"MaxIter",10),true);
            Vector<float> dm0(op.getDomain());
            string refname = valparse<std::string>(*pars,"ref0");
            if (refname.size()>0){
                AssignFilename dmfn(refname);
                Components<float> cdm0(dm0);
                cdm0[0].eval(dmfn);
            }
            else { dm0.zero(); }
            LinFitLS<float, CGNEPolicy<float>, CGNEPolicyData<float> > f(top,preop,td,dm0,pd,false,res);
            // compose with grid extension op
            FcnlOpComp<float> gf(f,gext);
            
            // lower, upper bds for csq
            Vector<float> lb(gf.getDomain());
            Vector<float> ub(gf.getDomain());
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
            
            FunctionalBd<float> fbd(gf,ultest);
            
            // if gradient requested, compute and store, along with other quantities
            // at reference model
            string gradname = valparse<std::string>(*pars,"grad","");
            if (gradname.size()>0) {
                
                FunctionalEvaluation<float> Fm(fbd,m);
                // cerr<<" \n  compute gradient \n";
                AssignFilename gradfn(gradname);
                Vector<ireal> grad(dom);
                Components<ireal> cgrad(grad);
                cgrad[0].eval(gradfn);
                //strgrad << "\n getgradient norm = " << (Fm.getGradient()).norm() << endl;
                grad.copy(Fm.getGradient());
                //strgrad << "\n grad norm = " << grad.norm() << endl;
                
                Vector<ireal> dltm(op.getDomain());
                AssignFilename dltmfn(valparse<std::string>(*pars,"reflectivity"));
                Components<ireal> cdltm(dltm);
                cdltm[0].eval(dltmfn);
                dltm.zero();
                
                FunctionalBd<float> const & f1 =
                dynamic_cast<FunctionalBd<float> const &>(Fm.getFunctional()); // current clone of fbd
                FcnlOpComp<float> const & f2 =
                dynamic_cast<FcnlOpComp<float> const &>(f1.getFunctional()); // function in fbd = gf = fcnaopcomp
                FunctionalEvaluation<float> const & fe2 = f2.getFcnlEval(); // current feval part of gf
                LinFitLS<float, CGNEPolicy<float>, CGNEPolicyData<float> > const & f3 =
                dynamic_cast<LinFitLS<float, CGNEPolicy<float>, CGNEPolicyData<float> > const & >
                (fe2.getFunctional()); // current clone of LSLinFit
                dltm.copy(f3.getLSSoln()); // copy dx from LSLinFit
                
                std::string dataest = valparse<std::string>(*pars,"dataest","");
                std::string datares = valparse<std::string>(*pars,"datares","");
                std::string normalres = valparse<std::string>(*pars,"normalres","");
                if (dataest.size()>0) {
                    OperatorEvaluation<float> gopeval(gext,m);
                    OperatorEvaluation<float> opeval(op,gopeval.getValue());
                    Vector<float> est(op.getRange());
                    AssignFilename estfn(dataest);
                    est.eval(estfn);
                    opeval.getDeriv().applyOp(dltm,est);
                    if (datares.size()>0) {
                        Vector<float> res(op.getRange());
                        AssignFilename resfn(datares);
                        res.eval(resfn);
                        res.copy(est);
                        res.linComb(-1.0f,mdd);
                        if (normalres.size()>0){
                            OperatorEvaluation<float> topeval(top,gopeval.getValue());
                            Vector<float> nres(op.getDomain());
                            AssignFilename nresfn(normalres);
                            nres.eval(nresfn);
                            // create RHS of block system
                            Vector<float> tres(top.getRange());
                            Components<float> ctres(tres);
                            ctres[0].copy(res);
                            ctres[1].zero();
                            topeval.getDeriv().applyAdjOp(tres,nres);
                        }
                    }
                }
            }
            
            // if model perturbation supplied, try scan, gradient tests
            string pertname = valparse<std::string>(*pars,"csq_d1","");
            if (pertname.size()>0) {
                
                Vector<ireal> dm(dom);
                AssignFilename mfn_d1(pertname);
                Components<ireal> cm_d1(dm);
                cm_d1[0].eval(mfn_d1);
                
                // perform gradient test if nhval > 0
                if(valparse<int>(*pars,"nhval",0)>0){
                    GradientTest(fbd,m,dm,strgrad,
                                 valparse<int>(*pars,"nhval",0),
                                 valparse<float>(*pars,"hmin",0.1f),
                                 valparse<float>(*pars,"hmax",1.0f));
                }
                
                // perform scan if nscan > -1
                if(valparse<int>(*pars,"nscan",-1)>-1){
                    Scan(fbd,m,dm,
                         valparse<int>(*pars,"nscan",-1),
                         valparse<float>(*pars,"hmin",-1.0f),
                         valparse<float>(*pars,"hmax",1.0f),
                         strgrad);	     
                }
            }
            if (retrieveRank() == 0) {
                std::string outfile = valparse<std::string>(*pars,"outfile","");
                if (outfile.size()>0) {
                    ofstream outf(outfile.c_str());
                    outf<<strgrad.str();
                    outf<<"\n ================================================= \n";
                    outf<<res.str();
                    outf.close();
                }
                else {
                    cout<<strgrad.str();
                    cout<<"\n ================================================= \n";
                    cout<<res.str();
                }
            }
            
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
