#include "gtest/gtest.h"
#include "acd_defn.hh"
#include "grid.h"
#include "istate.hh"
#include "adjsteptest.hh"

#define GTEST_VERBOSE

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq",    0, true,  true },
  {"data",   1, false, true },
  {"source", 1, true,  false},
  {"movie",  1, false, false},
  {"init",   1, true,  false},
  {"",       0, false, false}
};

namespace {

  using RVL::parse;
  using RVL::RVLException;
  using TSOpt::IWaveEnvironment;
  using TSOpt::IWaveTree;
  using TSOpt::IWaveSim;
  using TSOpt::TASK_RELN;
  using TSOpt::IOTask;

  void create_hfile(string hfile, string dfile, grid g, float val) {
    
    if (hfile.size()>120) {
      RVLException e;
      e<<"Error: create_hfile\n";
      e<<"filename "<<hfile<<" longer than 120 chars\n";
      throw e;
    }
    
    char * fname=(char *)malloc(128*sizeof(char));
    FILE * fp;

    // first try to open for read
    strcpy(fname,hfile.c_str());
    if (fp = iwave_fopen(&fname,"r",NULL,stderr)) {
      strcpy(fname,dfile.c_str());
      fp = iwave_fopen(&fname,"r",NULL,stderr);
    }
    
    // if either header or data file not present, create both
    
    if (!fp) {
      
      strcpy(fname,hfile.c_str());
      
      fp = iwave_fopen(&fname,"w",NULL,stderr);
      if (!fp) {
	RVLException e;
	e<<"Error: create_hfile\n";
	e<<"file testgrid.rsf not opened\n";
	throw e;
      }
      
      fprint_grid(fp,g);
      
      fprintf(fp,"data_format=native_float\n");

      fprintf(fp,"data_type = csq\n");
      fprintf(fp,"in=%s\n",dfile.c_str());

      fflush(fp);
      
      iwave_fclose(fp);
      
      strcpy(fname,dfile.c_str());
      float * buf = (float *)malloc(get_datasize_grid(g)*sizeof(float));
      for (int i=0;i<get_datasize_grid(g);i++) buf[i]=val;
      
      //      FILE * 
      fp = iwave_fopen(&fname,"w",NULL,stderr);
      if (!fp) {
	RVLException e;
	e<<"Error: create_hfile\n";
	e<<"file "<<fname<<" not opened\n";
	throw e;
      }
      for (int i=0;i<get_panelnum_grid(g);i++) 
	fwrite(buf,sizeof(float),get_datasize_grid(g),fp);
      
      free(buf);
      fflush(fp);

      iwave_fclose(fp);

    }
    
    free(fname);
  }

  ireal acd_rdom_inner(std::vector<RDOM *> updt, std::vector <RDOM*> old, bool fwd){
    try{
      if (updt.size() != old.size()){
	RVLException e;
	e<<"Error:acd_rdom_inner -- two input vector have different sizes.\n";
	throw e;
      }
      int n = updt.size();
        if (fwd==true) {
            ireal axy, ip;
            if (n == 2) {
                RARR *upd0_n = &(updt[1]->_s[D_UC]);
                RARR *ucd0_n = &(updt[1]->_s[D_UP]);
                RARR *upb_o = &(old[1]->_s[D_UP]);
                RARR *ucb_o = &(old[1]->_s[D_UC]);
                
                ra_a_inner(upd0_n, upb_o, &axy);
                ra_a_inner(ucd0_n, ucb_o, &ip);
                
                axy = axy + ip;
            }
            if (n == 4) {
                RARR *upd0_n = &(updt[2]->_s[D_UC]);
                RARR *ucd0_n = &(updt[2]->_s[D_UP]);
                RARR *upb_o = &(old[2]->_s[D_UP]);
                RARR *ucb_o = &(old[2]->_s[D_UC]);
                
                ra_a_inner(upd0_n, upb_o, &axy);
                ra_a_inner(ucd0_n, ucb_o, &ip);
                
                axy = axy + ip;
                
                RARR *updd_n = &(updt[n-1]->_s[D_UC]);
                RARR *ucdd_n = &(updt[n-1]->_s[D_UP]);
                RARR *updb_o = &(old[n-1]->_s[D_UP]);
                RARR *ucdb_o = &(old[n-1]->_s[D_UC]);
                
                ireal ip1;
                ra_a_inner(updd_n, updb_o, &ip1);
                ra_a_inner(ucdd_n, ucdb_o, &ip);
                
                axy = axy + ip + ip1;
            }
            return axy;
        }
        else if (fwd == false){
            ireal xaty, ip, ip1;
            if (n==2) {
                
                RARR *ucd0_o  = &(old[1]->_s[D_UC]);
                RARR *upd0_o  = &(old[1]->_s[D_UP]);
                RARR *csqd0_o = &(old[1]->_s[D_CSQ]);
                RARR *ucb_n   = &(updt[1]->_s[D_UC]);
                RARR *upb_n   = &(updt[1]->_s[D_UP]);
                RARR *csqb_n  = &(updt[1]->_s[D_CSQ]);
                
                ra_a_inner(upd0_o, upb_n, &xaty);
                ra_a_inner(csqd0_o, csqb_n, &ip);
                ra_a_inner(ucd0_o, ucb_n, &ip1);
                xaty = xaty + ip + ip1;
            }
            if (n == 4) {
                RARR *ucd0_o  = &(old[2]->_s[D_UC]);
                RARR *upd0_o  = &(old[2]->_s[D_UP]);
                RARR *csqd0_o = &(old[2]->_s[D_CSQ]);
                RARR *ucb_n   = &(updt[2]->_s[D_UC]);
                RARR *upb_n   = &(updt[2]->_s[D_UP]);
                RARR *csqb_n  = &(updt[2]->_s[D_CSQ]);
                
                ra_a_inner(upd0_o, upb_n, &xaty);
                ra_a_inner(csqd0_o, csqb_n, &ip);
                ra_a_inner(ucd0_o, ucb_n, &ip1);
                xaty = xaty + ip + ip1;
                
                RARR *updd_o  = &(old[n-1]->_s[D_UP]);
                RARR *ucdd_o  = &(old[n-1]->_s[D_UC]);
                RARR *ucdb_n  = &(updt[n-1]->_s[D_UC]);
                RARR *updb_n  = &(updt[n-1]->_s[D_UP]);
                
                ra_a_inner(ucdd_o, ucdb_n, &ip);
                ra_a_inner(updd_o, updb_n, &ip1);
                xaty = xaty + ip + ip1;
            }
            return xaty;
        }

/*      if (fwd==true) {
	RARR *upd_n = &(updt[n-1]->_s[D_UC]);
	RARR *ucd_n = &(updt[n-1]->_s[D_UP]);
	RARR *upb_o = &(old[n-1]->_s[D_UP]);
	RARR *ucb_o = &(old[n-1]->_s[D_UC]);
	ireal axy, ip;
	ra_a_inner(upd_n, upb_o, &axy);
	ra_a_inner(ucd_n, ucb_o, &ip);
	axy = axy + ip;
            
	return axy;
      }
      else if (fwd == false){
	RARR *updd_o  = &(old[n-1]->_s[D_UP]);
	RARR *ucdd_o  = &(old[n-1]->_s[D_UC]);
	RARR *csqdd_o = &(old[n-1]->_s[D_CSQ]);
	RARR *ucdb_n  = &(updt[n-1]->_s[D_UC]);
	RARR *updb_n  = &(updt[n-1]->_s[D_UP]);
	RARR *csqdb_n = &(updt[n-1]->_s[D_CSQ]);
	ireal xaty, ip, ip1;
	ra_a_inner(updd_o, updb_n, &xaty);
	ra_a_inner(csqdd_o, csqdb_n, &ip);
	ra_a_inner(ucdd_o, ucdb_n, &ip1);
	xaty = xaty + ip + ip1;
	if (n == 4) {
	  RARR *ucd0_o  = &(old[1]->_s[D_UC]);
	  RARR *csqd0_o = &(old[1]->_s[D_CSQ]);
	  RARR *ucb_n   = &(updt[1]->_s[D_UC]);
	  RARR *csqb_n  = &(updt[1]->_s[D_CSQ]);
	  ra_a_inner(ucd0_o, ucb_n, &ip);
	  ra_a_inner(csqd0_o, csqb_n, &ip1);
	  xaty = xaty + ip + ip1;
	}
	return xaty;
      }
 */   }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }
    ireal acd_rdom_inner(std::vector<RDOM *> updt, std::vector <RDOM> old, bool fwd){
        try{
            if (updt.size() != old.size()){
                RVLException e;
                e<<"Error:acd_rdom_inner -- two input vector have different sizes.\n";
                throw e;
            }
            int n = updt.size();
            if (fwd==true) {
                ireal axy, ip;
                if (n == 2) {
                RARR *upd0_n = &(updt[1]->_s[D_UC]);
                RARR *ucd0_n = &(updt[1]->_s[D_UP]);
                RARR *upb_o = &(old[1]._s[D_UP]);
                RARR *ucb_o = &(old[1]._s[D_UC]);
                
                ra_a_inner(upd0_n, upb_o, &axy);
                ra_a_inner(ucd0_n, ucb_o, &ip);
                
                axy = axy + ip;
                }
                if (n == 4) {
                    RARR *upd0_n = &(updt[2]->_s[D_UC]);
                    RARR *ucd0_n = &(updt[2]->_s[D_UP]);
                    RARR *upb_o = &(old[2]._s[D_UP]);
                    RARR *ucb_o = &(old[2]._s[D_UC]);
                    
                    ra_a_inner(upd0_n, upb_o, &axy);
                    ra_a_inner(ucd0_n, ucb_o, &ip);
                    
                    axy = axy + ip;

                    RARR *updd_n = &(updt[n-1]->_s[D_UC]);
                    RARR *ucdd_n = &(updt[n-1]->_s[D_UP]);
                    RARR *updb_o = &(old[n-1]._s[D_UP]);
                    RARR *ucdb_o = &(old[n-1]._s[D_UC]);
                    
                    ireal ip1;
                    ra_a_inner(updd_n, updb_o, &ip1);
                    ra_a_inner(ucdd_n, ucdb_o, &ip);
                    
                    axy = axy + ip + ip1;
                }
                return axy;
            }
            else if (fwd == false){
                ireal xaty, ip, ip1;
                if (n==2) {
                
                RARR *ucd0_o  = &(old[1]._s[D_UC]);
                RARR *upd0_o  = &(old[1]._s[D_UP]);
                RARR *csqd0_o = &(old[1]._s[D_CSQ]);
                RARR *ucb_n   = &(updt[1]->_s[D_UC]);
                RARR *upb_n   = &(updt[1]->_s[D_UP]);
                RARR *csqb_n  = &(updt[1]->_s[D_CSQ]);
                
                ra_a_inner(upd0_o, upb_n, &xaty);
                ra_a_inner(csqd0_o, csqb_n, &ip);
                ra_a_inner(ucd0_o, ucb_n, &ip1);
                xaty = xaty + ip + ip1;
                }
                if (n == 4) {
                    RARR *ucd0_o  = &(old[2]._s[D_UC]);
                    RARR *upd0_o  = &(old[2]._s[D_UP]);
                    RARR *csqd0_o = &(old[2]._s[D_CSQ]);
                    RARR *ucb_n   = &(updt[2]->_s[D_UC]);
                    RARR *upb_n   = &(updt[2]->_s[D_UP]);
                    RARR *csqb_n  = &(updt[2]->_s[D_CSQ]);
                    
                    ra_a_inner(upd0_o, upb_n, &xaty);
                    ra_a_inner(csqd0_o, csqb_n, &ip);
                    ra_a_inner(ucd0_o, ucb_n, &ip1);
                    xaty = xaty + ip + ip1;
                    
                    RARR *updd_o  = &(old[n-1]._s[D_UP]);
                    RARR *ucdd_o  = &(old[n-1]._s[D_UC]);
                    RARR *ucdb_n  = &(updt[n-1]->_s[D_UC]);
                    RARR *updb_n  = &(updt[n-1]->_s[D_UP]);

                    ra_a_inner(ucdd_o, ucdb_n, &ip);
                    ra_a_inner(updd_o, updb_n, &ip1);
                    xaty = xaty + ip + ip1;
                }
                return xaty;
            }
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }

  ireal acd_rdom_norm(std::vector<RDOM *> rd){
    int n=rd.size();
    RARR *upb = &(rd[n-1]->_s[D_UP]);
    RARR *ucb = &(rd[n-1]->_s[D_UC]);
        
    ireal norm, tmp;
    ra_a_inner(upb, upb, &norm);
    ra_a_inner(ucb, ucb, &tmp);
        
    norm = sqrt(norm + tmp);
    return norm;
  }
    
    ireal acd_rdom_inner_csq(RDOM *rd, int cflag){
        RARR *upb = &(rd->_s[D_UP]);
        RARR *ucb = &(rd->_s[D_UC]);

        ireal norm, tmp, tmp1;
        tmp1=0.0;
        ra_a_inner(upb, upb, &norm);
        ra_a_inner(ucb, ucb, &tmp);
        if (cflag==1) {
            RARR *csq = &(rd->_s[D_CSQ]);
            ra_a_inner(csq, csq, &tmp1);
        }
        norm = norm + tmp + tmp1;
        return norm;
    }
    ireal acd_rdom_inner_csq(RDOM rd, int cflag){
        RARR *upb = &(rd._s[D_UP]);
        RARR *ucb = &(rd._s[D_UC]);
        
        ireal norm, tmp, tmp1;
        tmp1=0.0;
        ra_a_inner(upb, upb, &norm);
        ra_a_inner(ucb, ucb, &tmp);
        if (cflag==1) {
            RARR *csq = &(rd._s[D_CSQ]);
            ra_a_inner(csq, csq, &tmp1);
        }
        norm = norm + tmp + tmp1;
        return norm;
    }
  ireal acd_rdom_norm(std::vector<RDOM> rd){
    try{
      int n=rd.size();
      RARR *upb = &(rd[n-1]._s[D_UP]);
      RARR *ucb = &(rd[n-1]._s[D_UC]);
            
      ireal norm, tmp;
      ra_a_inner(upb, upb, &norm);
      ra_a_inner(ucb, ucb, &tmp);
            
      norm = sqrt(norm + tmp);
      return norm;
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }
    
  class ACDStepTest : public ::testing::Test {
  public:

    string hfile;
    string dfile;
    IWaveInfo ic;
    ACDStepTest(): ic() {

      hfile = "csq_4layer.rsf";
      dfile = "csq_4layer.rsf@";

      // create simple par file
      ofstream ps("parfile");
      ps<<"INPUT DATA FOR IWAVE\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"FD:\n";
      ps<<"\n";
      ps<<"         order = 4           scheme half-order\n";
      ps<<"           cfl = 0.5        cfl number - frac of max stable\n";
      ps<<"          cmin = 1.0         min velocity - checked\n";
      ps<<"          cmax = 5.0         max velocity - checked\n";
      ps<<"      max_step = 0           1 = set adaptively, 0 = use standard cfl from cmax\n";
      ps<<"         fpeak = 0.010       nominal central frequency \n";
      ps<<"\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"Model info:\n";
      ps<<"\n";
      ps<<"           csq = csq_4layer.rsf\n";
      ps<<"        csq_d1 = csq_4layer.rsf\n";
      ps<<"        csq_d2 = csq_4layer.rsf\n";
      ps<<"        csq_b1 = csq_4layer.rsf\n";
      ps<<"        csq_b2 = csq_4layer.rsf\n";
      ps<<"\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"MPI info:\n";
      ps<<"\n";
      ps<<"       mpi_np1 = 1      n_doms along axis 1\n";
      ps<<"       mpi_np2 = 1      n_doms along axis 2\n";
      ps<<"       mpi_np3 = 1      n_doms along axis 3\n";
      ps<<"       partask = 1      task parallelization\n";
      ps<<"\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"Source info:\n";
      ps<<"\n";
      ps<<"        source = wavelet_fake.su\n";
      ps<<"       sampord = 1             sampling order\n";
      ps<<"\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"Trace info:\n";
      ps<<"\n";
      ps<<"            data = data_fake.su    output data file\n";
      ps<<"\"        movie1 = p\"\n";
      ps<<"\"     moviestep = 100\"\n";
      ps<<"\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"Output info:\n";
      ps<<"\n";
      ps<<"     printact = 1           per-time-step verbosity level\n";
      ps<<"                            0 - none\n";
      ps<<"                            1 - time step index\n";
      ps<<"                            2 - internal time step info\n";
      ps<<"                            > 5: dump everything\n";
      ps<<"      dump_pi = 1           dump parallel/dom. decomp info\n";
      ps<<"     dump_lda = 1           dump grid data for allocated arrays\n";
      ps<<"     dump_ldc = 1           dump grid data for computational arrays\n";
      ps<<"     dump_lds = 1           dump grid data for send arrays\n";
      ps<<"     dump_ldr = 1           dump grid data for receive arrays\n";
      ps<<"    dump_term = 0           dump terminator data\n";
      ps<<"    dump_pars = 0           print parameter table in IWaveOp\n";
      ps<<"   dump_steps = 0           print major steps in IWaveOp\n";
      ps.flush();
      ps.close();

      grid g;
      g.dim=2;
      g.gdim=2;
      g.axes[0].n=416;
      g.axes[1].n=800;
      g.axes[2].n=4;
      g.axes[0].d=25.0;
      g.axes[1].d=25.0;
      g.axes[2].d=1.0;
      g.axes[0].o=0.0;
      g.axes[1].o=0.0;
      g.axes[2].o=0.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 2;

      float val = 1.5;
      create_hfile(hfile, dfile, g, val);
    }
  };
    
  TEST_F(ACDStepTest, acd_tsf_order0test) {
    try {
      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;
            
      // build order zero IWaveTree, check
      int order=0;
      IWaveTree * fwd = new IWaveTree(*par,stream,ic,order);
      IWaveTree * adj = new IWaveTree(*par,stream,ic,order);
      std::vector<RDOM *> rdfwd, rdadj;
            
      IWAVE * w = (fwd->getStateArray()[0]);
            
      //      FD_MODEL *specs = (FD_MODEL *)(w->model.specs);
            
      rdfwd = fwd->getRDOMArray();
      rdadj = adj->getRDOMArray();
            
      for (int i=0; i<rdfwd.size(); i++) {
	rdom_rand(rdfwd[i]);
	rdom_rand(rdadj[i]);
      }
      //      acd_timestep(rdfwd, true,  0, specs->fdpars);
      acd_timestep(rdfwd, true,  0, w->model.specs);
      try {
	acd_timestep(rdadj, false, 0, w->model.specs);
      }
      catch (RVLException & e) {
	std::string ref = "Error: acd_timestep().  iw.size() = 1 has no adjoint!\n";
	std::stringstream err;
	e.write(err);
	EXPECT_EQ(ref,err.str());
      }

      delete fwd;
      delete adj;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
	  
      e.write(cerr);
      exit(1);
    }
  }    

  TEST_F(ACDStepTest, acd_tsf_deriv1_adjtest_iwave) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;
      srand(time(NULL));
      
      // build order one IWaveTree, check
      int order=1;
      IWaveTree * fwd = new IWaveTree(*par,stream,ic,order);
      IWaveTree * adj = new IWaveTree(*par,stream,ic,order);
      std::vector<RDOM *> rdfwd, rdadj;

      IWAVE * w = (fwd->getStateArray()[0]);

      // FD_MODEL *specs = (FD_MODEL *)(w->model.specs);
        
      rdfwd = fwd->getRDOMArray();
      rdadj = adj->getRDOMArray();
      for (int i=0; i<rdfwd.size(); i++) {
        //rdom_rand(rdfwd[i]);
        //rdom_rand(rdadj[i]);
          iwave_rdom_rand(fwd->getStateArray()[i]);
          iwave_rdom_rand(adj->getStateArray()[i]);
      }
      // initialize rdadj[1], D_CSQ to zero for adjoint operator
      //rdom_copy(rdfwd[0],rdadj[0]);
      iwave_rdom_copy(fwd->getStateArray()[0],adj->getStateArray()[0]);
      ra_a_zero(&(rdadj[1]->_s[D_CSQ]));
        
        IWaveTree * fwdcp = new IWaveTree(*par,stream,ic,order);
        IWaveTree * adjcp = new IWaveTree(*par,stream,ic,order);

      for (int i=0; i<rdfwd.size();i++) {
        iwave_rdom_copy(fwd->getStateArray()[i],fwdcp->getStateArray()[i]);
        iwave_rdom_copy(adj->getStateArray()[i],adjcp->getStateArray()[i]);
      }

        std::vector<RDOM *> rdfwdcp, rdadjcp;

        
        rdfwdcp = fwdcp->getRDOMArray();
        rdadjcp = adjcp->getRDOMArray();

      ra_a_swap(&(rdadj[1]->_s[D_UC]),&(rdadj[1]->_s[D_UP]));

      acd_timestep(rdfwd, true,  0, w->model.specs);
      acd_timestep(rdadj, false, 0, w->model.specs);
        
      ireal yn = acd_rdom_norm(rdadjcp);
      ireal axn = acd_rdom_norm(rdfwd);
      ireal axy = acd_rdom_inner(rdfwd, rdadjcp, true);
      ireal xaty = acd_rdom_inner(rdadj, rdfwdcp, false);
      ireal adjrl = abs(axy - xaty)/axn/yn;
        
#ifdef GTEST_VERBOSE
      cout << "        ||y|| = " << yn << endl;
      cout << "       ||Ax|| = " << axn << endl;
      cout <<  "< Ax,    y > = " << axy<< endl;
      cout << "<  x, A^Ty > = " << xaty<< endl;
      cout << "< Ax, y> - < x, A^Ty> \n";
      cout << "--------------------- = " << abs(axy - xaty)/axn/yn << endl;
      cout << "       |Ax||y|        \n";
#endif
        
      ireal eps = 200 * numeric_limits<float>::epsilon();
      EXPECT_LE(adjrl, eps);
        
        delete fwdcp;
        delete adjcp;
      delete fwd;
      delete adj;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }
    
    TEST_F(ACDStepTest, acd_tsf_deriv1_adjtest) {
        try {
            
            // fake command line environment
            int argc = 2;
            char * argvv = new char[128];
            char ** argv = new char*[2];
            argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
            strcpy(argv[1],"par=parfile");
            PARARRAY * par = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &par, &stream);
            delete [] argv;
            delete [] argvv;
            srand(time(NULL));
            
            // build order one IWaveTree, check
            int order=1;
            int err;
            IWaveTree * fwd = new IWaveTree(*par,stream,ic,order);
            IWaveTree * adj = new IWaveTree(*par,stream,ic,order);
            std::vector<RDOM *> rdfwd, rdadj;

      cerr<<"1\n";
        
            
            IWAVE * w = (fwd->getStateArray()[0]);
            
            // FD_MODEL *specs = (FD_MODEL *)(w->model.specs);
            
            rdfwd = fwd->getRDOMArray();
            rdadj = adj->getRDOMArray();
	    cerr<<"2\n";
        

            for (int i=0; i<rdfwd.size(); i++) {
                //rdom_rand(rdfwd[i]);
                //rdom_rand(rdadj[i]);
                iwave_rdom_rand(fwd->getStateArray()[i]);
                iwave_rdom_rand(adj->getStateArray()[i]);
            }
      cerr<<"3\n";
        

            // initialize rdadj[1], D_CSQ to zero for adjoint operator
            //rdom_copy(rdfwd[0],rdadj[0]);
            iwave_rdom_copy(fwd->getStateArray()[0],adj->getStateArray()[0]);
            ra_a_zero(&(rdadj[1]->_s[D_CSQ]));
      cerr<<"4\n";
        

            int ndim = (w->model).g.dim;
            
            IPNT dgs[RDOM_MAX_NARR], dge[RDOM_MAX_NARR];    /*< computational domain */
            IPNT dgsa[RDOM_MAX_NARR], dgea[RDOM_MAX_NARR];  /*< allocated domain */
            for (int i=0; i<RDOM_MAX_NARR; i++) {
	      IASN(dgs[i],IPNT_1); IASN(dgsa[i],IPNT_1);
	      IASN(dge[i],IPNT_0); IASN(dgea[i],IPNT_0);
	      rd_a_gse(&((w->model).ld_c),i,dgs[i],dge[i]);
              rd_gse(&((w->model).ld_c),i,dgsa[i],dgea[i]);
            }

            std::vector<RDOM > fwdcp(rdfwd.size());
            std::vector<RDOM > adjcp(rdadj.size());
            for (int i=0; i<rdfwd.size();i++) {
                /*-declare computational domain---------------------------------------------*/

                err = rd_a_create(&(fwdcp[i]), RDOM_MAX_NARR, dgs, dge);
                err = err || rd_a_create(&(adjcp[i]), RDOM_MAX_NARR, dgs, dge);
                
                if ( err ){
                    RVLException e;
                    e<<"Error: ACDStepTest - acdfwd\n";
                    e<<"  error from rd_a_create = "<<err<<"\n";
                    throw e;
                }
                
                rdom_copy(rdfwd[i],&fwdcp[i]);
                rdom_copy(rdadj[i],&adjcp[i]);
            }
            ra_a_swap(&(rdadj[1]->_s[D_UC]),&(rdadj[1]->_s[D_UP]));
            
            acd_timestep(rdfwd, true,  0, w->model.specs);
            acd_timestep(rdadj, false, 0, w->model.specs);
            
            
            ireal yn = acd_rdom_norm(adjcp);
            ireal axn = acd_rdom_norm(rdfwd);
            ireal axy = acd_rdom_inner(rdfwd, adjcp, true);
            ireal xaty = acd_rdom_inner(rdadj, fwdcp, false);
            ireal adjrl = abs(axy - xaty)/axn/yn;
            
#ifdef GTEST_VERBOSE
            cout << "        ||y|| = " << yn << endl;
            cout << "       ||Ax|| = " << axn << endl;
            cout <<  "< Ax,    y > = " << axy<< endl;
            cout << "<  x, A^Ty > = " << xaty<< endl;
            cout << "< Ax, y> - < x, A^Ty> \n";
            cout << "--------------------- = " << abs(axy - xaty)/axn/yn << endl;
            cout << "       |Ax||y|        \n";
#endif
            
            ireal eps = 200 * numeric_limits<float>::epsilon();
            EXPECT_LE(adjrl, eps);
            for (int i=0; i<rdfwd.size();i++) {
                rd_a_destroy(&(fwdcp[i]));
                rd_a_destroy(&(adjcp[i]));
            }
            
            delete fwd;
            delete adj;
            ps_delete(&par);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }    

    
  TEST_F(ACDStepTest, acd_tsf_deriv2_adjtest) {
    try {
            
      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;
        
      srand(time(NULL));
      // build order one IWaveTree, check
      int order=2;
      int err;
      IWaveTree * fwd = new IWaveTree(*par,stream,ic,order);
      IWaveTree * adj = new IWaveTree(*par,stream,ic,order);
      std::vector<RDOM *> rdfwd, rdadj;
            
      IWAVE * w = (fwd->getStateArray()[0]);
            
      // FD_MODEL *specs = (FD_MODEL *)(w->model.specs);
            
      rdfwd = fwd->getRDOMArray();
      rdadj = adj->getRDOMArray();
            
      for (int i=0; i<rdfwd.size(); i++) {
	    //rdom_rand(rdfwd[i]);
	    //rdom_rand(rdadj[i]);
          iwave_rdom_rand(fwd->getStateArray()[i]);
          iwave_rdom_rand(adj->getStateArray()[i]);
      }

        iwave_rdom_copy(fwd->getStateArray()[0],adj->getStateArray()[0]);
        iwave_rdom_copy(fwd->getStateArray()[1],adj->getStateArray()[1]);

      // initialize rdadj[1], D_CSQ to zero for adjoint operator
      ra_a_zero(&(rdadj[2]->_s[D_CSQ]));
      ra_a_zero(&(rdadj[2]->_s[D_UC]));
      ra_a_zero(&(rdadj[3]->_s[D_CSQ]));
            
      int ndim = (w->model).g.dim;
            
      IPNT dgs[RDOM_MAX_NARR], dge[RDOM_MAX_NARR];    /*< computational domain */
      IPNT dgsa[RDOM_MAX_NARR], dgea[RDOM_MAX_NARR];  /*< allocated domain */
      for (int i=0; i<RDOM_MAX_NARR; i++) {
	rd_a_gse(&((w->model).ld_a),i,dgs[i],dge[i]);
	rd_gse(&((w->model).ld_a),i,dgsa[i],dgea[i]);
      }
      std::vector<RDOM > fwdcp(rdfwd.size());
      std::vector<RDOM > adjcp(rdadj.size());
      for (int i=0; i<rdfwd.size();i++) {
	/*-declare computational domain---------------------------------------------*/
	err = rd_a_create(&(fwdcp[i]), RDOM_MAX_NARR, dgs, dge);
	err = err || rd_a_create(&(adjcp[i]), RDOM_MAX_NARR, dgs, dge);
                
	if ( err ){
	  RVLException e;
	  e<<"Error: ACDStepTest - acdfwd\n";
	  e<<"  error from rd_a_create = "<<err<<"\n";
	  throw e;
	}
	rdom_copy(rdfwd[i],&fwdcp[i]);
	rdom_copy(rdadj[i],&adjcp[i]);
      }
            
      ra_a_swap(&(rdadj[2]->_s[D_UC]),&(rdadj[2]->_s[D_UP]));
      ra_a_swap(&(rdadj[3]->_s[D_UC]),&(rdadj[3]->_s[D_UP]));

      acd_timestep(rdfwd, true,  0, w->model.specs);
      acd_timestep(rdadj, false, 0, w->model.specs);
            
       // ireal yn = acd_rdom_norm(adjcp);

      ireal yn = sqrt(acd_rdom_inner_csq(adjcp[2],0)+acd_rdom_inner_csq(adjcp[3],0));
      ireal axn = sqrt(acd_rdom_inner_csq(rdfwd[2],0)+acd_rdom_inner_csq(rdfwd[3],0));
      ireal axy = acd_rdom_inner(rdfwd, adjcp, true);
      ireal xaty = acd_rdom_inner(rdadj, fwdcp, false);
      ireal adjrl = abs(axy - xaty)/axn/yn;

#ifdef GTEST_VERBOSE
      cout << "        ||y|| = " << yn << endl;
      cout << "       ||Ax|| = " << axn << endl;
      cout <<  "< Ax,    y > = " << axy<< endl;
      cout << "<  x, A^Ty > = " << xaty<< endl;
      cout << "< Ax, y> - < x, A^Ty> \n";
      cout << "--------------------- = " << abs(axy - xaty)/axn/yn << endl;
      cout << "       |Ax||y|        \n";
#endif
            
      ireal eps = 200 * numeric_limits<float>::epsilon();
      EXPECT_LE(adjrl, eps);
      for (int i=0; i<rdfwd.size();i++) {
          rd_a_destroy(&(fwdcp[i]));
          rd_a_destroy(&(adjcp[i]));
      }
            
      delete fwd;
      delete adj;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

}

int xargc; 
char ** xargv;

int main(int argc, char **argv) {
  int ts=0;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
  try {
    ::testing::InitGoogleTest(&argc, argv);
    int res = RUN_ALL_TESTS();
    MPI_Finalize();
    return res;
  }
  catch (RVLException & e) {
    e.write(cerr);
    MPI_Abort(MPI_COMM_WORLD,0);
    MPI_Finalize();
    exit(1);
  }
}

