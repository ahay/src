#include "gtest/gtest.h"
#include "acd_defn.hh"
#include "grid.h"
#include "istate.hh"
#include "adjsteptest.hh"

//#define GTEST_VERBOSE
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
    
    void create_fixed_2D_data(string fnm,
                              int nt,
                              float dt,
                              float ot,
                              int nrec,
                              int ntr_in_rec,
                              float sx0,
                              float dsx,
                              float sz,
                              float rx0,
                              float drx,
                              float rz,
                              int scalel,
                              int scalco) {
        
        string CWPROOT = getenv("CWPROOT");
        stringstream cmd;
        cmd << CWPROOT + "/bin/sunull nt=" <<nt<<" ntr="<<ntr_in_rec*nrec<<" dt="<<dt;
        cmd << "| " + CWPROOT + "/bin/sushw key=delrt a="<<ot;
        cmd << "| " + CWPROOT + "/bin/sushw key=sx a="<<sx0<<" c="<<dsx<<" j="<<ntr_in_rec;
        cmd << "| " + CWPROOT + "/bin/sushw key=selev a=" << (-1.0f)*sz;
        cmd << "| " + CWPROOT + "/bin/sushw key=gx a="<<rx0<<" b="<<drx<<" j="<<ntr_in_rec;
        cmd <<" | " + CWPROOT + "/bin/sushw key=gelev a="<<(-1.0f)*rz;
        cmd << "| " + CWPROOT + "/bin/sushw key=scalel a="<<scalel<<" | sushw key=scalco a=";
        cmd << "| " + CWPROOT + "/bin/suchw key1=offset key2=gx key3=sx b=1 c=-1 > "<<fnm<<"\n";
        
        if (system(cmd.str().c_str())) {
            RVLException e;
            e<<"Error: create_fixed_2D_data\n";
            e<<"  failed system call on command \n";
            e<<cmd.str();
            throw e;
        }
        
    }

    void create_hfile(string hfile, string dfile, grid g, float val, bool var=false) {
        
        if (hfile.size()>120) {
            RVLException e;
            e<<"Error: create_hfile\n";
            e<<"filename "<<hfile<<" longer than 120 chars\n";
            throw e;
        }
        
        char * fname=(char *)malloc(128*sizeof(char));
        FILE * fp = NULL;
        
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
        if (var) {
            for (int i=0;i<get_datasize_grid(g);i++) buf[i]=val/((float)(i+1));
        }
        else {
            for (int i=0;i<get_datasize_grid(g);i++) buf[i]=val;
        }
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
        free(fname);
    }

    void create_gauss(grid g, float lam) {
        
        char * fname=(char *)malloc(128*sizeof(char));
        FILE * fp = NULL;
        
        strcpy(fname,"gauss.rsf");
        
        fp = iwave_fopen(&fname,"w",NULL,stderr);
        if (!fp) {
            RVLException e;
            e<<"Error: create_hfile\n";
            e<<"file gauss.rsf not opened\n";
            throw e;
        }
        
        fprint_grid(fp,g);
        
        fprintf(fp,"data_format=native_float\n");
        
        fprintf(fp,"data_type = csq\n");
        fprintf(fp,"in=%s\n","gauss.rsf@");
        
        fflush(fp);
        
        iwave_fclose(fp);
        
        strcpy(fname,"gauss.rsf@");
        
        // compute center
        RPNT ctr;
        float radsq=0.0f;
        float sgm[3];
        float sgmpi=1.0f;
        const double pi = 3.14159265;
        for (int i=0;i<g.dim;i++) {
            ctr[i] = g.axes[i].o + 0.5 * (g.axes[i].d * (g.axes[i].n-1)+1);
            radsq += 0.25 * (g.axes[i].d * g.axes[i].n)*(g.axes[i].d * g.axes[i].n);
            sgm[i] = g.axes[i].d * (g.axes[i].n-1);
            sgmpi = sgmpi * sgm[i];
        }
        sgmpi = sqrt(sgmpi);
        sgmpi = sgmpi * pi;

        float x;
        float y;
        float z;
        float rsq;
        //cout << "size = "<< get_datasize_grid(g) <<endl;
        //cout << g.axes[2].n <<" * " << g.axes[1].n <<" * "
        //<<g.axes[0].n << " = " <<g.axes[2].n*g.axes[1].n*g.axes[0].n <<endl;
        float * buf = (float *)malloc(get_datasize_grid(g)*sizeof(float));
        if (g.dim==2) {
            //	cerr<<"g.axes[0].n="<<g.axes[0].n<<" g.axes[1].n="<<g.axes[1].n<<endl;
            //	cerr<<"g.axes[0].o="<<g.axes[0].o<<" g.axes[1].o="<<g.axes[1].o<<endl;
            //	cerr<<"g.axes[0].d="<<g.axes[0].d<<" g.axes[1].d="<<g.axes[1].d<<endl;
            //	cerr<<"radsq="<<radsq<<endl;
            for (int j=0;j<g.axes[1].n;j++) {
                for (int i=0;i<g.axes[0].n;i++) {
                    z=g.axes[0].o + i*g.axes[0].d;
                    x=g.axes[1].o + j*g.axes[1].d;
                    rsq=(z-ctr[0])*(z-ctr[0])+(x-ctr[1])*(x-ctr[1]);
                    //	    cerr<<"rsq/radsq="<<rsq/radsq<<endl;
                    buf[i+j*g.axes[0].n]=exp(-lam*lam*rsq/radsq);
                }
            }
        }
        else if (g.dim==3) {
            for (int k=0;k<g.axes[2].n;k++) {
                for (int j=0;j<g.axes[1].n;j++) {
                    for (int i=0;i<g.axes[0].n;i++) {
                        z=g.axes[0].o + i*g.axes[0].d;
                        x=g.axes[1].o + j*g.axes[1].d;
                        y=g.axes[2].o + k*g.axes[2].d;
                        rsq=(z-ctr[0])*(z-ctr[0])/sgm[0]/sgm[0]+(x-ctr[1])*(x-ctr[1])/sgm[1]/sgm[1]+(y-ctr[2])*(y-ctr[2])/sgm[2]/sgm[2];
                        //rsq = rsq;
                        //cerr<< i+j*(g.axes[0].n+k*g.axes[1].n) << endl;
                        //buf[i+j*g.axes[0].n]=exp(-lam*lam*rsq/radsq);
                        buf[i+(j+k*g.axes[1].n)*g.axes[0].n]=sqrt(lam*lam*lam)*exp(-0.5f*lam*lam*rsq)/sgmpi;
                    }
                }
            }
        }
        else {
            RVLException e;
            e<<"Error: create_hfile\n";
            e<<"  wake up and smell the roses\n";
            throw e;
        }
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
        free(fname);
    }
    
  class ACD3DStepTest : public ::testing::Test {
  public:
      
      string hfile_csq3d;
      string dfile_csq3d;
      string hfilem;
      string dfilem;
      string hfile1;
      string dfile1;
      string hfile2;
      string dfile2;
      string hfile3;
      string dfile3;
      string hfiles;
      string dfiles;
      string hfiled;
      string dfiled;
      
      IWaveInfo ic;
      
      ACD3DStepTest(): ic() {
          
          // hfile etc. = simple 2D
          // hfile1 etc. = x, z reversed
          // hfile2      = gdim=3, n3=3 etc.
          // hfile3      = gdim=3, h=1, z=2, x=3
          hfile_csq3d = "csq3d.rsf";
          dfile_csq3d = "csq3d.rsf@";
          hfile1= "csq1.rsf";
          dfile1= "csq1.rsf@";
          hfile2= "csq2.rsf";
          dfile2= "csq2.rsf@";
          hfile3= "csq3.rsf";
          dfile3= "csq3.rsf@";
          hfiles= "csqsmall.rsf";
          dfiles= "csqsmall.rsf@";
          hfilem= "movie3d.rsf";
          dfilem= "movie3d.rsf@";
          //     hfilei= "init.rsf";
          //     dfilei= "init.rsf@";
          
          // create simple par file
          ofstream ps("parfile");
          ps<<"INPUT DATA FOR IWAVE\n";
          ps<<"------------------------------------------------------------------------\n";
          ps<<"FD:\n";
          ps<<"\n";
          ps<<"         order = 2           scheme half-order\n";
          ps<<"           cfl = 0.5        cfl number - frac of max stable\n";
          ps<<"          cmin = 1.0         min velocity - checked\n";
          ps<<"          cmax = 5.0         max velocity - checked\n";
          ps<<"      max_step = 0           1 = set adaptively, 0 = use standard cfl from cmax\n";
          ps<<"         fpeak = 0.010       nominal central frequency \n";
          ps<<"\n";
          ps<<"------------------------------------------------------------------------\n";
          ps<<"Model info:\n";
          ps<<"\n";
          ps<<"           csq = csq3d.rsf\n";
          ps<<"        csq_d1 = csq3d.rsf\n";
          ps<<"        csq_d2 = csq3d.rsf\n";
          ps<<"        csq_b1 = csq3d.rsf\n";
          ps<<"        csq_b2 = csq3d.rsf\n";
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
          
          ofstream qs("initfile");
          qs<<"INPUT DATA FOR IWAVE\n";
          qs<<"------------------------------------------------------------------------\n";
          qs<<"FD:\n";
          qs<<"\n";
          qs<<"         order = 2           scheme half-order\n";
          qs<<"           cfl = 0.5        cfl number - frac of max stable\n";
          qs<<"          cmin = 1.0         min velocity - checked\n";
          qs<<"          cmax = 2.0         max velocity - checked\n";
          qs<<"      max_step = 0           1 = set adaptively, 0 = use standard cfl from cmax\n";
          qs<<"         fpeak = 0.010       nominal central frequency \n";
          qs<<"\n";
          qs<<"------------------------------------------------------------------------\n";
          qs<<"Model info:\n";
          qs<<"\n";
          qs<<"           csq = csq3d.rsf\n";
          qs<<"          init = gauss.rsf\n";
          qs<<"         movie = movie3d.rsf\n";
          qs<<"\n";
          qs<<"------------------------------------------------------------------------\n";
          qs<<"MPI info:\n";
          qs<<"\n";
          qs<<"       mpi_np1 = 1      n_doms along axis 1\n";
          qs<<"       mpi_np2 = 1      n_doms along axis 2\n";
          qs<<"       mpi_np3 = 1      n_doms along axis 3\n";
          qs<<"       partask = 1      task parallelization\n";
          qs<<"\n";
          qs<<"------------------------------------------------------------------------\n";
          qs<<"Output info:\n";
          qs<<"\n";
          qs<<"     printact = 1           per-time-step verbosity level\n";
          qs<<"                            0 - none\n";
          qs<<"                            1 - time step index\n";
          qs<<"                            2 - internal time step info\n";
          qs<<"                            > 5: dump everything\n";
          qs<<"      dump_pi = 1           dump parallel/dom. decomp info\n";
          qs<<"     dump_lda = 1           dump grid data for allocated arrays\n";
          qs<<"     dump_ldc = 1           dump grid data for computational arrays\n";
          qs<<"     dump_lds = 1           dump grid data for send arrays\n";
          qs<<"     dump_ldr = 1           dump grid data for receive arrays\n";
          qs<<"    dump_term = 0           dump terminator data\n";
          qs<<"    dump_pars = 0           print parameter table in IWaveOp\n";
          qs<<"   dump_steqs = 0           print major steqs in IWaveOp\n";
          qs.flush();
          qs.close();
          
          grid g;
          g.dim=3;
          g.gdim=3;
          g.axes[0].n=122;
          g.axes[1].n=98;
          g.axes[2].n=90;
          //g.axes[3].n=3;
          g.axes[0].d=25.0;
          g.axes[1].d=25.0;
          g.axes[2].d=25.0;
          //g.axes[3].d=100.0;
          g.axes[0].o=0.0;
          g.axes[1].o=0.0;
          g.axes[2].o=0.0;
          //g.axes[3].o=-100.0;
          g.axes[0].id = 0;
          g.axes[1].id = 1;
          g.axes[2].id = 2;
          //g.axes[3].id = 4;
          
          float val = 1.;
          create_hfile(hfile_csq3d, dfile_csq3d, g, val);
          float lam = 10;
          create_gauss(g,lam);

          g.dim=3;
          g.gdim=4;
          g.axes[0].n=122;
          g.axes[1].n=98;
          g.axes[2].n=90;
          g.axes[3].n=3;
          g.axes[0].d=25.0;
          g.axes[1].d=25.0;
          g.axes[2].d=25.0;
          g.axes[3].d=200.0;
          g.axes[0].o=0.0;
          g.axes[1].o=0.0;
          g.axes[2].o=1200.0;
          g.axes[3].o=2.0;
          g.axes[0].id = 0;
          g.axes[1].id = 1;
          g.axes[2].id = 2;
          g.axes[3].id = 3;
          
          val = 0.0;
          create_hfile(hfilem, dfilem, g, val);
      }
  };
    
  TEST_F(ACD3DStepTest, dryrun3d_sim_fwd_ord0_init) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=initfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=0;
      bool fwd = true;
      int snaps=0;
      int printact=0;
      
      bool dryrun=true;
      ofstream drystr("dryrun3d_sim_fwd_ord0_init");

      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps,dryrun,drystr);
      sim->run();
      drystr.close();
      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACD3DStepTest, wetrun3d_sim_fwd_ord0_init) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=initfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=0;
      bool fwd = true;
      int snaps=0;
      int printact=0;
      bool dryrun=false;
      ofstream drystr("wetrun3d_sim_fwd_ord0_init");
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps,dryrun,drystr);
      sim->run();
      drystr.close();
      delete sim;
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
char **xargv;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  int err = RUN_ALL_TESTS();
  iwave_fdestroy();
  return err;
}
