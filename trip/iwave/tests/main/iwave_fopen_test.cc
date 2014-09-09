#include "gtest/gtest.h"
#include "acd_defn.hh"
#include "grid.h"
#include "istate.hh"
#include "iwop.hh"
#include "gridpp.hh"

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
  using RVL::Vector;
  using RVL::Operator;
  using RVL::OperatorEvaluation;
  using RVL::LinearOp;
  using RVL::SymmetricBilinearOp;
  using RVL::AssignFilename;
  using RVL::AssignParams;
  using TSOpt::IWaveEnvironment;
  using TSOpt::IWaveTree;
  using TSOpt::IWaveSampler;
  using TSOpt::IWaveSim;
  using TSOpt::TASK_RELN;
  using TSOpt::IOTask;
  using TSOpt::IWaveOp;
  using TSOpt::GridSpace;

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
    
    for (int i=0;i<g.dim;i++) {
      ctr[i] = g.axes[i].o + 0.5 * g.axes[i].d * g.axes[i].n;
      radsq += (0.5 * g.axes[i].d * g.axes[i].n)*(0.5 * g.axes[i].d * g.axes[i].n);
    }
    float x;
    float y; 
    float z;
    float rsq;
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
	    buf[i+j*g.axes[0].n]=exp(-lam*lam*((z-ctr[0])*(z-ctr[0])+(x-ctr[1])*(x-ctr[1])+(y-ctr[2])*(y-ctr[2])));
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

  class ACDSimTest : public ::testing::Test {
  public:

    string hfile;
    string dfile;
    string hfile1;
    string dfile1;
    string hfile2;
    string dfile2;
    string hfile3;
    string dfile3;
    string hfiles;
    string dfiles;
    string hfilem;
    string dfilem;
    string hfilei;
    string dfilei;

    IWaveInfo ic;

    ACDSimTest(): ic() {

      // hfile etc. = simple 2D
      // hfile1 etc. = x, z reversed
      // hfile2      = gdim=3, n3=3 etc.
      // hfile3      = gdim=3, h=1, z=2, x=3
      hfile = "csq_4layer.rsf";
      dfile = "csq_4layer.rsf@";
      hfile1= "csq1.rsf";
      dfile1= "csq1.rsf@";
      hfile2= "csq2.rsf";
      dfile2= "csq2.rsf@";
      hfile3= "csq3.rsf";
      dfile3= "csq3.rsf@";
      hfiles= "csqsmall.rsf";
      dfiles= "csqsmall.rsf@";
      hfilem= "movie.rsf";
      dfilem= "movie.rsf@";
      hfilei= "init.rsf";
      dfilei= "init.rsf@";

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
      qs<<"            dt = 2.0";
      qs<<"\n";
      qs<<"------------------------------------------------------------------------\n";
      qs<<"Model info:\n";
      qs<<"\n";
      qs<<"           csq = csq_4layer.rsf\n";
      qs<<"          init = gauss.rsf\n";
      qs<<"         movie = movie.rsf\n";
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
      g.dim=2;
      g.gdim=2;
      g.axes[0].n=416;
      g.axes[1].n=800;
      g.axes[2].n=3;
      g.axes[0].d=25.0;
      g.axes[1].d=25.0;
      g.axes[2].d=100.0;
      g.axes[0].o=0.0;
      g.axes[1].o=0.0;
      g.axes[2].o=-100.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 3;

      float val = 1.5;
      create_hfile(hfile, dfile, g, val);

      g.dim=2;
      g.gdim=2;
      g.axes[0].n=416;
      g.axes[1].n=800;
      g.axes[0].d=25.0;
      g.axes[1].d=25.0;
      g.axes[0].o=0.0;
      g.axes[1].o=0.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;

      val = 1.0;
      create_hfile(hfilei, dfilei, g, val);

      g.axes[0].n=800;
      g.axes[1].n=416;      
      g.axes[0].id = 1;
      g.axes[1].id = 0;

      val=1.5;
      create_hfile(hfile1, dfile1, g, val);

      g.gdim=3;
      g.axes[0].n=416;
      g.axes[1].n=800;      
      g.axes[2].n=3;
      g.axes[0].d=25.0;
      g.axes[1].d=25.0;
      g.axes[2].d=100;
      g.axes[0].o=0.0;
      g.axes[1].o=0.0;
      g.axes[2].o=-100.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 3;
      create_hfile(hfile2, dfile2, g, val);

      g.axes[0].n=3;
      g.axes[0].d=100.0;
      g.axes[0].o=-100.0;
      g.axes[1].n=416;
      g.axes[1].d=25.0;
      g.axes[1].o=0.0;
      g.axes[2].n=800;
      g.axes[2].d=25.0;
      g.axes[2].o=0.0;
      g.axes[0].id = 101;
      g.axes[1].id = 0;
      g.axes[2].id = 1;
      create_hfile(hfile3, dfile3, g, val);

      g.dim=2;
      g.gdim=2;
      g.axes[0].n=5;
      g.axes[1].n=10;
      g.axes[2].n=1;
      g.axes[0].d=25.0;
      g.axes[1].d=25.0;
      g.axes[2].d=1.0;
      g.axes[0].o=0.0;
      g.axes[1].o=0.0;
      g.axes[2].o=0.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 2;
      create_hfile(hfiles, dfiles, g, 1.0, true);

      g.dim=2;
      g.gdim=3;
      g.axes[0].n=416;
      g.axes[1].n=800;
      g.axes[2].n=1;
      g.axes[0].d=25.0;
      g.axes[1].d=25.0;
      g.axes[2].d=2.0;
      g.axes[0].o=0.0;
      g.axes[1].o=0.0;
      g.axes[2].o=2.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 2;

      val = 0.0;
      create_hfile(hfilem, dfilem, g, val);

      g.dim=2;
      g.gdim=3;
      g.axes[0].n=416;
      g.axes[1].n=800;
      g.axes[2].n=1;
      g.axes[0].d=25.0;
      g.axes[1].d=25.0;
      g.axes[2].d=2.0;
      g.axes[0].o=0.0;
      g.axes[1].o=0.0;
      g.axes[2].o=0.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 2;

      float lam = 10;
      create_gauss(g,lam);

    }
  };


  TEST_F(ACDSimTest, iwave_fopen_getproto) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;

      string filename;
      filename = "csq_4layer.rsf";
      char * cname;
      cname = new char[128];
      strcpy(cname,filename.c_str());
      FILE * fp1 = iwave_fopen(&cname,"r",NULL,stderr);
      strcpy(cname,"csq4.rsf");
      FILE * fp2 = iwave_fopen(&cname,"r",filename.c_str(),stderr);
      char * nname=NULL;
      FILE * fp3 = iwave_fopen(&nname,"w+",filename.c_str(),stderr);
      //      iwave_fprintall(stderr);
      const char* pname;
      //      cerr<<"file = "<<filename<<" proto = ";
      pname = iwave_getproto(filename.c_str());
      /*
      if (pname) cerr<<pname<<endl;
      else cerr<<"NULL\n";
      cerr<<"file = "<<cname<<" proto = ";
      */
      pname = iwave_getproto(cname);
      /*
      if (pname) cerr<<pname<<endl;
      else cerr<<"NULL\n";
      cerr<<"file = "<<nname<<" proto = ";
      */
      /*
      pname = iwave_getproto(nname);
      if (pname) cerr<<pname<<endl;
      else cerr<<"NULL\n";
      */
      if (fp1) fclose(fp1);
      if (fp2) fclose(fp2);
      if (fp3) fclose(fp3);
      ps_delete(&pars);
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
