#include "gtest/gtest.h"
#include "asg_defn.hh"
#include "grid.h"
#include "traceio.h"
#include "istate.hh"
#include "sgcoeffs.h" // igor's tabulated coeffs

//#define GTEST_VERBOSE

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"bulkmod",    0, true,  true },
  {"buoyancy",   1, true,  true },
  {"source",     2, true,  false},
  {"data_p",     2, false, true },
  {"data_v0",    5, false, true },
  {"data_v1",    6, false, true },
  {"movie_p",    2, false, false},
  {"movie_v0",   5, false, false},
  {"movie_v1",   6, false, false},
  {"",           0, false, false}
};

extern float * sgcoeffs(int k);

namespace {

  using RVL::valparse;
  using RVL::RVLException;
  using TSOpt::IWaveEnvironment;
  using TSOpt::IWaveTree;
  using TSOpt::IWaveSampler;
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
    cmd << "| " + CWPROOT + "/bin/sushw key=scalel a="<<scalel;
    cmd << "| " + CWPROOT + "/bin/sushw key=scalco a="<<scalco;
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

  class ASGSimTest : public ::testing::Test {
  public:

    string bulkh;
    string bulkd;
    string buoyh;
    string buoyd;
    string src;
    string data;
    grid gp;
    grid gd;

    // w/o pml
    int argc1;
    char ** argv1;

    // w pml
    int argc2;
    char ** argv2;

    IWaveInfo ic;

    ASGSimTest(): ic() {

      int nt;
      float dt;
      float ot;
      int nrec;
      int ntr_in_rec;
      float sx0;
      float dsx;
      float sz;
      float rx0;
      float drx;
      float rz;
      int scalel; 
      int scalco;

      // first segy - single gather
      nt = 101;
      dt = 0.004;
      ot = 0.0;
      nrec = 3;
      ntr_in_rec = 101;
      sx0=3300;
      dsx=100;
      sz=40;
      rx0=100;
      drx=20;
      rz=50;
      scalel=0;
      scalco=0;
	
      data = "data.su";
	
      create_fixed_2D_data(data,nt,dt,ot,
			   nrec,ntr_in_rec,
			   sx0,dsx,sz,
			   rx0,drx,rz,
			   scalel,scalco);

      // segcond segy - source pulses
      nt = 51;
      dt = 0.004;
      ot = -100.0;
      nrec = 3;
      ntr_in_rec = 1;
      sx0=3300;
      dsx=100;
      sz=40;
      rx0=0.0;
      drx=0.0;
      rz=20;
      scalel=0;
      scalco=0;
	
      src = "wavelet.su";
	
      create_fixed_2D_data(src,nt,dt,ot,
			   nrec,ntr_in_rec,
			   sx0,dsx,sz,
			   rx0,drx,rz,
			   scalel,scalco);

      init_default_grid(&gp);
      gp.dim=2;
      gp.gdim=2;
      gp.axes[0].n=416;
      gp.axes[1].n=800;
      gp.axes[0].d=25.0;
      gp.axes[1].d=25.0;
      gp.axes[0].o=0.0;
      gp.axes[1].o=0.0;
      gp.axes[0].id = 0;
      gp.axes[1].id = 1;

      init_default_grid(&gd);
      gd.dim=2;
      gd.gdim=2;
      gd.axes[0].n=418;
      gd.axes[1].n=802;
      gd.axes[0].d=25.0;
      gd.axes[1].d=25.0;
      gd.axes[0].o=-12.5;
      gd.axes[1].o=-12.50;
      //      gd.axes[0].o=-25.0;
      //      gd.axes[1].o=-25.0;
      
      gd.axes[0].id = 0;
      gd.axes[1].id = 1;
	
	
      // vp = 1.5 m/ms
      bulkh = "bulkmod.rsf";
      bulkd = "bulkmod.rsf@";
      float val = 4.5;
      create_hfile(bulkh, bulkd, gp, val);

      // density = 2 g/cc
      buoyh = "buoyancy.rsf";
      buoyd = "buoyancy.rsf@";
      val = 0.5;
      create_hfile(buoyh, buoyd, gd, val);
	
      argc1 = 19;
      argv1 = new char*[argc1];
      for (int i=0;i<argc1;i++) { argv1[i]=new char[128]; }
      strcpy(argv1[1],"order=2");
      strcpy(argv1[2],"cfl=0.5");
      strcpy(argv1[3],"cmin=1.0");
      strcpy(argv1[4],"cmax=2.0");
      strcpy(argv1[5],"fpeak=0.010");
      strcpy(argv1[6],"dt=2.0");
      strcpy(argv1[7],"bulkmod=bulkmod.rsf");
      strcpy(argv1[8],"buoyancy=buoyancy.rsf");
      strcpy(argv1[9],"source=wavelet.su");
      strcpy(argv1[10],"data_p=data.su");
      strcpy(argv1[11],"nl1=0.0");
      strcpy(argv1[12],"nl2=0.0");
      strcpy(argv1[13],"nr1=0.0");
      strcpy(argv1[14],"nr2=0.0");
      strcpy(argv1[15],"partask=1");
      strcpy(argv1[16],"dump_lda=1");
      strcpy(argv1[17],"dump_ldc=1");
      strcpy(argv1[18],"dump_term=1");

      argc2 = 19;
      argv2 = new char*[argc2];
      for (int i=0;i<argc2;i++) { argv2[i]=new char[128]; }
      strcpy(argv2[1],"order=2");
      strcpy(argv2[2],"cfl=0.5");
      strcpy(argv2[3],"cmin=1.0");
      strcpy(argv2[4],"cmax=2.0");
      strcpy(argv2[5],"fpeak=0.010");
      strcpy(argv2[6],"dt=2.0");
      strcpy(argv2[7],"bulkmod=bulkmod.rsf");
      strcpy(argv2[8],"buoyancy=buoyancy.rsf");
      strcpy(argv2[9],"source=wavelet.su");
      strcpy(argv2[10],"data_p=data.su");
      strcpy(argv2[11],"nl1=0.0");
      strcpy(argv2[12],"nl2=500.0");
      strcpy(argv2[13],"nr1=500.0");
      strcpy(argv2[14],"nr2=500.0");
      strcpy(argv2[15],"partask=1");
      strcpy(argv2[16],"dump_lda=1");
      strcpy(argv2[17],"dump_ldc=1");
      strcpy(argv2[18],"dump_term=1");
	
    }

    ~ASGSimTest() {
      for (int i=0;i<argc1; i++) delete [] argv1[i];
      delete [] argv1;
    }  
  };

  TEST_F(ASGSimTest, check_FD_coeffs) {
    for (int k=1;k<8;k++) {
      float * c = sgcoeffs(k);
      cerr<<"FD coeffs k="<<k<<endl;
      cerr<<"  computed:\n";
      for (int i=0;i<k;i++) cerr<<c[i]<<" ";
      cerr<<"\n  tabulated:\n";
      for (int i=0;i<k;i++) cerr<<SCHEME_COEFFS[k-1][i]<<" ";
      cerr<<"\n";
    }
  }

  TEST_F(ASGSimTest, setup_environment) {
    try {

      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc1, argv1, 0, &par, &stream);

#ifdef GTEST_VERBOSE
      ps_printall(*par,stderr);
#endif

      // check one of the main par entries
      std::string bulk=valparse<std::string>(*par,"bulkmod");
      EXPECT_EQ("bulkmod.rsf",bulk);
      std::string buoy=valparse<std::string>(*par,"buoyancy");
      EXPECT_EQ("buoyancy.rsf",buoy);
      
      ps_delete(&par);
      fclose(stream);

#ifndef GTEST_VERBOSE
      //      unlink("parfile");
#endif
     
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ASGSimTest, iwavetree_2D_serial_order0) {
    try {

      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc1, argv1, 0, &par, &stream);

      // build order zero IWaveTree, check 
      int order=0;
      IWaveTree * wt = new IWaveTree(*par,stream,ic,order);
      IWAVE & w = *(wt->getStateArray()[0]);
      
      iwave_printf(&w, par, stream);

      IPNT cgs0, cge0, gs0, ge0;
      ra_a_gse(&((w.model).ld_c._s[0]),gs0,ge0);
      ra_gse(&((w.model).ld_c._s[0]),cgs0,cge0);
      EXPECT_EQ(0,gs0[0]);
      EXPECT_EQ(415,ge0[0]);
      EXPECT_EQ(0,gs0[1]);
      EXPECT_EQ(799,ge0[1]);
      EXPECT_EQ(0,cgs0[0]);
      EXPECT_EQ(415,cge0[0]);
      EXPECT_EQ(0,cgs0[1]);
      EXPECT_EQ(799,cge0[1]);
      IPNT cgs1, cge1, gs1, ge1;
      ra_a_gse(&((w.model).ld_c._s[1]),gs1,ge1);
      ra_gse(&((w.model).ld_c._s[1]),cgs1,cge1);
      EXPECT_EQ(0,gs1[0]);
      EXPECT_EQ(415,ge1[0]);
      EXPECT_EQ(0,gs1[1]);
      EXPECT_EQ(799,ge1[1]);
      EXPECT_EQ(0,cgs1[0]);
      EXPECT_EQ(415,cge1[0]);
      EXPECT_EQ(0,cgs1[1]);
      EXPECT_EQ(799,cge1[1]);
      IPNT cgs2, cge2, gs2, ge2;
      ra_a_gse(&((w.model).ld_c._s[2]),gs2,ge2);
      ra_gse(&((w.model).ld_c._s[2]),cgs2,cge2);
      EXPECT_EQ(-1,gs2[0]);
      EXPECT_EQ(416,ge2[0]);
      EXPECT_EQ(1,gs2[1]);
      EXPECT_EQ(798,ge2[1]);
      EXPECT_EQ(1,cgs2[0]);
      EXPECT_EQ(414,cge2[0]);
      EXPECT_EQ(1,cgs2[1]);
      EXPECT_EQ(798,cge2[1]);
      IPNT cgs3, cge3, gs3, ge3;
      ra_a_gse(&((w.model).ld_c._s[3]),gs3,ge3);
      ra_gse(&((w.model).ld_c._s[3]),cgs3,cge3);
      EXPECT_EQ(1,gs3[0]);
      EXPECT_EQ(414,ge3[0]);
      EXPECT_EQ(-1,gs3[1]);
      EXPECT_EQ(800,ge3[1]);
      EXPECT_EQ(1,cgs3[0]);
      EXPECT_EQ(414,cge3[0]);
      EXPECT_EQ(1,cgs3[1]);
      EXPECT_EQ(798,cge3[1]);
      IPNT cgs5, cge5, gs5, ge5;
      ra_a_gse(&((w.model).ld_c._s[5]),gs5,ge5);
      ra_gse(&((w.model).ld_c._s[5]),cgs5,cge5);
      EXPECT_EQ(-1,gs5[0]);
      EXPECT_EQ(415,ge5[0]);
      EXPECT_EQ(1,gs5[1]);
      EXPECT_EQ(798,ge5[1]);
      EXPECT_EQ(0,cgs5[0]);
      EXPECT_EQ(414,cge5[0]);
      EXPECT_EQ(1,cgs5[1]);
      EXPECT_EQ(798,cge5[1]);
      IPNT cgs6, cge6, gs6, ge6;
      ra_a_gse(&((w.model).ld_c._s[6]),gs6,ge6);
      ra_gse(&((w.model).ld_c._s[6]),cgs6,cge6);
      EXPECT_EQ(1,gs6[0]);
      EXPECT_EQ(414,ge6[0]);
      EXPECT_EQ(-1,gs6[1]);
      EXPECT_EQ(799,ge6[1]);
      EXPECT_EQ(1,cgs6[0]);
      EXPECT_EQ(414,cge6[0]);
      EXPECT_EQ(0,cgs6[1]);
      EXPECT_EQ(798,cge6[1]);

      delete wt;
      ps_delete(&par);
      fclose(stream);
#ifndef GTEST_VERBOSE

#endif
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }    

  TEST_F(ASGSimTest, spacetime_grid_build_ext_before) {
    try {

      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc1, argv1, 0, &par, &stream);

      // build order zero IWaveTree, check 
      int order=0;
      bool fwd = true;

      // initialize grid
      grid g;
      copy_grid(&g,&gp);
      
      // iotask list
      std::vector<TASK_RELN *> t;
      // step 1: create list of i/o tasks

      IOTask(t,order,fwd,ic);

      IWaveTree * wt = new IWaveTree(*par,stream,ic,order);

      for (int it=0; it<t.size(); it++) {
	IWaveSampler * tmp = NULL;
	tmp = new IWaveSampler(wt->getStateArray()[0], t[it]->keyword, *par, stream);
	// note these are internal simulation axes, not archival, 
	// which will require some re-arranging in traceio
	// only test non-spatial axes per id
	for (int i=0;i<tmp->getNumAxes();i++) {
	  if (tmp->getAxis(i).id > g.dim -1) {
	    bool good=grid_union(&g,&(tmp->getAxis(i)));
	    if (!good) {
	      RVLException e;
	      e<<"Error: IWaveSim constructor from grid_union\n";
	      fprintf(stream,"Error: IWaveSim constructor from grid_union\n");
	      fprintf(stream,"  failed to add these axes:\n");
	      for (int j=0;j<tmp->getNumAxes(); j++) 
		fprint_axis(stream,tmp->getAxis(j));
	      fprintf(stream,"  to grid:\n");
	      fprint_grid(stream,g);
	      throw e;
	    }
	  }
	}
	// for this test no need to build sampler vector
	//	s.push_back(tmp);
	delete tmp;
      }

      fprint_grid(stream,g);

      EXPECT_EQ(4,g.gdim);
      EXPECT_EQ(416,g.axes[0].n);
      EXPECT_EQ(25.0,g.axes[0].d);
      EXPECT_EQ(0.0,g.axes[0].o);
      EXPECT_EQ(0,g.axes[0].id);
      EXPECT_EQ(800,g.axes[1].n);
      EXPECT_EQ(25.0,g.axes[1].d);
      EXPECT_EQ(0.0,g.axes[1].o);
      EXPECT_EQ(1,g.axes[1].id);
      EXPECT_EQ(251,g.axes[2].n);
      EXPECT_GT(1.e-6,fabs(g.axes[2].d-2.0));
      EXPECT_GT(1.e-4,fabs(-100.0-g.axes[2].o));
      EXPECT_EQ(2,g.axes[2].id);
      EXPECT_EQ(3,g.axes[3].n);
      EXPECT_EQ(1.0,g.axes[3].d);
      EXPECT_EQ(0.0,g.axes[3].o);
      EXPECT_EQ(3,g.axes[3].id);

      delete wt;
  
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ASGSimTest, iwavetree_2D_serial_pml_order0) {
    try {

      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc2, argv2, 0, &par, &stream);

      // build order zero IWaveTree, check 
      int order=0;
      IWaveTree * wt = new IWaveTree(*par,stream,ic,order);
      IWAVE & w = *(wt->getStateArray()[0]);
      fprintf(stream,"PML buffer lengths:\n");
      fprintf(stream,"nls[0]=%d\n",w.model.nls[0]);
      fprintf(stream,"nrs[0]=%d\n",w.model.nrs[0]);
      fprintf(stream,"nls[1]=%d\n",w.model.nls[1]);
      fprintf(stream,"nrs[1]=%d\n",w.model.nrs[1]);

      iwave_printf(&w, par, stream);

      IPNT cgs0, cge0, gs0, ge0;
      ra_a_gse(&((w.model).ld_c._s[0]),gs0,ge0);
      ra_gse(&((w.model).ld_c._s[0]),cgs0,cge0);
      EXPECT_EQ(0,gs0[0]);
      EXPECT_EQ(435,ge0[0]);
      EXPECT_EQ(-20,gs0[1]);
      EXPECT_EQ(819,ge0[1]);
      EXPECT_EQ(0,cgs0[0]);
      EXPECT_EQ(435,cge0[0]);
      EXPECT_EQ(-20,cgs0[1]);
      EXPECT_EQ(819,cge0[1]);
      IPNT cgs1, cge1, gs1, ge1;
      ra_a_gse(&((w.model).ld_c._s[1]),gs1,ge1);
      ra_gse(&((w.model).ld_c._s[1]),cgs1,cge1);
      EXPECT_EQ(0,gs1[0]);
      EXPECT_EQ(435,ge1[0]);
      EXPECT_EQ(-20,gs1[1]);
      EXPECT_EQ(819,ge1[1]);
      EXPECT_EQ(0,cgs1[0]);
      EXPECT_EQ(435,cge1[0]);
      EXPECT_EQ(-20,cgs1[1]);
      EXPECT_EQ(819,cge1[1]);
      IPNT cgs2, cge2, gs2, ge2;
      ra_a_gse(&((w.model).ld_c._s[2]),gs2,ge2);
      ra_gse(&((w.model).ld_c._s[2]),cgs2,cge2);
      EXPECT_EQ(-1,gs2[0]);
      EXPECT_EQ(436,ge2[0]);
      EXPECT_EQ(-19,gs2[1]);
      EXPECT_EQ(818,ge2[1]);
      EXPECT_EQ(1,cgs2[0]);
      EXPECT_EQ(434,cge2[0]);
      EXPECT_EQ(-19,cgs2[1]);
      EXPECT_EQ(818,cge2[1]);
      IPNT cgs3, cge3, gs3, ge3;
      ra_a_gse(&((w.model).ld_c._s[3]),gs3,ge3);
      ra_gse(&((w.model).ld_c._s[3]),cgs3,cge3);
      EXPECT_EQ(1,gs3[0]);
      EXPECT_EQ(434,ge3[0]);
      EXPECT_EQ(-21,gs3[1]);
      EXPECT_EQ(820,ge3[1]);
      EXPECT_EQ(1,cgs3[0]);
      EXPECT_EQ(434,cge3[0]);
      EXPECT_EQ(-19,cgs3[1]);
      EXPECT_EQ(818,cge3[1]);
      IPNT cgs5, cge5, gs5, ge5;
      ra_a_gse(&((w.model).ld_c._s[5]),gs5,ge5);
      ra_gse(&((w.model).ld_c._s[5]),cgs5,cge5);
      EXPECT_EQ(-1,gs5[0]);
      EXPECT_EQ(435,ge5[0]);
      EXPECT_EQ(-19,gs5[1]);
      EXPECT_EQ(818,ge5[1]);
      EXPECT_EQ(0,cgs5[0]);
      EXPECT_EQ(434,cge5[0]);
      EXPECT_EQ(-19,cgs5[1]);
      EXPECT_EQ(818,cge5[1]);
      IPNT cgs6, cge6, gs6, ge6;
      ra_a_gse(&((w.model).ld_c._s[6]),gs6,ge6);
      ra_gse(&((w.model).ld_c._s[6]),cgs6,cge6);
      EXPECT_EQ(1,gs6[0]);
      EXPECT_EQ(434,ge6[0]);
      EXPECT_EQ(-21,gs6[1]);
      EXPECT_EQ(819,ge6[1]);
      EXPECT_EQ(1,cgs6[0]);
      EXPECT_EQ(434,cge6[0]);
      EXPECT_EQ(-20,cgs6[1]);
      EXPECT_EQ(818,cge6[1]);

      delete wt;
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
#ifdef IWAVE_USE_MPI
  int ts=0;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif
  try {
    ::testing::InitGoogleTest(&argc, argv);
    int res = RUN_ALL_TESTS();
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
    return res;
  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
    MPI_Finalize();
#endif
    exit(1);
  }
}
