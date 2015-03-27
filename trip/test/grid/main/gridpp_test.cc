#include "gridpp.hh"
#include "gridops.hh"
#include "gtest/gtest.h"
#include "adjtest.hh"
#include "functions.hh"
#include "mpiserialfo.hh"

//#define VERBOSE

#ifdef IWAVE_USE_MPI
typedef TSOpt::MPIGridSpace myGridSpace;
#else
typedef TSOpt::GridSpace myGridSpace;
#endif

namespace {

  using TSOpt::GridDCF;
  using TSOpt::GridDC;
  using TSOpt::GridSpace;
  using TSOpt::GridWindowOp;
  using TSOpt::GridDerivOp;
  using TSOpt::GridExtendOp;
  using RVL::AssignFilename;
  using RVL::parse;
  using RVL::RVLException;
  using RVL::RVLAssignConst;
  using RVL::RVLMax;
  using RVL::RVLMin;
  using RVL::MPISerialFunctionObjectRedn;
  using RVL::OperatorEvaluation;
  using RVL::RVLRandomize;

  void create_hfile(string hfile, grid g) {
    cerr<<"begin create_hfile hfile="<<hfile<<"\n";
    if (hfile.size()>120) {
      RVLException e;
      e<<"Error: create_hfile\n";
      e<<"filename "<<hfile<<" longer than 120 chars\n";
      throw e;
    }
    
    cerr<<"-3\n";
    char * fname=(char *)malloc(128*sizeof(char));
    FILE * fp1 = NULL;
    FILE * fp2 = NULL;
    string dfile = hfile+"@";

    cerr<<"-2\n";
    // if either header or data file not present, create both
    
    strcpy(fname,hfile.c_str());
    
    cerr<<"-1 fname="<<fname<<"\n";

    //    fp = iwave_fopen(&fname,"w",NULL,stderr);
    fp1 = fopen(fname,"w");
    cerr<<"0\n";
    if (!fp1) {
      RVLException e;
      e<<"Error: create_hfile\n";
      e<<"file testgrid.rsf not opened\n";
      throw e;
    }
    
    cerr<<"1 fp1="<<fp1<<"\n";
    fprint_grid(fp1,g);
    cerr<<"1.5\n";
    fprintf(fp1,"data_format=native_float\n");
    
    fprintf(fp1,"in=%s\n",dfile.c_str());
    cerr<<"2\n";   
    fflush(fp1);
    fclose(fp1);   
    //    iwave_fclose(fp);
    
    strcpy(fname,dfile.c_str());
    float * buf = (float *)malloc(get_extended_datasize_grid(g)*sizeof(float));
    for (int i=0;i<get_extended_datasize_grid(g);i++) buf[i]=0.0;
    cerr<<"3\n";      
    //      FILE * 
    //    fp = iwave_fopen(&fname,"w",NULL,stderr);
    fp2 = fopen(fname,"w");
    if (!fp2) {
      RVLException e;
      e<<"Error: create_hfile\n";
      e<<"file "<<fname<<" not opened\n";
      throw e;
    }
    for (int i=0;i<get_panelnum_grid(g);i++) 
      fwrite(buf,sizeof(float),get_extended_datasize_grid(g),fp2);
    
    free(buf);
    cerr<<"4\n";
    //    iwave_fclose(fp);
    fflush(fp2);
    fclose(fp2);
    
    
    cerr<<"5\n";
    /*
    // reopen prototype files for read
    strcpy(fname,hfile.c_str());
    fp = iwave_fopen(&fname,"r",NULL,stderr);
    strcpy(fname,dfile.c_str());
    fp = iwave_fopen(&fname,"r",NULL,stderr);
    */
    fprintf(stderr,"end create_hfil\n");
    free(fname);    
  }
    
  /*
    float dot() {
    int n1=11;
    int n2=21;
    int n3=41;
    float d1=0.1;
    float d2=0.1;
    float d3=0.1;
    return n1*n2*n3*d1*d2*d3;
    }
  */

  class GridTest : public ::testing::Test {
  protected:
    
    // You can remove any or all of the following functions if its body
    // is empty.
    grid g;
    grid g0;
    grid g1;
    grid g2;
    grid g3;
    grid gshort;
    grid eg0;
    grid eg1;
    string fname;
  
    GridTest() {

      init_default_grid(&g);
      init_default_grid(&g0);
      init_default_grid(&g1);
      init_default_grid(&g2);
      init_default_grid(&g3);
      init_default_grid(&gshort);
      init_default_grid(&eg0);
      init_default_grid(&eg1);

      g.dim=3;
      g.gdim=3;
      g.axes[0].n=11;
      g.axes[1].n=21;
      g.axes[2].n=41;
      g.axes[0].d=0.1;
      g.axes[1].d=0.1;
      g.axes[2].d=1.0;
      g.axes[0].o=-0.5;
      g.axes[1].o=-1.0;
      g.axes[2].o=0.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 2;

      g0.gdim = 2;
      g0.dim  = 2;

      g0.axes[0].n=10;
      g0.axes[0].d=1.0;
      g0.axes[0].o=0.0;
      g0.axes[0].id=0;

      g0.axes[1].n=20;
      g0.axes[1].d=2.0;
      g0.axes[1].o=10.0;
      g0.axes[1].id = 1;

      // g1: internal extension, axis 3 = shot axis
      g1.gdim = 3;
      g1.dim  = 2;

      g1.axes[0].n=10;
      g1.axes[0].d=1.0;
      g1.axes[0].o=0.0;
      g1.axes[0].id=0;

      g1.axes[1].n=20;
      g1.axes[1].d=2.0;
      g1.axes[1].o=10.0;
      g1.axes[1].id = 1;

      g1.axes[2].n=30;
      g1.axes[2].d=3.0;
      g1.axes[2].o=100.0;
      g1.axes[2].id = 3;

      // gshort: internal extension, axis 3 = itty bitty shot axis
      gshort.gdim = 3;
      gshort.dim  = 2;

      gshort.axes[0].n=10;
      gshort.axes[0].d=1.0;
      gshort.axes[0].o=0.0;
      gshort.axes[0].id=0;

      gshort.axes[1].n=20;
      gshort.axes[1].d=2.0;
      gshort.axes[1].o=10.0;
      gshort.axes[1].id = 1;

      gshort.axes[2].n=2;
      gshort.axes[2].d=1.0;
      gshort.axes[2].o=0.0;
      gshort.axes[2].id = 3;

      // g2 same as g1

      g2.gdim = 3;
      g2.dim = 2;

      g2.axes[0].n=10;
      g2.axes[0].d=1.0;
      g2.axes[0].o=0.0;
      g2.axes[0].id=0;

      g2.axes[1].n=20;
      g2.axes[1].d=2.0;
      g2.axes[1].o=10.0;
      g2.axes[1].id=1;

      g2.axes[2].n=30;
      g2.axes[2].d=3.0;
      g2.axes[2].o=100.0;
      g2.axes[2].id=3;

      // window test grid - window in 
      // physical dimensions
      g3.gdim = 3;
      g3.dim = 2;

      g3.axes[0].n=6;
      g3.axes[0].d=1.0;
      g3.axes[0].o=2.0;
      g3.axes[0].id=0;

      g3.axes[1].n=14;
      g3.axes[1].d=2.0;
      g3.axes[1].o=6.0;
      g3.axes[1].id=1;

      g3.axes[2].n=30;
      g3.axes[2].d=3.0;
      g3.axes[2].o=100.0;
      g3.axes[2].id=3;

      // ext grid 0 - same as g1 but axis 3 is internal
      eg0.gdim = 3;
      eg0.dim = 2;

      eg0.axes[0].n=10;
      eg0.axes[0].d=1.0;
      eg0.axes[0].o=2.0;
      eg0.axes[0].id=0;

      eg0.axes[1].n=20;
      eg0.axes[1].d=2.0;
      eg0.axes[1].o=6.0;
      eg0.axes[1].id=1;

      eg0.axes[2].n=30;
      eg0.axes[2].d=3.0;
      eg0.axes[2].o=100.0;
      eg0.axes[2].id=100;

      // ext grid 1 - same as g3 (window in spatial axes) 
      // but axis 3 is internal
      eg1.gdim = 3;
      eg1.dim = 2;

      eg1.axes[0].n=6;
      eg1.axes[0].d=1.0;
      eg1.axes[0].o=2.0;
      eg1.axes[0].id=0;

      eg1.axes[1].n=14;
      eg1.axes[1].d=2.0;
      eg1.axes[1].o=6.0;
      eg1.axes[1].id=1;

      eg1.axes[2].n=30;
      eg1.axes[2].d=3.0;
      eg1.axes[2].o=100.0;
      eg1.axes[2].id=100;
    }

    virtual ~GridTest() {
      // You can do clean-up work that doesn't throw exceptions here.
      //      unlink(fname.c_str());
    }

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    virtual void SetUp() {
      // Code here will be called immediately after the constructor (right
      // before each test).
    }

    virtual void TearDown() {
      // Code here will be called immediately after each test (right
      // before the destructor).
    }

    // Objects declared here can be used by all tests in the test case for Grid.
  };


  TEST_F(GridTest, get_datasize_grid) {
    EXPECT_EQ(200, get_datasize_grid(g1));
  }

  TEST_F(GridTest, get_global_datasize_grid) {
    EXPECT_EQ(6000, get_global_datasize_grid(g1));
  }

  TEST_F(GridTest, compare_grid) {
    EXPECT_EQ(0,compare_grid(g1,g2));
  }

  TEST_F(GridTest, GridDCF_getGrid) {

    // create test file
    fname="test0.rsf";
    create_hfile(fname,g);
    //    cerr<<"GridDCF:::\n";
    //    fprint_grid(stderr,g);
    //    cerr<<"GridDCF:::fac.getGrid\n";
    GridDCF fac(fname);
    //    fprint_grid(stderr,fac.getGrid());
    EXPECT_EQ(0,compare_grid(g,fac.getGrid()));
    iwave_fdestroy();
    unlink(fname.c_str());
    string dname=fname+"@";
    unlink(dname.c_str());
  }

  TEST_F(GridTest, GridDCF_getCellVol) {
    // create test file
    fname="test0.rsf";
    create_hfile(fname,g);
    GridDCF fac(fname);
    EXPECT_GT(1.0e-06,fabs(0.01-fac.getCellVol()));
    iwave_fdestroy();
    unlink(fname.c_str());
    string dname=fname+"@";
    unlink(dname.c_str());
  }

  TEST_F(GridTest, GridDCF_getFilename) {
    // create test file
    fname="test0.rsf";
    create_hfile(fname,g);
    GridDCF fac(fname);
    EXPECT_EQ(fname,fac.getFilename());
    iwave_fdestroy();
    unlink(fname.c_str());
    string dname=fname+"@";
    unlink(dname.c_str());
  }

  TEST_F(GridTest, GridSpace_getGrid) {

    // create test file
    fname="test0.rsf";
    create_hfile(fname,g);
    myGridSpace sp(fname);
    EXPECT_EQ(0,compare_grid(g,sp.getGrid()));
    iwave_fdestroy();
    unlink(fname.c_str());
    string dname=fname+"@";
    unlink(dname.c_str());
  }

  TEST_F(GridTest, GridSpace_get_key) {
    // create test file
    fname="test0.rsf";
    create_hfile(fname,g);
    myGridSpace sp(fname,"bulkmod");
    string key="bulkmod";
    EXPECT_EQ(key,sp.get().key);
    iwave_fdestroy();
    unlink(fname.c_str());
    string dname=fname+"@";
    unlink(dname.c_str());
  }

  TEST_F(GridTest, GridSpace_Vector_norm) {
    // create test file
    fname="test0.rsf";
    create_hfile(fname,g);
    myGridSpace sp(fname,"bulkmod");
    Vector<float> x(sp);
    RVLAssignConst<float> ac(1.0);
    x.eval(ac);
    ireal lres = sqrt(11.0*21.0*41.0*0.01*0.01);
    ireal rres = x.norm();
    EXPECT_LT(1.0e-6,fabs(lres-rres));
    iwave_fdestroy();
    unlink(fname.c_str());
    string dname=fname+"@";
    unlink(dname.c_str());
  }

  TEST_F(GridTest, GridSpace_Vector_norm_incore) {
    // create test file
    fname="test0.rsf";
    create_hfile(fname,g);
    myGridSpace sp(fname,"bulkmod",true);
    Vector<float> x(sp);
    RVLAssignConst<float> ac(1.0);
    x.eval(ac);
    ireal lres = sqrt(11.0*21.0*41.0*0.01*0.01);
    ireal rres = x.norm();
    EXPECT_LT(1.0e-6,fabs(lres-rres));
    string dname=fname+"@";
    iwave_fdestroy();
    unlink(fname.c_str());
    unlink(dname.c_str());
  }

  TEST_F(GridTest, GridSpace_Assign_copy) {
    string bgname="test2.rsf";
    create_hfile(bgname,g2);
    myGridSpace rng(bgname,"fungus");
    Vector<float> bg(rng);
    AssignFilename bgaf("test1.rsf");
    bg.eval(bgaf);
    RVLAssignConst<float> bgac(1.0);
    bg.eval(bgac);
    Vector<float> tmp(rng);
    tmp.copy(bg);
    float a = sqrt(get_global_datasize_grid(g2)*get_global_cellvol_grid(g2));
    EXPECT_GT(1.0e-6,fabs(a-tmp.norm()));
    iwave_fdestroy();
    unlink(bgname.c_str());
    string dname=bgname+"@";
    unlink(dname.c_str());
  }

  TEST_F(GridTest, GridSpace_GridWindowOp_ZeroTaper) {

    string bgname="test2.rsf";
    create_hfile(bgname,g2);

    string wname="test3.rsf";
    create_hfile(wname,g3);

    myGridSpace rng(bgname,"fungus");
    myGridSpace dom(wname,"fungus");
    Vector<float> bg(rng);
    AssignFilename bgaf(bgname);
    bg.eval(bgaf);
    RVLAssignConst<float> bgac(1.0);
    bg.eval(bgac);
    Vector<float> w(dom);
    AssignFilename waf(wname);
    w.eval(waf);
    RVLAssignConst<float> wac(2.0);
    w.eval(wac);
    Vector<float> x(rng);
    string xname="wd.rsf";
    AssignFilename xaf(xname);
    x.eval(xaf);
    
    GridWindowOp op(dom,bg);
    OperatorEvaluation<float> opeval(op,w);
    x.copy(opeval.getValue());

    RVLMax<ireal> mx;
    MPISerialFunctionObjectRedn<ireal,ireal> mpimx(mx);
    x.eval(mpimx);
      
    RVLMin<ireal> mn;
    MPISerialFunctionObjectRedn<ireal,ireal> mpimn(mn);
    x.eval(mpimn);    


    EXPECT_FLOAT_EQ(mpimx.getValue(),3.0f);
    EXPECT_FLOAT_EQ(mpimn.getValue(),1.0f);
    iwave_fdestroy();

    string dname;
    dname=bgname+"@";
    unlink(bgname.c_str());
    unlink(dname.c_str());
    dname=wname+"@";
    unlink(wname.c_str());
    unlink(dname.c_str());
    dname=xname+"@";
    unlink(xname.c_str());
    unlink(dname.c_str());

  }

  TEST_F(GridTest, GridSpace_GridWindowOp_2Pt3PtTaper) {
    try {
      string bgname="test1.rsf";
      create_hfile(bgname,g2);
      
      string wname="test2.rsf";
      create_hfile(wname,g3);

      myGridSpace rng(bgname,"fungus");
      myGridSpace dom(wname,"fungus");
      Vector<float> bg(rng);
      AssignFilename bgaf(bgname);
      bg.eval(bgaf);
      RVLAssignConst<float> bgac(1.0);
      bg.eval(bgac);
      Vector<float> w(dom);
      AssignFilename waf(wname);
      w.eval(waf);
      RVLAssignConst<float> wac(2.0);
      w.eval(wac);
      Vector<float> x(rng);
      string xname="wd23.rsf";
      AssignFilename xaf(xname);
      x.eval(xaf);
      RPNT tap;
      for (int i=0;i<RARR_MAX_NDIM;i++) tap[i]=ireal(2*i+2)/10.0;
      GridWindowOp op(dom,bg,tap);
      OperatorEvaluation<float> opeval(op,w);
      x.copy(opeval.getValue());

      RVLMax<ireal> mx;
      MPISerialFunctionObjectRedn<ireal,ireal> mpimx(mx);
      x.eval(mpimx);
      
      RVLMin<ireal> mn;
      MPISerialFunctionObjectRedn<ireal,ireal> mpimn(mn);
      x.eval(mpimn);    

      EXPECT_FLOAT_EQ(mpimx.getValue(),3.0f);
      EXPECT_FLOAT_EQ(mpimn.getValue(),1.0f);

      iwave_fdestroy();

      string dname;
      dname=bgname+"@";
      unlink(bgname.c_str());
      unlink(dname.c_str());
      dname=wname+"@";
      unlink(wname.c_str());
      unlink(dname.c_str());
      dname=xname+"@";
      unlink(xname.c_str());
      unlink(dname.c_str());

    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, GridSpace_GridWindowOp_Adjtest) {
    try {
      string bgname="test1.rsf";
      create_hfile(bgname,g1);
      
      string wname="test2.rsf";
      create_hfile(wname,g2);

      myGridSpace rng(bgname,"fungus");
      myGridSpace dom(wname,"fungus");
      Vector<float> bg(rng);
      AssignFilename bgaf(bgname);
      bg.eval(bgaf);
      RVLAssignConst<float> bgac(1.0);
      bg.eval(bgac);
      Vector<float> w(dom);
      AssignFilename waf(wname);
      w.eval(waf);
      RVLAssignConst<float> wac(2.0);
      w.eval(wac);
      Vector<float> x(rng);
      string xname="wd23.rsf";
      AssignFilename xaf(xname);
      x.eval(xaf);
      RPNT tap;
      for (int i=0;i<RARR_MAX_NDIM;i++) tap[i]=ireal(2*i+2)/10.0;
      GridWindowOp op(dom,bg,tap);
      OperatorEvaluation<float> opeval(op,w);

      RVLRandomize<float> rnd(getpid(),-1.0,1.0);
      std::stringstream res;
      EXPECT_EQ(true,AdjointTest(opeval.getDeriv(),rnd,res));
#ifdef VERBOSE
      cout<<res.str();
#endif
      iwave_fdestroy();

      string dname;
      dname=bgname+"@";
      unlink(bgname.c_str());
      unlink(dname.c_str());
      dname=wname+"@";
      unlink(wname.c_str());
      unlink(dname.c_str());
      dname=xname+"@";
      unlink(xname.c_str());
      unlink(dname.c_str());

    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, GridSpace_GridDerivOp_dir0_Adjtest) {

    try {
      string inname="test1.rsf";
      create_hfile(inname,g1);
      myGridSpace dom(inname,"fungus");
    
      int dir=0;
      GridDerivOp op(dom,dir);
      
#ifdef IWAVE_USE_MPI
      RVLRandomize<float> srnd(getpid(),-1.0,1.0);
      RVL::MPISerialFunctionObject<ireal> rnd(srnd);
#else
      RVLRandomize<float> rnd(getpid(),-1.0,1.0);
#endif
      std::stringstream res;
      EXPECT_EQ(true,AdjointTest(op,rnd,res));

#ifdef VERBOSE
      cout<<res.str();
#endif      
      iwave_fdestroy();
      string dname=inname+"@";
      unlink(inname.c_str());
      unlink(dname.c_str());

    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, GridSpace_GridDerivOp_dir1_Adjtest) {

    try {
      string inname="test1.rsf";
      create_hfile(inname,g);
      myGridSpace dom(inname,"fungus");

      int dir=1;
      GridDerivOp op(dom,dir);

      RVLRandomize<float> rnd(getpid(),-1.0,1.0);
      std::stringstream res;
      EXPECT_EQ(true,AdjointTest(op,rnd,res));
#ifdef VERBOSE
      cout<<res.str();
#endif      
      iwave_fdestroy();
      string dname=inname+"@";
      unlink(inname.c_str());
      unlink(dname.c_str());

    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, GridSpace_GridDerivOp_dir2_Adjtest) {

    try {
      string inname="test1.rsf";
      create_hfile(inname,g1);
      myGridSpace dom(inname,"fungus",true);

      int dir=2;
      GridDerivOp op(dom,dir);

      RVLRandomize<float> rnd(getpid(),-1.0,1.0);
      std::stringstream res;
      EXPECT_EQ(true,AdjointTest(op,rnd,res));
#ifdef VERBOSE
      cout<<res.str();
#endif      
      iwave_fdestroy();
      string dname=inname+"@";
      unlink(inname.c_str());
      unlink(dname.c_str());

    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, GridSpace_GridDerivOp_dir2short_Adjtest) {

    try {
      string inname="testshort.rsf";
      create_hfile(inname,gshort);
      myGridSpace dom(inname,"fungus",true);

      int dir=2;
      GridDerivOp op(dom,dir);

      RVLRandomize<float> rnd(getpid(),-1.0,1.0);
      std::stringstream res;
      EXPECT_EQ(true,AdjointTest(op,rnd,res));
#ifdef VERBOSE
      cout<<res.str();
#endif      
      iwave_fdestroy();
      string dname=inname+"@";
      unlink(inname.c_str());
      unlink(dname.c_str());

    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, GridSpace_GridExtendOp_voltest_ext_d2r3) {

    try {
      string inname="testd2r3.rsf";
      create_hfile(inname,g0);
      myGridSpace dom(inname,"fungus");

      string outname="test1.rsf";
      create_hfile(outname,g1);
      myGridSpace rng(outname,"fungus",true);

      GridExtendOp op(dom,rng);

      RVLAssignConst<float> ac(1.0f);
      Vector<float> dvec(dom);
      Vector<float> rvec(rng);
      dvec.eval(ac);
      rvec.zero();
      op.applyOp(dvec,rvec);
      std::stringstream res;
      EXPECT_EQ(10*20*2.0f,dvec.normsq());
      EXPECT_EQ(10*20*30*6.0f,rvec.normsq());
#ifdef VERBOSE
      cout<<res.str();
#endif      
      iwave_fdestroy();
      string dname=inname+"@";
      unlink(inname.c_str());
      unlink(dname.c_str());

    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, GridSpace_GridExtendOp_adjtest_ext_d2r3) {

    try {
      string inname="testd2r3.rsf";
      create_hfile(inname,g0);
      myGridSpace dom(inname,"fungus");

      string outname="test1.rsf";
      create_hfile(outname,g1);
      myGridSpace rng(outname,"fungus",true);

      GridExtendOp op(dom,rng);

      RVLRandomize<float> rnd(getpid(),-1.0,1.0);
      std::stringstream res;
      EXPECT_EQ(true,AdjointTest(op,rnd,res));
#ifdef VERBOSE
      cout<<res.str();
#endif      

      iwave_fdestroy();
      string dname=inname+"@";
      unlink(inname.c_str());
      unlink(dname.c_str());

    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, GridUnion_one_axis) {
    try {
      axis a;
      a.n=10;
      a.d=2.0;
      a.o=-3.0;
      a.id=3;
      
      grid g;
      init_default_grid(&g);

      EXPECT_EQ(true,grid_union(&g, &a));
#ifdef VERBOSE
      fprint_grid(stderr,g);
#endif
    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }	

  TEST_F(GridTest, GridUnion_three_disjoint_axes) {
    try {
      axis a0;
      a0.n=10;
      a0.d=2.0;
      a0.o=-3.0;
      a0.id=3;
      
      axis a1;
      a1.n=100;
      a1.d=10.0;
      a1.o=500.0;
      a1.id=5;

      axis a2;
      a2.n=20;
      a2.d=1.0;
      a2.o=0.0;
      a2.id=0;
      grid g;
      init_default_grid(&g);

      EXPECT_EQ(true,
		(grid_union(&g, &a0) &&
		 grid_union(&g, &a1) &&
		 grid_union(&g, &a2) ));
#ifdef VERBOSE
      fprint_grid(stderr,g);
#endif
    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }	

  TEST_F(GridTest, GridUnion_two_axes_one_overlap) {
    try {
      axis a0;
      a0.n=10;
      a0.d=2.0;
      a0.o=-3.0;
      a0.id=3;
      
      axis a1;
      a1.n=12;
      a1.d=1.0;
      a1.o=15.0;
      a1.id=5;

      axis a2;
      a2.n=21;
      a2.d=1.0;
      a2.o=0.0;
      a2.id=5;
      grid g;
      init_default_grid(&g);

      EXPECT_EQ(true,(grid_union(&g, &a0) &&
		      grid_union(&g, &a1) &&
		      grid_union(&g, &a2) ));
#ifdef VERBOSE
      fprint_grid(stderr,g);
#endif
    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }	

  TEST_F(GridTest, GridUnion_two_axes_no_overlap) {
    try {
      axis a0;
      a0.n=10;
      a0.d=2.0;
      a0.o=-3.0;
      a0.id=3;
      
      axis a1;
      a1.n=12;
      a1.d=1.0;
      a1.o=15.0;
      a1.id=5;

      axis a2;
      a2.n=7;
      a2.d=1.0;
      a2.o=0.0;
      a2.id=5;

      grid g;
      init_default_grid(&g);

      EXPECT_EQ(true,(grid_union(&g, &a0) &&
		      grid_union(&g, &a1) &&
		      grid_union(&g, &a2) ) );
#ifdef VERBOSE
      fprint_grid(stderr,g);
#endif
    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, Grid_dims_and_vols) {  
    try {
      int sz_grid1 = get_datasize_grid(g1);
      int sz_egrid0 = get_datasize_grid(eg0);
      int esz_grid1 = get_extended_datasize_grid(g1);
      int esz_egrid0 = get_extended_datasize_grid(eg0);
      float v_grid1 = get_cellvol_grid(g1);
      float v_egrid0 = get_cellvol_grid(eg0);
      float ev_grid1 = get_extended_cellvol_grid(g1);
      float ev_egrid0 = get_extended_cellvol_grid(eg0);
      int p_grid1 = get_panelnum_grid(g1);
      int p_egrid0 = get_panelnum_grid(eg0);
      EXPECT_EQ(sz_grid1,200);
      EXPECT_EQ(sz_egrid0,200);
      EXPECT_EQ(esz_grid1,200);
      EXPECT_EQ(esz_egrid0,6000);
      EXPECT_EQ(v_grid1,2.0f);
      EXPECT_EQ(v_egrid0,2.0f);
      EXPECT_EQ(ev_grid1,2.0f);
      EXPECT_EQ(ev_egrid0,6.0f);
      EXPECT_EQ(p_grid1,30);
      EXPECT_EQ(p_egrid0,1);
    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, GridSpace_Vector_extd_norm) {
    try {
      // create test file
      fname="etest0.rsf";
      create_hfile(fname,eg0);
      myGridSpace sp(fname);
      Vector<float> x(sp);
      RVLAssignConst<float> ac(1.0);
      x.eval(ac);
      ireal lres = sqrt(get_extended_cellvol_grid(eg0)*get_extended_datasize_grid(eg0));
      ireal rres = x.norm();
      //      cerr<<"lres="<<lres<<" rres="<<rres<<endl;
      EXPECT_GT(1.0e-6,fabs(lres-rres));
      iwave_fdestroy();
      unlink(fname.c_str());
      string dname=fname+"@";
      unlink(dname.c_str());
    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, gridio_extd) {
    try {
      std::string file_in="etest_in.rsf";
      create_hfile(file_in,eg0);

      myGridSpace sp(file_in);
      AssignFilename afin(file_in);
      Vector<float> x(sp);
      x.eval(afin);
      RVLAssignConst<float> ac(1.0);
      x.eval(ac);
      cerr<<"norm of x = "<<x.norm()<<endl;

      float * a = (float *)usermalloc_(get_extended_datasize_grid(eg0));
      for (int i=0;i<get_extended_datasize_grid(eg0);i++) a[i]=0.0f;
      IPNT ran;  get_n(ran,eg0);
      IPNT rags; IASN(rags,IPNT_0);
      //      IPNT rags; get_gs(rags,eg0);
      float scale=1.0f;
      int panelindex=0;
      //      EXPECT_EQ(0, rsfread(a, rags, ran, file_in.c_str(), scale, stderr, panelindex));
      float ck = 0;
      cerr<<"size="<<get_extended_datasize_grid(eg0)<<endl;
      //      for (int i=0;i<get_extended_datasize_grid(eg0);i++) { 
      //	ck += a[i];
	//	fprintf(stderr,"i=%d a=%e ck=%e\n",i,a[i],ck);
      //      }

      //      fprintf(stderr,"ck=%e\n",ck);
      //      EXPECT_GT(1.e-6,ck-(float)get_extended_datasize_grid(eg0));
      //      fprintf(stderr,"output\n");

      std::string file_out="etest_out.rsf";
      create_hfile(file_out,eg0);
      fprintf(stderr,"before output vector\n");
      Vector<float> y(sp);
      AssignFilename afout(file_out);
      fprintf(stderr,"eval\n");
      y.eval(afout);

      /*


      fprintf(stderr,"rsfwrite\n");
      panelindex=0;
      //      EXPECT_EQ(0, rsfwrite(a, rags, ran, file_out.c_str(), scale, stderr, panelindex));
      fprintf(stderr,"after rsfwrite\n");      
      userfree_(a);

      fprintf(stderr,"exit\n");
      */
      //      y.linComb(-1.0f,x);

      //      EXPECT_GT(1.e-6,y.norm());
    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }


}  // namespace

int main(int argc, char **argv) {
#ifdef IWAVE_USE_MPI
  int ts=0;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif
  try {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef IWAVE_USE_MPI
    storeGlobalComm(MPI_COMM_WORLD);
#endif
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


