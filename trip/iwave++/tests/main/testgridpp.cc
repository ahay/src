#include "parserpp.hh"
#include "gridpp.hh"
#include "gridops.hh"
#include "create_hfile.hh"
#include "gtest/gtest.h"
#include "adjtest.hh"
#include "functions.hh"

namespace {

  using TSOpt::GridDCF;
  using TSOpt::GridDC;
  using TSOpt::GridSpace;
  using TSOpt::GridWindowOp;
  using TSOpt::GridDerivOp;
  using RVL::parse;

  // The fixture for testing class Foo.
  class GridTest : public ::testing::Test {
  protected:
    // You can remove any or all of the following functions if its body
    // is empty.
  
    grid g1;
    grid g2;
    string fname;
  
    GridTest() {
      g1.gdim = 3;
      g1.dim  = 3;

      g1.axes[0].n=10;
      g1.axes[0].d=1.0;
      g1.axes[0].o=0.0;

      g1.axes[1].n=20;
      g1.axes[1].d=2.0;
      g1.axes[1].o=10.0;

      g1.axes[2].n=30;
      g1.axes[2].d=3.0;
      g1.axes[2].o=100.0;

      g2.gdim = 3;
      g2.dim = 2;

      g2.axes[0].n=10;
      g2.axes[0].d=1.0;
      g2.axes[0].o=0.0;

      g2.axes[1].n=20;
      g2.axes[1].d=2.0;
      g2.axes[1].o=10.0;

      g2.axes[2].n=30;
      g2.axes[2].d=3.0;
      g2.axes[2].o=100.0;

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

  // Tests that the Grid::Bar() method does Abc.
  TEST_F(GridTest, get_datasize_grid) {
    EXPECT_EQ(6000, get_datasize_grid(g1));
  }

  TEST_F(GridTest, compare_grid) {
    EXPECT_LT(0,compare_grid(g1,g2));
  }

  TEST_F(GridTest, GridDCF_getGrid) {

    grid g;
    g.axes[0].n=11;
    g.axes[1].n=21;
    g.axes[2].n=41;
    g.axes[0].d=0.1;
    g.axes[1].d=0.1;
    g.axes[2].d=1.0;
    g.axes[0].o=-0.5;
    g.axes[1].o=-1.0;
    g.axes[2].o=0.0;
    g.axes[0].id=0;
    g.axes[1].id=1;
    g.axes[2].id=2;
    g.dim=2;
    g.gdim=3;

    // create test file
    fname="test3.rsf";
    create_hfile(fname,0);
    
    GridDCF fac(fname);

    EXPECT_EQ(0,compare_grid(g,fac.getGrid()));
    iwave_fdestroy();
    unlink(fname.c_str());
    string dname=fname+"@";
    unlink(dname.c_str());
  }

  TEST_F(GridTest, GridDCF_getCellVol) {
    // create test file
    fname="test4.rsf";
    create_hfile(fname,0);
    GridDCF fac(fname);
    EXPECT_GT(1.0e-06,fabs(0.01-fac.getCellVol()));
    iwave_fdestroy();
    unlink(fname.c_str());
    string dname=fname+"@";
    unlink(dname.c_str());
  }

  TEST_F(GridTest, GridDCF_getFilename) {
    // create test file
    fname="test5.rsf";
    create_hfile(fname,0);
    GridDCF fac(fname);
    EXPECT_EQ(fname,fac.getFilename());
    iwave_fdestroy();
    unlink(fname.c_str());
    string dname=fname+"@";
    unlink(dname.c_str());
  }

  TEST_F(GridTest, GridSpace_getGrid) {

    grid g;
    g.axes[0].n=11;
    g.axes[1].n=21;
    g.axes[2].n=41;
    g.axes[0].d=0.1;
    g.axes[1].d=0.1;
    g.axes[2].d=1.0;
    g.axes[0].o=-0.5;
    g.axes[1].o=-1.0;
    g.axes[2].o=0.0;
    g.axes[0].id=0;
    g.axes[1].id=1;
    g.axes[2].id=2;
    g.dim=2;
    g.gdim=3;

    // create test file
    fname="test6.rsf";
    create_hfile(fname,0);

    GridSpace sp(fname);

    EXPECT_EQ(0,compare_grid(g,sp.getGrid()));
    iwave_fdestroy();
    unlink(fname.c_str());
    string dname=fname+"@";
    unlink(dname.c_str());
  }

  TEST_F(GridTest, GridSpace_get_key) {
    // create test file
    fname="test7.rsf";
    create_hfile(fname,0);
    GridSpace sp(fname,"bulkmod",false,cerr);
    string key="bulkmod";
    EXPECT_EQ(key,sp.get().key);
    iwave_fdestroy();
    unlink(fname.c_str());
    string dname=fname+"@";
    unlink(dname.c_str());
  }

  TEST_F(GridTest, GridSpace_Vector_norm) {
    // create test file
    fname="test8.rsf";
    create_hfile(fname,0);
    GridSpace sp(fname,"bulkmod");
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
    fname="test9.rsf";
    create_hfile(fname,0);
    GridSpace sp(fname,"bulkmod",true);
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
    string bgname="test10.rsf";
    create_hfile(bgname,0);
    GridSpace rng(bgname,"fungus");
    Vector<float> bg(rng);
    AssignFilename bgaf("test10.rsf");
    bg.eval(bgaf);
    RVLAssignConst<float> bgac(1.0);
    bg.eval(bgac);
    Vector<float> tmp(rng);
    tmp.copy(bg);
    EXPECT_GT(1.0e-6,fabs(sqrt(1.21)-tmp.norm()));
    iwave_fdestroy();
    unlink(bgname.c_str());
  }

  TEST_F(GridTest, GridSpace_GridWindowOp_ZeroTaper) {

    string bgname="test10.rsf";
    create_hfile(bgname,0);

    string wname="test11.rsf";
    create_hfile(wname,0);

    GridSpace rng(bgname,"fungus");
    GridSpace dom(wname,"fungus");
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
    unlink(bgname.c_str());
    unlink(wname.c_str());
    unlink(xname.c_str());
  }

  TEST_F(GridTest, GridSpace_GridWindowOp_2Pt3PtTaper) {
    try {
      string bgname="test12.rsf";
      create_hfile(bgname,0);
      
    string wname="test13.rsf";
    create_hfile(wname,0);

    GridSpace rng(bgname,"fungus");
    GridSpace dom(wname,"fungus");
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
    unlink(bgname.c_str());
    unlink(wname.c_str());
    unlink(xname.c_str());
    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, GridSpace_GridWindowOp_Adjtest) {
    try {
      string bgname="test12.rsf";
      create_hfile(bgname,0);
      
      string wname="test13.rsf";
      create_hfile(wname,0);

      GridSpace rng(bgname,"fungus");
      GridSpace dom(wname,"fungus");
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
      EXPECT_EQ(true,AdjointTest(opeval.getDeriv(),rnd,cout));
      
      iwave_fdestroy();
      unlink(bgname.c_str());
      unlink(wname.c_str());
      unlink(xname.c_str());
    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, GridSpace_GridDerivOp_dir0_Adjtest) {

    try {
      string inname="test15.rsf";
      create_hfile(inname,0);
      GridSpace dom(inname,"fungus");

      int dir=0;
      GridDerivOp op(dom,dir);

      RVLRandomize<float> rnd(getpid(),-1.0,1.0);
      EXPECT_EQ(true,AdjointTest(op,rnd,cout));
      
      iwave_fdestroy();
      unlink(inname.c_str());
    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, GridSpace_GridDerivOp_dir1_Adjtest) {

    try {
      string inname="test15.rsf";
      create_hfile(inname,0);
      GridSpace dom(inname,"fungus");

      int dir=1;
      GridDerivOp op(dom,dir);

      RVLRandomize<float> rnd(getpid(),-1.0,1.0);
      EXPECT_EQ(true,AdjointTest(op,rnd,cout));
      
      iwave_fdestroy();
      unlink(inname.c_str());
    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }

  TEST_F(GridTest, GridSpace_GridDerivOp_dir2_Adjtest) {

    try {
      string inname="test15.rsf";
      create_hfile(inname,0);
      GridSpace dom(inname,"fungus",true);

      int dir=2;
      GridDerivOp op(dom,dir);

      RVLRandomize<float> rnd(getpid(),-1.0,1.0);
      EXPECT_EQ(true,AdjointTest(op,rnd,cout));
      
      iwave_fdestroy();
      unlink(inname.c_str());
    }
    catch (RVLException & e) {
      e.write(cerr);
    }
  }
}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
