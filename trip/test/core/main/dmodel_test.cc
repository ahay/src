#include "gtest/gtest.h"
#include "acd_defn.hh"

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

  using TSOpt::IOTask;
  using TSOpt::TASK_RELN;
  using TSOpt::IOTaskWriter;
  using RVL::RVLException;

  class ACD_ModelDefn_Test: public ::testing::Test {
  public:
    IWaveInfo ic;
    ofstream details;
    ACD_ModelDefn_Test(): ic(), details("details.txt",ios_base::app) {}
    ~ACD_ModelDefn_Test() { details.close(); }

    virtual void SetUp() {
      // Code here will be called immediately after the constructor (right
      // before each test).
    }

    virtual void TearDown() {
      // Code here will be called immediately after each test (right
      // before the destructor).
    }

  };

  TEST_F(ACD_ModelDefn_Test, model_defn_test) {
    std::ostringstream s;
    string t = "IO Definition: name = acd\n  keyword[0]=csq index=0 input=1 active=1\n  keyword[1]=data index=1 input=0 active=1\n  keyword[2]=source index=1 input=1 active=0\n  keyword[3]=movie index=1 input=0 active=0\n  keyword[4]=init index=1 input=1 active=0\n";
    ic.write_iwave_iokeys(s);
#ifdef GTEST_VERBOSE
    details<<"======================================================\n";
    details<<"ACD_ModelDefn_Test, model_defn_test\n";
    ic.write_iwave_iokeys(details);
    ic.write_iwave_fields(details);
#endif
    EXPECT_EQ(t,s.str());
  }

  TEST_F(ACD_ModelDefn_Test, model_task_order0_fwd) {
    std::ostringstream s;
    std::vector<TASK_RELN *> tr;
    IOTask(tr,0,true,ic);
    IOTaskWriter(tr,s);
    for (int i=0;i<tr.size();i++) delete tr[i];
    std::string t = "index=0 keyword=csq rarrindex=0 input=1\nindex=0 keyword=data rarrindex=1 input=0\nindex=0 keyword=source rarrindex=1 input=1\nindex=0 keyword=movie rarrindex=1 input=0\nindex=0 keyword=init rarrindex=1 input=1\n";
    EXPECT_EQ(t,s.str());
#ifdef GTEST_VERBOSE
    details<<"======================================================\n";
    details<<"ACD_ModelDefn_Test, model_task_order0_fwd\n";
    details<<s.str();
#endif
  }

  TEST_F(ACD_ModelDefn_Test, model_task_order1_fwd) {
    std::vector<TASK_RELN *> tr;
    IOTask(tr,1,true,ic);
    std::ostringstream s;
    IOTaskWriter(tr,s);
    for (int i=0;i<tr.size();i++) delete tr[i];
    std::string t = "index=0 keyword=csq rarrindex=0 input=1\nindex=0 keyword=source rarrindex=1 input=1\nindex=0 keyword=movie rarrindex=1 input=0\nindex=0 keyword=init rarrindex=1 input=1\nindex=1 keyword=csq_d1 rarrindex=0 input=1\nindex=1 keyword=data rarrindex=1 input=0\n";
    EXPECT_EQ(t,s.str());
#ifdef GTEST_VERBOSE
    details<<"======================================================\n";
    details<<"ACD_ModelDefn_Test, model_task_order1_fwd\n";
    details<<s.str();
#endif
  }

  TEST_F(ACD_ModelDefn_Test, model_task_order2_fwd) {
    std::vector<TASK_RELN *> tr;
    IOTask(tr,2,true,ic);
    std::ostringstream s;
    IOTaskWriter(tr,s);
    for (int i=0;i<tr.size();i++) delete tr[i];
    std::string t = "index=0 keyword=csq rarrindex=0 input=1\nindex=0 keyword=source rarrindex=1 input=1\nindex=0 keyword=movie rarrindex=1 input=0\nindex=0 keyword=init rarrindex=1 input=1\nindex=1 keyword=csq_d1 rarrindex=0 input=1\nindex=2 keyword=csq_d2 rarrindex=0 input=1\nindex=3 keyword=data rarrindex=1 input=0\n";
    EXPECT_EQ(t,s.str());
#ifdef GTEST_VERBOSE
    details<<"======================================================\n";
    details<<"ACD_ModelDefn_Test, model_task_order2_fwd\n";
    details<<s.str();
#endif
  }

  TEST_F(ACD_ModelDefn_Test, model_task_order3_fwd) {
    std::vector<TASK_RELN *> tr;
    IOTask(tr,3,true,ic);
    std::ostringstream s;
    IOTaskWriter(tr,s);
    for (int i=0;i<tr.size();i++) delete tr[i];
    std::string t = "index=0 keyword=csq rarrindex=0 input=1\nindex=0 keyword=source rarrindex=1 input=1\nindex=0 keyword=movie rarrindex=1 input=0\nindex=0 keyword=init rarrindex=1 input=1\nindex=1 keyword=csq_d1 rarrindex=0 input=1\nindex=2 keyword=csq_d2 rarrindex=0 input=1\nindex=4 keyword=csq_d3 rarrindex=0 input=1\nindex=7 keyword=data rarrindex=1 input=0\n";
    EXPECT_EQ(t,s.str());
#ifdef GTEST_VERBOSE
    details<<"======================================================\n";
    details<<"ACD_ModelDefn_Test, model_task_order3_fwd\n";
    details<<s.str();
#endif
  }

  TEST_F(ACD_ModelDefn_Test, model_task_order0_adj) {
    try {
      std::vector<TASK_RELN *> tr;
      IOTask(tr,0,false,ic);
      std::ostringstream s;
      IOTaskWriter(tr,s);
      std::string t = "";
      for (int i=0;i<tr.size();i++) delete tr[i];
      EXPECT_EQ(t,s.str());
    }
    catch (RVLException & e) {
      std::ostringstream s;
      e.write(s);
      std::string t = "Error: IOTask constructor\n  order 0 invalid for adjoint branch\n\n";
      EXPECT_EQ(t,s.str());
#ifdef GTEST_VERBOSE
      details<<"======================================================\n";
      details<<"ACD_ModelDefn_Test, model_task_order0_adj\n";
      details<<s.str();
#endif
    }
  }

  TEST_F(ACD_ModelDefn_Test, model_task_order1_adj) {
    std::vector<TASK_RELN *> tr;
    IOTask(tr,1,false,ic);
    std::ostringstream s;
    IOTaskWriter(tr,s);
    for (int i=0;i<tr.size();i++) delete tr[i];
    std::string t = "index=0 keyword=csq rarrindex=0 input=1\nindex=0 keyword=source rarrindex=1 input=1\nindex=0 keyword=movie rarrindex=1 input=0\nindex=0 keyword=init rarrindex=1 input=1\nindex=1 keyword=csq_b1 rarrindex=0 input=0\nindex=1 keyword=data rarrindex=1 input=1\n";
    EXPECT_EQ(t,s.str());
#ifdef GTEST_VERBOSE
    details<<"======================================================\n";
    details<<"ACD_ModelDefn_Test, model_task_order1_adj\n";
    details<<s.str();
#endif

  }

  TEST_F(ACD_ModelDefn_Test, model_task_order2_adj) {
    std::vector<TASK_RELN *> tr;
    IOTask(tr,2,false,ic);
    std::ostringstream s;
    IOTaskWriter(tr,s);
    for (int i=0;i<tr.size();i++) delete tr[i];
    std::string t = "index=0 keyword=csq rarrindex=0 input=1\nindex=0 keyword=source rarrindex=1 input=1\nindex=0 keyword=movie rarrindex=1 input=0\nindex=0 keyword=init rarrindex=1 input=1\nindex=1 keyword=csq_d1 rarrindex=0 input=1\nindex=2 keyword=csq_b2 rarrindex=0 input=0\nindex=3 keyword=data rarrindex=1 input=1\n";
    EXPECT_EQ(t,s.str());
#ifdef GTEST_VERBOSE
    details<<"======================================================\n";
    details<<"ACD_ModelDefn_Test, model_task_order2_adj\n";
    details<<s.str();
#endif
  }

  TEST_F(ACD_ModelDefn_Test, model_task_order3_adj) {
    std::vector<TASK_RELN *> tr;
    IOTask(tr,3,false,ic);
    std::ostringstream s;
    IOTaskWriter(tr,s);
    for (int i=0;i<tr.size();i++) delete tr[i];
    std::string t = "index=0 keyword=csq rarrindex=0 input=1\nindex=0 keyword=source rarrindex=1 input=1\nindex=0 keyword=movie rarrindex=1 input=0\nindex=0 keyword=init rarrindex=1 input=1\nindex=1 keyword=csq_d1 rarrindex=0 input=1\nindex=2 keyword=csq_d2 rarrindex=0 input=1\nindex=4 keyword=csq_b3 rarrindex=0 input=0\nindex=7 keyword=data rarrindex=1 input=1\n";
    EXPECT_EQ(t,s.str());
#ifdef GTEST_VERBOSE
    details<<"======================================================\n";
    details<<"ACD_ModelDefn_Test, model_task_order3_adj\n";
    details<<s.str();
#endif
  }
}

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


