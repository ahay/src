#include "gtest/gtest.h"
#include "treeutils.hh"

namespace {

  using TSOpt::level;
  using TSOpt::b2d;
  using TSOpt::d2b;
  using TSOpt::index_incr;
  using TSOpt::print_xcode;
  using TSOpt::print_idx;
  using TSOpt::TreeTimeStep;
  using TSOpt::TestState0;
  using TSOpt::TestState0Fun;
  using TSOpt::TestTreeState0;
  using TSOpt::TestState0Data;

  class TreeTest: public ::testing::Test {
  public:

    TreeTest() {}
    ~TreeTest() {}

    virtual void SetUp() {
      // Code here will be called immediately after the constructor (right
      // before each test).
    }

    virtual void TearDown() {
      // Code here will be called immediately after each test (right
      // before the destructor).
    }

  };

  TEST_F(TreeTest, xcode_level_test) {
    deque<bool> xcode;
    xcode.push_back(1);
    cerr<<"XCodeLevelTest: xcode = ";
    print_xcode(xcode,cerr);
    cerr<<"XCodeLevelTest: lvl = "<<level(xcode)<<endl;
    cerr<<"XCodeLevelTest: val = "<<b2d(xcode)<<endl;
    xcode.push_back(0);
    cerr<<"XCodeLevelTest: xcode = ";
    print_xcode(xcode,cerr);
    cerr<<"XCodeLevelTest: lvl = "<<level(xcode)<<endl;
    cerr<<"XCodeLevelTest: val = "<<b2d(xcode)<<endl;
    xcode.push_back(1);
    cerr<<"XCodeLevelTest: xcode = ";
    print_xcode(xcode,cerr);
    cerr<<"XCodeLevelTest: lvl = "<<level(xcode)<<endl;
    cerr<<"XCodeLevelTest: val = "<<b2d(xcode)<<endl;
    xcode.push_back(1);
    cerr<<"XCodeLevelTest: xcode = ";
    print_xcode(xcode,cerr);
    cerr<<"XCodeLevelTest: lvl = "<<level(xcode)<<endl;
    cerr<<"XCodeLevelTest: val = "<<b2d(xcode)<<endl;
  }
  
  TEST_F(TreeTest, xcode_index_test) {
    vector<int> idx;
    deque<bool> xcode;

    for (int i=15; i>=0; i--) {
      cerr<<"-------------------------------------------------------------"<<endl;
      xcode = d2b(i);
      cerr<<"XCodeLevelTest: xcode = ";
      print_xcode(xcode,cerr);
      cerr<<"XCodeLevelTest: xcode value = "<<b2d(xcode)<<endl;
      idx.clear();
      index_incr(xcode,idx);
      cerr<<"XCodeLevelTest: index = ";
      print_idx(idx,cerr);
    }
  }

  TEST_F(TreeTest, treetimestep_teststate0_lvl1ord0) { 
    // level = 1; 
    TestState0Data d;
    d.nlvl = 1;
    // order = 0
    int order = 0;
    TestTreeState0 state(d,order);
    TreeTimeStep<TestState0> step(state,TestState0Fun,true);
    step.run();
    step.run();
  }

  TEST_F(TreeTest, treetimestep_teststate0_lvl1ord1) {  
    // level = 1
    TestState0Data d;
    d.nlvl = 1;
    // order = 1
    int order = 1;
    TestTreeState0 state(d,order);
    TreeTimeStep<TestState0> step(state,TestState0Fun,true);
    step.run();
    step.run();
  }

  TEST_F(TreeTest, treetimestep_teststate0_lvl1ord2) {  
    // level = 1
    TestState0Data d;
    d.nlvl = 1;
    // order = 2
    int order = 2;
    TestTreeState0 state(d,order);
    TreeTimeStep<TestState0> step(state,TestState0Fun,true);
    step.run();
    step.run();
  }

  TEST_F(TreeTest, treetimestep_teststate0_lvl2ord0) { 
    // level = 1; 
    TestState0Data d;
    d.nlvl = 2;
    // order = 0
    int order = 0;
    TestTreeState0 state(d,order);
    TreeTimeStep<TestState0> step(state,TestState0Fun,true);
    step.run();
    step.run();
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
