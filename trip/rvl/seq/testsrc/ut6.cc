#include "seqspace.hh"
#include "multop.hh"
#include "adjtest.hh"
#include "terminator.hh"


using namespace RVL;
using namespace RVLAlg;

int main() {

  // freeze random seed for regr tests
  //  PlantSeeds(getpid());
  PlantSeeds(19490615);

  // length of trial vectors
  int n = 25;
  // length of factor
  int m = 10;

  // coeff list for op
  list<double> flist;
  for (int i=0;i<m;i++)
    flist.push_back(-0.5+ rand()/(RAND_MAX+1.0));

  // op of multiplication by 1+x
  PolyMultOp op(flist);

  // randomizer
  RandList rnd(n);
 
  // well did you screw up or not
  ofstream str("testsrc/ut6.aux");
  cout<<"AdjointTest return = "<<AdjointTest(op,rnd,str)<<"; details in ut6.rpt"<<endl;
  str.close();
}
