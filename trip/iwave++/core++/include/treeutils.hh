#include "except.hh"
#include <deque>
#include "sim.hh"
#include "iwavetime.hh"

namespace TSOpt {

  using RVL::RVLException;

  size_t pow2(int);
  deque<bool> d2b(int);
  int b2d(deque<bool>);
  int level(deque<bool>);
  bool index_incr(deque<bool>, vector<int> &);
  void print_xcode(deque<bool> xcode, ostream & str);
  void print_idx(vector<int> idx, ostream & str);

  template<typename T> 
  void TreeStep(std::vector<T *> & u, void (*F)(std::vector<T *>), bool verb=false) {
    // first update 2nd half given first half
    F(u);

    for (size_t i=u.size(); i>0; i--) {
      cerr<<"  component "<<i-1<<endl;
      // extract index vector
      deque<bool> xcode = d2b(i-1);
      vector<int> idx;
      idx.clear();
      index_incr(xcode,idx);
      // build an argument vector 
      vector<T *> w(idx.size());
      for (size_t j=0; j<idx.size(); j++) 
	w.at(j) = u.at(idx.at(j));
      // call timestep function interface
      F(w);
      if (verb) {
	cerr<<"----------------------------------------------------------"<<endl;
	cerr<<"TreeStep: updated field "<<idx[0];
	if (idx.size()>1) {
	  cerr<<" using fields ";
	  for (size_t j=1;j<idx.size()-1;j++) cerr<<idx[j]<<", ";
	  cerr<<idx[idx.size()-1];
	}
        cerr<<endl;
	cerr<<"  field states:\n";
	for (size_t j=0;j<idx.size();j++) {
	  cerr<<"  *** "<<idx[j]<<":"<<endl;
	  w[j]->write(cerr);
        }
	cerr<<"----------------------------------------------------------"<<endl;
      }
    }
    cerr<<"exit TreeStep\n";
  }

  template<typename T>
  class TreeState {
  public:
    virtual void setTime(Time const & t) = 0;
    virtual Time & getTime() = 0;
    virtual Time const & getTime() const = 0;
    virtual std::vector<T *> & getStateArray() = 0;
    virtual std::vector<T *> const & getStateArray() const = 0;
  };

  //template<typename R, typename T> 
  template<typename T>
  class TreeTimeStep: public TimeStep<TreeState<T> > {

    typedef void (*TSF)(std::vector<T *>);

  private:
    // R & mystate;
    TreeState<T> & mystate;
    TSF f;
    TreeTimeStep();
    //    TreeTimeStep(TreeTimeStep<R,T> const &);
    TreeTimeStep(TreeTimeStep<T> const &);
    bool verb;

  public:
    TreeTimeStep(TreeState<T> & _mystate, TSF _f, bool _verb=false): mystate(_mystate), f(_f), verb(_verb) {cerr<<"TreeTimeStep constructor\n";}
    TreeState<T> & getState() { return mystate; }
    TreeState<T> const & getState() const { return mystate; }
    void run()  {
      cerr<<"call TreeStep\n";
      TreeStep(this->getState().getStateArray(),f,verb);
      cerr<<"return TreeStep\n";
      return;
    }

    Time & getNextTime() const { 
      RVLException e;
      e<<"ERROR - DONT CALL THIS\n";
      throw e;
    }
    ostream & write(ostream & str) const {
      return str;
    }
  };

  // T = state type, D = instance data for T 
  // n = order of derivative tree
  template<typename T, typename D>
  class IWaveTreeState: public TreeState<T> {
  private:
    bool fwd;
    std::vector<T *> sa;
    size_t get_time_index() const {
      if (fwd) return 0;
      return sa.size()-1;
    }
    IWaveTreeState();
    IWaveTreeState(IWaveTreeState const &);
  public:
    IWaveTreeState(D const & sdata, int n=0, bool _fwd=true): fwd(_fwd), sa(pow2(n)) {
      cerr<<"IWaveTreeState Constructor:";
      for (size_t i=0; i<sa.size(); i++) {
	sa[i]=new T(sdata);
      }
    }
    ~IWaveTreeState() {
      for (size_t i=0; i<sa.size(); i++) delete sa[i];
    }

    void setTime(Time const & t) {
      try { sa[get_time_index()]->setTime(t);
      }
      catch (RVLException & e) {
	e<<"\ncalled from IWaveTreeState::setTime\n";
	throw e;
      }
    }
    Time const & getTime() const { return sa[get_time_index()]->getTime(); }
    Time & getTime() { return sa[get_time_index()]->getTime(); }
    Time & getNextTime() const { sa[get_time_index()]->getTime(); }

    virtual std::vector<T *> & getStateArray() { return sa; }
    virtual std::vector<T *> const & getStateArray() const { return sa; }

    ostream & write(ostream & str) const {
      str<<"IWaveTreeState, length "<<sa.size()<<"\n";
      for (size_t i=0; i<sa.size(); i++) {
	str<<"** component "<<i<<":\n";
	sa[i].write(str);
      }
      return str;
    }
  };

    


  /////////////////////// Test Structures ///////////////////////

  typedef struct s_TestState0Data {
    int nlvl;
  } TestState0Data;

  class TestState0 {
  private:
    int nlvl;
    mutable IWaveTime ts;
    mutable TIMESTEPINDEX lts;
  public:
    TestState0(TestState0Data const & d): nlvl(d.nlvl), ts(lts) {lts.it=0; lts.iv=0; lts.dt=1.0; lts.niv=nlvl;}
    TestState0(TestState0 const & s): nlvl(s.nlvl), ts(s.ts) {}
    void setTime(Time const & t) {
      try {
	ts = t;
      }
      catch (RVLException & e) {
	e<<"\ncalled from TestState0::setTime\n";
	throw e;
      }
    }
    Time & getTime() { return ts; }
    Time const & getTime() const { return ts; }

    void next_step() { 
      TIMESTEPINDEX & dts = ts.getCstruct();
      dts.iv++;
      if (dts.iv >= nlvl) {
	dts.iv=0;
	dts.it++;
      }
    }
    ostream & write(ostream & str) const {
      str<<"TestState0: IWaveTime = \n";
      ts.write(str);
      return str;
    }
  };

  typedef IWaveTreeState<TestState0, TestState0Data> TestTreeState0;

  void TestState0Fun(std::vector<TestState0 *> w);

  
}
