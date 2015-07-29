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
}
