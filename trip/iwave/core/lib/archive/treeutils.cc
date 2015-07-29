#include "treeutils.hh"

//#define VERBOSE
// "little endian" bool deques: most significant bit is first

namespace TSOpt {

  deque<bool> d2b(int n) {
    deque<bool> xcode;
    xcode.clear();
    if (n <= 0) {
      xcode.push_front(false);
      return xcode;
    }
    while (n>0) {
      if (n % 2) xcode.push_front(true);
      else xcode.push_front(false);
      n = n/2;
    }
    return xcode;
  }

  int level(deque<bool> xcode) {
    int lvl = 0;
    for (int i=0; i<xcode.size(); i++) 
      if (xcode[i]) lvl++;
    return lvl;
  }

  int b2d(deque<bool> xcode) {
    int val = 0; 
    if (xcode[0]) val=1;
    for (int i=1; i<xcode.size(); i++)  {
      val *= 2;
      if (xcode[i]) val++;
    }
    return val;
  }

  // index arrays are also "little-endian": highest value 
  // index is last
  bool index_incr(deque<bool> xcode, vector<int> & idx) {
#ifdef VERBOSE
    cerr<<"index_incr: before pop: xcode = ";
    print_xcode(xcode,cerr);
    cerr<<"  idx = ";
    print_idx(idx,cerr);
#endif
    while (xcode.size()>1 && !xcode[0]) {
      xcode.pop_front();
#ifdef VERBOSE
      cerr<<"pop: xcode = ";
      print_xcode(xcode,cerr);
#endif
    }
    if (xcode.size()>0) {
#ifdef VERBOSE
      cerr<<"index_incr: xcode after stripping leading 0's: ";
      print_xcode(xcode,cerr);
#endif
      int val = (int)pow(2.0f,(int)(xcode.size()-1));
#ifdef VERBOSE
      cerr<<"index_incr: value for this level = "<<val<<endl;
#endif
      if (xcode.size()>1) {
	xcode.pop_front();
#ifdef VERBOSE
	cerr<<"after pop: ";
	print_xcode(xcode,cerr);
#endif
	vector<int> tmp;
	tmp.clear();
#ifdef VERBOSE
	cerr<<"index_incr: recursive call\n";
#endif
	if (index_incr(xcode,tmp)) {
	  // add truncated index array to idx, both incremented and unincremented
#ifdef VERBOSE
	  cerr<<"idx returnd from recursive call: ";
	  print_idx(tmp,cerr);
	  cerr<<"  for xcode = ";
	  print_xcode(xcode,cerr);
#endif
	  for (int i=0;i<tmp.size();i++) {
#ifdef VERBOSE
	    cerr<<"  pushing back "<<val + tmp[i]<<endl;
#endif
	    idx.push_back(val + tmp[i]);
	  }
#ifdef VERBOSE
	  cerr<<"after first push-back: ";
	  print_idx(idx,cerr);
#endif
	  for (int i=0;i<tmp.size();i++) {
#ifdef VERBOSE
	    cerr<<"  pushing back "<< tmp[i] <<endl;
#endif
	    idx.push_back(tmp[i]);
	  }
#ifdef VERBOSE
	  cerr<<"after second push-back: ";
	  print_idx(idx,cerr);
#endif
	  return true;
	}
      }
      else {
	if (xcode[0]) idx.push_back(1);
	idx.push_back(0);
#ifdef VERBOSE
	cerr<<"idx after level 1: ";
	print_idx(idx,cerr);
#endif
	return true;
      }
    }
    return false;
  }

  void print_idx(vector<int> idx, ostream & str) {
    for (int i=0;i<idx.size();i++)
      str<<" "<<idx[i];
    str<<endl;
  }

  void print_xcode(deque<bool> xcode, ostream & str) {
    for (int i=0;i<xcode.size();i++) {
      if (xcode[i]) str<<" 1";
      else str<<" 0";
    }
    str<<endl;
  }

}
