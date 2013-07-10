#include "tsind.hh"

namespace TSOpt {

  TSIndex & TSIndex::operator=(TSIndex const & t) {
    ts.it=t.ts.it;
    ts.iv=t.ts.iv;
    ts.dt=t.ts.dt;
    synclts();
    return *this;
  }

  Time & TSIndex::operator=(Time const &t) {
    try {
      TSIndex const & ts1 = dynamic_cast< TSIndex const &>(t);
      return *this=ts1;
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: TSIndex::operator=\n";
      e<<"input Time object not TSIndex\n";
      throw e;
    }
  }

  TSIndex & TSIndex::operator=(int it){
    ts.it = it;
    ts.iv = 0;
    synclts();
    return *this;
  }

  bool TSIndex::operator<(Time const & t) const {
    try {
      TSIndex const & ts1 = dynamic_cast< TSIndex const &>(t);
      return less_than(this->ts,ts1.ts);
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: TSIndex::operator<\n";
      e<<"input Time object not TSIndex\n";
      throw e;
    }
  }	

  bool TSIndex::operator>(Time const & t) const {
    try {
      TSIndex const & ts1 = dynamic_cast< TSIndex const &>(t);
      return less_than(ts1.ts,this->ts);
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: TSIndex::operator>\n";
      e<<"input Time object not TSIndex\n";
      throw e;
    }
  }	  

  TSIndex & TSIndex::operator++() { 
    this->ts.it++;
    this->ts.iv = 0;
    synclts();
    return *this;
  }

  TSIndex & TSIndex::operator--() { 
    this->ts.it--; 
    this->ts.iv = 0;
    synclts();
    return *this;
  }
   
  ostream & TSIndex::write(ostream & str) const {
    str<<"TSIndex: it="<<ts.it<<" iv="<<ts.iv<<" t="<<get_time(this->ts)<<"\n";
    return str;
  }
}
