#include "iwavetime.hh"

namespace TSOpt {

  IWaveTime & IWaveTime::operator=(IWaveTime const & t) {
    ts.it=t.ts.it;
    ts.iv=t.ts.iv;
    ts.dt=t.ts.dt;
    return *this;
  }

  Time & IWaveTime::operator=(Time const &t) {
    try {
      IWaveTime const & ts1 = dynamic_cast< IWaveTime const &>(t);
      return *this=ts1;
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: IWaveTime::operator=\n";
      e<<"input Time object not IWaveTime\n";
      throw e;
    }
  }

  IWaveTime & IWaveTime::operator=(int it){
    ts.it = it;
    ts.iv = 0;
    return *this;
  }

  bool IWaveTime::operator<(Time const & t) const {
    try {
      IWaveTime const & ts1 = dynamic_cast< IWaveTime const &>(t);
      return less_than(this->ts,ts1.ts);
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: IWaveTime::operator<\n";
      e<<"input Time object not IWaveTime\n";
      throw e;
    }
  }	

  bool IWaveTime::operator>(Time const & t) const {
    try {
      IWaveTime const & ts1 = dynamic_cast< IWaveTime const &>(t);
      return less_than(ts1.ts,this->ts);
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: IWaveTime::operator>\n";
      e<<"input Time object not IWaveTime\n";
      throw e;
    }
  }	  

  IWaveTime & IWaveTime::operator++() { 
    this->ts.it++;
    this->ts.iv = 0;
    return *this;
  }

  IWaveTime & IWaveTime::operator--() { 
    this->ts.it--; 
    this->ts.iv = 0;
    return *this;
  }
   
  ostream & IWaveTime::write(ostream & str) const {
    str<<"IWaveTime: it="<<(this->ts).it<<" iv="<<(this->ts).iv<<" dt="<<(this->ts).dt<<" t="<<get_time(this->ts)<<"\n";
    return str;
  }
}
