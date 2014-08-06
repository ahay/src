#ifndef __SEQ_SPACE
#define __SEQ_SPACE

#include "space.hh"
#include "seqfos.hh"

class SeqSpace: public Space<double> {

public:

  SeqSpace() {}
  SeqSpace(const SeqSpace & sp) {}
  ~SeqSpace(){}
  
  DataContainer * buildDataContainer() const { return new SeqDC; }

  bool operator ==(const Space<double> & sp) const;
  bool isCompatible(DataContainer const & dc) const;

  double inner(DataContainer const & x, 
	       DataContainer const & y) const;
  void zero(DataContainer & x) const;
  void linComb(double a, DataContainer const & x,
	       double b, DataContainer & y) const;

  ostream & write(ostream & str) const;
  
};


#endif
