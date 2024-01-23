#include "seqspace.hh"

bool SeqSpace::operator==(const Space<double> & sp) const {
  // only the type need be tested, as there is - up to trivial
  // isomorphism - only one finite sequence space
  SeqSpace const * t = NULL;
  if ((t = dynamic_cast<SeqSpace const *>(&sp))) return true;
  return false;
}

bool SeqSpace::isCompatible(DataContainer const & dc) const {
  // similarly any SeqDC is compatible with any SeqSpace
  SeqDC const * t = NULL;
  if ((t = dynamic_cast<SeqDC const *>(&dc))) return true;
  return false;
}

double SeqSpace::inner(DataContainer const & x,
		       DataContainer const & y) const {
  try {
    SeqDC const & sx = dynamic_cast<SeqDC const &>(x);
    SeqDC const & sy = dynamic_cast<SeqDC const &>(y);
    InnerProdList ip;
    ip(sx,sy);
    return ip.getValue();
  }
  catch (bad_cast const&) {
    RVLException e;
    e<<"Error: SeqSpace::inner\n";
    e<<"one or the other DC inputs not a SeqDC\n";
    throw e;
  }
}

void SeqSpace::zero(DataContainer & x) const {
  try {
    SeqDC & sx = dynamic_cast<SeqDC &>(x);
    ZeroList zip;
    zip(sx);
  }
  catch (bad_cast const&) {
    RVLException e;
    e<<"Error: SeqSpace::zero\n";
    e<<"DC input not a SeqDC\n";
    throw e;
  }
}

void SeqSpace::linComb(double a, DataContainer const & x,
		       double b, DataContainer & y) const {
  try {
    SeqDC const & sx = dynamic_cast<SeqDC const &>(x);
    SeqDC & sy = dynamic_cast<SeqDC &>(y);
    LinCombList lc(a,b);
    lc(sy,sx);
  }
  catch (bad_cast const&) {
    RVLException e;
    e<<"Error: SeqSpace::linComb\n";
    e<<"one or the other DC inputs not a SeqDC\n";
    throw e;
  }
}

ostream & SeqSpace::write(ostream & str) const {
  str<<"SeqSpace: preHilbert space of finite sequences\n";
  return str;
}
