#include "iwinfo.hh"

size_t pow2(int n) {
  size_t m = 1;
  while (n>0) {
    m *=2;
    n--;
  }
  return m;
}
    
int IWaveInfo::get_num_fields() const {
  int num=0;
  while ((get_iwave_fields()[num].field != "") && num<IWAVEMAXDATA) num++;
  if (num >= IWAVEMAXDATA) {
    RVLException e;
    e<<"Error: get_num_fields\n";
    e<<"  over limit for number - probably left off last entry with field=\"\"\n";
    throw e;
  }
  return num;
}
  
int IWaveInfo::get_num_iokeys() const {
  int num=0;
  while ((get_iwave_iokeys()[num].keyword != "") && num<IWAVEMAXDATA) num++;
  if (num >= IWAVEMAXDATA) {
    RVLException e;
    e<<"Error: get_num_iokeys\n";
    e<<"  over limit for number - probably left off last entry with keyword=\"\"\n";
    throw e;
  }
  return num;
}
  
ostream & IWaveInfo::write_iwave_fields(ostream & str) const {
  str <<"Field Definition: name = "<<get_iwave_model()<<"\n";
  for (int i=0;i<get_num_fields();i++) {
    str<<"  field["<<i<<"]="<<get_iwave_fields()[i].field
       <<" dynamic="<<get_iwave_fields()[i].dynamic
       <<" substep="<<get_iwave_fields()[i].substep
       <<" gtype=[";
    for (int j=0; j<RARR_MAX_NDIM-1; j++) 
      str<<get_iwave_fields()[i].gtype[j]<<",";
    str<<get_iwave_fields()[i].gtype[RARR_MAX_NDIM-1]<<"]\n";
  }
  return str;
}
  
ostream & IWaveInfo::write_iwave_iokeys(ostream & str) const {
  str <<"IO Definition: name = "<<get_iwave_model()<<"\n";
  for (int i=0;i<get_num_iokeys();i++) {
    str<<"  keyword["<<i<<"]="<<get_iwave_iokeys()[i].keyword
       <<" index="<<get_iwave_iokeys()[i].rarrindex
       <<" input="<<get_iwave_iokeys()[i].input
       <<" active="<<get_iwave_iokeys()[i].active
       <<"\n";
  }
  return str;
}

namespace TSOpt {

  void IOTask(std::vector<TASK_RELN *> & tr,int order, bool fwd, IWaveInfo const & ic) {
    // loop over 0,...order, adj flag
    if (fwd) {
      for (int n = 0; n<=order; n++) {
	if (n==0) {
	  // all inputs kept for n=0, also all passive outputs, and active if order==0
	  for (int i=0;i<ic.get_num_iokeys();i++) {
	    if (ic.get_iwave_iokeys()[i].input || 
		(!(ic.get_iwave_iokeys()[i].input) && !(ic.get_iwave_iokeys()[i].active))) {
	      TASK_RELN * p = new TASK_RELN;
	      p->iwaveindex = 0;
	      p->keyword = ic.get_iwave_iokeys()[i].keyword;
	      p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	      p->input = ic.get_iwave_iokeys()[i].input;
	      tr.push_back(p);
	    }
	    if (order==0 && !(ic.get_iwave_iokeys()[i].input) && ic.get_iwave_iokeys()[i].active) {
	      TASK_RELN * p = new TASK_RELN;
	      p->iwaveindex = n;
	      p->keyword = ic.get_iwave_iokeys()[i].keyword;
	      p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	      p->input = ic.get_iwave_iokeys()[i].input;
	      tr.push_back(p);
	    }
	  }
	}
	else if (0<n && n<order) {
	  // only active inputs; keywords altered
	  for (int i=0;i<ic.get_num_iokeys();i++) {
	    if (ic.get_iwave_iokeys()[i].input && ic.get_iwave_iokeys()[i].active) {
	      TASK_RELN * p = new TASK_RELN;
	      p->iwaveindex = n;
	      std::ostringstream t;
	      t<<n;
	      p->keyword = ic.get_iwave_iokeys()[i].keyword + "_d"+t.str();
	      p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	      p->input = ic.get_iwave_iokeys()[i].input;
	      tr.push_back(p);
	    }
	  }
	}
	// active output only in comp 2^order-1
	else {
	  for (int i=0;i<ic.get_num_iokeys();i++) {
	    if (ic.get_iwave_iokeys()[i].input && ic.get_iwave_iokeys()[i].active) {
	      TASK_RELN * p = new TASK_RELN;
	      p->iwaveindex = n;
	      std::stringstream t;
	      t<<n;
	      p->keyword = ic.get_iwave_iokeys()[i].keyword + "_d"+t.str();
	      p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	      p->input = ic.get_iwave_iokeys()[i].input;
	      tr.push_back(p);
	    }
	    if (!(ic.get_iwave_iokeys()[i].input) && ic.get_iwave_iokeys()[i].active) {
	      TASK_RELN * p = new TASK_RELN;
	      p->iwaveindex = pow2(n)-1;
	      p->keyword = ic.get_iwave_iokeys()[i].keyword;
	      p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	      p->input = ic.get_iwave_iokeys()[i].input;
	      tr.push_back(p);
	    }
	  }
	}
      }
    }
    // adjoint branch 
    else {
      // order 0 doesn't make sense for adjoint
      if (order==0) {
	RVLException e;
	e<<"Error: IOTask constructor\n";
	e<<"  order 0 invalid for adjoint branch\n";
	throw e;
      }
      for (int n = 0; n<=order; n++) {
	if (n==0) {
	  // all inputs kept for n=0, also all passive outputs
	  for (int i=0;i<ic.get_num_iokeys();i++) {
	    if (ic.get_iwave_iokeys()[i].input || 
		(!(ic.get_iwave_iokeys()[i].input) && !(ic.get_iwave_iokeys()[i].active))) {
	      TASK_RELN * p = new TASK_RELN;
	      p->iwaveindex = 0;
	      p->keyword = ic.get_iwave_iokeys()[i].keyword;
	      p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	      p->input = ic.get_iwave_iokeys()[i].input;
	      tr.push_back(p);
	    }
	  }
	}
	else if (0<n && n<order) {
	  // only active inputs; keywords altered
	  for (int i=0;i<ic.get_num_iokeys();i++) {
	    if (ic.get_iwave_iokeys()[i].input && ic.get_iwave_iokeys()[i].active) {
	      TASK_RELN * p = new TASK_RELN;
	      p->iwaveindex = n;
	      std::ostringstream t;
	      t<<n;
	      p->keyword = ic.get_iwave_iokeys()[i].keyword + "_d"+t.str();
	      p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	      p->input = ic.get_iwave_iokeys()[i].input;
	      tr.push_back(p);
	    }
	  }
	}

	// active output only in comp = order
	else {
	  // only active inputs; keywords altered
	  for (int i=0;i<ic.get_num_iokeys();i++) {
	    if (ic.get_iwave_iokeys()[i].input && ic.get_iwave_iokeys()[i].active) {
	      TASK_RELN * p = new TASK_RELN;
	      p->iwaveindex = n;
	      std::ostringstream t;
	      t<<n;
	      p->keyword = ic.get_iwave_iokeys()[i].keyword + "_b"+t.str();
	      p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	      p->input = !ic.get_iwave_iokeys()[i].input;
	      tr.push_back(p);
	    }
	    if (!(ic.get_iwave_iokeys()[i].input) && ic.get_iwave_iokeys()[i].active) {
	      TASK_RELN * p = new TASK_RELN;
	      p->iwaveindex = pow2(n)-1;
	      p->keyword = ic.get_iwave_iokeys()[i].keyword;
	      p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	      p->input = !ic.get_iwave_iokeys()[i].input;
	      tr.push_back(p);	  
	    }
	  }
	}
      }
    }
  }

  void IOTaskWriter(std::vector<TASK_RELN *> const & tr, ostream & str) {
    for (int i=0;i<tr.size();i++) {
      str<<"index="<<tr[i]->iwaveindex<<" keyword="<<tr[i]->keyword<<" rarrindex="<<tr[i]->rarrindex<<" input=" << tr[i]->input<<"\n";
    }
  }
}

