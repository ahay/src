#include "rkstream.hh"

namespace RVL {

  void makeRankStream(ofstream & outfile, int rk, string prefix) {
    /* create output file - all processes */
    stringstream nm;
    if (prefix.size() > 0) {
      if (prefix[prefix.size()-1] == '/') nm<<prefix<<"cout"<<rk<<".txt";
      else nm<<prefix<<"/cout"<<rk<<".txt";
    }
    else nm << "./cout"<<rk<<".txt";
    string outnm;
    nm>>outnm;
    outfile.open(outnm.c_str());
    outfile<<"output file opened on "<<outnm<<endl;
  }

}
