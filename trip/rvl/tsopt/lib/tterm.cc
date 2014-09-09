#include "tterm.hh"

namespace TSOpt {

  ostream & OrTimeTerm::write(ostream & fp) const {
    fp<<"--------- beg  OrTimeTerm::write ---------------\n";    
    fp<<" first term in OR (TimeTerm):\n";
    first.write(fp);
    fp<<"--------- end  OrTimeTerm::write ---------------\n";    
    return fp;
  }

  ostream & AndTimeTerm::write(ostream & fp) const {
    fp<<"--------- beg AndTimeTerm::write ---------------\n";    
    fp<<" first term in AND (TimeTerm):\n";
    first.write(fp);
    fp<<"--------- end AndTimeTerm::write ---------------\n";    
    return fp;
  }
}
