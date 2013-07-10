#ifndef __IWAVE_SRC_COORDS
#define __IWAVE_SRC_COORDS

#include "seamx_headers.hh"
#include "std_cpp_includes.hh"

namespace TSOpt {

  /** interface for coord access */
  class Coords {
  public:
    /** get current value of source acquisition par*/
    virtual void getCoords(RPNT & c) const = 0; 
    /** get current index of record/model-panel */
    virtual int getPanelInd() const = 0;
    /** get the number of model-panels/records*/
    virtual int getPanelNum() const = 0;
  };

}

#endif
