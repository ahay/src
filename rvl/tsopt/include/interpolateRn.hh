#ifndef INTERP_RN_HH
#define INTERP_RN_HH

namespace TSOpt {

  /** An Interpolator class, providing a static member function capable of 
      interpolating RnStates' internal data. Currently a test class performing
      piecewise constants
  */
  class RnStateInterp {
  public:
   static void getInterpState(RnState & r, const Time & tt, std::deque<RnState> & stateHist) {


      int nu = stateHist.front().getStateDim();
      int nc = stateHist.front().getControlDim();
      // int it = static_cast<const StdDiscreteTime&>(t).getint();

      // make an array of size numNodes to hold interp data
      std::vector<float> data;
      std::vector<float> control;
      std::vector<StdDiscreteTime> timeList;

      // find out what data elements should be...
      // seriously may consider threading here to parallelize code
      for (int i=0; i<nu; ++i) {

	int count = 0;
	std::deque<RnState > :: iterator j;
	
	// <ME?> have to do this for the controls as well!!
	for(j=stateHist.begin(); j != stateHist.end(); j++) {     
	  data.push_back((*j).getStateElement(i));
	  control.push_back((*j).getControlElement(i));
	  timeList.push_back(static_cast<const StdDiscreteTime &>((*j).getTime()));
	  ++count;
	}

	// initialize interpolator
	// TEST: do piecewise constants
	
	int k;  
	for (k=0; k<count-1; ++k) {
	  if ( ( (tt > (timeList.at(k))) || (tt == (timeList.at(k))) ) && 
	       (tt < (timeList.at(k+1)))  )
	    {
	      break;
	    }
	}


	// set proper part of state
	r.setStateElement(data[k],i); 
	r.setControlElement(control[k],i);

      }


   }


 };

}

#endif
