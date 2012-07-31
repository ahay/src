#ifndef SIGNAL_UNWRAPPER_2D_H
#define SIGNAL_UNWRAPPER_2D_H

#include "support.hh"

class SignalUnwrapper {

public:

	       SignalUnwrapper ();
	      ~SignalUnwrapper ();

			void setParams (VolumeParams* dp, float* data) {dp_ = dp; ptrToData_ = data;}

			void unwrap (float* zo);

private:

		VolumeParams*    dp_;
		float*           ptrToData_;
}; 

#endif
