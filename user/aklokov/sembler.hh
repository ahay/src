#ifndef SEMBLER_H
#define SEMBLER_H

class Sembler {

public:

    Sembler ();
   ~Sembler ();

	// semblance for whole trace
    static void getSemblanceForTrace (int tracesNum, float* data, float* dataSq, int zNum, int sembWindow, float* semb, int k = 1); 

};
#endif
