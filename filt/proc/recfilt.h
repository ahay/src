#ifndef _recfilt_h
#define _recfilt_h

#include <rsf.h>

/*
  Recfilt
  -------
  Recursive convolution (polynomial division). 
  nd     - data size
  nb     - filter size
  bb[nb] - filter
*/
void recfilt_init( int nd, int nb, float* bb);
void recfilt_lop( bool adj, bool add, int nx, int ny, float* xx, float*yy);
void recfilt_close (void);

#endif
