#ifndef _xcorr_h
#define _xcorr_h

void xcorr_init (int nx, int n2, int maxshift);
void xcorr_close (void);
float xcorr (const float *x1, const float *x2);

#endif
