#ifndef _dijkstra_h
#define _dijkstra_h

void dijkstra_init(int m1, int m2);
void dijkstra_close(void);
void dijkstra(int s1, int s2, float **ud, float **lr);
void dijkstra_start(int s1, int s2);
bool dijkstra_next(int *ud, int *lr);

#endif
