#ifndef _divide1_h
#define _divide1_h

void divide1_init(int n, float eps, float lam);
void divide1_close (void);
void divide1 (float* num, float* den, float* rat0, float* rat1, float* rat);
void divide2 (int n2, float** num, float** den,  float** rat);
void divide3 (int n3, int n2, float*** num, float*** den,  float*** rat);

#endif
