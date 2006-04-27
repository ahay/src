#ifndef _divide1_h
#define _divide1_h

/* abstract data type */
typedef struct Div1 *div1;

div1 divide1_init(int n, float eps, float lam);
void divide1_close (div1 div);
void divide1 (div1 div, float* num, float* den, 
	      float* rat0, float* rat1, float* rat);
void divide2 (div1 div, int n2, float** num, float** den,  float** rat);
void divide3 (div1 div, int n3, int n2, 
	      float*** num, float*** den,  float*** rat);

#endif

/* 	$Id$	 */
