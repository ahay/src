#ifndef _divide2_h
#define _divide2_h

/* abstract data type */
typedef struct Div2 *div2;

div2 divide2_init(int n1, int n2, float eps1, float eps2);
void divide2_close (div2 div);
void divide2 (div2 div, float** num, float** den,  float** rat);

#endif

/* 	$Id$	 */
