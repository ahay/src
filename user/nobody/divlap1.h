#ifndef _divlap1_h
#define _divlap1_h

/* abstract data type */
typedef struct Divlap *divlap;

divlap divlap1_init(int n, float eps, float lam);
void divlap1_close (divlap div);
void divlap1 (divlap div, float* num, float* den, float* ref, 
	      float* rat0, float* rat1, float* rat);
void divlap2 (divlap div, int n2, 
	      float** num, float** den, float** ref, float** rat);
void divlap3 (divlap div, int n3, int n2, 
	      float*** num, float*** den, float*** ref, float*** rat);

#endif

/* 	$Id: divlap1.h,v 1.1 2004/05/25 00:46:12 fomels Exp $	 */
