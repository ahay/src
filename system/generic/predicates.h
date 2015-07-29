/*
 * File: predicates.h
 * ------------------
 * Provides an interface for J.Shewchuk's geometric predicates.
 */
#ifndef _predicates_h
#define _predicates_h

/*****************************************************************************/
/*  exactinit()   Initialize the variables used for exact arithmetic.        */
/*****************************************************************************/
void exactinit(void);

/*****************************************************************************/
/*  orient2dfast()   Approximate 2D orientation test.  Nonrobust.            */
/*****************************************************************************/
double orient2dfast(double* pa, double* pb, double* pc);

/*****************************************************************************/
/*  orient2d()   Adaptive exact 2D orientation test.  Robust.                */
/*****************************************************************************/
double orient2d(double* pa, double* pb, double* pc);

/*****************************************************************************/
/*  incirclefast()   Approximate 2D incircle test.  Nonrobust.               */
/*****************************************************************************/
double incirclefast(double* pa, double* pb, double* pc, double* pd);

/*****************************************************************************/
/*  incircle()   Adaptive exact 2D incircle test.  Robust.                   */
/*****************************************************************************/
double incircle(double* pa, double* pb, double* pc, double* pd);

#endif
