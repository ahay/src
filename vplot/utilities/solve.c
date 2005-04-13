/*
 * Copyright 1987 the Board of Trustees of the Leland Stanford Junior
 * University. Official permission to use this software is included in
 * the documentation. It authorizes you to use this file for any
 * non-commercial purpose, provided that this copyright notice is not
 * removed and that any modifications made to this file are commented
 * and dated in the style of my example below.
 */

/*
 *
 *  source file:   ./filters/utilities/solve.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

int solve (int pnot, int p1, int q1, int p2, int q2)
/*< compute intersection - floating point version >*/
{
    double          invslope;
    register int    qnot;
    register float  tempf;

    if (pnot == p1) return q1;
    if (pnot == p2) return q2;
    if (q1 == q2)   return q1;

    invslope = (q1 - q2) / ((double) (p1 - p2));
    tempf = (pnot - p1) * invslope + (double) q1;
    qnot = (tempf > 0) ? tempf + 0.5 : tempf - 0.5;
    return qnot;
}
