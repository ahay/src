#ifndef ALMOST_H
#define ALMOST_H

#include <limits.h>
#include <float.h>
#include <math.h>

class Almost
{
public:
  /* return 1 if x is between x1 and x2, including almost equality, else 0 */
  static int between (double x, double x1, double x2)
  {
    return (Almost::cmp (x, x1) * Almost::cmp (x, x2) <= 0) ? 1 : 0;
  }

  /* return 0 if between x1 and x2, -1 if outside and closer to x1,
   * 1 if outside and closer to x2 */
  static int outside (double x, double x1, double x2)
  {
    /* return 0 if between x1 and x2, -1 if outside and closer to x1,
     * 1 if outside and closer to x2 */
    int i = 0;

    if (Almost::between (x, x1, x2))
      i = 0;
    else if (Almost::between (x1, x, x2))
      i = -1;
    else if (Almost::between (x2, x, x1))
      i = 1;
    return i;
  }

  /* return 1 if looks like a legitimate float, else 0 */
  static int is_a_float (float r)
  {
    return ((fabs (r) < FLT_MAX)
            && (fabs (r) > FLT_MIN || Almost::zero (r))) ? 1 : 0;
  }

  /* return 1 if zero or very close to zero, for float precision  */
  static int zero (double r)
  {
    return (fabs (r) < 10. * FLT_MIN) ? 1 : 0;
  }

  /* return 1 if these are close to being equal for float precision, else 0 */
  static int equal (double r1, double r2)
  {
    return (!Almost::cmp (r1, r2)) ? 1 : 0;
  }

  /* return 1 if r1>r2, -1 if r1<r2, 0 if r1==r2, within float precision
   * always assume r1 == r2, if within precision
   * Test "if (Almost::cmp(r1,r2)>0)" for a strict r1>r2 */
  static int cmp (double r1, double r2)
  {
    if (Almost::zero (r1) && Almost::zero (r2))
      return 0;
    if (r1 - (10. * FLT_EPSILON) * fabs (r1) >
        r2 + (10. * FLT_EPSILON) * fabs (r2))
      return 1;
    if (r1 + (10. * FLT_EPSILON) * fabs (r1) <
        r2 - (10. * FLT_EPSILON) * fabs (r2))
      return -1;
    return 0;

  }

  /* take reciprocal and return large number if dividing by zero */
  static double reciprocal (double a)
  {
    if (!Almost::zero (a))
      a = 1. / a;
    else
      a = (a >= 0.) ? 0.1 * FLT_MAX : -0.1 * FLT_MAX;
    return a;
  }

  /* divide top by bottom, with clipping of overflow; 0/0 = 1 */
  static double divide (double top, double bottom)
  {
    double sign;

    if ((top > 0. && bottom < 0.) || (top < 0. && bottom > 0.))
      sign = -1.;
    else
      sign = 1.;
    top = fabs (top);           // remove signs to simplify if's
    bottom = fabs (bottom);
    if (bottom >= 1.            // Might underflow, but don't care.
        || top < bottom * 0.1 * FLT_MAX)        // (bottom<1) and won't overflow
      return sign * top / bottom;       // safe
    else if (Almost::equal (top, bottom)) {
      if (Almost::zero (top))   /* define 0/0 = 1 */
        return 1.;
      else
        return sign;            /* ratio is within precision of unity */
    } else                      /* now top >= bottom*0.1*FLT_MAX, clip overflow */
      return sign * 0.1 * FLT_MAX;
  }

  /* pseudo random number, idum is a seed */
  static double ran0 (long *idum)
  {
    long k;
    double ans;

    static long IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836, MASK =
      123459876;
    static double AM = (1.0 / (IM));

    *idum ^= MASK;
    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if (*idum < 0)
      *idum += IM;
    ans = AM * ((double) (*idum));
    *idum ^= MASK;
    return ans;
  }

  /* pseudo random number, keeps own static seed */
  static double ran (void)
  {
    static long iseed = 294257;

    return Almost::ran0 (&iseed);
  }

  /* one plus a small number, detectable within float precision */
  static double one_plus_epsilon (void)
  {
    return (1. + 15. * FLT_EPSILON);
  }

  /* a small number, detectable within float precision */
  static double epsilon (void)
  {
    return (15. * FLT_MIN);
  }

  /* return 1 if corresponding entries are almost equal, else 0 */
  static int equal_arrays (float *a1, float *a2, long n)
  {
    int same;
    long i;

    for (i = 0, same = 1; i < n && same; i++) {
      same = (same && Almost::equal (a1[i], a2[i])) ? 1 : 0;
    }
    return same;
  }

  /* return 1 if all appear to be floats */
  static int float_arrays (float *a, long n)
  {
    int good;
    long i;

    for (i = 0, good = 1; i < n && good; i++) {
      good = (good && Almost::is_a_float (a[i])) ? 1 : 0;
    }
    return good;
  }

};

#endif
