#define FLT_EPSILON 0.000001
 #include <math.h> 

/* ** */
/*  * Special-purpose eigensolvers for digital signal processing. */
/*  * Methods of this class solve small eigen-problems efficiently. */
/*  * @author Dave Hale, Colorado School of Mines */
/*  * @version 2008.09.04 */


  
/* /\** */
/*    * Computes eigenvalues and eigenvectors for a symmetric 2x2 matrix A. */
/*    * If the eigenvectors are placed in columns in a matrix V, and the  */
/*    * eigenvalues are placed in corresponding columns of a diagonal  */
/*    * matrix D, then AV = VD. */
/*    * @param a the symmetric matrix A. */
/*    * @param v the array of eigenvectors v[0] and v[1]. */
/*    * @param d the array of eigenvalues d[0] and d[1]. */

void solveSymmetric22(float a[2][2], float v[2][2], float d[2]) 
/*< test >*/
{

    // Copy matrix to local variables.
    float a00 = a[0][0];
    float a01 = a[0][1],  a11 = a[1][1];

    // Initial eigenvectors. 
    float v00 = 1.0,     v01 = 0.0;
    float v10 = 0.0,     v11 = 1.0;

    // If off-diagonal element is non-zero, zero it with a Jacobi rotation.
    if (a01!=0.0f) {
      float tiny = 0.1*sqrt(FLT_EPSILON); // avoid overflow in r*r below
      float c,r,s,t,u,vpr,vqr;
      u = a11-a00;
      if (fabs(a01)<tiny*fabs(u)) {
        t = a01/u;
      } else {
        r = 0.5*u/a01;
        t = (r>=0.0)?1.0/(r+sqrt(1.0+r*r)):1.0/(r-sqrt(1.0+r*r));
      }
      c = 1.0f/sqrt(1.0+t*t);
      s = t*c;
      u = s/(1.0+c);
      r = t*a01;
      a00 -= r;
      a11 += r;
      //a01 = 0.0f;
      vpr = v00;
      vqr = v10;
      v00 = vpr-s*(vqr+vpr*u);
      v10 = vqr+s*(vpr-vqr*u);
      vpr = v01;
      vqr = v11;
      v01 = vpr-s*(vqr+vpr*u);
      v11 = vqr+s*(vpr-vqr*u);
    }

    // Copy eigenvalues and eigenvectors to output arrays.
    d[0] = a00;
    d[1] = a11;
    v[0][0] = v00;  v[0][1] = v01;
    v[1][0] = v10;  v[1][1] = v11;

    // Sort eigenvalues (and eigenvectors) in descending order.
    if (d[0]<d[1]) {
      float dt = d[1];
      d[1] = d[0];
      d[0] = dt;
      float vt[2];
      vt[0]= v[1][0];vt[1]=v[1][1];
      v[1][0] = v[0][0];v[1][1] = v[0][1];
      v[0][0] = vt[0];v[0][1]=vt[1];
    }
  }
