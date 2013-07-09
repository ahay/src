#ifndef __RVL_ALG_LINEAR_SOLVER_H_
#define __RVL_ALG_LINEAR_SOLVER_H_

#include "linop.hh"
#include "alg.hh"

namespace RVLUmin {

  using namespace RVL;

  /**  Compute the first k eigenpairs of a linear operator A.   
       Unfortunately, the storage must be fully allocated on construction,
       so a fixed number of eigenpairs must be specified up front.
       Attempting to iterate too many times results in algorithmic failures.
       
       The linear operator A can be any form.  The resulting V will be a product vector
       in the kth cartesian power of the domain of A.  The upper hessenberg matrix
       H is stored in compressed column form in a valarray.
  */  
  template<class Scalar>
  class ArnoldiSolverStep : public Algorithm {
  protected: 
    LinearOp<Scalar> & A;
    int k;
    Vector<Scalar> v;
    Vector<Scalar> f;
    CartesianPowerSpace<Scalar> vspc;
    ProductVector<Scalar> V; // length k.  Will push_back new column DC each iter.
    valarray<Scalar> H; // arranged columnwise.  Will add k+1 elements each iter.
    int hcount;
    int iter;

    ArnoldiSolverStep();

  public:
    ArnoldiSolverStep(LinearOp<Scalar> & _A, int num_eigs, const Vector<Scalar> & init_guess)
      :A(_A), k(num_eigs), v(init_guess), f(A.getRange()), 
       vspc(k, A.getDomain()), V(vspc), H(k*(k+1)/2+k), hcount(0), iter(0) {
      Vector<Scalar> w(A.getRange());
      v.scale(Scalar(1.0/v.norm()));
      A.applyOp(v,w); // w <- Av
      Scalar alpha = v.inner(w); // alpha <- v'*Av
      f.copy(w); f.linComb(-alpha, v); // f<- Av - alpha v
      Scalar c = v.inner(f); // c <- v'*(Av-alpha v)
      f.linComb(-c, v); // f<- f-c v
      alpha+=c;
      V[iter].copy(v);
      H[0] = alpha;
      hcount++;
    }

    void run() {
      // note: no-op if you're done - needs coupling to terminator in driver
      if( iter >= k-1 ) return;
      int i = 0;
      Vector<Scalar> w(A.getRange());
      iter++;
      Scalar beta = f.norm();
      H[hcount] = beta; // fill in subdiagonal entry on previous column;
      hcount++;
      v.scale(Scalar(1.0)/beta, f);
      V[iter].copy(v);
      A.applyOp(v,w);
      valarray<Scalar> h(iter+1), c(iter+1);
      // h = V'*w
      for(i=0; i <= iter; i++)
	h[i] = V[i].inner(w);
      // f = w-V*h
      f.copy(w);
      for(i=0; i<= iter; i++)
	f.linComb(-h[i], V[i]);
      // c = V'*f
      for(i=0; i <= iter; i++)
	c[i] = V[i].inner(f);
      // f = f-V*c
      for(i=0; i<= iter; i++)
	f.linComb(-c[i], V[i]);
 
      // fill in upper triangle
      for(i = 0; i < c.size(); i++, hcount++) {
	H[hcount] = c[i]+h[i];
      }
    }

    /** Access the upper Hessenberg matrix H, stored in a valarray in column-major order,
	skipping all the zero entries.  
	H[0] = H(1,1), H[1] = H(2,1), 
	H[3] = H(1,2), H[4] = H(2,2), H[5] = H(3,2),
	H[6] = H(1,3), ....
    */
    const valarray<Scalar> & getHessenberg() { return H; }
  
    /** Return the current number of entries set in H. */
    int getHNumEntries() { return hcount; }

    /** Access the orthogonal vectors V.  Note that only the 
	vectors 0,...,iter have been set.
    */
    const ProductVector<Scalar> & getOrthogonalVecs() { return V; }

    int getIterationCount() { return iter; }
    int getMaxIter() { return k; }
  };

  template<typename Scalar>
  class Arnoldi: public Algorithm {



}
#endif
