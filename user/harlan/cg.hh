/**

  These C++ classes allow you to optimize quadratic objective
  functions, including least-squares inversion with model
  damping and preconditioning.  Nonlinear inversion is
  supported by iterative linearization.

  To minimize the simplest quadratic form 0.5 x'Hx + b'x, you
  need to implement a model vector from CgVector, with a few
  essential operations such as copy, scale, and add.  Then
  implement a CgQuadInterface to multiply the Hessian H times
  any arbitrary model x. Create a CgQuad solver with the
  interface and a model b as arguments to the constructor.
  Call solve() for a solution after a specified number of
  conjugate-gradient iterations

  You can also solve a more convenient objective function for
  inverse problems: (Fx-e)'N(Fx-e) + (x+m)'M(x+m), where F is
  a forward transform, e is the data error, m is a starting
  model, and N and M are inverse covariance matrices for the
  noise and model.  Implement a data vector from CgVector.
  For the solver, you implement CgTransLinearInterface, and
  define the forward and adjoint versions of the transform.
  (The Hessian and b vector are implicitly derived from your
  new operations.)  Provide this interface, the starting model
  vector, and a data vector to the constructor of a CgTrans
  solver.  The inverse covariance matrices, if not identities,
  can be defined in the vector classes.  You can also define
  approximate inverse Hessians as preconditioning to speed
  your convergence.

  If you implement the more general nonlinear CgTransInterface
  for a nonlinear forward transform, then you can use CgTrans
  to solve the non-quadratic objective function by a Gauss
  Newton method.

  See the prototypes below for more details.
  To run the example code, set "#define TESTCG".

  Author, 1997:
  William S. Harlan              harlan[at]advance.com; www.lgc.com
  Landmark                       harlan[at]sep.stanford.edu;
  7409 S. Alton Court              sepwww.stanford.edu/oldsep/harlan
  Englewood, CO 80112, USA

  This code is being distributed to universities under the conditions
  of the GNU FSF library public license.

*/

#ifndef CG_H
#define CG_H 1

class ProgressMonitor;

/** Models and data should both implement this interface */
class CgVector
{
public:
  /** create new exact copy of this vector, all new memory */
  virtual CgVector * clone () const = 0;

  /** scale this vector by scalar factor */
  virtual void scale (double factor) = 0;

  /** add vector to this one */
  virtual void add (const CgVector & toBeAdded) = 0;

  /** return dot product of this vector with another vector */
  virtual double dot (const CgVector & toBeDotted) const = 0;

  /** redefine the destruction of this object, unless you want a memory leak */
    virtual ~ CgVector ()
  {
  };

  /** Optionally multiply a vector by the inverse covariance matrix.
    Used by CgTrans but not CgQuad.
    The default no-op is equivalent to assuming identity. */
  virtual void multiply_inverse_covariance ()
  {
  };

  /** Optionally apply a hard constraint (such as inequality) to current vector
     This is used only by the nonlinear optimization.  Default does nothing.*/
  virtual void apply_constraint ()
  {
  };

  /* Override these composite functions if you want something more efficient
     for the vector operations.  Used on CgTrans model but not on the data. */
  virtual void copy (const CgVector & input)
  {
    this->scale (0.);
    this->add (input);
  }
  virtual CgVector *clone_zero () const
  {
    CgVector *out = this->clone ();
      out->scale (0.);
      return out;
  }
  virtual void scale_add (double scale, const CgVector & other,
                          CgVector & scratch)
  {
    scratch.copy (other);
    scratch.scale (scale);
    this->add (scratch);
  }

  // printout self for debugging (default prints out magnitude)
  virtual void printme (const char *label) const;

  CgVector () {
  };                            // some compilers want an explicit default

  // Test the vector properties of an instance of a CgVector.
  static void test(const CgVector& sample) ;
};

/** This interface must be implemented to use the CgQuad solver */
class CgQuadInterface
{
public:
  /** Define this to multiply model x by the hessian H. */
  virtual void hessian (CgVector & x) const = 0;

  /** To speed convergence can multiply model by approximate inverse Hessian
    This no-op default will work too. */
  virtual void inverse_hessian (CgVector & model) const
  {
  };

};

/** This class optimizes quadratic objective functions with
    conjugate gradients */
class CgQuad
{
public:
  /** Method 1: Implement the CgQuadInterface and pass to this constructor.
    Pass a vector b to the constructor.   CgQuad clones its own copy.

    Hx + b = 0 are the Normal Equations of your least-squares problem.
    Minimizes quadratic function 0.5 x'Hx + b'x .  Return new model x.

    This solver calls CgVector and CgQuadInterface methods.
    Does not call CgVector::multiply_inverse_covariance().  */

  /** instantiate solver with this. */
  CgQuad (const CgQuadInterface & cqi, const CgVector & b):_cqi (&cqi)
  {
    _b = b.clone ();
  };

  /** return a new solution after the number of conjugate grad iterations
      you delete with delete. */
  CgVector *solve (long numberIterations) const;

  /** optionally change the definition of b after creation */
  void set_b (const CgVector & b)
  {
    delete _b;
    _b = b.clone ();
  };

  /** If true, then run extra code to test internal consistency of operators */
  static void set_debug (int debug)
  {
    _debug = debug;
  }

  /** If true, then run extra code to test internal consistency of operators */
  static int is_debug ()
  {
    return _debug;
  }

  /** dtor */
  virtual ~ CgQuad () {
    delete _b;
  };
private:
  CgVector * _b;
  const CgQuadInterface *_cqi;

  CgQuad () {
  };                            // don't use default ctor
  CgQuad & operator= (const CgQuad &);  // don't copy
  CgQuad (const CgQuad &);      // don't copy

  static int _debug;            // run extra code for internal consistency
};

/** This interface must be implemented to use the CgTrans solver in non-linear
   form.  The CgTransLinearInterface implements the simpler linear version.*/
class CgTransInterface
{
public:
  /** Nonlinear forward transform: data = f(model) */
  virtual void forward_nonlinear (CgVector & data,
                                  const CgVector & model) const = 0;

  /** Linearized forward transform for a given reference model
    data = F model ; data ~= f(model+modelReference) - f(modelReference) */
  virtual void forward_linearized (CgVector & data, const CgVector & model,
                                   const CgVector & modelReference) const = 0;

  /** The adjoint of the linearized forward transform: model = F' data
      Output will be initialized to 0.
   */
  virtual void adjoint_linearized (const CgVector & data, CgVector & model,
                                   const CgVector & modelReference) const = 0;

  /** To speed convergence can multiply model by approximate inverse Hessian
    This no-op default will work too. */
  virtual void inverse_hessian (CgVector & model, CgVector & modelReference) const
  {
  }
};

/** Implement this interface as a linear simplification of the more
   general CgTransInterface for the forward and adjoint transforms. */
class CgTransLinearInterface:public CgTransInterface
{
public:
  /** Define linear forward transform: data = F model
    You can assume data is initialized to zero */
  virtual void forward (CgVector & data, const CgVector & model) const = 0;

  /** Adjoint transform: model = F' data
    You can assume model is initialized to zero going into adjoint(). */
  virtual void adjoint (const CgVector & data, CgVector & model) const = 0;

  /** To speed convergence can multiply model by approximate inverse Hessian
    This no-op default will work too. */
  virtual void inverse_hessian (CgVector & model) const
  {
  };

  /* DO NOT CHANGE THE FOLLOWING, which implement the nonlinear interface
     by ignoring the reference model */
  virtual void forward_nonlinear (CgVector & data, const CgVector & model) const
  {
    forward (data, model);
  }
  virtual void forward_linearized (CgVector & data, const CgVector & model,
                                   const CgVector & modelReference) const
  {
    forward (data, model);
  }
  virtual void adjoint_linearized (const CgVector & data, CgVector & model,
                                   const CgVector & modelReference) const
  {
    adjoint (data, model);
  }
  virtual void inverse_hessian (CgVector & model, CgVector & modelReference) const
  {
    inverse_hessian (model);
  };
};

/** Implement this interface to transform a model to data */
class CgTrans:private CgQuadInterface
{
public:
  /** Method 2: Define the linear transforms in CgTransLinearInterface
    and create this solver with the appropriate data and starting model.
    Solve damped, overdetermined system of equations: F(m+x) = d,
    and return m+x.  x is the perturbation to the starting model m.
    Equivalently, minimize [F(m+x)-data]'N[F(m+x)-data] + (x+m)'M(x+m), where
    The Hessian H = (F'NF + M) and b = -F'N(data-Fm) + Mm, where
    Hx + b = 0 are the Normal equations of your least-squares problem.
    N and M are the INVERSE covariance matrices for noise and the model.
    This function returns a new copy of m+x.
    To damp x'Mx instead of (x+m)'M(x+m), change dampPerturb to 1;
    To return x instead of m+x, change addReferenceModel to 0.
    Calls CgVector methods, including Cgvector::multiply_inverse_covariance().
    */
  CgTrans (const CgTransInterface & cti, const CgVector & data0,
           const CgVector & model0, int addReferenceModel = 1,
           int dampPerturb = 0);

  /** return a new solution after the specified number of conj. grad. iterations
     This is the only solver you need for the linear solution. */
  CgVector *solve (long numberIterations) const;

  /** Method 3: Solve nonlinear system of equations:  f(m) = d.
    Define the CgTransInterface methods for the nonlinear transform
    and for the linearized forward and adjoint transforms.
    Solve nonquadratic objective function with Gauss Newton iterations:
    [f(m+x)-data]'N[f(m+x)-data] + (m+x)'M(m+x)
    Always returns full solution m+x, not a perturbation x;
    Iterative linearization of f(m+x) ~= f(m) + Fx makes the objective
    function quadratic in x: [f(m)+Fx-data]'N[f(m)+Fx-data] + (m+x)'M(m+x)
    x is solved with the specified number of conjugate gradient iterations.
    This perturbation is then scaled after searching the nonquadratic
    objective function with the specified number of line search iterations.
    The scaled perturbation x is added to the previous reference model m
    to update the new reference model m.  Relinearization is repeated for
    the specified number of linearization iterations. Cost is proprotional to
    linearizationIterations*( 2* conjugateGradIterations +
                               lineSearchIterations );
    Hard constraints, if any, will be applied during line searches, and
    to the final result.
    "Line search error" is an acceptable fraction of imprecision
    in the scale factor for the line search.   The default will
    be good enough that the maximum number of linearization iterations
    will be used.
    */
  CgVector *solve (int conjugateGradIterations,
                   int lineSearchIterations,
                   int linearizationIterations,
                   double lineSearchError = 0.00001,
                   ProgressMonitor * pm = 0);

  // You may use the following for customized nonlinear optimization.
  /** evaluate objective function for a given model
    if (is_a perturbation) then the input x is added to the reference
    model m before evaluation.  Else, the input x is interpreted as m+x.
     [f(m+x) - d]'N [f(m+x) -d]  + (m+x)'N(m+x)
     if (dampPerturb) then return
     [f(m+x) - d]'N [f(m+x) -d]  + x'Nx
     */
  double objective_function (const CgVector & x, int is_a_perturbation) const;

  /** Returns true if linearized operations pass definition of adjoint.
     dot(d,Fm) == dot(F'd,m) */
  int adjoint_is_good () const;

  /** Minimize objective_function(alpha*x) and return alpha.
    alphamin < alpha < alphamax.  Perform a maximum of nter iterations,
    or when the error (dalpha) in alpha satisfies
    dalpha/(alphamax-alphamin) <= okerror or
    dalpha < okfraction*(alpha-dalpha).  Returns optimum scale factor alpha.
    Performs a minimum of 10 iterations. I recommend at least 20. */
  double line_search (const CgVector & x, double alphamin, double alphamax,
                      double okerror, double okfraction, long nter) const;

  // Implementation of CgQuadInterface
  virtual void hessian (CgVector & x) const;    // H = F'NF + M
  // Implementation of CgQuadInterface
  virtual void inverse_hessian (CgVector & model) const
  {
    _cti->inverse_hessian (model, *_model0);
  };

  // dtor
    virtual ~ CgTrans ()
  {
    delete _cq;
    delete _model0;
  };

private:
  CgQuad * _cq;                 // corresponding quadratic solver
  CgVector *_model0;            // reference model m, for linearization
  const CgTransInterface *_cti; // passed by constructor
  const CgVector *_data0;       // data to be inverted
  // configuration of this solver
  int _dampPerturb;             // damp perturbations only
  int _addReferenceModel;       // solvers add reference model
  // create internal
  CgVector *create_b () const;  // b = -F' N [data - f(m)] + Mm
  // avoid some defaults
  CgTrans & operator= (const CgTrans &);        // don't copy
  CgTrans (const CgTrans &);    // don't copy
  CgTrans () {
  };                            // default constructor
  // internal convenience functions
  static void index (int *ix, const double *rx, long n);        // slothsort
  void line_sample (double *xmin, double *xerr, double *xnew, long *iter, double fnew) const;   // sample line search intellegently
  static int parabola_min (double *xmin,        // find the minimum of best parabola
                           double x1, double f1, double x2, double f2,
                           double x3, double f3);
  static int almost_equal (double r1, double r2);       // ignore roundoff
  static int almost_cmp (double r1, double r2); // cmp(a,b) w/ roundoff
  static int almost_zero (double r);    // define fabs(zero) < 10*FLT_MIN
};

/** to minimize a custom objective function, define data as a small scalar: */
class CgScalar:public CgVector
{
  double _d;
public:
  CgScalar (double d = 0.):_d (d) {;
  }                             // ctor
  // can use default copy ctor, copy assign, and dtor
  //CgScalar(const CgScalar&);CgScalar& operator=(const CgScalar&);~CgScalar();
  // set value
  CgScalar & operator= (double d)
  {
    _d = d;
    return *this;
  }                             //set value
  void set (double d)
  {
    _d = d;
  }                             // set same value
  // get value value
  double operator* () const
  {
    return _d;
  }                             // dereference returns value
  double get () const
  {
    return _d;
  }                             // get same value
  operator  double () const
  {
    return _d;
  }                             // (conversion) get value
  // CgVector methods
  virtual CgVector *clone () const
  {
    CgScalar *sf = new CgScalar (*this);
      return (CgVector *) sf;
  }
  virtual void scale (double factor)
  {
    _d *= factor;
  }
  virtual void add (const CgVector & rhs)
  {
    _d += ((const CgScalar &) rhs).get ();
  }
  virtual double dot (const CgVector & rhs) const
  {
    return (_d * ((const CgScalar &) rhs).get ());
  }
};

#endif
