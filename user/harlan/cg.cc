/* See Cg.h for documentation of this class */

#include "almost.hh"
#include "cg.hh"
#include "progress_monitor.hh"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static const int DEBUG = 0;

int CgQuad::_debug = 0;

  // Test the vector properties of an instance of a CgVector.
void 
CgVector::test(const CgVector& sample) 
{
  CgVector *t = sample.clone ();
  CgVector *tt = sample.clone_zero ();
  t->scale_add (-3., sample, *tt);
  double ss = sample.dot (sample);
  if (Almost::zero(ss)) {
    // fprintf (stdout, "Cannot test vector with zero magnitude\n");
    return;
  }
  double check = t->dot (sample) / ss;
  if (!(-check > 1.999 && -check < 2.001)) {
    fprintf (stderr, "Numerical problem: dot product should be -2 = %g\n",
	     check);
    assert (0);
  }
  t -> multiply_inverse_covariance ();  // just make sure it can be called.
  tt -> apply_constraint ();  // just make sure it can be called.
  assert( t->dot(*t) > 0);
  assert( tt->dot(*tt) > 0);
  delete t;
  delete tt;
}

void
CgVector::printme (const char *label) const
{
  printf ("%s: magnitude=%g\n", label, sqrt (this->dot (*this)));
  fflush (stdout);
}

CgVector *
CgQuad::solve (long nter) const
{
  /* minimize quadratic function 0.5 x'Hx + b'x , and return vector x.
     b is the gradient for x=0, and H is the hessian
     algorithm:
     x = p = u = 0;
     beta = 0;
     g = b;
     a = A g;
     q = H a;
     do {
     p = -a + beta p
     u = Hp = -q + beta u
     alpha = -p'g/p'u
     x = x + alpha p
     composite functions  if (done) return x
     g = H x + b = g + alpha u
     a = A g
     q = H a
     beta = p'H A g / p'H p = p'q/p'u
     }
   */

  double beta, alpha, pu, pg;
  long iter;
  if (DEBUG)
    printf ("Starting conj grad solution with %ld iterations\n", nter);
  double bb = _b->dot (*_b);
  if (bb < 10. * FLT_MIN) {
    if (DEBUG)
      fprintf (stdout,
               "Right hand side of conjugate gradient equation is negligible.  "
               "Not solving.\n");
    fflush (stdout);
    return _b->clone_zero ();
  }
  if (CgQuad::is_debug()) {                 
    CgVector::test(*_b);
  }
  CgVector *t = _b->clone ();   // scratch
  CgVector *x = _b->clone_zero ();
  CgVector *p = _b->clone_zero ();
  CgVector *u = _b->clone_zero ();
  beta = 0.;
  CgVector *g = _b->clone ();
  CgVector *a = g->clone ();
  _cqi->inverse_hessian (*a);
  CgVector *q = a->clone ();
  _cqi->hessian (*q);
  for (iter = 0; iter < nter; iter++) {
    // fprintf(stdout, "beta = %g; ",beta); fflush(stdout);
    p->scale (beta);
    p->scale_add (-1., *a, *t);
    u->scale (beta);
    u->scale_add (-1., *q, *t);
    pg = p->dot (*g);
    pu = p->dot (*u);
    // fprintf(stdout, "pg = %g pu = %g;",pg,pu); fflush(stdout);
    if (fabs (pg) < 20. * FLT_MIN || fabs (pu) < 20. * FLT_MIN)
      goto done;
    alpha = -pg / pu;
    // fprintf(stdout, "alpha = %g\n",alpha); fflush(stdout);
    x->scale_add (alpha, *p, *t);
    if (iter == nter - 1)
      goto done;
    g->scale_add (alpha, *u, *t);
    a->copy (*g);
    _cqi->inverse_hessian (*a);
    q->copy (*a);
    _cqi->hessian (*q);
    beta = p->dot (*q) / pu;
    if (beta > 5.)
      beta = 5.;
    if (beta < -5.)
      beta = -5.;
  }
done:
  // fprintf(stdout,"\n"); fflush(stdout);
  delete p;
  delete u;
  delete g;
  delete a;
  delete q;
  delete t;
  return x;
}

/** unique constructor */
CgTrans::CgTrans (const CgTransInterface & cti, const CgVector & data0, 
		  const CgVector & model0, int addReferenceModel, 
		  int dampPerturb):
_cti (&cti), _data0 (&data0), _dampPerturb (dampPerturb),
_addReferenceModel (addReferenceModel)
{
  _model0 = model0.clone ();
  _model0->apply_constraint ();
  CgVector *b = this->create_b ();      // requires _model0 first
  _cq = new CgQuad (*this, *b); // implements the CgQuadInterface
  delete b;
}

/** solve quadratic problem */
CgVector *
CgTrans::solve (long numberIterations) const
{
  if (CgQuad::is_debug ()) {
    adjoint_is_good ();
  }
  CgVector *answer = _cq->solve (numberIterations);
  if (_addReferenceModel) {
    answer->add (*_model0);
  }
  return answer;
}

/** Solve nonquadratic objective function with Gauss Newton iterations */
CgVector *
CgTrans::solve (int conjugateGradIterations, int lineSearchIterations,
                int linearizationIterations, double lineSearchError,
                ProgressMonitor * pm)
{
  if (CgQuad::is_debug ()) {
    adjoint_is_good ();
  }
  if (DEBUG)
    printf ("Begin nonlinear optimization\n");
  int iter;
  double alphamin, alphamax, okerror, okfraction;
  // require good accuracy so that specified iterations are used.
  alphamin = 0.;
  alphamax = 1.1;
  okerror = lineSearchError;
  okfraction = lineSearchError;
  if (pm) {
    if (pm->set (this, 0.)) {
      goto done;
    }
  }
  for (iter = 0; iter < linearizationIterations; iter++) {
    // use parent solver to avoid adding reference model
    if (pm) {
      if (pm->set (this, double (iter) / double (linearizationIterations)))
      {
        goto done;
      }
    }
    CgVector *perturbation = _cq->solve (conjugateGradIterations);
    double pp;
    pp = perturbation->dot (*perturbation);
    if (DEBUG)
      perturbation->printme ("linear perturbation");
    if (almost_zero (pp)) {
      if (DEBUG) {
        printf ("done early with linearization, iteration %d\n", iter + 1);
      }
      delete perturbation;
      goto done;
    }
    double scale;
    if (DEBUG)
      printf ("Beginning line search\n");
    if (pm) {
      if (pm->
          set (this, (double (iter) + 0.5) /double (linearizationIterations)))
      {
        goto done;
      }
    }
    scale = this->line_search (*perturbation, alphamin, alphamax, okerror,
                               okfraction, lineSearchIterations);
    if (DEBUG)
      printf ("Scale factor = %g\n", scale);
    if (almost_zero (scale)) {
      if (DEBUG) {
        printf ("finished line search early, iteration %d\n", iter + 1);
      }
      delete perturbation;
      goto done;
    }
    perturbation->scale (scale);
    if (DEBUG)
      perturbation->printme ("scaled perturbation");
    perturbation->add (*_model0);
    if (DEBUG)
      perturbation->printme ("scaled perturb'n plus reference model");
    perturbation->apply_constraint ();
    if (DEBUG)
      perturbation->printme ("updated model with constraint");

    delete _model0;
    _model0 = perturbation;
    perturbation = 0;
    CgVector *b = this->create_b ();
    _cq->set_b (*b);
    delete b;
  }
done:
  if (pm)
    pm->set (this, 1.);
  return _model0->clone ();
}

CgVector *
CgTrans::create_b () const
{
  if (DEBUG)
    _model0->printme ("reference model to begin calculating b");
  CgVector *data = _data0->clone_zero ();       // data is data with zeros
  CgVector *model = _model0->clone ();  // model is m
  CgVector *outmodel = _model0->clone_zero ();  // outmodel is m with zeros
  _cti->forward_nonlinear (*data, *model);      // data is f(m) (defaults to Fm)
  if (DEBUG)
    data->printme (" current modeled data: f(model)");
  data->scale (-1);             // data is -f(m)
  data->add (*_data0);          // data is e = d - f(m)
  if (DEBUG)
    data->printme (" current error: data-f(model)");
  data->multiply_inverse_covariance (); // data is N e
  _cti->adjoint_linearized (*data, *outmodel, *_model0);        // outmodel is F'Ne
  outmodel->scale (-1.);        // outmodel is -F'Ne
  if (!_dampPerturb) {
    model->multiply_inverse_covariance ();      // model is M m
    outmodel->add (*model);     // outmodel is -F'Ne + Mm
  }
  delete data;
  delete model;
  if (DEBUG)
    printf ("Done with calculating b\n");
  return outmodel;
}

void
CgTrans::hessian (CgVector & x) const
{
  if (DEBUG)
    x.printme ("x for hessian");
  CgVector *data = _data0->clone_zero ();       // data is full of zeros
  CgVector *model = x.clone_zero ();    // model is full of zeros
  _cti->forward_linearized (*data, x, *_model0);        // data is Fx
  data->multiply_inverse_covariance (); // data is NFx
  _cti->adjoint_linearized (*data, *model, *_model0);   // model is F'NFx
  delete data;
  x.multiply_inverse_covariance ();     // x is Mx
  x.add (*model);               // x is (F'NF + M)x
  delete model;
  if (DEBUG)
    printf ("Done with hessian\n");
}

/** evaluate objective function for nonquadratic linear search*/
double
CgTrans::objective_function (const CgVector & x, int is_a_perturbation) const
{
  CgVector *model = x.clone (); // model is x or m+x
  CgVector *data = _data0->clone_zero ();       // data is full of zeros
  if (is_a_perturbation)
    model->add (*_model0);      // model is m+x
  model->apply_constraint ();   // model is constrained(m+x)
  _cti->forward_nonlinear (*data, *model);      // data is f(m+x)
  data->scale (-1);             // data is -f(m+x)
  data->add (*_data0);          // data is e = d - f(m+x)
  CgVector *datap = data->clone ();     // datap is e
  data->multiply_inverse_covariance (); // data is Ne
  double eNe = datap->dot (*data);      // eNe is e'Ne
  delete data;
  delete datap;
  if (_dampPerturb) {           // if requested, replace model by x
    model->scale (-1);          // model is -(m+x)
    model->add (*_model0);      // model is -x
    model->scale (-1);          // model is x (after constraint)
  }                             // model is m+x or x
  CgVector *modelp = model->clone ();   // modelp is m+x or x
  model->multiply_inverse_covariance ();        // model is M(m+x) or Mx
  double xMx = modelp->dot (*model);    // xMx is (m+x)'M(m+x) or x'Mx
  delete model;
  delete modelp;
  return (eNe + xMx);
}

int
CgTrans::adjoint_is_good () const
{
  if (CgQuad::is_debug ()) {
    if (DEBUG)
      printf ("Testing dot product\n");
    fflush (stdout);
  }
  CgVector *b = this->create_b ();      // use as an arbitrary non-zero model
  CgVector *d = _data0->clone ();       // use as arbitrary non-zero data
  CgVector::test(*b);
  CgVector::test(*d);
  CgVector *Fb = d->clone_zero ();      // for forward transform of b
  CgVector *Ad = b->clone_zero ();      // for adjoint of data
  _cti->forward_linearized (*Fb, *b, *_model0);
  _cti->adjoint_linearized (*d, *Ad, *_model0);
  double dFb = d->dot (*Fb);
  double Adb = Ad->dot (*b);
  if (b->dot(*b) > 0 && d->dot(*d) > 0) {
    assert (dFb != 0);
    assert (Adb != 0);
  }
  // see if equality is very good.
  int matches = almost_equal (dFb, Adb);        
  // allow rough equality (magnitudes most important)
  if (!matches) {               
    matches = ((dFb * Adb >= 0.) ? 1 : 0);      // assert have the same sign
    double r1 = (dFb > 0.) ? dFb : -dFb;        // abs
    double r2 = (Adb > 0.) ? Adb : -Adb;        // abs
    double oneplus = 1.2;       // may still fail for degenerate b and data
    matches = matches && (almost_cmp (oneplus * r1, r2) >= 0);
    matches = matches && (almost_cmp (oneplus * r2, r1) >= 0);
  }
  if (CgQuad::is_debug ()) {
    if (!matches) {
      printf ("\ndot(d,Fb)=%g dot(F'd,b)=%g\n", dFb, Adb);
      printf ("Failed dot-product test.\n");
    }
    fflush (stdout);
  }
  delete b;
  delete d;
  delete Fb;
  delete Ad;
  return matches;
}

// search scaled perturbation to minimize objective function
double
CgTrans::line_search (const CgVector & x, double alphamin,
                      double alphamax, double okerror,
                      double okfraction, long nter) const
{

  CgVector *xx = 0;
  long iter;
  double fnew, xmin, xerr, xnew, alphabest = 0., alphanew;

  // allocate space
  iter = 0;
  nter = (nter < 10) ? 10 : nter;
  fnew = 0.;
  line_sample (&xmin, &xerr, &xnew, &iter, fnew);
  while (iter <= nter &&
         (iter < 5 || xerr > okerror || xerr >= okfraction * xmin)) {
    alphanew = (alphamax - alphamin) * xnew + alphamin;
    xx = x.clone ();
    xx->scale (alphanew);
    fnew = this->objective_function (*xx, 1);
    delete xx;
    line_sample (&xmin, &xerr, &xnew, &iter, fnew);
    alphabest = (alphamax - alphamin) * xmin + alphamin;
  }
  // remove temp space */
  return alphabest;
}

/** search an interval (0,1) for a minimum, parabolas and golden section
  long iter = 0, nter = 20;
  float okerr = 0.001, factor = 1., xerr=1.;
  line_sample(&xmin,&xerr,&xnew,&iter,fnew);
  while (xerr>okerr*factor && iter<=nter) {
    fnew = function(xnew)
    line_sample(&xmin,&xerr,&xnew,&iter,fnew);
    factor=xmin  / * optional * /
  }
  No need to initialize xmin,xerr,xnew.  fnew not used
  in first iteration.  iter must be 0 in first iteration, and fnew
  must be updated for the returned xnew on later iterations.

  Originally patterned after a line search in 1st ed. of Numerical Recipes
  in 1988. This Ctran code could be cleaned up, but it has seen a lot of data
  and is well debugged already.  This is NOT REENTRANT.  Notice the
  internal static variables.
  */
void
CgTrans::line_sample (double *xmin, double *xerr, double *xnew, long *iter,
                      double fnew) const
{
  static double x[4] = { 0, 0, 0, 0 }, fx[4] = {
  0, 0, 0, 0}, dx[3] = {
  1., 1., 1.};
  static double r1, r2, rgold = 0;
  static int imin = 0, inew = 0, ix[4], ifx[4];
  int ip1, ip2, ip3;

  if (DEBUG) {
    printf ("beginning iteration %ld\n", *iter);
  }
  if (*iter <= 3) {
    if (DEBUG)
      printf ("first three iter\n");
    if (*iter == 0) {
      int i;
      if (DEBUG)
        printf ("first iter\n");
      // first iteration initialization
      imin = 0;
      inew = 0;
      rgold = (sqrt (5.) - 1.) * .5;
      // starting samping points
      x[0] = 0;
      x[1] = 1. - rgold;
      x[2] = rgold;
      x[3] = 1.;
      // dx are errors
      for (i = 0; i < 3; i++)
        dx[i] = 1.;
      // fx are function values at x's
      for (i = 0; i < 4; i++)
        fx[i] = 0.;
      // return
      goto done;
    }
    // 2nd and 3rd iteration, fill in fx
    // fill in function values
    fx[inew] = fnew;
    ++inew;
    // return
    goto done;
  }
  // four and more iterations, f's are full
  if (DEBUG)
    printf ("four and more\n");
  fx[inew] = fnew;              // store previously requested value of f
  // index by increasing x and by increasing fx
  index (ix, x, 4);
  index (ifx, fx, 4);
  if (DEBUG) {
    for (int i = 0; i < 4; i++) {
      printf ("x[%d]=%g fx[%d]=%g ", ix[i], x[ix[i]], ix[i], fx[ix[i]]);
      printf ("ix[%d]=%d ifx[%d]=%d\n", i, ix[i], i, ifx[i]);
    }
  }
  imin = ifx[0];                // minimum value of function
  inew = ix[0];                 // x(inew) should be replaced
  if (imin == ix[0] || imin == ix[1]) {
    inew = ix[3];
  }
  // set flag if smallest x was original minimum of interval
  int ltenth;
  ltenth = almost_zero (x[imin]);
  // if we have a smallest or largest x, then we do not have enough
  //   points for a parabolic method.  set flag for golden section
  int lgold;
  lgold = ((imin == ix[3]) || (imin == ix[0]) || ltenth) ? 1 : 0;
  // dx's are the successive errors, dx(2) last, dx(3) before last,
  //   dx(1) to be for this one
  dx[2] = dx[1];
  dx[1] = dx[0];
  // look at three points not to be replaced
  if (inew == ix[0] || inew == ix[1]) {
    ip1 = ix[1];
    ip2 = ix[2];
    ip3 = ix[3];
  } else {
    ip1 = ix[0];
    ip2 = ix[1];
    ip3 = ix[2];
  }
  // ip1, ip2, and ip3 are indices for three points to be saved, in
  //    ascending order of x
  // calculate distance between two points that span minimum = error
  dx[0] = x[ip3] - x[ip1];
  // use golden section if errors have not been decreasing rapidly enough
  //   (error must be less than half after two previous iterations)
  lgold = (lgold || dx[0] > dx[2] * .5) ? 1 : 0;
  // r1 is the span of left interval, r1 of the right.
  r1 = x[ip2] - x[ip1];
  r2 = x[ip3] - x[ip2];
  // go to golden section search, if flag set
  if (DEBUG)
    printf ("ip1=%d ip2=%d ip3=%d\n", ip1, ip2, ip3);
  if (!lgold) {
    if (DEBUG)
      printf ("Use parabolic method\n");
    // here use parabolic method, if have no previous objections
    int igood = parabola_min (&x[inew], x[ip1], fx[ip1], x[ip2],
                              fx[ip2], x[ip3], fx[ip3]);
    if (DEBUG) {
      printf ("x[ip1]=%g fx[ip1]=%g ", x[ip1], fx[ip1]);
      printf ("x[ip2]=%g fx[ip2]=%g ", x[ip2], fx[ip2]);
      printf ("x[ip3]=%g fx[ip3]=%g\n", x[ip3], fx[ip3]);
      printf ("parabolic x[inew]=%g\n", x[inew]);
    }
    /* if parabolic search was numerically stable, then return.
       (Straight lines have no interior minimum, for instance.)
       If the new point to be searched is already the minimum point,
       then might have already converged.  In case not, go on and
       do golden section to reduce the width of the interval searched. u
     */
    if (igood && x[inew] != x[imin]) {
      goto done;
    }
    // if parabolic search failed, do golden search
  }
  if (DEBUG)
    printf ("Use golden search\n");
  // Do golden search.  Divide larger section.
  int ll;
  ll = (r1 > r2) ? 1 : 0;       // ll==1 if left span is largest.
  if (ll)                       // find new x in larger left interval
    x[inew] = x[ip1] + rgold * (x[ip2] - x[ip1]);
  if (!ll)                      // find new x in larger right interval
    x[inew] = x[ip3] - rgold * (x[ip3] - x[ip2]);
  // if smallest x has the minimum, then scale down drastically
  //  we may have orders of magnitude to go.
  if (ltenth) {
    if (DEBUG)
      printf ("Scale by tenth\n");
    x[inew] = x[ip1] + (x[ip2] - x[ip1]) * .1;
  }
  // the following is the only valid exit from this routine
done:
  ++(*iter);                    // update iteration number
  *xmin = x[imin];              // note present best value of x
  *xerr = dx[0];                // note size of interval spanning minimum
  // give the value of x that needs to be searched, user will return fnew
  *xnew = x[inew];
  if (DEBUG) {
    printf ("imin=%d inew=%d\n", imin, inew);
    printf ("current minimum x is %g with f(x)=%g\n", x[imin], fx[imin]);
  }
  return;
}

// slothsort
void
CgTrans::index (int *ix, const double *rx, long n)
{
  assert (n < 20);
  long ixx, i, j;
  for (i = 0; i < n; i++) {
    ix[i] = i;
  }
  for (i = 0; i < n - 1; i++) {
    for (j = 0; j < n - 1; j++) {
      if (rx[ix[j]] > rx[ix[j + 1]]) {
        ixx = ix[j];
        ix[j] = ix[j + 1];
        ix[j + 1] = ixx;
      }
    }
  }
}

/** find minimum x of parabola with values f(x1) = f1, f(x2) = f2, f(x3) = f3
  f(xmin) <= f(x) for all x.
  return 1 if successful; return 0 and minimum of 0.5 if unsucessful
  x1 < x2 < x3; f(x2) < f(x1), f(x2) < f(x3)
  */
int
CgTrans::parabola_min (double *xmin, double x1, double f1,
                       double x2, double f2, double x3, double f3)
{
  double xm, a, b;

  if (almost_equal (x3, x1))
    goto middle;
  xm = (x2 - x1) / (x3 - x1);
  if (xm < 0.001 || xm > 0.999)
    goto middle;
  a = ((f3 - f1) * xm - (f2 - f1)) / (xm - xm * xm);
  b = f3 - f1 - a;
  if (a * b >= 0 || 0.5 * fabs (b) > fabs (a))
    goto middle;
  xm = -0.5 * b / a;
  if (xm < 0. || xm > 1.)
    goto middle;
  *xmin = xm * (x3 - x1) + x1;
  return 1;
middle:
  xm = 0.5;
  *xmin = xm * (x3 - x1) + x1;
  return 0;
}

int
CgTrans::almost_equal (double r1, double r2)
{
  return (!almost_cmp (r1, r2));
}

int
CgTrans::almost_cmp (double r1, double r2)
  /* return 1 if r1>r2, -1 if r1<r2, 0 if r1==r2, within 10*precision
   * always assume r1 == r2, if within precision
   * Test "if (almost_cmp(r1,r2)>0)" for a strict r1>r2
   */
{
  if (almost_zero (r1) && almost_zero (r2))
    return 0;
  if (r1 - (10. * FLT_EPSILON) * fabs (r1) >
      r2 + (10. * FLT_EPSILON) * fabs (r2))
    return 1;
  if (r1 + (10. * FLT_EPSILON) * fabs (r1) <
      r2 - (10. * FLT_EPSILON) * fabs (r2))
    return -1;
  return 0;
}

int
CgTrans::almost_zero (double r)
{
  if (fabs (r) < 10. * FLT_MIN)
    return 1;
  else
    return 0;
}

#ifdef TESTCG
/*
   Test method 1 to
   Minimize  0.5 x'Hx + b'x = 0, where H = |2 4 | and b = (2 1)'
                                           |4 11|
   solution is x = (-3 1)'
   Override default hessian and b, ignore data space and transform.
   The examples are explained in "Introduction to Applied Mathematics"
   by G. Strang.)
*/
class SampleModel:public CgVector
{
public:
  virtual ~ SampleModel ()
  {
    delete[]m;
  }
  SampleModel (float m0, float m1)
  {
    this->m = new float[2];
    this->m[0] = m0;
    this->m[1] = m1;
  }
  virtual CgVector *clone () const
  {
    SampleModel *out = new SampleModel (this->m[0], this->m[1]);
      return out;
  }
  virtual void scale (double factor)
  {
    this->m[0] *= factor;
    this->m[1] *= factor;
  }
  virtual void add (const CgVector & toBeAdded)
  {
    //const SampleModel *in = dynamic_cast<const SampleModel*>(&toBeAdded);
    const SampleModel *in = (const SampleModel *) &toBeAdded;
    this->m[0] += in->m[0];
    this->m[1] += in->m[1];
  }
  virtual double dot (const CgVector & toBeDotted) const
  {
    const SampleModel *in = (const SampleModel *) &toBeDotted;
      return (in->m[0] * this->m[0] + in->m[1] * this->m[1]);
  }
  float *m;
  virtual char *str ()
  {
    char *result = new char[80];
    // should be snprintf, but not supported on Windows(TM)
    sprintf (result, "[%g %g]", m[0], m[1]);
    return result;
  }
};

/**
ostream& operator<<(ostream&output, const SampleModel& x) {
  return output << "[" << x.m[0] << " " << x.m[1] << "]";
}
*/

class SampleQuad:public CgQuadInterface
{
public:
  virtual void hessian (CgVector & cm) const
  {
    float m0, m1;
    SampleModel *sm = (SampleModel *) & cm;
      m0 = 2. * sm->m[0] + 4. * sm->m[1];
      m1 = 4. * sm->m[0] + 11. * sm->m[1];
      sm->m[0] = m0;
      sm->m[1] = m1;
  }
  // This approximate inverse Hessian can be omitted entirely.
  virtual void inverse_hessian (CgVector & cm) const
  {
    SampleModel *sm = (SampleModel *) & cm;
    float m0, m1;
      m0 = 5. * sm->m[0] - 2. * sm->m[1];
      m1 = -2. * sm->m[0] + 1. * sm->m[1];
      sm->m[0] = m0;
      sm->m[1] = m1;
  }
};

/* Fit straight line to points (0,0) (1,8) (3,8) (4,20) */
class SampleData:public CgVector
{
public:
  virtual ~ SampleData ()
  {
    delete[]m;
  }
  SampleData (float m0, float m1, float m2, float m3)
  {
    this->m = new float[4];
    this->m[0] = m0;
    this->m[1] = m1;
    this->m[2] = m2;
    this->m[3] = m3;
  }
  virtual CgVector *clone () const
  {
    SampleData *out = new SampleData (this->m[0], this->m[1],
                                      this->m[2], this->m[3]);
      return out;
  }
  virtual void scale (double factor)
  {
    this->m[0] *= factor;
    this->m[1] *= factor;
    this->m[2] *= factor;
    this->m[3] *= factor;
  }
  virtual void add (const CgVector & toBeAdded)
  {
    //const SampleData *in = dynamic_cast<const SampleData*>(&toBeAdded);
    const SampleData *in = (const SampleData *) &toBeAdded;
    this->m[0] += in->m[0];
    this->m[1] += in->m[1];
    this->m[2] += in->m[2];
    this->m[3] += in->m[3];
  }
  virtual double dot (const CgVector & toBeDotted) const
  {
    const SampleData *in = (const SampleData *) &toBeDotted;
      return (in->m[0] * this->m[0] + in->m[1] * this->m[1] +
              in->m[2] * this->m[2] + in->m[3] * this->m[3]);
  }
  virtual void multiply_inverse_covariance ()
  {
    this->m[0] *= 10000.;
    this->m[1] *= 10000.;       // emphasize data error
    this->m[2] *= 10000.;
    this->m[3] *= 10000.;       // rather than model size
  }
  float *m;
};

class SampleTrans:public CgTransLinearInterface
{
  virtual void forward (CgVector & data, const CgVector & model) const
  {
    SampleData *d = (SampleData *) & data;
    const SampleModel *m = (const SampleModel *) &model;
      d->m[0] = m->m[0] + 0. * m->m[1];
      d->m[1] = m->m[0] + 1. * m->m[1];
      d->m[2] = m->m[0] + 3. * m->m[1];
      d->m[3] = m->m[0] + 4. * m->m[1];
  }
  virtual void adjoint (const CgVector & data, CgVector & model) const
  {
    const SampleData *d = (const SampleData *) &data;
    SampleModel *m = (SampleModel *) & model;
      m->m[0] = d->m[0] + d->m[1] + d->m[2] + d->m[3];
      m->m[1] = 0. * d->m[0] + 1. * d->m[1] + 3. * d->m[2] + 4. * d->m[3];
  }
};

int
main (int argc, char **argv)
{
  int verbose = (argc == 1) ? 1 : 0;
  SampleModel *x = 0;
  if (verbose) {
    fprintf (stdout, "---------------------------------------------------\n");
  }
  {                             // each of the following solutions is independent
    SampleQuad st;
    SampleModel b (2., 1.);
    CgQuad t (st, b);
    if (verbose) {
      fprintf (stdout, "x should converge to [-3 1 ] in two iterations\n");
    }
    for (int num_iterations = 3; num_iterations >= 0; num_iterations--) {
      x = (SampleModel *) t.solve (num_iterations);
      if (verbose)
        fprintf (stdout, "%d iteration(s) x=%s\n", num_iterations, x->str ());
      if (num_iterations > 1) {
        assert (Almost::equal (x->m[0], -3.));
        assert (Almost::equal (x->m[1], 1.));
      }
      delete x;
    }
  }
  {
    SampleTrans st;
    SampleData d (0., 8., 8., 20.);
    SampleModel m (0., 0.);
    CgTrans t (st, d, m);
    if (verbose)
      fprintf (stdout, "\nx should converge to [1 4] in two iterations "
               "(damp solution)\n");
    for (int num_iterations = 3; num_iterations >= 0; num_iterations--) {
      x = (SampleModel *) t.solve (num_iterations);
      if (verbose) {
        fprintf (stdout, "%d iteration(s) x=%s", num_iterations, x->str ());
        fprintf (stdout, " objective function = %g\n",
                 t.objective_function (*x, 0));
      }
      if (num_iterations > 1) {
        assert (x->m[0] > 1. && x->m[0] < 1.0001);
        assert (x->m[1] > 3.9999 && x->m[1] < 4.);
      }
      delete x;
    }
  }
  {
    SampleTrans st;
    SampleData d (0., 8., 8., 20.);
    SampleModel m (0.9, 3.9);
    CgTrans t (st, d, m, 0, 1); //don't add reference model, damp perturbations
    if (verbose) {
      fprintf (stdout, "\nx converges to [0.1 0.1] (damp perturbations)\n");
    }
    for (int num_iterations = 3; num_iterations >= 0; num_iterations--) {
      x = (SampleModel *) t.solve (num_iterations);
      if (verbose) {
        fprintf (stdout, "%d iteration(s) x=%s", num_iterations, x->str ());
        fprintf (stdout, " objective function = %g\n",
                 t.objective_function (*x, 1));
      }
      if (num_iterations > 1) {
        assert (x->m[0] >= 0.09999 && x->m[0] <= 0.1);
        assert (x->m[1] >= 0.1 && x->m[1] <= 0.10001);
      }
      delete x;
    }
  }
  {
    if (verbose) {
      fprintf (stdout,
               "\nLine search solution for [1 4], "
               "each with 1 conjugate grad\n");
    }
    for (int num_iterations = 4; num_iterations >= 0; num_iterations--) {
      SampleTrans st;
      SampleData d (0., 8., 8., 20.);
      SampleModel m (0., 0.);
      CgTrans t (st, d, m);
      x = (SampleModel *) t.solve (1, 10, num_iterations);
      if (verbose) {
        fprintf (stdout, "%d iteration(s) x=%s", num_iterations, x->str ());
        fprintf (stdout, " objective function = %g\n",
                 t.objective_function (*x, 0));
      }
      if (num_iterations > 2) {
        assert (x->m[0] > 0.999 && x->m[0] < 1.001);
        assert (x->m[1] > 3.999 && x->m[1] < 4.001);
      }
      delete x;
    }
  }
  return EXIT_SUCCESS;
}

/* Output
x should converge to [-3 1 ] in two iterations
3 iteration(s) x=[-3 1]
2 iteration(s) x=[-3 1]
1 iteration(s) x=[-2.97143 1.11429]
0 iteration(s) x=[0 0]

x should converge to [1 4] in two iterations (damp solution)
3 iteration(s) x=[1.00002 3.99998] objective function = 440017
2 iteration(s) x=[1.00002 3.99998] objective function = 440017
1 iteration(s) x=[1.25869 3.91591] objective function = 441052
0 iteration(s) x=[0 0] objective function = 5.28e+06

x converges to [0.1 0.1] (damp perturbations)
3 iteration(s) x=[0.0999957 0.100001] objective function = 440000
2 iteration(s) x=[0.0999957 0.100001] objective function = 440000
1 iteration(s) x=[0.0419805 0.118945] objective function = 440052
0 iteration(s) x=[-0 -0] objective function = 444600

Line search solution for [1 4], each with 1 conjugate grad
4 iteration(s) x=[1 3.99993] objective function = 440017
3 iteration(s) x=[1 3.99993] objective function = 440017
2 iteration(s) x=[0.999751 3.99915] objective function = 440017
1 iteration(s) x=[1.25869 3.91591] objective function = 441052
0 iteration(s) x=[0 0] objective function = 5.28e+06
*/

#endif
