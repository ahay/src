/* 
Fourier domain cross-correlation
Paul Sava
Copyright (C) 2022 Colorado School of Mines
*/
#include <rsf.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

/*------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  bool verb;

  /* I/O files */
  sf_file Fa = NULL;      /*  input: x - f */
  sf_file Fb = NULL;      /* input:  x - f */
  sf_file Fc = NULL;      /* output: x - t */

  int           nt;
  float      ot,dt;
  sf_axis ax,af,at;
  int     jx,jf,jt;
  float       f, t;

  sf_complex ** aa = NULL;
  sf_complex ** bb = NULL;
  float      ** cc = NULL;

  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc,argv);

  /*------------------------------------------------------------*/
  /* OMP init */
  #ifdef _OPENMP
    omp_init();
  #endif

  /* default behavior */
  if (!sf_getbool("verb",&verb)) verb=false; /* verbosity  */

  if (!sf_getint  ("nt",&nt)) nt = 1;
  if (!sf_getfloat("ot",&ot)) ot = 0.0;
  if (!sf_getfloat("dt",&dt)) dt = 1.0;
  at = sf_maxa(nt,ot,dt);
  sf_setlabel(at,"t");
  sf_setunit(at,"s");


  /* setup i/o and auxiliary file */
  Fa = sf_input ( "in");
  Fb = sf_input ("ref");
  Fc = sf_output("out");

  ax = sf_iaxa(Fa,1);
  af = sf_iaxa(Fa,2);

  if(verb) sf_raxa(ax);
  if(verb) sf_raxa(af);
  if(verb) sf_raxa(at);

  sf_oaxa(Fc, ax, 1);
  sf_oaxa(Fc, at, 2);
  sf_settype(Fc,SF_FLOAT);

  /*------------------------------------------------------------*/
  /* allocate arrays */
  aa = sf_complexalloc2( sf_n(ax),sf_n(af) );
  bb = sf_complexalloc2( sf_n(ax),sf_n(af) );
  cc = sf_floatalloc2  ( sf_n(ax),sf_n(at) );

  /*------------------------------------------------------------*/
  sf_complexread(aa[0], sf_n(ax) * sf_n(af) , Fa);
  sf_complexread(bb[0], sf_n(ax) * sf_n(af) , Fb);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) \
    private(jt,jx,jf, t,f) \
    shared( at,ax,af, aa,bb,cc)
  #endif
  for(jt = 0; jt < sf_n(at); jt++) {
    t = sf_o(at) + jt * sf_d(at);

    for(jx = 0; jx < sf_n(ax); jx++) {
      cc[jt][jx] = 0;

      for(jf = 0; jf < sf_n(af); jf++) {
        f = sf_o(af) + jf * sf_d(af);

        cc[jt][jx] += aa[jf][jx] * conj(bb[jf][jx]) * cexpf(2*SF_PI*f *I*t);
      }
    }
  }

  sf_floatwrite (cc[0], sf_n(ax) * sf_n(at) , Fc);

  /*------------------------------------------------------------*/
  /* deallocate arrays */
  free(aa);
  free(bb);
  free(cc);

}
