/* 3D CLoud DATA merge */
#include <rsf.h>
#include "hash.h"

#define NCO 9         // x,y,z, nx,ny,nz, vx,vy,vz

int main(int argc, char *argv[])
{
  bool verb, isreal, norm;

  const char **fname; // array of file names
  int nFILES;         //     number of input files
  int nCLOUD;         //     number of point clouds
  int nopen;          // max number of open files

  /* I/O files */
  sf_file Fdwin = NULL; /* win  data */
  sf_file Fcwin = NULL; /* win cloud */
  sf_file Fdall = NULL; /* all  data */
  sf_file Fcall = NULL; /* all cloud */

  sf_axis aa,aw,af;
  int     ja,jw,jf;

  pt3d       * aco  = NULL; /* all cloud coordinates */
  pt3d       * wco  = NULL; /* win cloud coordinates */
  float      * allR = NULL; /* all data real */
  float      * winR = NULL; /* win data real */
  sf_complex * allC = NULL; /* all data complex */
  sf_complex * winC = NULL; /* win data complex */

  float      * fold = NULL;
  int ihash;

  float      * jnk = NULL;
  jnk = sf_floatalloc(NCO);

  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc,argv);
  if (!sf_getbool("verb", &verb)) verb = false; /* verbosity */
  if (!sf_getbool("norm", &norm)) norm = true;  /* fold normalization */

  /*------------------------------------------------------------*/
  nopen = sysconf(_SC_OPEN_MAX);
  if(verb) sf_warning("max open files: %d",nopen);

  // array of file names
  fname = (const char**) sf_alloc( (size_t) argc, sizeof(char*) );

  // count the clouds
  nFILES = 0;
  for(int i = 1; i < argc; i++) { /* collect inputs */
    if (NULL != strchr(argv[i],'=')) {
      continue; /* not a file */
    } else {
      fname[nFILES] = argv[i];
      nFILES++;
    }
  }
  if( nFILES == 0 ) sf_error("no input");
  nCLOUD = nFILES / 2; // number of point clouds

  if(verb)
    for(int jCLOUD = 0; jCLOUD < nCLOUD; jCLOUD++)
      sf_warning("%d %s %s",jCLOUD,fname[jCLOUD],fname[jCLOUD+nCLOUD]);

  /*------------------------------------------------------------*/
  Fdwin = sf_input(fname[0+nCLOUD]); // open win  data
  af = sf_iaxa(Fdwin, 2);            // get time/freq axis

  // check file type
  if (SF_FLOAT == sf_gettype(Fdwin)) isreal = true;
  else                               isreal = false;

  sf_fileclose(Fdwin);               // close win  data
  if(verb) sf_raxa(af);

  /*------------------------------------------------------------*/
  Fcall = sf_input ("all"); // open all cloud
  Fdall = sf_output("out"); // open all  data
  aa = sf_iaxa(Fcall, 2);

  // shape all data file
  sf_oaxa(Fdall, aa, 1);
  sf_oaxa(Fdall, af, 2);

  // set output type
  if(isreal) sf_settype (Fdall,SF_FLOAT);
  else       sf_settype (Fdall,SF_COMPLEX);

  // read all cloud
  pt3d o; // reference point
  aco = (pt3d *) sf_alloc(sf_n(aa), sizeof(*aco));
  for(ja = 0; ja < sf_n(aa); ja++) {
      sf_floatread( jnk,NCO,Fcall);
      aco[ja].x  = jnk[0]; o.x += aco[ja].x;
      aco[ja].y  = jnk[1]; o.y += aco[ja].y;
      aco[ja].z  = jnk[2]; o.z += aco[ja].z;
  }
  o.x /= sf_n(aa);
  o.y /= sf_n(aa);
  o.z /= sf_n(aa);

  // allocate all data storage
  if(isreal) allR = sf_floatalloc  ( sf_n(aa) );
  else       allC = sf_complexalloc( sf_n(aa) );

  if(norm) fold = sf_floatalloc( sf_n(aa) );

  /*------------------------------------------------------------*/
  /* define the hash table */
  unsigned int nhash = 2 * sf_n(aa);
  htInit( nhash );

  if(verb) fprintf(stderr,"Generate hash table: ");
  for(ja = 0; ja < sf_n(aa); ja++) {
      if(verb) fprintf(stderr,"%8d\b\b\b\b\b\b\b\b",sf_n(aa)-ja-1);
      htInsert( nhash, &aco[ja], &o, ja );
  }
  if(verb) fprintf(stderr,"\n");

  /*------------------------------------------------------------*/
  /* main loop */
  /*------------------------------------------------------------*/
  for(jf = 0; jf < sf_n(af); jf++) {
    if(verb) fprintf(stderr,"%8d\b\b\b\b\b\b\b\b",sf_n(af)-jf-1);

    // reset all data
    if(isreal) {
      #ifdef _OPENMP
      #pragma omp parallel for schedule(dynamic) \
        private(ja) shared(aa, allR, fold)
      #endif
      for(ja = 0; ja < sf_n(aa); ja++) {
        allR[ ja ] = 0.0;
        if(norm) fold[ ja ] = 0;
      }

    } else {
      #ifdef _OPENMP
      #pragma omp parallel for schedule(dynamic) \
        private(ja) shared(aa, allC, fold)
      #endif
      for(ja = 0; ja < sf_n(aa); ja++) {
        allC[ ja ] = 0.0;
        if(norm) fold[ ja ] = 0;
      }
    }

    // loop over clouds
    for(int jCLOUD = 0; jCLOUD < nCLOUD; jCLOUD++) {

      // open cloud/data files
      Fcwin = sf_input(fname[jCLOUD]       );
      Fdwin = sf_input(fname[jCLOUD+nCLOUD]);
      aw = sf_iaxa(Fcwin, 2);

      //  read win data
      if(isreal) {
        winR =   sf_floatalloc( sf_n(aw) );
        sf_floatread  (winR, sf_n(aw), Fdwin);
      } else {
        winC = sf_complexalloc( sf_n(aw) );
        sf_complexread(winC, sf_n(aw), Fdwin);
      }

      // read win cloud
      wco = (pt3d *) sf_alloc(sf_n(aw), sizeof(*wco));
      for(jw = 0; jw < sf_n(aw); jw++) {
          sf_floatread( jnk,NCO,Fcwin);
          wco[jw].x  = jnk[0];
          wco[jw].y  = jnk[1];
          wco[jw].z  = jnk[2];
      }

      // merge win clouds
      if(isreal) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) \
          private(jw,ihash) shared(aw, wco, winR, allR, o, fold)
        #endif
        for(jw = 0; jw < sf_n(aw); jw++) {
            ihash = htLookup( nhash, &wco[jw], &o);
            allR[ ihash ] += winR[ jw ];
            if(norm) fold[ ihash ]++;
        }
      } else {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) \
          private(jw,ihash) shared(aw, wco, winC, allC, o, fold)
        #endif
        for(jw = 0; jw < sf_n(aw); jw++) {
            ihash = htLookup( nhash, &wco[jw], &o);
            allC[ ihash ] += winC[ jw ];
            if(norm) fold[ ihash ]++;
        }
      }

      // clear win cloud
      free(wco);

      // clear win data
      if(isreal) free(winR);
      else       free(winC);

      // close cloud/data files
      sf_fileclose(Fcwin);
      sf_fileclose(Fdwin);

    } // end loop over clouds

    // normalize by fold
    if(norm) {
      for(ja = 0; ja < sf_n(aa); ja++) {
        fold[ja] = SF_MAX(fold[ja],1);
      }
      if(isreal) for(ja = 0; ja < sf_n(aa); ja++) { allR[ja] /= fold[ja]; }
      else       for(ja = 0; ja < sf_n(aa); ja++) { allC[ja] /= fold[ja]; }
    }

    // output all data
    if(isreal) sf_floatwrite  (allR, sf_n(aa), Fdall); // write all data
    else       sf_complexwrite(allC, sf_n(aa), Fdall);

  } // end loop over frequencies
  if(verb) fprintf(stderr,"\n");

  /*------------------------------------------------------------*/
  htClose();

  /*------------------------------------------------------------*/
  free(aco);
  if(isreal) free(allR);
  else       free(allC);
  sf_fileclose(Fcall);
  sf_fileclose(Fdall);

  if(norm) free(fold);
  free(jnk);
}
