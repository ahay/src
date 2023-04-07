/*
3D CLoud DATA merge
Paul Sava
Copyright (C) 2022 Colorado School of Mines
*/
#include <rsf.h>
#include "hash.h"

#define NCO 9         /* x,y,z, nx,ny,nz, vx,vy,vz */

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  bool verb, isreal, norm;

  bool allopen = false;
  const char **fname; /* array of file names */
  int nFILES;         /*     number of input files */
  int nCLOUD;         /*     number of point clouds */
  int nopen;          /* max number of open files */

  /* I/O files */
  sf_file * Fdwin = NULL; /* win  data */
  sf_file * Fcwin = NULL; /* win cloud */
  sf_file   Fdall = NULL; /* all  data */
  sf_file   Fcall = NULL; /* all cloud */

  sf_axis            aa,aw,af;
  unsigned long long ja,jw,jf;

  pt3d       * aco  = NULL; /* all cloud coordinates */
  pt3d       * wco  = NULL; /* win cloud coordinates */
  float      * allR = NULL; /* all data real */
  float      * winR = NULL; /* win data real */
  sf_complex * allC = NULL; /* all data complex */
  sf_complex * winC = NULL; /* win data complex */

  float      * fold = NULL;

  unsigned long long nhash;
  unsigned long long ihash;
  float hashscale;

  float      * jnk = NULL;
  jnk = sf_floatalloc(NCO);

  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc,argv);

  // default behavior
  if (!sf_getbool("verb", &verb)) verb = false; /* verbosity */
  if (!sf_getbool("norm", &norm)) norm = true;  /* fold normalization */

  /*------------------------------------------------------------*/
  nopen = sysconf(_SC_OPEN_MAX); /* max number of open files (system) */
  if(verb) sf_warning("max open files: %d",nopen);
  if(verb) sf_warning("total in files: %d",argc);
  if(argc < nopen) {
    allopen = true;
    if(verb) sf_warning("opening all files");
  }

  // array of file names
  fname = (const char**) sf_alloc( (size_t) argc, sizeof(char*) );

  // count the clouds
  nFILES = 0;
  for(int i = 1; i < argc; i++) {
    if (NULL != strchr(argv[i],'=')) {
      continue; /* not a file */
    } else {
      fname[nFILES] = argv[i];
      nFILES++;
    }
  }
  if( nFILES == 0 ) sf_error("no input");
  nCLOUD = nFILES / 2; // number of point clouds

  // list the clouds
  if(verb)
    for(int jCLOUD = 0; jCLOUD < nCLOUD; jCLOUD++)
      sf_warning(">>> %4d %s %s",jCLOUD,fname[jCLOUD],fname[jCLOUD+nCLOUD]);

  /*------------------------------------------------------------*/
  // allocate files
  if(allopen) {
    Fdwin = (sf_file*) sf_alloc( (size_t) nCLOUD, sizeof(sf_file) );
    Fcwin = (sf_file*) sf_alloc( (size_t) nCLOUD, sizeof(sf_file) );

    // open cloud/data files
    for(int jCLOUD = 0; jCLOUD < nCLOUD; jCLOUD++) {
      Fcwin[jCLOUD] = sf_input( fname[jCLOUD]        );
      Fdwin[jCLOUD] = sf_input( fname[jCLOUD+nCLOUD] );
    }

  } else {
    Fdwin = (sf_file*) sf_alloc( (size_t)    1, sizeof(sf_file) );
    Fcwin = (sf_file*) sf_alloc( (size_t)    1, sizeof(sf_file) );
  }

  /*------------------------------------------------------------*/
  // check file type
  if(!allopen) Fdwin[0] = sf_input(fname[0+nCLOUD]); // open win  data

  af = sf_iaxa(Fdwin[0], 2); // get time/freq axis
  if(verb) sf_raxa(af);

  // check file type
  if (SF_FLOAT == sf_gettype(Fdwin[0])) isreal = true;
  else                                  isreal = false;
  if(verb) sf_warning("file is real:%d",isreal);

  if(!allopen) sf_fileclose(Fdwin[0]);               // close win  data

  /*------------------------------------------------------------*/
  Fcall = sf_input ("all"); // open all cloud
  Fdall = sf_output("out"); // open all  data
  aa = sf_iaxa(Fcall, 2);

  if(    !sf_getfloat("hashscale",  &hashscale)) hashscale = 2.0;
  if(verb) sf_warning("hashscale=%f",hashscale);
  nhash = hashscale * sf_n(aa);
  htInit( nhash ); /* initialize the hash table */

  // shape all data file
  sf_oaxa(Fdall, aa, 1);
  sf_oaxa(Fdall, af, 2);

  // set output type
  if(isreal) sf_settype (Fdall,SF_FLOAT);
  else       sf_settype (Fdall,SF_COMPLEX);

  // read all cloud
  pt3d o; // reference point
  o.x = o.y = o.z = 0.0;
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

  /*------------------------------------------------------------*/
  /* find smallest point separation */
  double dd, ddMIN = DBL_MAX;
  for(ja = 1; ja < sf_n(aa); ja++) {

    dd = sqrtf(pow( aco[ja].x - aco[ja-1].x, 2) + 
               pow( aco[ja].y - aco[ja-1].y, 2) + 
               pow( aco[ja].z - aco[ja-1].z, 2) );
    ddMIN = SF_MIN(ddMIN, dd);
  }

  /*------------------------------------------------------------*/
  // allocate all data storage
  if(isreal) allR = sf_floatalloc  ( sf_n(aa) );
  else       allC = sf_complexalloc( sf_n(aa) );

  if(norm)   fold = sf_floatalloc  ( sf_n(aa) );

  /*------------------------------------------------------------*/
  /* define the hash table */
  if(verb) fprintf(stderr,"make hash table: ");
  for(ja = 0; ja < sf_n(aa); ja++) {
      if(verb && (sf_n(aa)-ja-1)%1000==0) fprintf(stderr,"%12llu\b\b\b\b\b\b\b\b\b\b\b\b",sf_n(aa)-ja-1);
      htInsert( nhash, &aco[ja], &o, ja );
  }
  if(verb) fprintf(stderr,"            \n");


  /*------------------------------------------------------------*/
  /* main loop */
  /*------------------------------------------------------------*/
  for(jf = 0; jf < sf_n(af); jf++) {
    if(verb) fprintf(stderr,"%12llu\b\b\b\b\b\b\b\b\b\b\b\b",sf_n(af)-jf-1);

    // reset out data & fold
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) \
      private(ja) shared(aa, allR,allC, fold, isreal)
    #endif
    for(ja = 0; ja < sf_n(aa); ja++) {
      if(isreal) allR[ ja ] = 0.0;
      else       allC[ ja ] = 0.0;
      if(norm && jf==0) fold[ ja ] = 0;
    }

    // loop over clouds
    for(int jCLOUD = 0; jCLOUD < nCLOUD; jCLOUD++) {

      // open cloud/data files
      if( ! allopen ) {
        Fcwin[0] = sf_input( fname[jCLOUD]          );
        Fdwin[0] = sf_input( fname[jCLOUD + nCLOUD] );
      }

      // get the size of the cloud
      if( allopen ) aw = sf_iaxa(Fcwin[jCLOUD], 2);
      else          aw = sf_iaxa(Fcwin[     0], 2);

      // skip empty clouds
      if(sf_n(aw) == 0) continue;

      // read win cloud
      wco = (pt3d *) sf_alloc(sf_n(aw), sizeof(*wco));
      if( allopen ) sf_seek(Fcwin[jCLOUD],0,SEEK_SET);
      for(jw = 0; jw < sf_n(aw); jw++) {
        if( allopen ) sf_floatread( jnk,NCO,Fcwin[jCLOUD]);
        else          sf_floatread( jnk,NCO,Fcwin[     0]);
        wco[jw].x = jnk[0];
        wco[jw].y = jnk[1];
        wco[jw].z = jnk[2];
      }

      // read win data
      if(isreal) {
        winR = sf_floatalloc( sf_n(aw) );
        if( allopen ){
          sf_floatread  (winR, sf_n(aw), Fdwin[jCLOUD]);
        } else {
          sf_seek(Fdwin[0], jf*sf_n(aw)*sizeof(float),SEEK_SET);
          sf_floatread  (winR, sf_n(aw), Fdwin[0]);
        }
      } else {
        winC = sf_complexalloc( sf_n(aw) );
        if( allopen ){
          sf_complexread(winC, sf_n(aw), Fdwin[jCLOUD]);
        } else {
          sf_seek(Fdwin[0], jf*sf_n(aw)*sizeof(sf_complex),SEEK_SET);
          sf_complexread(winC, sf_n(aw), Fdwin[0]);
        }
      }

      // merge win clouds
      #ifdef _OPENMP
      #pragma omp parallel for schedule(dynamic) \
        private(jw,ihash) \
        shared(aw, wco, winR,winC, allR,allC, o, fold,isreal, ddMIN)
      #endif
      for(jw = 0; jw < sf_n(aw); jw++) {
          ja = htLookup( nhash, &wco[jw], &o, ddMIN);
          if(isreal) allR[ ja ] += winR[ jw ];
          else       allC[ ja ] += winC[ jw ];
          if(norm && jf==0) fold[ ja ]++;
      }

      free(wco);             // clear win cloud
      if(isreal) free(winR); // clear win data
      else       free(winC);

      if( ! allopen ) {
        sf_fileclose(Fcwin[0]); // close cloud files
        sf_fileclose(Fdwin[0]); // close  data files
      }

    } // end loop over clouds

    // normalize by fold
    if(norm) {
      for(ja = 0; ja < sf_n(aa); ja++) {
        fold[ ja ] = SF_MAX(fold[ja],1);
      }
      if(isreal) for(ja = 0; ja < sf_n(aa); ja++) { allR[ja] /= fold[ja]; }
      else       for(ja = 0; ja < sf_n(aa); ja++) { allC[ja] /= fold[ja]; }
    }

    // output all data
    if(isreal) sf_floatwrite  (allR, sf_n(aa), Fdall); // write all data
    else       sf_complexwrite(allC, sf_n(aa), Fdall);

  } // end loop over frequencies
  if(verb) fprintf(stderr,"            \n");

  /*------------------------------------------------------------*/
  htClose();

  /*------------------------------------------------------------*/
  // close cloud/data files
  if(allopen) {
    for(int jCLOUD = 0; jCLOUD < nCLOUD; jCLOUD++) {
      sf_fileclose( Fcwin[jCLOUD] );
      sf_fileclose( Fdwin[jCLOUD] );
    }
  }
  free(Fcwin);
  free(Fdwin);

  sf_fileclose(Fcall);
  sf_fileclose(Fdall);

  free(aco);
  if(isreal) free(allR);
  else       free(allC);

  if(norm) free(fold);
  free(jnk);
}
