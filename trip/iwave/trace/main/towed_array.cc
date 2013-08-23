#include <su.h>
#include <segy.h>
#include <header.h>
#include <cubic.h>
#include <parser.h>

#define DT_TOL 0.001

char * sdoc[] = { 
  "Usage: towed_array.x data= src= towed_src=",
  "",
  "Purpose: create towed streamer source traces, by transferring source",
  "coordinates to both source and receiver coordinates of source trace",
  "copies. That is, if a source trace has source and receiver coordinates",
  "src_sx, src_sy, src_selev, src_gx, src_gy, and src_gelev, and an input",
  "data trace has source coordinates data_sx, data_sy, and data_selev, then",
  "a trace with coordinates ",
  "  src_sx+data_sx,",
  "  src_sy+data_sy,",
  "  src_selev+data_selev,",
  "  src_gx+data_sx,",
  "  src_gy+data_sy,",
  "  src_gelev+data_selev,",
  "and the other headers and data samples of the source trace is written to",
  "the output towed source file. This algorithm produces a towed source data ",
  "file, representing the input source data translated by the source position",
  "read from the data file. Note that the array source option for the various",
  "IWAVE simulators expects to read the source position from the receiver",
  "coordinates of the source data file.",
  "",
  "Required parameters:",
  "  data           = filename for first modeled trace data",
  "  src            = source traces to be translated, source and receiver",
  "                   coordinates will be overwritten",
"  towed_src      = output file - input source data for IWAVE",
  "",
  "Optional Parameters:",
  "  tol [float]   = relative error threshhold for sample-by-sample",
  "                  comparison.",
  "",
  NULL};

int main(int argc, char ** argv) {
  
/* ARGUMENTS */
  char * data;     /* trace data file name */
  char * src;      /* source file name */
  char * towed;    /* towed source file name */

  /* INTERNAL VARIABLES */
  PARARRAY * par;  /* param array */
  FILE * fp1;      /* input file pointer */
  FILE * fp2;      /* input file pointer */
  FILE * fp3;      /* output file pointer */ 
  segy tr1;        /* input 1 trace workspace */
  segy * tr2;      /* input 2 trace workspace */
  int sx=0;        /* workspace */
  int sy=0;        /* workspace */
  int selev=0;     /* workspace */
  int tsx;         /* workspace */
  int tsy;         /* workspace */
  int tselev;      /* workspace */
  int tgx;         /* workspace */
  int tgy;         /* workspace */
  int tgelev;      /* workspace */
/*  int err;         global error flag */
  int itr;         /* trace counter */
  int jtr;         /* trace counter */
  int ktr;         /* trace counter */
  int ltr;         /* trace counter */
  int ntr;         /* trace counter */

  xargc=argc; xargv=argv;
  requestdoc(1);

  /* err=0; */

  /* extract input parameters */
  par=ps_new();
  if ( ps_createargs(par, argc - 1, argv + 1) ) {
    printf("Error parsing input data. ABORT.\n");
    exit(1);
  }
             
  if (ps_flcstring(*par,"data",&data)) {
    printf("TOWED_ARRAY reading data. ABORT.\n");
    exit(1);
  }

  if (ps_flcstring(*par,"src",&src)) {
    printf("TOWED_ARRAY reading src. ABORT.\n");
    exit(1);
  }

  if (ps_flcstring(*par,"towed",&towed)) {
    printf("TOWED_ARRAY reading towed. ABORT.\n");
    exit(1);
  }

  /* open data files */
  if (!(fp1=fopen(data,"r"))) {
    printf("COMP: failed to open 1st input file = %s. ABORT.\n",data);
    exit(1);
  }

  if (!(fp2=fopen(src,"r"))) {
    printf("COMP: failed to open 2nd input file = %s. ABORT.\n",src);
    exit(1);
  }

  if (!(fp3=fopen(towed,"w"))) {
    printf("COMP: failed to open output file = %s. ABORT.\n",towed);
    exit(1);
  }

  /* read loop */
  ntr=0;
  while (fgettr(fp2,&tr1)) ntr++;
  tr2=(segy *)malloc(ntr*sizeof(segy));
  fseek(fp2,0L,SEEK_SET);
  jtr=0;
  while (fgettr(fp2,&(tr2[jtr]))) jtr++;
  itr=0;
  ktr=0;
  ltr=0;
  while (fgettr(fp1,&tr1)) {
    if (itr==0) {
      sx=tr1.sx;
      sy=tr1.sy;
      selev=tr1.selev;
    }
    
    if ((itr==0) ||
	(sx != tr1.sx) ||
	(sy != tr1.sy) ||
	(selev != tr1.selev)) {
     
      for (jtr=0;jtr<ntr;jtr++) {
	tsx=tr2[jtr].sx;
	tsy=tr2[jtr].sy;
	tselev=tr2[jtr].selev;
	tgx=tr2[jtr].gx;
	tgy=tr2[jtr].gy;
	tgelev=tr2[jtr].gelev;
	tr2[jtr].sx += tr1.sx;
	tr2[jtr].sy += tr1.sy;
	tr2[jtr].selev += tr1.selev;
	tr2[jtr].gx += tr1.sx;
	tr2[jtr].gy += tr1.sy;
	tr2[jtr].gelev += tr1.selev;
	tr2[jtr].tracl = ktr;
	tr2[jtr].tracr = ktr;
	tr2[jtr].tracf = jtr;
	tr2[jtr].fldr  = ltr;
	fputtr(fp3,&(tr2[jtr]));
	tr2[jtr].sx = tsx;
	tr2[jtr].sy = tsy;
	tr2[jtr].selev = tselev;
	tr2[jtr].gx = tgx;
	tr2[jtr].gy = tgy;
	tr2[jtr].gelev = tgelev;
	ktr++;
      }
      ltr++;
      sx=tr1.sx;
      sy=tr1.sy;
      selev=tr1.selev;
    }
    itr++;
  }
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  ps_delete(&par);

  if (!itr) exit(1);
  exit(0);
}
  
