#include <rsf.h>
#include <rsfplot.h>

int main(int argc, char* argv[])
{
    int n1, n2, n3, gainstep, panel;
    float o1, o2, o3, d1, d2, d3, tpow, gpow, clip, pclip, phalf, bias;
    bool transp, yreverse, gain=false, allpos, coltab, polarity;
    char *gainpanel, *color;
    sf_file in;

    sf_init(argc,argv);
    in = sf_input("in");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;
    if (!sf_histfloat(in,"o3",&o3)) o3=0.;

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_histfloat(in,"d3",&d3)) d3=1.;
    
    if (!sf_getbool("transp",&transp)) transp=true;
    if (!sf_getbool("yreverse",&yreverse)) yreverse=true;

    if (!sf_getfloat("tpow",&tpow)) tpow=0.;
    if (!sf_getfloat("gpow",&gpow)) {
	gpow=1.;
    } else if (gpow <= 0.) {
	gpow=0.;
	if (!sf_getfloat("phalf",&phalf)) {
	    phalf=85.;
	} else if (phalf <=0. || phalf > 100.) {
	    sf_error("phalf=%g should be > 0 and <= 100",phalf);
	}
	gain = true;
    }

    if (!sf_getfloat("clip",&clip)) {
	if (!sf_getfloat("pclip",&pclip)) {
	    pclip=99.;
	} else if (pclip <=0. || pclip > 100.) {
	    sf_error("pclip=%g should be > 0 and <= 100",pclip);
	}
	gain = true;
    }

    if (gain) {
	if (!sf_getint("gainstep",&gainstep)) gainstep=0.5+n1/256.;
	gainpanel = sf_getstring("gainpanel");
	switch (gainpanel[0]) {
	    case 'a': /* all */
		panel=-1;
		break;
	    case 'e': /* each */
		panel=1;
		break;
	    default:
		if (!sf_getint("gainpanel",&panel)) 
		    sf_error("Wrong gainpanel=");
		break;
	}
	free (gainpanel);
    }

    if (!sf_getbool("allpos",&allpos)) allpos=false;
    if (!sf_getfloat("bias",&bias)) bias=0.; 
    if (!sf_getbool("polarity",&polarity)) polarity=false;

    if (!sf_getbool("coltab",&coltab)) coltab=true;
    if (NULL == (color = sf_getstring("color"))) color="I";

    vp_stdplot_init (o1-0.5*d1, o1+(n1-1)*d1+0.5*d1, 
		     o2-0.5*d2, o2+(n2-1)*d2+0.5*d2,
		     transp, false, yreverse, false);

    exit (0);
}

#ifdef jhvkjhvb

if(datain. esize==4) make_unpipe("in");

if (!getch ("nreserve", "d", &nreserve)) nreserve = 8;

   numorient = getch ("orient", "d", &orient);
   numinvert = getch ("invert", "1", &invert);
   blast = 1;
   getch ("hurry", "1", &blast);
	return(0);
}

/*------------------------------------------------------------------------*/
/*----------    TIME TO DO SOME INITIALIZING STUFF    --------------------*/
/*------------------------------------------------------------------------*/

  if(datain.esize==4){
    data = (float*) alloc( datain.n1[0]* datain.n2* sizeof(float) );
    tgain = (float*) alloc( datain.n1[0]*sizeof(float) );
    buf = ( char*) alloc( nslow*nfast*sizeof(float) );
    data2 = ( unsigned char*) alloc( nslow*nfast*sizeof(float) );
    tbl = ( char*) alloc( TSIZE +1); te=tbl;
	
    ierr=init_esize_4(data,tgain,tbl,buf,garg,te);
  }
  else {
		buf = ( char *) alloc (datain.n1[0] * datain.n2 * datain.esize);
		data2 = ( unsigned  char *) alloc (datain.n1[0] * datain.n2 * datain.esize);
	}

  if(dataout.esize==1) ierr=rite_eout_1_head ();
  else  ierr=rite_vplot_head();

  if(dataout.esize !=1){  
    if(wantscalebar) bdata = (unsigned char *) malloc (datain.esize * 256);
    ierr=init_vplot_out();
  }





/*------------------------------------------------------------------------*/
/*-------------    LOOP OVER PLANES DOING WORK    ------------------------*/
/*------------------------------------------------------------------------*/


  for (i3=0; i3<datain.n3; i3++) {

  /* convert from esize =4 to esize =1 */
    if(datain.esize==4) ierr=convert_4_to_1(data,tgain,tbl,buf,i3,te);
    else ierr=read_in_byte(buf,i3);
    
    if(dataout.esize ==1) ierr=rite_byte_frame(buf);
    else ierr=rite_vplot_frame(buf,data2,bdata,i3);
    }  
    if(dataout.esize==1) ierr=finish_byte_plot();
    else ierr=finish_vplot_plot();
	return(0);
}














/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*----------------  Begining of old Taplot subroutines   --------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*  Function initializes the conversion table, the gain vector,
  and sets fastinc and slowinc.
*/
#if defined(__STDC__) || defined(__stdc__)
tbinit (float gpow,float tpow,float clp,register char *tbl,float *tgain,
float t0,int nt,int *fastinc,int *slowinc,int nfast,float dt)
#else
tbinit (gpow,tpow,clp,tbl,tgain,t0,nt,fastinc,slowinc,nfast,dt)
register char *tbl;
float gpow,tpow,clp,*tgain,t0,dt;
int nt,*fastinc,*slowinc,nfast;
#endif
{
int it;
double t;

/* initialize the conversion table */
if ( gpow != 1.)
  {
  if      ( allpos == NO ) 
  for (it=1; it<=TSIZE/2; it++)
    {
    tbl[TSIZE-it] =
    (252*(pow((double)((TSIZE-2.*it)/TSIZE),gpow)+1.))/2.+3.;
    tbl[it] = 255 - tbl[TSIZE-it] + 2.;
    }
  else if ( allpos == YES )
  for( it=1; it < TSIZE ; it++ )
    tbl[it] = 256*( (double) (it-1)/TSIZE );
  }
else
  {
  if      ( allpos == NO ) 
  for (it=1; it<=TSIZE/2; it++)

    {
    tbl[TSIZE-it] =
    (252*((double)((TSIZE-2.*it)/TSIZE)+1.))/2.+3.;
    tbl[it] = 255 - tbl[TSIZE-it] + 2.;
    }
  else if ( allpos == YES ){
  for( it=1; it < TSIZE ; it++ ){
    tbl[it] = 256*( (double) (it-1)/TSIZE );
	}
	}
  }
tbl[0] = tbl[1];
tbl[TSIZE] = tbl[TSIZE-1];

/* initialize the gain vector */
if ( tpow != 0.)
  {
  if( allpos == NO )
  for (it=0; it<nt; it++)
    {
    t=(it+1)*dt+t0;
  tgain[it] = ((t == 0.0) && (tpow == 0.0)) ? 1.0 : pow (fabs(t),tpow);
    tgain[it] *= TSIZE/(2.*clp);
    }

  else if ( allpos == YES )
       for (it=0; it<nt; it++)
              {
              t=(it+1)*dt+t0;
       tgain[it] = ((t == 0.0) && (tpow == 0.0)) ? 1.0 : pow (fabs(t),tpow);
              tgain[it] *= TSIZE/(clp);
              }
  } 
else {
  if( allpos == NO ) tgain[0] = TSIZE/(2.*clp);
  else         tgain[0] = TSIZE/(clp);
      }

/*Put here to duplicate Taplot behavior for backward compatibility*/
/*if (*transp!='y' && dataout.esize==1)*/
/*if (*transp=='y')*/
 /*the only time the fastinc will not be 1 is when
   we are doing old sytle Taplot behavior which did
   a transpose on the fly */
if (*transp=='y' && dataout.esize==1 )
  {
  *fastinc = nfast;
  *slowinc = 1;
  }
else
  {
  *fastinc = 1;
  *slowinc = nfast;
  }
		return(0);
}


#if defined(__STDC__) || defined(__stdc__)
int cnvrt(int i3,float bias,float qclp,int np,int nx,int nt,float *data,char 
*buf,register char *tbl,int sfil,int nfast,int nslow,int slowinc,int fastinc,char *te,
float *tgain,float pbias,float tpow)
#else
int cnvrt(i3,bias,qclp,np,nx,nt,data,buf,tbl,sfil,nfast,nslow,slowinc,fastinc,te,
tgain,pbias,tpow)
int i3,np,nx,nt,nfast,nslow,sfil,slowinc,fastinc; 
float qclp,*data,*tgain,pbias,tpow; /*register*/ float bias;
 char *buf; 
register char *tbl;
char *te;
#endif

{
int ix,it; float *dp,*de; register float *sp,ierr;
register char *bp; char *bs;
register char *tp;
float tmp;
char *min,*max;


 if(i3 != 0 || qclp == 0.) sreed("in",data,nt*nx*4);
/*	minval[i3]=*data; maxval[i3]=*data;*/

/* convert each panel */
min=te;
max=tbl;

if (tpow != 0.)
  {
  for (ix=0, bs=buf; ix<nx; ix++,bs += slowinc)
    {
    sp = tgain; bp = bs;
    for (dp=data+ix*nt, de=dp+nt; dp<de; dp++,sp++,bp += fastinc)
      {
      tp = tbl + (int)((*dp - pbias)*(*sp)+bias);
      if (tp>te) tp = te;
      else if (tp<tbl) tp = tbl;
      *bp = *tp;
      }
    }
   } 
else 
  {
  tmp = tgain[0];
  for (ix=0, bs=buf; ix<nx; ix++,bs += slowinc)
    {
    bp = bs;
    for (dp=data+ix*nt, de=dp+nt; dp<de; dp++,bp += fastinc)
      {
      tp = tbl + (int)((*dp - pbias)*tmp+bias);
      if (tp>te) tp = te;
      else if (tp<tbl) tp = tbl;
      *bp = *tp;
			if(tp < min) min=tp;
			else if(tp > max) max=tp;
      }
    }
  } 
	minval[i3]=((float)(min-tbl)-bias)/tmp + pbias;
	maxval[i3]=((float)(max-tbl)-bias)/tmp + pbias;

		return(0);
}

#if defined(__STDC__) || defined(__stdc__)
int rite_eout_1_head ()
#else
int rite_eout_1_head ()
#endif
{
  char label1[128], label2[128], label3[128];
  float r; int i;
  struct dimensions {
    int ni;   /* input dimensions */
    float oi;
    float di;
    char labeli[100];
    } dm[4],*p0,*p1,*p2,*p3;
  for(i=0; i<=3; i++) {
    dm[i].ni = 1; dm[i].oi = 0.; dm[i].di = 1.;
    sprintf(dm[i].labeli,"");
    }

  dm[1].ni= datain.n1[0]; dm[2].ni= datain.n2; dm[3].ni= datain.n3;
  dm[1].oi= datain.o1[0]; dm[2].oi= datain.o2; dm[3].oi= datain.o3;
  dm[1].di= datain.d1[0]; dm[2].di= datain.d2; dm[3].di= datain.d3;

  fetch("label1","s",dm[1].labeli);
  fetch("label2","s",dm[2].labeli);
  fetch("label3","s",dm[3].labeli);

  p0 = dm; p1 = dm + 1; p2 = dm + 2; p3 = dm + 3;
  if (*transp=='y') { p0 = p1 ; p1 = p2 ; p2 = p0; }
  
  puthead ("  transp=%s \n",transp);
  puthead ("  gainpanel=%s   gainstep=%d\n",gainpanel,gainstep);
  puthead ("  n1=%d   n2=%d   \n",datain.n1[0],datain.n2);
  puthead ("    tpow=%g   bias=%g\n",tpow,pbias);
  puthead ("    clip=%g   pclip=%f\n",clp,qclp);
  puthead ("    gpow=%g   phalf=%g\n",gpow,phalf);
  puthead("  esize=1   n1=%d   n2=%d   \n",p1->ni,p2->ni);
  puthead("  o1=%g   o2=%g   o3=%g\n",p1->oi,p2->oi,p3->oi);
  puthead("  d1=%g   d2=%g   d3=%g\n",p1->di,p2->di,p3->di);
  puthead("  label1=\"%s\"  label2=\"%s\"  label3=\"%s\"\n",
    p1->labeli,p2->labeli,p3->labeli);
  puthead("        maxval=%f \n", maxval[1]);
  puthead("        minval=%f \n", minval[1]);
  hclose();
		return(0);
}

#if defined(__STDC__) || defined(__stdc__)
int init_esize_4(float *data,float *tgain, register char *tbl,char *buf,int garg,char *te)
#else
int init_esize_4(data,tgain,tbl,buf,garg,te)
float *data,*tgain;
register char *tbl;
 char *buf;
char *te;
int garg;
#endif
{
int it;

    /* gain and clip parameters, initialize the gain table, tgain, */
    /*  fastinc and slowinc */
    if (garg) {
      /*   gainpar is called if user does not specify clip or gpow. */
  
  
  sgainpar ("in",data,datain.n1,&hbytes,&datain.n2,&gainstep,&tpow,
            datain.o1,&qclp,&phalf,&clp,&gpow,&pbias,datain.d1,&datain.n3,&gainip);

  if ((gainpanel[0] != 'e') && (clp == 0)) {
      if (clp==0) {
      }
  }
  if ((gainpanel[0] == 'e') && (clp == 0)) { clp = 1.e-10; }

    }
    sseek("in",0,0);
    sreed("in",(char *)data,datain.n1[0]*datain.n2*4);
    te = tbl + TSIZE;
    tbinit (gpow,tpow,clp,tbl,tgain,datain.d1[0],datain.n1[0],&fastinc,&slowinc,nfast,datain.d1[0]);
    if( allpos == NO )  bias = TSIZE/2.;
    else          bias = 0;

    /* print parameters */
    
    /* set the correct format in the output */
    set_format("out","native_byte");
/*		putch("data_format","s","native_byte");*/

    /* Compute maximum and minimum of data (first panel only) */
    max = min = data[0];
    for (it=1; it<datain.n1[0]*datain.n2; it++) {
  if(max < data[it]) max=data[it];
  if(min > data[it]) min=data[it];
    }
    /* Compute maximum and minimum value that will be displayed on plot */
    /* Values are used by Ta2vplot, when wantscalebar option is on */
  minval[0]=maxval[0]=0.;
    if(max <=0. && min < 0.) {
  /* Case all data values are negative */
  maxval[1] = MIN (max, clp + pbias);
  minval[1] = MAX (-clp + pbias, min);
    } else if(min >=0. && max > 0.) {
  /* Case all data values are positive */
  maxval[1] = MIN (max, clp + pbias);
  minval[1] = MAX (pbias, min);
    } else {
  /* Both positive and negative data values */
  maxval[1] = MIN (max, clp + pbias);
  minval[1] = MAX (-clp + pbias, min);
    }
/*		minval[1]=-clp+pbias; maxval[1]=clp+pbias;*/
    /* done processing header */
		return(0);
}

#if defined(__STDC__) || defined(__stdc__)
int convert_4_to_1(float *data,float *tgain,register char *tbl,char *buf,int i3,char *tm)
#else
int convert_4_to_1(data,tgain,tbl,buf,i3,tm)
float *data,*tgain;
  register char *tbl;
   char *buf;
  int i3;
  char *tm;
#endif
{
/*  register char *te;*/
  char *te;
      te = tbl + TSIZE;
  if (i3!=0 && gainpanel[0] == 'e') {
      if (garg) {
    gpow *= igpow;
    clp *= iclip;
    gainip = i3 + 1;
    sseek("in",0,0);

    sgainpar ("in",data,datain.n1,&hbytes,&datain.n2,&gainstep,&tpow,
          datain.o1,&qclp,&phalf,&clp,&gpow,&pbias,datain.d1,&datain.n3,&gainip);


    sseek("in",i3*(datain.n2*hbytes+datain.n1[0]*datain.n2*4),0);
    if (clp == 0.) {  clp = 1.e-10;  }
      }
      te = tbl + TSIZE;
      tbinit (gpow,tpow,clp,tbl,tgain,datain.o1[0],datain.n1[0],&fastinc,&slowinc,nfast,datain.d1[0]);
      if( allpos == NO ) bias = TSIZE/2.;
      else         bias = 0;
  }



  /* convert panel */
  cnvrt(i3,bias,qclp,datain.n3,datain.n2,datain.n1[0],data,buf,tbl,outfd,                                      nfast,nslow,slowinc,fastinc,te,tgain,pbias,tpow);
  
		return(0);
}


/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*-------------------  BEGINING OF OLD TA2VPLOT PARMS  -------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/



#if defined(__STDC__) || defined(__stdc__)
int get_to_vplot_params()
#else
int get_to_vplot_params()
#endif


#if defined(__STDC__) || defined(__stdc__)
int rite_vplot_head()
#else
int rite_vplot_head()
#endif
{
  putch ("movish", "d", &movish);
  puthead ("\tn1=-1\n");
  putch ("polarity", "d", &polarity);
  set_output_data_format("vplot");
  hclose();
	return(0);
}


#if defined(__STDC__) || defined(__stdc__)
void setcoordinate (int orient, int invert, struct coordinfo *coordinate)
#else
void setcoordinate (orient, invert, coordinate)
  int             orient;
  int             invert;
  struct coordinfo *coordinate;
#endif
{
  if (orient == 0 && invert == 1) {
    coordinate->transp = 0;
    coordinate->yreverse = 0;
    coordinate->xreverse = 0;
  }
  if (orient == 1 && invert == 1) {
    coordinate->transp = 1;
    coordinate->yreverse = 1;
    coordinate->xreverse = 0;
  }
  if (orient == 2 && invert == 1) {
    coordinate->transp = 0;
    coordinate->yreverse = 1;
    coordinate->xreverse = 1;
  }
  if (orient == 3 && invert == 1) {
    coordinate->transp = 1;
    coordinate->yreverse = 0;
    coordinate->xreverse = 1;
  }
  if (orient == 0 && invert == 0) {
    coordinate->transp = 0;
    coordinate->yreverse = 1;
    coordinate->xreverse = 0;
  }
  if (orient == 1 && invert == 0) {
    coordinate->transp = 1;
     coordinate->yreverse = 1;
     coordinate->xreverse = 1;
  }
  if (orient == 2 && invert == 0) {
    coordinate->transp = 0;
    coordinate->yreverse = 0;
    coordinate->xreverse = 1;
  }
  if (orient == 3 && invert == 0) {
    coordinate->transp = 1;
    coordinate->yreverse = 0;
    coordinate->xreverse = 0;
  }

}


#if defined(__STDC__) || defined(__stdc__)
void setorient (int *orient, int *invert, struct coordinfo coordinate)
#else
void setorient (orient, invert, coordinate)
  int            *orient;
  int            *invert;
  struct coordinfo coordinate;
#endif
  {
  if (coordinate.transp == 0 && coordinate.yreverse == 0 & 
    coordinate.xreverse== 0) {
    *orient = 0; *invert = 1;
  }
  if (coordinate.transp && coordinate.yreverse && coordinate.xreverse == 0) {
    *orient = 1; *invert = 1;
  }
  if (coordinate.transp == 0 && coordinate.yreverse && coordinate.xreverse) {
    *orient = 2; *invert = 1;
  }
  if (coordinate.transp && coordinate.yreverse == 0 && coordinate.xreverse) {
    *orient = 3; *invert = 1;
  }
  if (coordinate.transp == 0 && coordinate.xreverse == 0 && 
    coordinate.yreverse) {
    *orient = 0; *invert = 0;
  }
  if (coordinate.transp && coordinate.yreverse && coordinate.xreverse) {
    *orient = 1; *invert = 0;
  }
  if (coordinate.transp == 0 && coordinate.xreverse && 
    coordinate.yreverse == 0) {
    *orient = 2; *invert = 0;
  }
  if (coordinate.transp && coordinate.yreverse == 0 &&coordinate.xreverse== 0){
    *orient = 3; *invert = 0;
  }

}

#if defined(__STDC__) || defined(__stdc__)
int init_vplot_out()
#else
int init_vplot_out()
#endif
{
  int i3;

  if(wantscalebar) gl_barint (&position, &axis1, &barposit, &baraxis, minval, 
    maxval, bartype, &barreverse, dataout.n3,0);

  if (numorient || numinvert) setcoordinate (orient, invert, &coordinate);
  else setorient (&orient, &invert, coordinate);
 vp_erase ();
 vp_color (axis1.col[0]);
 if (coltab) {
  if (color[0] >= '0' && color[0] <= '9') {
    redbit = color[0] - '0';
    greenbit = color[1] - '0';
    bluebit = color[2] - '0';
    if (redbit + greenbit + bluebit != 8)
      seperr ("You must use exactly 8 bits!\n");

    redoff = 0;
    greenoff = redbit;
    blueoff = redbit + greenbit;

    for (i3 = 0; i3 < 256; i3++) {
      ii = ~(~0 << redbit);
      if (ii > 0) red[i3] = (float) ((i3 >> redoff) & ii) / (float) (ii);
      else red[i3] = 0.;
      ii = ~(~0 << greenbit);

       if (ii > 0) green[i3] = (float) ((i3 >> greenoff) & ii) / (float) (ii);
      else green[i3] = 0.;

      ii = ~(~0 << bluebit);
      if (ii > 0)
         blue[i3] = (float) ((i3 >> blueoff) & ii) / (float) (ii);
      else blue[i3] = 0.;
    }
    for (jj = 0; jj < 256; jj++) {
      ii = 0;
      greenbit2 = greenbit;
      bluebit2 = bluebit;
      redbit2 = redbit;
      kk = 0;
      while (kk < 8) {
        greenbit2--;
        if (greenbit2 >= 0) {
          if (jj & (1 << (greenbit2 + greenoff))) ii |= 1 << kk;
          kk++;
        }
        redbit2--;
        if (redbit2 >= 0) {
          if (jj & (1 << (redbit2 + redoff))) ii |= 1 << kk;
          kk++;
        }
        bluebit2--;
        if (bluebit2 >= 0) {
          if (jj & (1 << (bluebit2 + blueoff))) ii |= 1 << kk;
          kk++;
        }
      }
      map[ii] = jj;
    }
    for (i3 = nreserve; i3 < 256; i3++) {
      jj = i3 - nreserve;
      vp_coltab (i3, red[map[jj]], green[map[jj]], blue[map[jj]]);
    }
  }
  else { 
		vp_rascoltab (nreserve, color);}
 }

 /* Set the coordinate transformation */
 gl_vplotint (&position, &coordinate, &axis1, &axis2);
 gl_plotpram (&colorin, &coordinate);
	 multi_t=getch("titles","s",titles);
	 if(0==fetch("title","s",title_temp)) sprintf(title_temp,"%s"," ");
/* fastplt = fastplot ();*/
  fastplt=0;
 return(0);
}



int read_in_byte(data,i3)
char *data;
int i3;
{
if(datain.n1[0] * datain.esize * datain.n2 != sreed ("in", data, 
  datain.n1[0] * datain.n2 * datain.esize)) seperr("trouble reading in data \n");

return(0);
}






#if defined(__STDC__) || defined(__stdc__)
int rite_vplot_frame(char *data,unsigned char *data2,unsigned char *bdata,int i3)
#else
int rite_vplot_frame(data,data2,bdata,i3)
char *data;
unsigned char *bdata,*data2;
int i3;
#endif
{
int esize,i1;
char title_out[128],temp;
int bad=0;

if(strcmp(temp_ch,"bad")==0) bad=1;

if(datain.esize==3) esize=3;
else esize=1;

  if (i3 != 0) { 
	vp_purge(); 
	vp_erase (); 
	vp_color (axis1.col[i3]); 
	}

   for (ii = 0; ii < dataout.n1[0] * dataout.n2 * esize; ii++) {
			data2[ii]=(unsigned char)data[ii];
	}
	

	if(bad==1){
   for (ii = 0; ii < dataout.n1[0] * dataout.n2 * esize; ii++) {
		if((int)data2[ii]>253) data2[ii]=(unsigned char)128;
		if((int)data2[ii]<6) data2[ii]=(unsigned char)128;
	}
	}
		

/*CHANGE*/
  if (polarity < 0)
    for (ii = 0; ii < dataout.n1[0] * dataout.n2 * esize; ii++) {
      data2[ii] = (unsigned char) 255  - data2[ii];
    }



  /*
   * If esize=3, then map the RGB triples onto the closest available
   * color.
   */
  if (datain.esize == 3) {
    if (color[0] >= '0' && color[0] <= '9') {
      for (ii = 0; ii < dataout.n1[0] * dataout.n2; ii++) {
        ired_lev = data2[datain.esize * ii];
        igreen_lev =  data2[datain.esize * ii + 1];
        iblue_lev =  data2[datain.esize * ii + 2];
        win = 0;
        win |= ((ired_lev >> (8 - redbit)) & ~(~0 << redbit)) << redoff;
        win |= ((igreen_lev >> (8 - greenbit)) & ~(~0 << greenbit)) << greenoff;
        win |= ((iblue_lev >> (8 - bluebit)) & ~(~0 << bluebit)) << blueoff;
        data2[ii] = win;
      }
    }
    else {
      for (ii = 0; ii < dataout.n1[0] * dataout.n2; ii++) {
        red_lev = data2[datain.esize * ii] / 255.;
        green_lev = data2[datain.esize * ii + 1] / 255.;
        blue_lev = data2[datain.esize * ii + 2] / 255.;
        error_win = 8.;
        for (jj = 0; jj < 256; jj++) {
          error = 2. * SQUARE (red_lev - red[jj]) + 4. * 
            SQUARE (green_lev - green[jj]) + SQUARE (blue_lev - blue[jj]);
          if (error < error_win) {
            error_win = error;
            win = jj;
            if (error == 0.) break;
          }
        }
        data2[ii] = win;
      }
    }
  }


  /*
   * Only offset the colors if we have defined a color table.
   * Otherwise, leave them alone.
   */
  if (coltab) offset = 256;
  else offset = 0;

  /*  Set up coordinate transform for raster plot */
  xpix = dataout.n1[0]; ypix = dataout.n2;
  bit = 0; ppi = 0;

  if (coordinate.yreverse) {
      new.yll = position.yll + (coordinate.max2 - coordinatec.max2) * 
        (position.yur - position.yll) / (coordinate.max2 - coordinate.min2);
      new.yur = position.yur + (coordinate.min2 - coordinatec.min2) * 
        (position.yur - position.yll) / (coordinate.max2 - coordinate.min2);
   }
  else {
    new.yll = position.yll - (coordinate.min2 - coordinatec.min2) * 
      (position.yur - position.yll) / (coordinate.max2 - coordinate.min2);
    new.yur = position.yur - (coordinate.max2 - coordinatec.max2) * 
      (position.yur - position.yll) / (coordinate.max2 - coordinate.min2);
  }
  if (coordinate.xreverse) {
    new.xll = position.xll + (coordinate.max1 - coordinatec.max1) * 
      (position.xur - position.xll) / (coordinate.max1 - coordinate.min1);
    new.xur = position.xur + (coordinate.min1 - coordinatec.min1) * 
      (position.xur - position.xll) / (coordinate.max1 - coordinate.min1);
  }
  else {
    new.xll = position.xll - (coordinate.min1 - coordinatec.min1) * 
      (position.xur - position.xll) / (coordinate.max1 - coordinate.min1);
     new.xur = position.xur - (coordinate.max1 - coordinatec.max1) * 
      (position.xur - position.xll) / (coordinate.max1 - coordinate.min1);
  }

  /*
   * OK, now we know where to tell vp_raster to put it so it lines up
   * correctly with where the axes are to go. So tell vplot to clip it
   * for us, and then plot the entire raster!
   */
  gl_vplotint (&position, &coordinate, &axis1, &axis2);
  vp_raster (data2, blast, bit, offset, xpix, ypix, new.xll, new.yll, ppi, 
    &new.xur, &new.yur, orient, invert);

  /* * Now do the axes (first tells vplot to turn the clipping off).  */

	if(multi_t ==1) {
    if(0==search_title(i3,title_out)) sprintf(title.title,"%s",title_out);
		else sprintf(title.title,"%s",title_temp);
  }
  counter = i3;
  gl_stdplot (&dataout, &coordinate, &axis1, &axis2, &grid, &title, counter, 
    fastplt, wantframe,wantframenum);

  /* * Draw a scale bar with axis */

  if(wantscalebar) {
    barmn = data2[0];
    barmx = data2[0];
    for (ii = 0; ii < dataout.n1[0] * dataout.n2 * esize; ii++) {
      if(barmx < data2[ii]) barmx = data2[ii];
      if(barmn > data2[ii]) barmn = data2[ii];
    }
     numbcol = barmx - barmn;
		/* hack  for when all data2 values are the same*/
		if(numbcol==0) numbcol=1;
    for (ii = 0; ii < numbcol; ii++)  bdata[ii] = barmn + ii;
    /* Figure out orientation of scale */
    if (bartype[0] == 'h') {
      if (barreverse == 0) {
        barorient = 0;
        if (polarity < 0) barorient = 2;
       } 
      else {
        barorient = 2;
        if (polarity < 0) barorient = 0;
      } 
    } 
    else {
      if (barreverse == 0) {
        barorient = 3;
        if (polarity < 0) barorient = 1;
      } 
      else {
        barorient = 1;
        if (polarity < 0) barorient = 3;
      } 
    }

    /* Make scale bar */
    vp_raster(bdata,blast,bit,offset,numbcol,1,barposit.xll,barposit.yll,ppi,
      &barposit.xur,&barposit.yur,barorient,1);
    /* Plot axis for scale bar */
    gl_barplot(&barposit, &baraxis, minval, maxval, bartype, barreverse, 
      counter);
  }
return(0);

}





#if defined(__STDC__) || defined(__stdc__)
int rite_byte_frame( char *buf)
#else
int rite_byte_frame(buf)
char *buf;
#endif
{

    if(datain.n1[0] * datain.n2 != srite("out",buf,datain.n1[0] * datain.n2))
      seperr("trouble writing out byte data \n");

return(0);
}

#if defined(__STDC__) || defined(__stdc__)
int finish_vplot_plot()
#else
int finish_vplot_plot()
#endif
{
vp_purge();
return(0);
}

#if defined(__STDC__) || defined(__stdc__)
int finish_byte_plot()
#else
int finish_byte_plot()
#endif
{
return(0);
}



#if defined(__STDC__) || defined(__stdc__)
int search_title(int i2,char *title_out)
#else
int search_title(i2,title_out)
/*  find the i3rd membr of labels=first:second:third*/
int i2;
char *title_out;
#endif
{
  char *ptr;
  int i, colon,junk;
  colon = 0;
  title_out[0] = '\0';
	junk=1;
  for( ptr=titles; *ptr!='\0'; ptr++ ) {
    if(*ptr == ':') {
      colon++;
      }
    else if( colon == i2 ) {
      for( i=0; *ptr!='\0' && *ptr !=':'; ptr++) {
        title_out[i++] = *ptr;
        }
      title_out[i] = '\0';
			junk=0;
      break;
      }
    }
	return(junk);
}

#endif
