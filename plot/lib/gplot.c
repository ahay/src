/*
gl_arrow


*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_arrow(float x0,float y0,float x,float y ,float r)
#else
int gl_arrow (x0, y0, x, y, r)
    float           x0, y0, x, y, r;
#endif

{
	vp_arrow (x0, y0, x, y, r);
	return 0;
}
/*
*Reads input to initialize the length of axis and the number of tics between
*label tics
*
*end self doc
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_axisint (struct axisinfo *axis1, struct axisinfo *axis2, 
	struct coordinfo *coord, struct plotposition *position)
#else
int gl_axisint (axis1, axis2, coord, position)
    struct axisinfo *axis1;
    struct axisinfo *axis2;
    struct coordinfo *coord;
    struct plotposition *position; /* isn't used ?? */
#endif
{
int             wantaxis, i, kk, counter, n1tic, n2tic;
float           dnum, labelsz, tempmin;
int             axisfat[NPMAX + 1], axiscol[NPMAX + 1], labelfat;
char            wheretics[10], label[280], wherelabel1[1], wherelabel2[1];
int             numdxnum, numdynum; 
int             numo1, numo2, numaxisor1, numaxisor2;
float           tempnum;
/*
 * initializing axisfat and axiscol
 */

    for (i = 0; i < NPMAX + 1; i++)
    {
	axisfat[i] = 0;
	axiscol[i] = 7;
    }
	labelsz = 8;
	labelfat = 0;
    /*
     * fetching axisfat , if fetched then find # fetched 
     */
	kk = 0;
    if (getch ("axisfat", "d", axisfat))
    {
	/* calculating # of axisfat fetched */
	while (axisfat[kk] != 0)
	    kk++;
	/*
	 * checking to see if # of axisfats fetched exceed NPMAX if exceed
	 * NPMAX will exit  
	 */
	if (kk > NPMAX)
	{
	    seperr ("Too many values for axisfat were entered, exceeding NPMAX\n");
	    exit (-1);
	}
    }
    /* assigning axisfat to the structures */
    for (i = 0; i < kk; i++)
    {
	axis1->fat[i] = axisfat[i];
	axis2->fat[i] = axisfat[i];
    }
    /*
     * assigning axisfat to the structures after the initial count has been
     * done.  Axisfat will rotate the fatness 
     */
    for (i = kk; i < NPMAX + 1; i++)
    {
	if (kk == 0)
	    kk++;
	axis1->fat[i] = axisfat[i % kk];
	axis2->fat[i] = axisfat[i % kk];
    }

    /*
     * fetching axiscol , if fetched then find # fetched 
     */
    if (getch ("axiscol", "d", axiscol))
    {
	kk = 0;
	/* calculating # of axiscol fetched */
	while (axiscol[kk] != 7)
	    kk++;
	/* checking to see if # of axiscol fetched exceed NPMAX */
	/* if exceed NPMAX will exit  */
	if (kk > NPMAX)
	{
	    seperr ("Too many values for axiscol were entered, exceeding NPMAX\n");
	    exit (-1);
	}
    }
    /* assigning axiscol to the structures */
    for (i = 0; i < kk; i++)
    {
	axis1->col[i] = axiscol[i];
	axis2->col[i] = axiscol[i];
    }
    /*
     * assigning axiscol to the structures after the initial count has been
     * done.  Axiscol will rotate the color of the axis 
     */
    for (i = kk; i < NPMAX + 1; i++)
    {
	if (kk == 0)
	    kk++;
	axis1->col[i] = axiscol[i % kk];
	axis2->col[i] = axiscol[i % kk];
    }

    /*
     * fetching labelfat , if fetched then find # fetched 
     */
    getch ("labelfat", "d", &labelfat);
    axis1->labelfat = labelfat;    
    axis2->labelfat = labelfat;    
    /*
     * fetching labelsz , if fetched then find # fetched 
     */
    getch ("labelsz", "f", &labelsz);
    axis1->labelsz = labelsz;    
    axis2->labelsz = labelsz;    
    /* initialize the axis1.label to be blank */
    strcpy(axis1->label, " ");
    /* fetching the axis1.label */
    fetch ("label1", "s", axis1->label);
    /*
     * The defaults are not set in here for wherelabel a value needs to be
     * for it  your main program should set the defaults, will fetch the
     * override values 
     */
    if (getch ("wherexlabel", "s", wherelabel1))
       switch (*wherelabel1)
       {
       case 't':
       case 'b':
       strcpy ( axis1->wherelabel, wherelabel1);
       break;
       default:
           fprintf (stderr, "wherexlabel should be either t or b not %s \n", wherelabel1);
       }
    if (getch ("whereylabel", "s", wherelabel2))
       switch (*wherelabel2)
       {
       case 'r':
       case 'l':
       strcpy ( axis2->wherelabel, wherelabel2);
       break;
       default:
           fprintf (stderr, "whereylabel should be either l or r not %s \n", wherelabel2);
       }
    /* checking to see if wantaxis is fetched */
    if (getch ("wantaxis", "1", &wantaxis) == 0) 
    {
	/* setting the default to put the axes on the plot */
        wantaxis = 1;
	axis1->wantaxis = 1;
	axis2->wantaxis = 1;
    }
    /* if no axis is wanted then set the structures */
    if (wantaxis == 0)
    {
	axis1->wantaxis=0;
	axis2->wantaxis=0;
    }
    /* else check for each individual axis */
    else
    {
	if (getch ("wantaxis1", "1", &axis1->wantaxis) == 0)
	    axis1->wantaxis = 1;
	if (getch ("wantaxis2", "1", &axis2->wantaxis) == 0)
	   axis2->wantaxis = 1;
    }
    if(coord->transp)
    {
    wantaxis = axis1->wantaxis;
    axis1->wantaxis = axis2->wantaxis;
    axis2->wantaxis = wantaxis;
    }

    /* check to see were to put the tics for the axes */
    if (getch ("wheretics", "s", wheretics) == 0)
	strcpy(wheretics," ");	/* if not set by user default to frame */
    /* set to structures */
    strcpy(axis1->wheretics, wheretics);
    strcpy(axis2->wheretics, wheretics);
    /* fetch position of axisor1 */
    numaxisor1 = getch ("axisor1", "f", &axis1->axisor);
       if (numaxisor1 == 0 )
       {
      if (*axis1->wherelabel == 'b')
        axis1->axisor = coord->min2;  /* set default position to min2 */
       else
       if (*axis1->wherelabel == 't')
	axis1->axisor = coord->max2;  /* set default position to max2 */
       }
    numaxisor2 = getch ("axisor2", "f", &axis2->axisor); 
     if ( numaxisor2 == 0)
     {
       if (*axis2->wherelabel == 'l')
        axis2->axisor = coord->min1;
       else
       if (*axis2->wherelabel == 'r')
        axis2->axisor = coord->max1;
     }
     if (coord->transp)
      {
     if ( numaxisor2 != 0 || numaxisor1 != 0)
     {
       tempnum = axis1->axisor;
       axis1->axisor = axis2->axisor;
       axis2->axisor = tempnum;
      }
     } 
    if (getch ("n1tic", "d", &n1tic) == 0)
	n1tic = 1;
    if (getch ("n2tic", "d", &n2tic) == 0)
	n2tic = 1;
     if (coord->transp)
     {
     axis1->ntic = n2tic;
     axis2->ntic = n1tic;
     }
     else
     { 
     axis1->ntic = n1tic;
     axis2->ntic = n2tic;
     } 
    numdxnum = getch ("d1num", "f", &axis1->dnum);
    numdynum = getch ("d2num", "f", &axis2->dnum);
    if (coord->transp)
    {
    if (numdxnum != 0 || numdynum != 0 )
    {
     tempnum = axis2->dnum; 
     axis2->dnum = axis1->dnum; 
     axis1->dnum = tempnum; 
     tempnum = numdxnum;
     numdxnum = numdynum;
     numdynum = tempnum;
    }   
    }
	if (numdxnum == 0 )
        gl_getscl (coord, axis1);
	if (numdynum == 0 )
	gl_getscl (coord, axis2);
 
    if (axis1->dnum != 0.)
    {
       if (coord->min1 < coord->max1)
	tempmin = coord->min1;  /* set default position to min2 */
       else
	tempmin = coord->max1;  /* set default position to min2 */
	for (axis1->num0 = (int) (tempmin / axis1->dnum) * axis1->dnum - axis1->dnum; axis1->num0 < tempmin; axis1->num0 += axis1->dnum);
    }


    strcpy(axis2->label,  " ");
    fetch ("label2", "s", axis2->label);
    if (axis2->dnum != 0.)
    {
       if (coord->min2 < coord->max2)
	tempmin = coord->min2;  /* set default position to min2 */
       else
	tempmin = coord->max2;  /* set default position to min2 */
	for (axis2->num0 = (int) (tempmin / axis2->dnum) * axis2->dnum - axis2->dnum; axis2->num0 < tempmin; axis2->num0 += axis2->dnum);
    }
    if (coord->transp)
    {
    numo1 = getch ("o1num", "f", &axis2->num0);
    numo2 = getch ("o2num", "f", &axis1->num0);
    }
    else
    {
    numo1 = getch ("o1num", "f", &axis1->num0);
    numo2 = getch ("o2num", "f", &axis2->num0);
    }
    if ( axis1->dnum !=0.)
    {
      n1tic = axis1->ntic;
     if (axis1->ntic == 0)
    {
	n1tic = 1;
    }
    axis1->dtic = axis1->dnum / n1tic ;
       if (coord->min1 < coord->max1)
	tempmin = coord->min1;  /* set default position to min2 */
       else
        tempmin=coord->max1;
    for (axis1->tic0 = axis1->num0 - axis1->ntic * axis1->dtic; axis1->tic0 < tempmin; axis1->tic0 += axis1->dtic);
    }
    if (axis2->dnum !=0.)
    {
       if (coord->min2 < coord->max2)
	tempmin = coord->min2;  /* set default position to min2 */
       else
	tempmin = coord->max2;  /* set default position to min2 */

      n2tic = axis2->ntic;
    if (axis2->ntic == 0)
    {
	n2tic = 1;
    }
    axis2->dtic = axis2->dnum / n2tic;
    for (axis2->tic0 = axis2->num0 - axis2->ntic * axis2->dtic; axis2->tic0 < tempmin; axis2->tic0 += axis2->dtic);
    }
    
       
    if (coord->transp)
    {
       strcpy(label, axis1->label);
       strcpy(axis1->label,axis2->label);
       strcpy(axis2->label, label);
    }
    if( coord->yreverse )
        gl_rotate1(&axis1->axisor, coord->min2, coord->max2);
    if( coord->xreverse )
        gl_rotate1(&axis2->axisor, coord->min1, coord->max1);
	return 0;
}
/*
* This routine initializes scale bar device.
* 
* Author - Hector Urdaneta (SEP)   Jan 3 95
*
*Edit History
*	
*  Bob - bad hack to allow minval and maxval to be passed in rather
*         requiring from other source (getch, file)
*/
#include <glplot.h>
#ifdef USE_PROTO
int gl_barint (struct plotposition *position, struct axisinfo *axis1,                         struct plotposition *barposit, struct axisinfo *baraxis,                       float *minval, float *maxval, char *bartype,int *barreverse,                   int nplots,int cubeplot)
#else
int gl_barint (position, axis1, barposit, baraxis, minval, maxval, bartype,                   barreverse, nplots,cubeplot)
struct plotposition *position, *barposit;
struct axisinfo *axis1, *baraxis;
float *minval, *maxval;
char *bartype;
int *barreverse, nplots,cubeplot;
#endif
{

    int kk, ii;
    char wherebarlabel[2], wherebartics[10], label3[280];
    int esize, numminval, nummaxval;
    int wantbaraxis,nf,tempi;
		int flat,ierr_bar;
    float dbarnum,point1,point2;
		float mypct;

    /* A vertical scale bar is the default */
		if(cubeplot==1) bartype[0] = 'h';
		else bartype[0] = 'v';
    ierr_bar=getch("bartype","s",bartype);
    if(bartype[0] != 'h' && bartype[0] != 'v')
        seperr("bartype option not implemented\n");

    /* scale bar width in inches */
    if(!getch("barwidth", "f", &barposit->screenwd)) barposit->screenwd = .36;

    /* Read in maxvals and minvals. First asume that the user
     * is passing the values in the command line as: 
     * maxval=..,..,..,...
     * minval=..,..,..,...
     * If not look in the history file  */
    kk = 0;
    if(kk = fetch("minval", "f", minval)) {
	/* calculate # of minvals and exit if greater than NPMAX */   

	if (kk > NPMAX){
	    seperr ("Too many values for minval were entered, exceeding NPMAX\n");
	    exit(-1);
	/* look to see if the # of minvals is less than the # of plots.
	 * This could be the case in which the user may want all the
	 * scale bars to have the same minvals */
	} else if (kk < nplots) {
	    for (ii = kk; ii < nplots; ii++) {
		minval[ii] = minval[kk-1];
	    }
	}
    /* Finally check if the minvals are specified in a file */
    } else if(auxin("mintmp")) {
	auxpar("esize","d",&esize,"mintmp");
	auxpar("n1","d",&numminval,"mintmp");
	/* The # of minvals need to be equal to the # of plots */
	if(numminval > nplots) numminval = nplots;

	sreed("mintmp", minval, numminval*esize);

	if (numminval < nplots) {
	    for (ii = numminval; ii < nplots; ii++) { minval[ii] = minval[kk-1]; }
	}
    }
	/* this is a terrible hack, but without changing interface this seems
    to be to only way to support minval and maxval not to be calculated
    in the calling program */
	else if(minval[0] ==0 && maxval[0] ==0 && (minval[1] != 0 || maxval[1] !=0)){
			kk=1; minval[0]=minval[1]; maxval[0]=maxval[1];
	    for (ii = kk; ii < nplots; ii++) { minval[ii] = minval[kk-1];  maxval[ii]=maxval[0];}
	}
	else{
	   seperr ("In order to have a scale bar you must specify its minimum value");
	    exit(-1);
	}

    /* Same for the maxvals */
    kk = 0;
    if(kk = fetch("maxval","f",maxval)) {
	/* calculate # of maxvals and exit if greater than NPMAX */   

	if (kk > NPMAX) {
	    seperr ("Too many values for maxval were entered, exceeding NPMAX\n");
	    exit(-1);
	} else if (kk < nplots) {
	    for (ii = kk; ii < nplots; ii++) {
			maxval[ii] = maxval[kk-1];
	    }
	}
    } else if(auxin("maxtmp")) {
	auxpar("esize","d",&esize,"maxtmp");
	auxpar("n1","d",&nummaxval,"maxtmp");
	if(nummaxval > nplots) nummaxval = nplots;

	sreed("maxtmp", maxval, nummaxval*esize);

	if (nummaxval < nplots) {
	    for (ii = nummaxval; ii < nplots; ii++) {
		maxval[ii] = maxval[kk-1];
	    }
	}
    }

    /* Bar Axis info */
		if( 1==getch("xll","d",&tempi)|| 1==getch("xur","d",&tempi)||
      1==getch("yll","d",&tempi) || 1==getch("yur","d",&tempi)){
				nf=1;
				 if(0==getch("bar.xll","d",&tempi)) 
						seperr("If specifyin xll, yll, xur, or yur must specify bar.xll\n");
				 else barposit->xll=tempi;
				 if(0==getch("bar.yll","d",&tempi)) 
						seperr("If specifyin xll, yll, xur, or yur must specify bar.yll\n");
				 else barposit->yll=tempi;
				 if(0==getch("bar.xur","d",&tempi)) 
						seperr("If specifyin xll, yll, xur, or yur must specify bar.xur\n");
				 else barposit->xur=tempi;
				 if(0==getch("bar.yur","d",&tempi)) 
						seperr("If specifyin xll, yll, xur, or yur must specify bar.yur\n");
				 else barposit->yur=tempi;
			}
			else nf=0;


    /* For a horiz. bar && barreverse == FALSE then bar goes
       from min (left) to max (right) value, while for a vert. 
       bar the max is located at the top and the min at the
       bottom. */
    if(!getch("barreverse", "1", barreverse)) *barreverse = 0;

    if(!getch("barlabelsz", "f", &baraxis->labelsz)) {
	baraxis->labelsz = axis1->labelsz; }
    if(!getch("barlabelfat", "d", &baraxis->labelfat)) {
	baraxis->labelfat = axis1->labelfat; }

    /* Initialize barlabel to be label3. */
		strcpy(label3, " ");
    fetch ("label3", "s", label3);
    strcpy(baraxis->label, label3);
    getch ("barlabel", "s", baraxis->label);

    for (kk = 0; kk < NPMAX; kk++) {
	baraxis->fat[kk] = axis1->fat[kk];
	baraxis->col[kk] = axis1->col[kk];
    }

    /* check if user wants a scale bar axis */
    baraxis->wantaxis = 1;
    if (getch ("wantbaraxis", "1", &wantbaraxis)) {
	if(!wantbaraxis) baraxis->wantaxis=0;
    }

    /* check to see where to put the tics for the scale bar axis */
    if (! getch ("wherebartics", "s", wherebartics)) strcpy(wherebartics," ");
    strcpy(baraxis->wheretics, wherebartics);

    /* Spacing of tics for bar scale */
    baraxis->dnum = 0.;
    if (getch ("dbarnum", "f", &dbarnum) != 0)
	baraxis->dnum = dbarnum;

      
if(bartype[0] == 'h') {
		if(nf==0){
	barposit->xll = position->xll;
	barposit->xur = position->xur;
	barposit->yll = .12 * position->screenht;
	barposit->yur = barposit->yll + barposit->screenwd;
	

	position->yur += .4*barposit->screenwd;
	if(cubeplot==1)
	position->yll += .12*position->screenht + barposit->screenwd;
	else
	position->yll += .04*position->screenht + barposit->screenwd;
	}

	baraxis->inch = barposit->xur - barposit->xll;

	if (getch ("wherebarlabel","s",wherebarlabel)) {
	    switch (*wherebarlabel) {
	        case 't':
	        case 'b':
		strcpy ( baraxis->wherelabel, wherebarlabel );
		break;
	        default:
		seperr("wherebarlabel should be either t or b for a horizontal scale bar\n");
	    }
	}
    
    } else {    

	if(nf==0){
	barposit->xur = .88 * position->screenwd;
	barposit->xll = barposit->xur - barposit->screenwd;
	barposit->yll = position->yll;
	barposit->yur = position->yur;
	

	position->xur -= (.07*position->screenwd + barposit->screenwd);
	position->xll -= .5*barposit->screenwd;
	}
	baraxis->inch = barposit->yur - barposit->yll;

	if (getch ("wherebarlabel","s",wherebarlabel)) {
	    switch (*wherebarlabel) {
	        case 'r':
	        case 'l':
		strcpy ( baraxis->wherelabel, wherebarlabel );
		break;
	        default:
		seperr("wherebarlabel should be either r or l for a vertical scale bar\n");
	    }
	}
    }








	return 0;
}
/*
* This routine draws a frame and an axis around a rectangular bar.
* Part of the code was taken from Steve's gl_simpleaxis.c
*
* Author - Hector Urdaneta (SEP)   Jan 3 95
*/
#include <glplot.h>
#define UNIT      1.000001
#ifdef USE_PROTO
int gl_barplot (struct plotposition *posit, struct axisinfo *axis,                       float *minval, float *maxval, char *type, int reverse, int counter)
#else
int gl_barplot (posit, axis, minval, maxval, type, reverse, counter)
struct plotposition *posit;
struct axisinfo *axis;
float *minval, *maxval;
char *type;
int reverse;
int counter;
#endif
{
    float pad, ch, x1, x2, y1, y2;
    float loc, num, ltic, final;
    float xpath, ypath, xup, yup;
    char string[10];
    float costh, sinth, xpos, ypos;
    float dxtic, dytic;
		int logic;
	  float num0,loc0;

    /* plot frame */    

		  vp_fat (axis->labelfat);

    gl_color(axis->col[counter]);
    gl_move(posit->xll, posit->yll);
    gl_draw(posit->xll, posit->yur);
    gl_draw(posit->xur, posit->yur);
    gl_draw(posit->xur, posit->yll);
    gl_draw(posit->xll, posit->yll);

    if(axis->wantaxis) {

	/* pad is the space skipped between tic and number,
	   and number and label. ltic is the longitude of the tic */
	ltic = .1;
	pad = 0.15;
	ch = axis->labelsz / 33.;
 
	if(type[0] == 'h') {

	    x1 = posit->xll;
	    y1 = posit->yll;
	    x2 = posit->xur;
	    y2 = posit->yll;
	
	} else { /* type[0] == 'v' */
	    x2 = posit->xur;
	    y2 = posit->yur;
	    x1 = posit->xur;
	    y1 = posit->yll;
	}

	gl_move(x1,y1);
	gl_draw(x2,y2);

	/* x and y sizes of tic marks */
	costh = (x2 - x1) / axis->inch;
	sinth = (y2 - y1) / axis->inch;
	dxtic = ltic * sinth;
	dytic = -ltic * costh;

	/* 
	 * Figure out which quadrant we are in.
	 * This determines which side of the line the tic
	 * marks will go on, and how the text will be aligned.
	 * For each quadrant we compute the vplot text path
	 * and up vectors.
	 */
	if (x1 <= x2 && y1 <= y2) {
	    vp_tjust(TH_CENTER,TV_TOP);
	    xpath = ch * costh;
	    ypath = ch * sinth;
	    xup = -ch * sinth;
	    yup = ch * costh;
	} else if (x1 > x2 && y1 <= y2) {	
	    vp_tjust(TH_CENTER,TV_BOTTOM);
	    xpath = -ch * costh;
	    ypath = -ch * sinth;
	    xup = ch * sinth;
	    yup = -ch * costh;
	} else if (x1 > x2 && y1 > y2) {	
	    vp_tjust(TH_CENTER,TV_BOTTOM);
	    xpath = -ch * costh;
	    ypath = -ch * sinth;
	    xup = ch * sinth;
	    yup = -ch * costh;
	} else if (x1 <= x2 && y1 > y2) {	
	    vp_tjust(TH_CENTER,TV_TOP);
	    xpath = ch * costh;
	    ypath = ch * sinth;
	    xup = -ch * sinth;
	    yup = ch * costh;
	}

	if(axis->dnum==0.) logic=0;
	else logic=1.;
	/* find the optimum tic mark interval */
	if(!reverse) {
	    axis->num0 = minval[counter];
	    if(logic==0) {
		gl_opttic (minval[counter], maxval[counter], axis->inch, axis->num0, &axis->dnum, axis->labelsz);	}
	    if(maxval[counter] < minval[counter]) axis->dnum *= -1;
	} else {
	    axis->num0 = minval[counter];
	    if(logic==0) {
		gl_opttic (maxval[counter], minval[counter], axis->inch, axis->num0, &axis->dnum, axis->labelsz);       }
	    if(maxval[counter] > minval[counter]) axis->dnum *= -1;
	    axis->num0 = maxval[counter];
	}

	/* compute tic mark spacing */
	axis->dtic = axis->dnum / (maxval[counter]-minval[counter]) * axis->inch;
	if(axis->dtic < 0.) axis->dtic *= -1;


	if(axis->dnum>0.){
		if(axis->num0 < 0.){
			num0=axis->dnum*((int)(axis->num0/axis->dnum));
		}
		else {
			num0=axis->dnum*((int)(axis->num0/axis->dnum+.999999));
		}
			loc0=(num0-axis->num0)*axis->dtic/axis->dnum;
	}
	else{
		num0=axis->dnum*((int)(axis->num0/axis->dnum));
		loc0=(num0-axis->num0)*axis->dtic/axis->dnum;
	}
	



	
	/*
	 * move to each tic mark location, draw the tic and the number
	 */
	for (loc=loc0, num=num0; loc <= axis->inch; loc+=axis->dtic, num+=axis->dnum) { 
	    sprintf(string,"%1.3g",num);
	    xpos = x1 + loc * costh;
	    ypos = y1 + loc * sinth;
	    vp_move(xpos,ypos);
	    vp_draw(xpos+dxtic,ypos+dytic);
	    vp_gtext(xpos+dxtic+pad*sinth,ypos+dytic-pad*costh,xpath,ypath,xup,yup,string);
	}
	/* now the axis label */
	xpos = x1 + loc/2. * costh + dxtic + ch * sinth + 2. * pad * sinth;
	ypos = y1 + loc/2. * sinth + dytic - ch * costh - 2. * pad * costh;
  xpos=  x1 + axis->inch * costh * .5 + dxtic + ch * sinth + 2. * pad * sinth;
  ypos=  y1 + axis->inch * sinth * .5 + dytic - ch * costh - 2. * pad * costh;

	vp_gtext(xpos,ypos,xpath,ypath,xup*1.1,yup*1.1,axis->label);

		if(logic==0) axis->dnum=0.;
    }
	
	return 0;
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_clip(float pos1, float pos2, float pos3, float pos4)
#else
int gl_clip(pos1, pos2, pos3, pos4)
float pos1, pos2, pos3, pos4;
#endif
{
vp_clip(pos1, pos2, pos3, pos4);
return 0;
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_color(int color)
#else
int gl_color(color)
int color;
#endif
{
vp_color(color);
return 0;
}
/*
*
*
* The routine will fetch and initialize plot and axis parameters.
*
* 
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_colorint (struct colorinfo *color)
#else
int gl_colorint (color)
    struct colorinfo *color;
#endif
{
int             ii;
    for (ii = 0; ii < 3; ii++)
    {
	color->backcol[ii] = 0.;
    }
    getch ("backcol", "f", color->backcol);
    for (ii = 0; ii < 3; ii++)
    {
	color->fillcol[ii] = color->backcol[ii];
    }
    getch ("fillcol", "f", color->fillcol);
	return 0;
}
/*
*
*  This subroutine will getch in the variables needed to initialize
*  the plot window. This routine will fetch xll, yll, *  xur, yur.  
*  If the user does not set  any of these variables then
*  a default value will be used.  The structures and SCREENWD and
*  SCREENHT  are defined in glplot.h
*
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_coordint( struct plotposition *pos, struct coordinfo *coord, struct axisinfo
*axis1, struct axisinfo *axis2)
#else
int gl_coordint (pos, coord, axis1, axis2)
    struct plotposition *pos;
    struct coordinfo *coord;
    struct axisinfo *axis1;
    struct axisinfo *axis2;
#endif
{
int             temp1, temp2;
float		crowd, crowd1, crowd2;
    /*
     * initialization of lower left corner counter 
     */
    temp1 = 0;
    /*
     * initialization of upper right corner counter 
     */
    temp2 = 0;
    /*
     * if user does not specify a screenwd a default one defined in glplot.h
     * will be used 
     */
    if (getch ("screenratio", "f", &pos->screenratio) == 0)
	pos->screenratio = SCREEN_RATIO;
    if (getch ("screenht", "f", &pos->screenht) == 0)
	pos->screenht =SCREENHT;
    if (getch ("screenwd", "f", &pos->screenwd) == 0)
	pos->screenwd = pos->screenht / pos->screenratio;

/* want to crowd off the axes?   .75 < crowd  < 1.00
   want to crowd the  1-axis ?   .75 < crowd1 < 1.00
   want to crowd the  2-axis ?   .75 < crowd2 < 1.00
   */
     if (!getch ("crowd" , "f", &crowd ))  crowd  = .75;
     if (!getch ("crowd1", "f", &crowd1))  crowd1 = crowd;
     if (!getch ("crowd2", "f", &crowd2))  crowd2 = crowd;

/* calculate inches */
     if (!getch ("xinch", "f", &axis1->inch))  axis1->inch = pos->screenwd * crowd1;
     if (!getch ("yinch", "f", &axis2->inch))  axis2->inch = pos->screenht * crowd2;
    /*

     * initialization of xll, yll, xur, yur to 0. 
     */

    pos->xll = 0;
    pos->xur = pos->xll;
    pos->yll = 0;
    pos->yur = pos->yll;

    /*
     * if xll is input the temp1 is 1 
     */
    if (getch ("xll", "f", &pos->xll))
	temp1 = 1;
    /*
     * if yll is input the temp1 is temp1 + 1 
     */
    if (getch ("yll", "f", &pos->yll))
	temp1 = temp1 + 2;
    switch (temp1)
    {
	/*
	 * if temp1 is 1 then yll was not set 
	 */
    case 1:
	seperr ("yll was not set\n");
	/*
	 * reseting xll back to 0 will be using inch1 and inch2 to determine
	 * plot position and size 
	 */
	pos->xll = 0;
	temp1 = 0;
	break;
	/*
	 * if temp1 is 2 then xll was not set 
	 */
    case 2:
	seperr ("xll was not set \n");
	pos->yll = 0;
	temp1 = 0;
	/*
	 * reseting yll back to 0 will be using inch1 and inch2 to determine
	 * plot position and size 
	 */
	break;
    }
    if (getch ("xur", "f", &pos->xur))
	temp2 = 1;
    if (getch ("yur", "f", &pos->yur))
	temp2 = temp2 + 2;
    switch (temp2)
    {
    case 1:
	seperr ("yur was not set\n");
	pos->xur = pos->xll;
	temp2 = 0;
	/*
	 * reseting xur back to 0 will be using inch1 and inch2 to determine
	 * plot position and size 
	 */
	break;
    case 2:
	seperr ("xur was not set\n ");
	pos->yur = pos->yll;
	temp2 = 0;
	/*
	 * reseting yur back to 0 will be using inch1 and inch2 to determine
	 * plot position and size 
	 */
	break;
    }

    if (temp1 == 0 && temp2 != 0)
    {
	pos->xll = pos->xur - axis1->inch;
	pos->yll = pos->yur - axis2->inch;
    }
    if (temp1 == 3 && temp2 == 0)
    {
	pos->xur = pos->xll + axis1->inch;
	pos->yur = pos->yll + axis2->inch;
    }
    if (temp1 == 3 && temp2 == 3)
    {
       axis1->inch = pos->xur - pos->xll;
       axis2->inch = pos->yur - pos->yll;
    }
    if (temp1 == 0 && temp2 == 0)
    {
     float marg1, marg2;
     marg1 =  (pos->screenwd - axis1->inch);
     marg2 =  (pos->screenht - axis2->inch);
     pos->xll = marg1 * 2./3.;
     pos->yll = marg2 * 1./2.;
     pos->xur = pos->screenwd - marg1 * 1./3.;
     pos->yur = pos->screenht - marg2 * 1./2.;
    }
    /*
     * fetching transp,  transp default needs to be set in the main
     * program
     */
      getch ("transp","1", &coord->transp);
      getch ("xreverse","1", &coord->xreverse);
      getch ("yreverse","1", &coord->yreverse);
      if (getch ("labelrot", "1", &coord->labelrot) == 0 )
           coord->labelrot = -1 ;
      else
         if(coord->labelrot == 1)
            coord->labelrot = -1;  
         if(coord->labelrot == 0)
            coord->labelrot = 1;  
	return 0;
}
/*
*
*
*This routine will set the the device parameters
* namely color of the background, plot color and ploot fat
*
* Modifications:
* 05-14-91	W. Bauske
*		Make it use a pointer instead of a structure on the stack
*
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_dash (struct dashinfo *dash)
#else
int gl_dash (dash)
    struct dashinfo *dash;
#endif
{
    vp_setdash (dash->dash, dash->gap, 2);
	return 0;
}
/*
*  gl_dashfig(dash, n3_loop)
*
*  Input is dash->dashtype 0-9
*	0 continuous   DEFAULT
*	1 fine dash
*	2 fine dot
*	3 dash
*	4 large dash
*	5 dot dash
*	6 large dash small dash
*	7 double dot
*	8 double dash
*	9 loose dash  The part after the decimal point determines 
*                     the pattern repetition interval
*
*     	Output is dash[0], gap[0], dash[1], gap[1]	
*       determining the type and size of the dashing
*
*
*	END OF SELF DOC
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_dashfig (struct dashinfo *dash, int n3_loop)
#else
int gl_dashfig (dash, n3_loop)
    struct dashinfo *dash;
    int             n3_loop;
#endif
{
float           patsize;
    switch ( (int) dash->dashtype[n3_loop])	/* determining type of dash */
    {
    case 1:
	patsize = .3;
	dash->dash[0] = dash->dash[1] = 1.;
	dash->gap[0] = dash->gap[1] = 2.;
	break;
    case 2:
	patsize = .2;
	dash->dash[0] = dash->dash[1] = 1.;
	dash->gap[0] = dash->gap[1] = 6.;
	break;
    case 3:
	patsize = .4;
	dash->dash[0] = dash->dash[1] = 4.;
	dash->gap[0] = dash->gap[1] = 2.;
	break;
    case 4:
	patsize = .6;
	dash->dash[0] = dash->dash[1] = 3.;
	dash->gap[0] = dash->gap[1] = 2.;
	break;
    case 5:
	patsize = .5;
	dash->dash[0] = (.3);
	dash->dash[1] = 3.;
	dash->gap[0] = dash->gap[1] = 1.;
	break;
    case 6:
	patsize = .6;
	dash->dash[0] = 4.;
	dash->dash[1] = 2.;
	dash->gap[0] = dash->gap[1] = 1.;
	break;
    case 7:
	patsize = .4;
	dash->dash[0] = dash->dash[1] = 1.;
	dash->gap[0] = 2.;
	dash->gap[1] = 4.;
	break;
    case 8:
	patsize = .8;
	dash->dash[0] = dash->dash[1] = 5.;
	dash->gap[0] = 2.;
	dash->gap[1] = 4.;
	break;
    case 9:
	patsize = .6;
	dash->dash[0] = dash->dash[1] = 1.;
	dash->gap[0] = dash->gap[1] = 1.;
	break;
    case 0:
	dash->dash[0] = 0.;
	dash->gap[0] = 0.;
	dash->dash[1] = 0.;
	dash->gap[1] = 0.;
	break;
    default:
	dash->dash[0] = 0.;
	dash->gap[0] = 0.;
	dash->dash[1] = 0.;
	dash->gap[1] = 0.;
	break;
    }
    if (dash->dash[0] + dash->dash[1] + dash->gap[0] + dash->gap[1] != 0.)
    {
	/*
	 * checking to see if default case 
	 */
	/*If not default case then find the decimal part of dash->dashtype*/

	dash->dashtype[n3_loop] -= (float) (int ) (dash->dashtype[n3_loop]);
/* If dash->dashtypeis 0 then set dash->dashtypeto .4 */

	if (dash->dashtype[n3_loop] == 0.)
	    dash->dashtype[n3_loop] = (.4);
/* Computing the pattern size for scaling  */
	patsize *= dash->dashtype[n3_loop] / (.4 * (dash->dash[0] + dash->dash[1] + dash->gap[0] + dash->gap[1]));
/* scaling by the pattern size */
	dash->dash[0] *= patsize;
	dash->dash[1] *= patsize;
	dash->gap[0] *= patsize;
	dash->gap[1] *= patsize;
    }
	return 0;
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_draw (float position1, float position2)
#else
int gl_draw (position1, position2)
    float           position1, position2;
#endif
{
    vp_draw (position1, position2);
	return 0;
}
#include "glplot.h"
int gl_erase ()
{
return (vp_erase());

}
#include "glplot.h"
#ifdef USE_PROTO
int gl_fat(int fat)
#else
int gl_fat(fat)
int fat;
#endif
{
return(vp_fat(fat));
}
/*
*
*
*This routine will fill in side the frame around the picture
*
*
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_fillin (struct coordinfo *coordinate, struct colorinfo *color)
#else
int gl_fillin (coordinate, color)
    struct coordinfo *coordinate;
    struct colorinfo *color;
#endif
{
float           xp[4], yp[4];
int             lp, fat, xmask, ymask;
    vp_coltab (8, color->fillcol[0], color->fillcol[1], color->fillcol[2]);
    vp_color (8);
    xp[0] = coordinate->min1;
    xp[1] = coordinate->max1;
    xp[2] = coordinate->max1;
    xp[3] = coordinate->min1;
    yp[0] = coordinate->min2;
    yp[1] = coordinate->min2;
    yp[2] = coordinate->max2;
    yp[3] = coordinate->max2;
    lp = 4;
    fat = 0;
    ymask = 1;
    xmask = 1;
    vp_uarea (xp, yp, lp, fat, xmask, ymask);
	return 0;
} 
/*
*	gl_framenumber(n3_loop, d3, 03, xmin,ymin, labelsz)
*
*       Plots the frame number when in the movie option 
*
* modified - S. Cole - 24Nov92
*   changed default format to 3 decimal places. also if the
*   value passed is an integer, this routine figures it out
*   and plots it as such.
*
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_framenum (int n3_loop, float d3, float o3, float xmin, float ymin, 
	     float labelsz )
#else
int gl_framenum (n3_loop, d3, o3, xmin, ymin, labelsz)
    int             n3_loop;
    float           d3, o3, xmin, ymin, labelsz;
#endif
{

float           ch, vs, xc, yc, x, y, xpath, ypath, xup, yup, temp;
char            string[80];
    vp_fat (0);
    temp = n3_loop * d3 + o3;
    sprintf (string, "%.3f", temp);
    if (temp-((int) temp) == 0.) sprintf (string, "%.0f", temp);
    vp_umove (xmin, ymin);
    vp_where (&xc, &yc);
    ch = labelsz / 33.;
    vs = ch * 5. / 10.;
    gl_tjust ("b");
    x = xc;
    y = yc -  3 * ch - 2 * vs ;
    ypath = 0.;
    xpath = ch;
    yup = ch;
    xup = 0.;
    vp_gtext (x, y, xpath, ypath, xup, yup, string);
	return 0;
}
/*
 * Modifications:
 * W. Bauske IBM 04-01-91
 *	RS/6000 uses macros for log10() and pow()
 */
#include <stdio.h>
#include <math.h>
static float    sqr[3] = {1.414214, 3.162278, 7.071068};
static float    vint[4] = {1., 2., 5., 10.};
#include "glplot.h"
#ifdef USE_PROTO
int gl_getscl(struct coordinfo *coord, struct axisinfo *axis)
#else
int gl_getscl(coord, axis)
struct coordinfo *coord;
struct axisinfo *axis;
#endif
{
float           temp, temp1, temp2, temp3, temp4, temp5, tempdnum, num;
float           inch, min, max, num0, dnum;
double          div, div1, div3, a, b ;
int             div2, i, length, counter;
char            string[100];
    if (*axis->wherelabel == 't' || *axis->wherelabel == 'b' )
    {
     min = coord->min1;
     max = coord->max1;
     inch = axis->inch;
     num0 = axis->num0;
     
    }
    if (*axis->wherelabel == 'l' || *axis->wherelabel == 'r' )
    {
     min = coord->min2;
     max = coord->max2;
     inch = axis->inch;
     num0 = axis->num0;
    }
    if (inch <= 0)
    {
	fprintf (stderr, "non positive line length in get_scl");
	exit (-1);
    }
    if (fabs(max - min) < .000001)
    {
	axis->dnum = .5;
	return;
    }
    if (max > min)
	temp = max - min;
    else
	temp = min - max;
    if ((fabs(temp)/inch) < .0000001)
    {
    if (max > min)
    {
	temp = max - min;
	axis->dnum = .5;
    }
    else
    {
	temp = min - max;
	axis->dnum = -.5;
    }
	return;
    }
    div = temp / (inch);
    div1 = log10 (div);
    div2 = div1;
    if (div < 1.)
	div2--;
    b = div / pow (10., (double) div2);
    for (i = 0; i < 3 && b >= sqr[i]; i++);
    tempdnum = vint[i] * pow (10., (double) div2);
    temp4 = tempdnum;
    length = 0;
    counter = (int) (fabs((max - min)) / tempdnum);
    for (i= 0; i < counter; i++)
    {
        num = num0 + (i * tempdnum);
        temp5 = fabs(num);
        if (temp5 < (fabs((max - min)) / 10000))
             num = 0.;
	sprintf (string, "%1.5g", num);
	if (strlen (string) > length)
	    length = strlen (string);
    }
    temp1 = axis->labelsz / 33.;
    temp2 = inch / temp1;
    temp3 = (length + 1.5) * ((fabs((max - min)) / tempdnum));
    if (temp2 < temp3)
	tempdnum = tempdnum * 2.;
    if (tempdnum != temp4)
    {
	while (temp2 < temp3)
	{
        length = 0;
        counter = (int) (fabs((max - min)) / tempdnum);
        for (i= 0; i < counter; i++)
        {
            num = num0 + (i * tempdnum);
        temp5 = fabs(num);
        if (temp5 < (fabs((max - min)) / 10000))
             num = 0.;
	    sprintf (string, "%1.5g", num);
	    if (strlen (string) > length)
	       length = strlen (string);
        }
	    temp1 = axis->labelsz / 33.;
	    temp2 = inch / temp1;
	    temp3 = (length + 1) * ((fabs((max - min)) / tempdnum));
	    if (temp2 < temp3)
		tempdnum = tempdnum * 2.;
	}
    }
	axis->dnum = tempdnum;

	return 0;
}
/*
*	gl_gridint(xgrid, ygrid, grid, g1num, g2num, d1num, d2num, fastplt)
*
*        Input and initialization of variables needed to draw a grid 
*
*
* 	Biondo: 12/1996: Getch gridfat from command line to initialize fat
*
* end self-doc
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_gridint(struct gridinfo *gridin, struct coordinfo *coordinate, 
	   struct axisinfo *axis1, struct axisinfo *axis2)
#else
int gl_gridint(gridin, coordinate, axis1, axis2)
        struct coordinfo *coordinate;
	struct gridinfo *gridin;
	struct axisinfo *axis1;
	struct axisinfo *axis2;
#endif
{
    int grid1, i;
    float tempnum;
    int tempgrid;
    int grid2;
    int grid0;
    int gridfat;
    int gridcol;
    float g1num, g2num;
    grid0 = 0;
    grid1 = 0;
    grid2 = 0;
		if(coordinate->transp==0){
    g1num = axis1->dnum;
    g2num = axis2->dnum;
		}
		else{
    g1num = axis2->dnum;
    g2num = axis1->dnum;
		}
    tempgrid = getch ("grid", "1", &grid0);
    if (tempgrid == 0) 
    {
	getch ("grid1", "1", &grid1);
        gridin->grid1 = grid1;
	getch ("grid2", "1", &grid2);
        gridin->grid2 = grid2;
    }
   else
   {
    if (grid0 == 1)
    {
	gridin->grid1 = 1;
	gridin->grid2 = 1;
        grid1=1;
        grid2=1;
    }
    else
    if (grid0 == 0)
    {
	gridin->grid1 = 0;
	gridin->grid2 = 0;
        grid1=0;
        grid2=0;
    }
   }
       getch ("g1num", "f", &g1num);
       gridin->g1num=g1num;
       getch ("g2num", "f", &g2num);
       gridin->g2num=g2num;
   if (g1num == 0. )
	if (grid1)
		seperr("g1num is 0 needs to be set");
   if (g2num == 0. )
	if (grid2)
		seperr("g2num is 0 needs to be set");
   if (getch("gridcol", "d", &gridcol) == 0 )
     for ( i = 0; i <= NPMAX; i++ )
        gridin->col[i]=axis1->col[i];
   else
     for ( i = 0; i <= NPMAX; i++ )
        gridin->col[i]=gridcol;
   if (getch("gridfat", "d", &gridfat) == 0 )
     gridin->fat=1;
   else
     gridin->fat=gridfat;
   if (coordinate->transp)
   {
    tempgrid = gridin->grid1; 
    gridin->grid1 = gridin->grid2;
    gridin->grid2 = tempgrid;
   }


	return 0;
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_gtext (float x, float y, float xpath, float ypath, float xup, float yup, char
* string, char* labelpos)
#else
int gl_gtext (x, y, xpath, ypath, xup, yup, string, labelpos)
    float           x, y, xpath, ypath, xup, yup;
    char            *string, *labelpos;
#endif
{
float           temp, temp1, temp2;
        gl_tjust(labelpos);
	return(vp_gtext (x, y, xpath, ypath, xup, yup, string));
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_invmassage (float* min, float* max, float mid, float dev)
#else
int gl_invmassage (min, max, mid, dev)
    float           *min, *max,  mid, dev;
#endif
{
double          temp1num, temp2num, temp3num, tempnum, tempmin, tempmax;
    tempmin = *min;
    tempmax = *max;
    temp1num = (tempmin * dev) + mid;
    temp2num = (tempmax * dev) + mid;
    if (temp1num == temp2num) 
     {
        temp1num = temp1num - temp1num ;
        temp2num = temp2num + temp2num ;
     }
    *min = temp1num;
    *max = temp2num;
	return 0;
    
}
/*
*
*
*This routine will plot a label in what position it is told to
*
*
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_labelaxis(struct coordinfo *coord, struct axisinfo *axis)
#else
int gl_labelaxis(coord, axis)
    struct coordinfo *coord;
    struct axisinfo *axis;
#endif
{
float           vs, ch, xc, yc, position1, position2;
float           x, y, xup, yup, xpath, ypath;

    ch = axis->labelsz / 33.;
    vs = 1.5 * ch * 5. / 10.;
	switch (*(axis->wherelabel))
	{	
	case 't':
	    position1 = (coord->max1 + coord->min1) / 2;
	    position2 = coord->max2;
	    gl_umove (position1, position2);
	    gl_where (&xc, &yc);
	    y = yc + (ch) + (2 *  vs);
	    x = xc;
	    xpath = ch ;
	    ypath = 0.;
	    yup = ch;
	    xup = 0.;
	    break;
	case 'b':
	    position1 = (coord->max1 + coord->min1) / 2;
	    position2 = coord->min2;
	    gl_umove (position1, position2);
	    gl_where (&xc, &yc);
	    y = yc - ch - (2 * vs);  
	    x = xc;
	    xpath = ch ;
	    yup = ch;
	    ypath = 0.;
	    xup = 0.;
	    break;
	case 'l':
	    position1 = coord->min1;
	    position2 = (coord->min2 + coord->max2) / 2.;
	    gl_umove (position1, position2);
	    gl_where (&xc, &yc);
	    x = xc - ch - (2 *  vs);
	    y = yc;
            if (coord->labelrot == -1)
             {
               x = x -(ch  + vs) ; 
             }
	    xpath = 0.;
	    ypath = coord->labelrot * ch;
	    xup = coord->labelrot * -ch;
	    yup = 0.;
	    break;
	case 'r':
	    position1 = coord->max1;
	    position2 = (coord->min2 + coord->max2) / 2.;
	    gl_umove (position1, position2);
	    gl_where (&xc, &yc);
	    x = xc + ((1 * ch) + (2 * vs));
	    y = yc;
            if (coord->labelrot == -1)
             {
               x = x + (ch  + vs) ; 
             }
	    xpath = 0.;
	    ypath = coord->labelrot * ch ;
	    xup = coord->labelrot * -ch;
	    yup = 0.;
	    break;
	}
        vp_fat (axis->labelfat);
        gl_gtext(x,y,xpath,ypath,xup,yup,axis->label,axis->wherelabel );
	return 0;
}
/*
*	gl_labeltic(num0, ymin, ymax, xmin, xmax, dnum, wherenumber, ch, vs, 
*	            axisorig, labelsz, fastplt, whichaxis)
*
*	This routine will label the tics
*
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_labeltic (struct coordinfo *coord, struct axisinfo *axis)
#else
int gl_labeltic (coord, axis)
    struct coordinfo *coord;
    struct axisinfo *axis;
#endif
{
char            string[10];
int             temp1; 
float           nummax, num, xc, yc, min1, min2, max1, max2, point1, point2, xmax;
float           x, y, xup, yup, xpath, ypath, temp, num1, ch, vs, num0, dnum;

    ch = axis->labelsz / 33.;
    vs = 1.5 * (ch * 5. / 10.);
    vp_fat (axis->labelfat);
    num0 = axis->num0;
    dnum = axis->dnum;
    xmax = coord->max1;
    if (*axis->wherelabel == 'l' || *axis->wherelabel == 'r')
    {
	min1 = coord->min2;
	max1 = coord->max2;
	min2 = coord->min1;
	max2 = coord->max1;
    }
    if (*axis->wherelabel == 't' || *axis->wherelabel == 'b')
    {
	max1 = coord->max1;
	min1 = coord->min1;
	max2 = coord->max2;
	min2 = coord->min2;
    }
    if (max1 < min1)
    {
	temp = min1;
	min1 = max1;
	max1 = temp;
    }
    if (num0 >= max1)
	num0 = min1 + (max1 - min1) / 2.;


    if ((coord->xreverse) && (*axis->wherelabel == 't' || *axis->wherelabel == 'b'))
    {
        nummax = num0;
        num0 = min1 + max1 - num0;
	for (num = num0; num >= min1; num -= dnum)
         {
	    gl_makelabel (num, nummax, coord, axis); 
            nummax = nummax +  dnum ;  
         }

    }
    else 
    if ((coord->yreverse) && (*axis->wherelabel == 'l' || *axis->wherelabel == 'r'))
    {
        nummax = num0;
        num0 = min1 + max1 - num0;
        temp1 = 0; 
	for (num = num0; num >= min1; num -= dnum)
        {
	    gl_makelabel (num, nummax, coord, axis);
            nummax = nummax +  dnum ;  
    }
         }
    else
    {
	for (num = num0; num <= max1; num += dnum)
         {
	    gl_makelabel (num,  num, coord, axis);
          }
    }
	return 0;
}


int gl_makelabel (num, num2, coord, axis)
    float num, num2;
    struct coordinfo *coord;
    struct axisinfo *axis;
{
char            string[10];
float           xc, yc, min1, min2, max1, max2, point1, point2, xmax;
float           x, y, xup, yup, xpath, ypath, temp, num1, ch, vs ;
    ch = axis->labelsz / 33.;
    vs = 1.5 * (ch * 5. / 10.);
    if (*axis->wherelabel == 'l' || *axis->wherelabel == 'r')
    {
	min1 = coord->min2;
	max1 = coord->max2;
	min2 = coord->min1;
	max2 = coord->max1;
    }
    if (*axis->wherelabel == 't' || *axis->wherelabel == 'b')
    {
	max1 = coord->max1;
	min1 = coord->min1;
	max2 = coord->max2;
	min2 = coord->min2;
    }


    if (max1 < min1)
    {
	temp = min1;
	min1 = max1;
	max1 = temp;
    }



    if ((coord->xreverse)  && (*axis->wherelabel == 't' || *axis->wherelabel == 'b'))
    {
	temp = min1;
	min1 = max1;
	max1 = temp;

    }
    else 
    if ((coord->yreverse) && (*axis->wherelabel == 'l' || *axis->wherelabel == 'r'))
    {
	temp = min1;
	min1 = max1;
	max1 = temp;
    }



    if (fabs (num) < ((max1 - min1) / 10000))
    {
	num = 0.0;
    }
    if (*axis->wheretics != 'a')
    {
	switch (*axis->wherelabel)
	{
	case 't':
	    point1 = num;
	    point2 = max2;
	    break;
	case 'b':
	    point1 = num;
	    point2 = min2;
	    break;
	case 'l':
	    point1 = min2;
	    point2 = num;
	    break;
	case 'r':
	    point1 = max2;
	    point2 = num;
	    break;
	}
	gl_umove (point1, point2);
	gl_where (&xc, &yc);

    }
    else
    {
	if (*axis->wherelabel == 't' || *axis->wherelabel == 'b')
	{
	    point1 = num;
	    point2 = axis->axisor;
	}
	if (*axis->wherelabel == 'l' || *axis->wherelabel == 'r')
	{
	    point1 = axis->axisor;
	    point2 = num;
	}
	gl_umove (point1, point2);
	gl_where (&xc, &yc);
    }
	num1 = num2;
    if (num1 < 0)
    {
	num1 = -1 * num1;
	if (num1 < .0000001)
	{
	    num1 = 0.0;
	}
	else
	    num1 = num2;
    }
    else
    {
	if (num1 < .0000001)
	{
	    num1 = 0.0;
	}
	else
	    num1 = num2;
    }
    sprintf (string, "%1.5g", num1);
    switch (*axis->wherelabel)
    {
    case 'l':
	y = yc;
	x = xc - vs;
        if (coord->labelrot == -1)
        {
            x = x - ch - vs;
        }
	xpath = 0.;
	ypath =  coord->labelrot * ch;
	xup =  coord->labelrot * -ch;
	yup = 0;
	break;
    case 'r':
	y = yc;
	x = xc + vs;
        if (coord->labelrot == -1)
        {
            x = x + ch + vs;
        }
	xpath = 0.;
	ypath = coord->labelrot * ch;
	xup = coord->labelrot * -ch;
	yup = 0;
	break;
    case 't':
	y = yc + vs;
	x = xc;
	yup = ch;
	xpath = ch;
	ypath = 0.;
	xup = 0.;
	break;
    case 'b':
	y = yc - vs;
	x = xc;
	yup = ch;
	xup = 0.;
	xpath = ch;
	ypath = 0.;
	break;
    }
    gl_gtext (x, y, xpath, ypath, xup, yup, string, axis->wherelabel);

	return 0;

}
/*
 * Modifications:
 * W. Bauske IBM 04-01-91
 *	Removed re-declare of malloc()
 */
#include "glplot.h"
#include <stdio.h>
#if defined(HAVE_STDLIB_H)
#include <stdlib.h>
#else
extern char    *malloc ();
#endif
#ifdef USE_PROTO
int gl_massage (float *min, float *max, float *mid, float *dev)
#else
int gl_massage (min, max, mid, dev)
    float          *min, *max, *mid, *dev;
#endif

{
double          xtemp[2], mintemp, maxtemp, midtemp, nextnum;
double          devtemp;
int             ii, jj, kk;

/* first calculate max and min  of original data*/
    mintemp = *min;
    maxtemp = *max;
    midtemp = (mintemp + maxtemp) / 2.0;
    devtemp = 0.;
    xtemp[0] = mintemp - midtemp;
    xtemp[1] = maxtemp - midtemp;
    for (ii = 0; ii < 2; ii++)
    {
	devtemp = (fabs (xtemp[ii]) > devtemp) ? fabs (xtemp[ii]) : devtemp;
    }
    if (devtemp != 0.)
    {
	*min = (mintemp - midtemp) / devtemp;
	*max = (maxtemp - midtemp) / devtemp;
	*mid = midtemp;
	*dev = devtemp;
    }
    else
    {
	*min = mintemp - midtemp ;
	*max = mintemp - midtemp ;
	*mid = midtemp;
	*dev = devtemp;
    }

	return 0;
}
/*  fetch min and maxes   will also take care of transposing*/ 
#include "glplot.h"
#ifdef USE_PROTO
int gl_minmax (struct coordinfo *coordinate)
#else
int gl_minmax (coordinate)
    struct         coordinfo *coordinate;
#endif
{
int   temp1;
float temp;
    coordinate->fmin1 = 0; 
    coordinate->fmax1 = 0; 
    coordinate->fmin2 = 0; 
    coordinate->fmax2 = 0; 
    if (getch ("min1", "f", &coordinate->min1))
    coordinate->fmin1 = 1; 
    if (getch ("max1", "f", &coordinate->max1))
    coordinate->fmax1 = 1; 
    if (getch ("min2", "f", &coordinate->min2))
    coordinate->fmin2 = 1; 
    if (getch ("max2", "f", &coordinate->max2))
    coordinate->fmax2 = 1; 


    if (coordinate->transp)
    {
      temp = coordinate->min1;
      coordinate->min1 = coordinate->min2;
      coordinate->min2 = temp;
      temp = coordinate->max1;
      coordinate->max1 = coordinate->max2;
      coordinate->max2 = temp;
      temp1 = coordinate->fmin1;
      coordinate->fmin1 = coordinate->fmin2;
      coordinate->fmin2 = temp1;
      temp1 = coordinate->fmax1;
      coordinate->fmax1 = coordinate->fmax2;
      coordinate->fmax2 = temp;
    }

	return 0;
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_move(float point1, float point2)
#else
int gl_move(point1, point2)
float point1, point2;
#endif
{
	return(vp_move(point1, point2));
}

/*
 * gl_opttic(min,max,inch,num0,dnum,dnumtic,labelsz)
 * determines the optimum tic mark interval
 *
 * inputs:
 * outputs:
 *
 * Author - S. Cole   11 Dec92
 *  Code taken from gl_getscl and made into standalone routine.
 */
#include <stdio.h>
#include <math.h>
static float    sqr[3] = {1.414214, 3.162278, 7.071068};
static float    vint[4] = {1., 2., 5., 10.};
#include "glplot.h"
#ifdef USE_PROTO
int gl_opttic(float min, float max, float inch, float num0, float *dnumtic, float labelsz)
#else
int gl_opttic(min, max, inch, num0, dnumtic, labelsz)
float min, max, inch, num0;
float *dnumtic;
float labelsz;
#endif
{
float           temp, temp1, temp2, temp3, temp4, temp5, tempdnum, num;
double          div, div1, div3, a, b ;
int             div2, i, length, counter;
char            string[100];

    if (fabs(max - min) < .000001)
    {
	*dnumtic = .5;
	return;
    }
    if (max > min)
	temp = max - min;
    else
	temp = min - max;
    if ((fabs(temp)/inch) < .0000001)
    {
    if (max > min)
    {
	temp = max - min;
	*dnumtic = .5;
    }
    else
    {
	temp = min - max;
	*dnumtic = -.5;
    }
	return;
    }
    div = temp / (inch);
    div1 = log10 (div);
    div2 = div1;
    if (div < 1.)
	div2--;
    b = div / pow (10., (double) div2);
    for (i = 0; i < 3 && b >= sqr[i]; i++);
    tempdnum = vint[i] * pow (10., (double) div2);
    temp4 = tempdnum;
    length = 0;
    counter = (int) (fabs((max - min)) / tempdnum);
    for (i= 0; i < counter; i++)
    {
        num = num0 + (i * tempdnum);
        temp5 = fabs(num);
        if (temp5 < (fabs((max - min)) / 10000))
             num = 0.;
	sprintf (string, "%1.5g", num);
	if (strlen (string) > length)
	    length = strlen (string);
    }
    temp1 = labelsz / 33.;
    temp2 = inch / temp1;
    temp3 = (length + 1.5) * ((fabs((max - min)) / tempdnum));
    if (temp2 < temp3)
	tempdnum = tempdnum * 2.;
    if (tempdnum != temp4)
    {
	while (temp2 < temp3)
	{
        length = 0;
        counter = (int) (fabs((max - min)) / tempdnum);
        for (i= 0; i < counter; i++)
        {
            num = num0 + (i * tempdnum);
        temp5 = fabs(num);
        if (temp5 < (fabs((max - min)) / 10000))
             num = 0.;
	    sprintf (string, "%1.5g", num);
	    if (strlen (string) > length)
	       length = strlen (string);
        }
	    temp1 = labelsz / 33.;
	    temp2 = inch / temp1;
	    temp3 = (length + 1) * ((fabs((max - min)) / tempdnum));
	    if (temp2 < temp3)
		tempdnum = tempdnum * 2.;
	}
    }
	*dnumtic = tempdnum;

	return 0;
}
/*
* This routine will fecth the min,and max values  and pad. 
*/
/* 
3-6-90  fixed padding problem
        now do calculations on normalized data
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_padint(struct coordinfo *coordinate)
#else
int gl_padint(coordinate)
struct coordinfo *coordinate;
#endif
{
float           mid1, mid2, dev1, dev2;
int            pad, npad;
    if (coordinate->pad) 
    {
        gl_massage(&coordinate->min1, &coordinate->max1, &mid1, &dev1) ;
        if (coordinate->fmin1)
        {
        if (coordinate->npad)
        coordinate->min1 = 1.04 * coordinate->min1;
        }
        else
        coordinate->min1 = 1.04 * coordinate->min1;
        if (coordinate->fmax1)
        {
        if (coordinate->npad)
        coordinate->max1 = 1.04 * coordinate->max1;
        }
        else
        coordinate->max1 = 1.04 * coordinate->max1;
        gl_invmassage(&coordinate->min1, &coordinate->max1,  mid1, dev1) ;
        gl_massage(&coordinate->min2, &coordinate->max2, &mid2, &dev2) ;
        if (coordinate->fmin2)
        {
        if (coordinate->npad)
        coordinate->min2 = 1.04 * coordinate->min2;
        }
        else
        coordinate->min2 = 1.04 * coordinate->min2;
        if (coordinate->fmax2)
        {
        if (coordinate->npad)
        coordinate->max2 = 1.04 * coordinate->max2;
        }
        else
        coordinate->max2 = 1.04 * coordinate->max2;
        gl_invmassage(&coordinate->min2, &coordinate->max2, mid2, dev2) ;
    }
    if (coordinate->min1 == 0 && coordinate->max1 == 0)
    {
    coordinate->min1 = -1.;
    coordinate->max1 = 1.;
    }
    if (coordinate->min2 == 0 && coordinate->max2 == 0)
    {
    coordinate->min2 = -1.;
    coordinate->max2 = 1.;
    }
    if (coordinate->min1 == coordinate->max1)
    {
        coordinate->max1 = coordinate->max1 * 1.04;
        coordinate->min1 = coordinate->min1 * .96;
    }
    if (coordinate->min2 == coordinate->max2)
    {
        coordinate->max2 = coordinate->max2 * 1.04;
        coordinate->min2 = coordinate->min2 * .96;
    }

	return 0;
}
#include "glplot.h"
int gl_penup()
{
return(vp_penup());
}
/*
*
*
*    This routine will plot the axes
*
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_plotaxis(struct axisinfo *axis, struct coordinfo *coord, int counter )
#else
int  gl_plotaxis(axis, coord, counter )
    struct axisinfo *axis;
    struct coordinfo *coord;
#endif
{
        gl_fat(axis->fat[counter]);
    	if (*axis->wherelabel == 't' || *axis->wherelabel == 'b')
    	{
		gl_umove(coord->min1, axis->axisor);
		gl_udraw(coord->max1, axis->axisor);
    	}

    	if (*axis->wherelabel == 'r' || *axis->wherelabel == 'l')
    	{
		gl_umove(axis->axisor, coord->min2);
		gl_udraw(axis->axisor, coord->max2);
	}
	return 0;
} 
/*
*
*
*This routine will plot the frame around the picture
*
*
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_plotframe(struct coordinfo *coord, int col)
#else
int gl_plotframe(coord, col)
struct coordinfo *coord;
int col;
#endif
{
    gl_color(col);
    gl_umove(coord->min1, coord->min2);
    gl_udraw(coord->min1, coord->max2);
    gl_udraw(coord->max1, coord->max2);
    gl_udraw(coord->max1, coord->min2);
    gl_udraw(coord->min1, coord->min2);
	return 0;
}
/*
*
*
*   This routine will plot a grid for one axis depending on which one is 
*	specified
*
*
* 	Biondo: 12/1996: use fat
*
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_plotgrid(struct coordinfo *coordinate, struct axisinfo *axis, struct gridinfo
 *grid, int counter)
#else
int gl_plotgrid(coordinate, axis, grid, counter)
struct coordinfo *coordinate;
struct axisinfo *axis;
struct gridinfo *grid;
int counter;
#endif
{
float           pos1, pos2, pos3, pos4, num, xmax; 
float           gnum, min1, min2, max1, max2, num0;
       min1 = coordinate->min1;
       max1 = coordinate->max1;
       min2 = coordinate->min2;
       max2 = coordinate->max2; 
       num0 = axis->num0;
    if (*axis->wherelabel == 't' || *axis->wherelabel == 'b')
    { 
       min1 = coordinate->min1;
       max1 = coordinate->max1;
       min2 = coordinate->min2;
       max2 = coordinate->max2; 
       num0 = axis->num0;
       gnum = grid->g1num;
			 if(coordinate->transp==0) gnum = grid->g1num;
			 else gnum=grid->g2num;
    }
    if (*axis->wherelabel == 'l' || *axis->wherelabel == 'r')
    {
       min1 = coordinate->min2;
       max1 = coordinate->max2;
       min2 = coordinate->min1;
       max2 = coordinate->max1;
       num0 = axis->num0; 
			 if(coordinate->transp==0) gnum = grid->g2num;
			 else gnum=grid->g1num;
    }
   gl_fat(grid->fat);
    if (( grid->col[counter] >= 0 ) && ( grid->col[counter] <= 7 ))
    gl_color(grid->col[counter]);
    else 
    gl_color(axis->col[counter]);
   if ( num0 >= max1 )
    num0 = min1 + (max1 - min1 ) / 2.;
    for (num = num0; num <= max1; num += gnum)
    {
	if (fabs(num) < ((max1 - min1) / 10000))
	{
	    num = 0.0;
	}
	if (*axis->wherelabel == 't' || *axis->wherelabel == 'b')
	{
	    pos1 = num;
	    pos2 = min2;
	    pos3 = num;
	    pos4 = max2;
      xmax = max1;
			if(coordinate->xreverse) gl_rotate1(&pos1,min1,max1);
	    pos3 = pos1;
	}
	if (*axis->wherelabel == 'l' || *axis->wherelabel == 'r')
	{
	    pos1 = min2;
	    pos2 = num;
	    pos3 = max2;
      xmax = max2;
			/* time to rotate */
			if(coordinate->yreverse) gl_rotate1(&pos2,min1,max1);
	    pos4 = pos2;
	}
	gl_umove(pos1, pos2);
	gl_udraw(pos3, pos4);
    }
	return 0;
}
/*
*
*
* The routine will fetch and initialize plot and axis parameters.
*
* 
*/
/*
 * Edited Nov 12 1990 by Joe Dellinger
 * 	plotcol=0 is perfectly legitimate, not an error, so allow it.
 */
#include "glplot.h"
#define NOT_SET_YET	-1
#ifdef USE_PROTO
int gl_plotint(struct plotinfo *plot, struct dashinfo *dash)
#else
int gl_plotint(plot, dash)
struct plotinfo *plot;
struct dashinfo *dash;
#endif
{
int             colnum, i, j, length, k, nplotfat, nsymbolsz, ndashtype;
char            symbol[NPMAX];
    k = 0;
    for (i = 0; i < NPMAX; i++)
    {
	plot->symbol[i] = ' ';
	plot->symbolsz[i] = 2;
	plot->col[i] = NOT_SET_YET;
	plot->fat[i] = 0;
	dash->dashtype[i] = 0;
    }
    ndashtype = getch ("dash", "f", dash->dashtype); 
    if (ndashtype)
    {
	for (i = ndashtype; i < NPMAX; i++)
	{
	    dash->dashtype[i] = dash->dashtype[(i % ndashtype)];
	}
    }
    getch ("symbol", "s", plot->symbol);
    nplotfat = getch ("plotfat", "d", plot->fat);
    if(nplotfat) 
    {
	for (i = nplotfat; i < NPMAX; i++)
	{
	    plot->fat[i] = plot->fat[(i % nplotfat)];
	}
    }
    nsymbolsz = getch ("symbolsz", "d", plot->symbolsz);
    if (nsymbolsz)
    {
	for (i = nsymbolsz; i < NPMAX; i++)
	{
	    plot->symbolsz[i] = plot->symbolsz[(i % nsymbolsz)];
	}
    }
    colnum =  getch ("plotcol", "d", plot->col);
    if (colnum == 0)
    {
	for (i = 0; i < NPMAX; i++)
	{
	    plot->col[i] = 6 - (i % 6);
	}
    }
    else
    {
        if (colnum > NPMAX)
           seperr("ENTERED TOO MANY VALUES FOR PLOTCOL\n");
	while (plot->col[k] != NOT_SET_YET)
	    k++;
	for (i = k; i < NPMAX; i++)
	{
	    plot->col[i] = plot->col[(i % k)];
	}
    }
    length = strlen (plot->symbol);
    for (i = length; i < NPMAX; i = i + length)
    {
	for (j = 0; j <= length; j++)
	{
	    plot->symbol[i + j] = plot->symbol[j];
	}
    } 
	return 0;
}
/*
*
*
*This routine will set the the device parameters
* namely color of the background, plot color and plot fat
*
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_plotpram(struct colorinfo *color,struct coordinfo * coordinate)
#else
int gl_plotpram(color, coordinate)
struct colorinfo *color;
struct coordinfo *coordinate;
#endif
{
        if (color->fillcol[0] !=color->backcol[0] || color->fillcol[1] != color->backcol[1] || color->fillcol[2] != color->backcol[2] )
            gl_fillin(coordinate, color);
            vp_coltab(0, color->backcol[0], color->backcol[1], color->backcol[2]);
	return 0;
}
/*
*
*
*This routine will plot the tics for the axes.  
*
*
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_plottic ( struct coordinfo *coord,struct axisinfo * axis, int counter)
#else
int gl_plottic (coord, axis, counter)
    struct coordinfo *coord;
    struct axisinfo *axis;
    int             counter;
#endif
{
float           xc, yc, position1, position2, position3, position4, xmax, min1, min2, max1, max2, num, tic, num0, dtic, axisor, dnum, vs;
int             length, i, temp2;
float           temp1, ch;
char            string[10];

    ch = (axis->labelsz / 33.);
    vs = ch * 5. / 10.;
    num0 = axis->num0;
    dnum = axis->dnum;
    vp_fat (axis->fat[counter]);
    if (*axis->wherelabel == 'l' || *axis->wherelabel == 'r')
    {
	min1 = coord->min2;
	max1 = coord->max2;
	min2 = coord->min1;
	max2 = coord->max1;
	xmax = coord->max1;
    }
    if (*axis->wherelabel == 't' || *axis->wherelabel == 'b')
    {
	min1 = coord->min1;
	max1 = coord->max1;
	min2 = coord->min2;
	max2 = coord->max2;
	xmax = coord->max1;
    }

    if (max1 < min1)
    {
	temp1 = min1;
	min1 = max1;
	max1 = temp1;
    } 
    if (num0 >= max1)
	num0 = min1 + (max1 - min1) / 2.; 

    if ((coord->xreverse) && (*axis->wherelabel == 't' || *axis->wherelabel == 'b'))
    {
	temp1 = min1;
	min1 = max1;
	max1 = temp1;
	dnum = dnum;
        num0 = min1 + max1 - num0;

    }
    if ((coord->yreverse) && (*axis->wherelabel == 'l' || *axis->wherelabel == 'r'))
    {
        num0 = min1 + max1 - num0;
	temp1 = min1;
	min1 = max1;
	max1 = temp1;
	dnum = dnum;
    }
    if ((coord->yreverse) && (*axis->wherelabel == 'l' || *axis->wherelabel == 'r'))
    {
	for (num = num0; num >= max1; num -= dnum)
	    gl_maketic (num, coord, axis);
    }
    else
    if ((coord->xreverse) && (*axis->wherelabel == 't' || *axis->wherelabel == 'b'))
    {
	for (num = num0; num >= max1; num -= dnum)
	    gl_maketic (num, coord, axis);
    }
    else
    {
	for (num = num0; num <= max1; num += dnum)
	    gl_maketic (num, coord, axis);
    }
	return 0;
}


/**
*
*This routine will plot the tics for the axes.  
*
*
*/
gl_maketic (num, coord, axis)
    float           num;
    struct axisinfo *axis;
    struct coordinfo *coord;
{
float           xc, yc, position1, position2, position3, position4, xmax, min1, max1, min2, max2, tic, num0, dtic, axisor, dnum, vs;
int             length, i;
float           temp1, ch;
char            string[10];

    ch = (axis->labelsz / 33.);
    vs = ch * 5. / 10.;
    num0 = axis->num0;
    dnum = axis->dnum;
    if (*axis->wherelabel == 'l' || *axis->wherelabel == 'r')
    {
	min1 = coord->min2;
	min2 = coord->min1;
	max1 = coord->max2;
	max2 = coord->max1;
	xmax = coord->max1;
    }
    if (*axis->wherelabel == 't' || *axis->wherelabel == 'b')
    {
	min1 = coord->min1;
	min2 = coord->min2;
	max1 = coord->max1;
	max2 = coord->max2;
	xmax = coord->max1;
    }
    if (fabs (num) < ((max1 - min1) / 10000))
    {
	num = 0.0;
    }
    if (*axis->wheretics != 'a')
    {
	switch (*axis->wherelabel)
	{
	case 't':
	    position1 = num;
	    position2 = max2;
	    break;
	case 'b':
	    position1 = num;
	    position2 = min2;
	    break;
	case 'r':
	    position1 = max2;
	    position2 = num;
	    break;
	case 'l':
	    position1 = min2;
	    position2 = num;
	    break;
	}
	gl_umove (position1, position2);
	gl_where (&xc, &yc);
	switch (*axis->wherelabel)
	{
	case 'b':
	    position3 = xc;
	    position4 = yc - vs;
	    break;
	case 't':
	    position3 = xc;
	    position4 = yc + vs;
	    break;
	case 'l':
	    position3 = xc - vs;
	    position4 = yc;
	    break;
	case 'r':
	    position3 = xc + vs;
	    position4 = yc;
	    break;
	}
	gl_draw (position3, position4);
    }
    else
    {
	switch (*axis->wherelabel)
	{
	case 'b':
	    position1 = num;
	    position2 = axis->axisor;
	    break;
	case 't':
	    position1 = num;
	    position2 = axis->axisor;
	    break;
	case 'l':
	    position1 = axis->axisor;
	    position2 = num;
	    break;
	case 'r':
	    position1 = axis->axisor;
	    position2 = num;
	    break;
	}
	gl_umove (position1, position2);
	gl_where (&xc, &yc);
	switch (*axis->wherelabel)
	{
	case 'b':
	    position3 = xc;
	    position4 = yc - vs;
	    break;
	case 't':
	    position3 = xc;
	    position4 = yc + vs;
	    break;
	case 'l':
	    position3 = xc - vs;
	    position4 = yc;
	    break;
	case 'r':
	    position3 = xc + vs;
	    position4 = yc;
	    break;
	}
	gl_draw (position3, position4);
    }
    if (axis->dtic != dnum)
    {
	for (tic = axis->num0; tic <= max1; tic += axis->dtic)
	{
	    switch (*axis->wherelabel)
	    {
	    case 'b':
		position1 = tic;
		position2 = min2;
		break;
	    case 't':
		position1 = tic;
		position2 = max2;
		break;
	    case 'l':
		position1 = min2;
		position2 = tic;
		break;
	    case 'r':
		position1 = max2;
		position2 = tic;
		break;
	    }
	    gl_umove (position1, position2);
	    gl_where (&xc, &yc);
	    if (axis->axisor == min2)
	    {
		switch (*axis->wherelabel)
		{
		case 'b':
		    position3 = xc;
		    position4 = yc - vs / 2.;
		    break;
		case 't':
		    position3 = xc;
		    position4 = yc + vs / 2.;
		    break;
		case 'l':
		    position3 = xc - vs / 2.;
		    position4 = yc;
		    break;
		case 'r':
		    position3 = xc - vs / 2.;
		    position4 = yc;
		    break;
		}
	    }
	    else
	    {
		switch (*axis->wherelabel)
		{
		case 'b':
		    position3 = xc;
		    position4 = yc - vs / 4.;
		    break;
		case 't':
		    position3 = xc;
		    position4 = yc + vs / 4.;
		    break;
		case 'l':
		    position3 = xc - vs / 4.;
		    position4 = yc;
		    break;
		case 'r':
		    position3 = xc - vs / 4.;
		    position4 = yc;
		    break;
		}

	    }
	    gl_draw (position3, position4);
	}
    }
}
/*
*
*
*This routine will plot the plot's title
*
*
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_plottitle(struct coordinfo *coord, struct titleinfo *title, struct axisinfo *
axis, int counter)
#else
int gl_plottitle(coord, title, axis, counter)
struct coordinfo *coord;
struct axisinfo *axis;
struct titleinfo *title;
int  counter;
#endif
{
int             temp;
float           labelvs, labelch, ch, vs, xc, yc, position1, position2;
float           x, y, xup, yup, xpath, ypath;
char            titletemp[280];
    ch = title->titlesz / 33.;
    vs = ch * 6. / 10.;
    temp = title->titlefat;
    gl_fat (temp);
    gl_color(axis->col[counter]);
    labelch = 0.;
    labelvs = 0.; 
    
    if ((*title->wheretitle == 't')  && ( *axis->wherelabel == 't' ))
     {
       labelch = axis->labelsz / 33.;
       labelvs = 1.5 * labelch * 5. / 10.;  
     }
    if ((*title->wheretitle == 'b')  && ( *axis->wherelabel == 'b' ))
     {
       labelch = axis->labelsz / 33.;
       labelvs = 1.5 *labelch * 5. / 10.;  
     }
    if ((*title->wheretitle == 'l')  && ( *axis->wherelabel == 'l' ))
     {
       labelch = axis->labelsz / 33.;
       labelvs = 1.5 *labelch * 5. / 10.;  
     }
    if ((*title->wheretitle == 'r')  && ( *axis->wherelabel == 'r' ))
     {
       labelch = axis->labelsz / 33.;
       labelvs = 1.5 *labelch * 5. / 10.;  
     }
     
	switch (*title->wheretitle)
	{
	case 't':
	    position1 = (coord->max1 + coord->min1) / 2;
	    position2 = coord->max2;
	    gl_umove (position1, position2);
	    gl_where (&xc, &yc);
	    x = xc;
	    y = yc  + vs + (1 * labelch)  +  (3 * labelvs);
	    ypath = 0.;
	    xpath = title->titlesz / 33;
	    yup = ch;
	    xup = 0.;
	    break;
	case 'b':
	    position1 = (coord->max1 + coord->min1) / 2;
	    position2 = coord->min2;
	    gl_umove (position1, position2);
	    gl_where (&xc, &yc);
            y = yc - vs - ((1 * labelch)  + (3 * labelvs));
	    x = xc;
            ypath = 0.;
	    xpath = ch;
	    yup = ch;
	    xup = 0.;
	    break;
	case 'l':
	    position1 = coord->min1;
	    position2 = (coord->min2 + coord->max2) / 2.;
	    gl_umove (position1, position2);
	    gl_where (&xc, &yc );
	    x = xc - vs - (1 * labelch) - (3 * labelvs);
	    y = yc;
            if (coord->labelrot == -1)
            {
              x = x - ch - vs - labelvs;
            }
              xpath = 0.;
	    ypath = coord->labelrot * ch;
	    xup = coord->labelrot * -ch;
	    yup = 0.;
	    break;
	case 'r':
	    position1 = coord->max1;
	    position2 = (coord->min2 + coord->max2) / 2.;
	    gl_umove (position1, position2);
	    gl_where (&xc, &yc);
	    x = xc + vs + ((1 * labelch) + ( 3 * labelvs));
	    y = yc;
            if (coord->labelrot == -1)
            {
              x = x  + ch + vs + labelvs;
            }
            xpath = 0.;
	    ypath = coord->labelrot * ch;
	    xup = coord->labelrot * -ch;
	    yup = 0.;
	    break;
	}
    gl_gtext(x, y, xpath, ypath, xup, yup ,title->title, title->wheretitle);
	return 0;
}
#include "glplot.h"
int gl_purge()
{
return(vp_purge());
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_rotate(float* x,float  min,float  max,struct datainfo * data)
#else
int gl_rotate(x, min, max, data)
float *x, min, max;
struct datainfo *data;
#endif
{
int i, j;
j = 0;
for ( i = 0; i < data->n2; i++ )
    j =  j + data->n1[i];
for (i  = 0; i < j; i++)
{
 x[i] = (min + max) - x[i];
} 
return 0;
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_rotate1(float *x, float min, float max)
#else
int gl_rotate1(x, min, max)
float *x, min, max;
#endif
{
float temp;
temp = *x;
*x = (min + max) - temp;
return 0;
}
/*
* This routine draws and annotates an axis. Any orientation.
* Axis goes from (x1,y1) to (x2,y2). Tic marks are perpendicular
* to it, on the right side as you go from point 1 to point 2.
*
* Author - S. Cole (SEP)   11 Dec 92
*/
#include "glplot.h"

#ifdef USE_PROTO
int gl_simpleaxis ( float x1, float y1, float x2, float y2, float num1, float num2,
  float dnum, float dnumticin, float ltic, char* label, float labelsz )
#else
int gl_simpleaxis ( x1, y1, x2, y2, num1, num2, dnum, dnumticin, ltic, label, labelsz )
float x1, x2, y1, y2, num1, num2, dnum, dnumticin, ltic;
char *label;
float labelsz;
#endif
{
float		ch, xpath, ypath, xup, yup;
char            string[10];
float		dist, costh, sinth, dtic, xpos, ypos;
float		dxtic, dytic, loc, num;
float		pad;
float		dnumtic;

    /* pad is the space skipped between tic and number, and number and label */
    pad = 0.15;
    ch = labelsz / 33.;

    gl_color (7);
    gl_move(x1,y1);
    gl_draw(x2,y2);

    /* compute sines and cosines of axis angle */
    dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    costh = (x2 - x1)/dist;
    sinth = (y2 - y1)/dist;
    /* x and y sizes of tic marks */
    dxtic = ltic * sinth;
    dytic = -ltic * costh;

    /* 
     * Figure out which quadrant we are in.
     * This determines which side of the line the tic
     * marks will go on, and how the text will be aligned.
     * For each quadrant we compute the vplot text path
     * and up vectors.
     */
    if (x1 <= x2 && y1 <= y2) {
        vp_tjust(TH_CENTER,TV_TOP);
        xpath = ch * costh;
        ypath = ch * sinth;
        xup = -ch * sinth;
        yup = ch * costh;
    } else if (x1 > x2 && y1 <= y2) {	
        vp_tjust(TH_CENTER,TV_BOTTOM);
        xpath = -ch * costh;
        ypath = -ch * sinth;
        xup = ch * sinth;
        yup = -ch * costh;
    } else if (x1 > x2 && y1 > y2) {	
        vp_tjust(TH_CENTER,TV_BOTTOM);
        xpath = -ch * costh;
        ypath = -ch * sinth;
        xup = ch * sinth;
        yup = -ch * costh;
    } else if (x1 <= x2 && y1 > y2) {	
        vp_tjust(TH_CENTER,TV_TOP);
        xpath = ch * costh;
        ypath = ch * sinth;
        xup = -ch * sinth;
        yup = ch * costh;
    }

    /* call gl_opttic to find the optimum tic mark interval */
    dnumtic = dnumticin;
    if (dnumtic == 0.) gl_opttic(num1,num2,dist,num1,&dnumtic,labelsz);
    if (num1 > num2) dnumtic *= -1.;

    /* figure out the tic mark spacing */
    dtic = dnumtic/(num2-num1) * dist; 
    if (dtic < 0.) dtic *= -1.;

    /*
     * move to each tic mark location, draw the tic and the number
     */
    for (loc=0., num=num1; loc < dist; loc+=dtic, num+=dnumtic) {
        sprintf(string,"%1.5g",num);
	xpos = x1 + loc * costh;
	ypos = y1 + loc * sinth;
        vp_move(xpos,ypos);
	vp_draw(xpos+dxtic,ypos+dytic);
	vp_gtext(xpos+dxtic+pad*sinth,ypos+dytic-pad*costh,xpath,ypath,xup,yup,string);
    }
    /* now the axis label */
    xpos = x1 + loc/2. * costh + dxtic + ch * sinth + 2. * pad * sinth;
    ypos = y1 + loc/2. * sinth + dytic - ch * costh - 2. * pad * costh;
    vp_gtext(xpos,ypos,xpath,ypath,xup*1.5,yup*1.5,label);

	return 0;
}
#include <glplot.h>
#ifdef USE_PROTO
int gl_stdplot (struct datainfo *data, struct coordinfo *coordinate,
	    struct axisinfo * axis1, struct axisinfo *axis2, 
	    struct gridinfo *grid, struct titleinfo *title,
 	    int counter, int fastplt, int wantframe, int wantframenum)
#else
int gl_stdplot (data, coordinate, axis1, axis2, grid, title, counter, fastplt, wantframe, wantframenum)
    struct datainfo *data;
    struct coordinfo *coordinate;
    struct axisinfo *axis1;
    struct axisinfo *axis2;
    struct gridinfo *grid;
    struct titleinfo *title;
    int             counter, fastplt;
    int             wantframe, wantframenum;
#endif
{

    /*
     * declaration of internal variables for Graph 
     */

    if (fastplt < 20)
    {

	/* initialization for plotting axes */

	gl_clip (-VP_MAX, -VP_MAX, VP_MAX, VP_MAX);
	gl_color (axis1->col[counter]);
	gl_fat (axis1->fat[counter]);

	/* If a frame is wanted draw it */
	if (wantframe)
	    if (fastplt < 16)
	    {
		gl_plotframe (coordinate, axis1->col[counter]);
	    }
	if (fastplt < 14)
	{

	    /* If  y axis is wanted draw it */

	    if (((wantframe == 0) || (axis1->axisor != coordinate->min1)) && (axis1->wantaxis))
		gl_plotaxis (axis1, coordinate, counter);
	    if ((axis2->wantaxis) && fastplt < 12 )
	    {
             if (axis2->dnum != 0.)
		{
                    if (axis2->ntic != 0)
		         gl_plottic (coordinate, axis2, counter);
		    gl_labeltic (coordinate, axis2);
                 }
		/* Label y axis */
		gl_labelaxis (coordinate, axis2);
	    }
	    /* If only x axis is wanted draw it */
	    if (((wantframe == 0) || (axis2->axisor != coordinate->min2)) && (axis2->wantaxis))
		gl_plotaxis (axis2, coordinate, counter);
	    if ((axis1->wantaxis) && (fastplt < 12))
	    {
                if ((axis1->dnum != 0.))
                {
		if (axis1->ntic != 0)
		    gl_plottic (coordinate, axis1, counter);
		gl_labeltic (coordinate, axis1);
		}
                /* Label x axis */
		gl_labelaxis (coordinate, axis1 );

	    }
	    if (fastplt < 8)
	    {
		if (grid->grid1)
		{
		    gl_plotgrid (coordinate, axis1, grid, counter);
		}
		if (grid->grid2)
		{
		    gl_plotgrid (coordinate, axis2, grid, counter);
		}
	    }

	}			/* end of fastplot 14 */
	/* label title */
	    if ( title->wanttitle  && (fastplt < 3))
        {
	if (*title->wheretitle == 't' || *title->wheretitle == 'b')
	    gl_plottitle (coordinate, title, axis1, counter);
	if (*title->wheretitle == 'l' || *title->wheretitle == 'r')
	    gl_plottitle (coordinate, title, axis2, counter);
        }
	/* label movie frame if needed */
	if (wantframenum > 0 && data->n3 > 1 && fastplt < 13)
	{
	    gl_framenum (counter, data->d3, data->o3, coordinate->min1, coordinate->min2, axis1->labelsz);

	}
    }				/* end of fastplot 20 */
	return 0;
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_tfont(int font, int prec, int ovly)
#else
int gl_tfont(font, prec, ovly)
int font, prec, ovly;
#endif
{
return(vp_tfont(font, prec, ovly));
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_titleint(struct titleinfo *title)
#else
int gl_titleint(title)
struct titleinfo *title;
#endif
{
getch ("wheretitle", "s", title->wheretitle);
title->wanttitle = 1;
getch ("wanttitle", "1", &title->wanttitle);
title->titlefat = 0;
fetch ("titlefat", "d", &title->titlefat);
title->titlesz = 10;
fetch ("titlesz", "f", &title->titlesz);
strcpy (title->title, " ");
if (!fetch ("title", "s", title->title))
   fetch ("in", "s", title->title);
	return 0;
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_tjust( char* where)
#else
int gl_tjust(where)
    char            *where;
#endif
{
    switch (*where)
    {
    case 'l':
	vp_tjust (TH_CENTER, TV_BOTTOM);
	break;
    case 'r':
	vp_tjust (TH_CENTER, TV_TOP);
	break;
    case 't':
	vp_tjust (TH_CENTER, TV_BOTTOM);
	break;
    case 'b':
	vp_tjust (TH_CENTER, TV_TOP);
	break;
    case 's':
	vp_tjust (TH_SYMBOL, TV_SYMBOL);
	break;
    }
	return 0;
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_transp (float *x, float* y, struct datainfo *data)
#else
int gl_transp (x, y, data)
    float          *x, *y;
    struct datainfo *data;
#endif
{
float          *xyexch;
int             i, tempalloc, tempval;

    tempalloc = 0;
    for (i = 0; i < data->n2; i++)
    {
	tempalloc = tempalloc + data->n1[i];
    }
    xyexch = (float *) calloc ((data->esize / 2) * (tempalloc), sizeof (float));
	tempval = 0;
    for (i = 0; i < data->n2; i++)
    {
	tempval = tempval + data->n1[i];
    }

    for (i = 0; i < tempval; i++)
    {
	xyexch[i] = x[i];
	x[i] = y[i];
	y[i] = xyexch[i];
    }

return 0;
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_uarea( float *px, float *py, int iii, int fatp, int ymask, int xmask, int fun
)
#else
int gl_uarea( px, py, iii, fatp, ymask, xmask, fun)
float *px, *py;
int iii, fatp, xmask, ymask, fun;
#endif
{
int j;
char string[10];
return(vp_uarea(px, py, iii, fatp, ymask, xmask));
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_uarrow (float x0, float y0, float x, float y, float r)
#else
int gl_uarrow (x0, y0, x, y, r)
    float           x0, y0, x, y, r;
#endif
{
return(	vp_uarrow (x0, y0, x, y, r));
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_uclip(float pos1, float pos2, float pos3, float pos4)
#else
int gl_uclip(pos1, pos2, pos3, pos4)
float pos1, pos2, pos3, pos4;
#endif
{
return(vp_uclip(pos1, pos2, pos3, pos4));
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_udraw(float point1, float point2)
#else
int gl_udraw(point1, point2)
float point1, point2;
#endif
{
return(vp_udraw(point1, point2));
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_umove(float point1, float point2)
#else
int gl_umove(point1, point2)
float point1, point2;
#endif
{
return(	vp_umove(point1, point2));
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_upendn(float point1, float point2)
#else
int gl_upendn(point1, point2)
float point1, point2;
#endif
{
return(	vp_upendn(point1, point2));
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_upmark(int npts, int mtype, int msize,float xp, float yp)
#else
int gl_upmark(npts, mtype, msize,xp, yp)
int npts, mtype, msize;
float xp, yp;
#endif
{
  return(  vp_upmark(npts, mtype, msize, &xp, &yp));
}
/*
*
*This routine initializes the device 
*
*/
#include "glplot.h"
#ifdef USE_PROTO
int gl_vplotint (struct plotposition *position, struct coordinfo *coord, 
	     struct axisinfo *axis1,struct axisinfo * axis2)
#else
int gl_vplotint (position, coord, axis1, axis2)
struct plotposition *position;
struct coordinfo *coord;
struct axisinfo *axis1;
struct axisinfo *axis2;
#endif
{
    float scale1, scale2, orig1, orig2, uorig1, uorig2;
    vp_style (STANDARD);
    if (coord->min2 == coord->max2)
	seperr ("min2 is equal to max2 change one of them\n");
    if (coord->min1 == coord->max1)
	seperr ("min1 is equal to max1 change one of them\n");
    if (position->xll != position->xur)
    {
	axis1->inch = position->xur - position->xll;
	orig1 = position->xll + (position->xur - position->xll) / 2.0;
    }
    else
	orig1 = (position->screenwd) / 2 ;
    if (position->yll != position->yur)
    {
	axis2->inch = position->yur - position->yll;
	orig2 = position->yll + (position->yur - position->yll) / 2.0;
    }
    else
	orig2 = (position->screenht) / 2 ;
    scale1 = axis1->inch / (coord->max1 - coord->min1);
    scale2 = axis2->inch / (coord->max2 - coord->min2);
    uorig2 = (coord->min2 + coord->max2) / 2;
    uorig1 = (coord->min1 + coord->max1) / 2;
    vp_scale (scale1, scale2);
    vp_orig (orig1, orig2);
    vp_uorig (uorig1, uorig2);
    vp_uclip (coord->min1, coord->min2, coord->max1, coord->max2); 
	return 0;
}
#include "glplot.h"
#ifdef USE_PROTO
int gl_where(float *xc, float *yc)
#else
int gl_where(xc, yc)
float *xc, *yc;
#endif
{
 	return(vp_where(&*xc, &*yc));
}
