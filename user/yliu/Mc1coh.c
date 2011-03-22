/* C1 coherency algorithm. */
/*
  Copyright (C) 2011 Jilin University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>	

int main(int argc, char* argv[])
{
    bool verb;
    int n1, i1, i2, n2, i3, n3, i4, n4, lag1, lag2, lag, ntw, win, w;
    float done, dtwo, corr;
    float *data, *indat, *crossdat;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) n3=1;
    n4 = sf_leftsize(in,3);

    if (!sf_getint("ntw",&ntw)) ntw=3;
    /* Temporal length of the correlation window (default=3) */

    if (!sf_getbool("verb",&verb)) verb = true;
    /* verbosity flag */

    if (0== ntw%2) ntw++;
    win = (ntw-1)/2;
    
    if (!sf_getint("lag1",&lag1)) lag1=3;
    /* Inline time lag (default=3) */
 
    if (!sf_getint("lag2",&lag2)) lag2=3;
    /* Crossline time lag (default=3) */

    data     = sf_floatalloc(n1*n2*n3);
    indat    = sf_floatalloc(n1*n2*n3);
    crossdat = sf_floatalloc(n1*n2*n3);

    for (i4=0; i4 < n4; i4++) {
	sf_floatread(data,n1*n2*n3,in);

	for (i3=0; i3 < n3; i3++) {
	    if (verb) sf_warning("slice %d of %d;",i3+1,n3);
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    /* inline l-lag crosscorrelation */
		    indat[i3*n2*n1+i2*n1+i1] = 0.;
		    for (lag=-lag1; lag <= lag1; lag++) {
			done = 0.;
			dtwo = 0.;
			corr = 0.;
			for (w=-win; w <= win; w++) {
			    if (((i2+1) < n2) && ((i1+w) >=0) && ((i1+w) < n1) && ((i1+w+lag) >=0) && ((i1+w+lag) < n1)) {
				done += data[i3*n2*n1+i2*n1+(i1+w)]*data[i3*n2*n1+i2*n1+(i1+w)];
				dtwo += data[i3*n2*n1+(i2+1)*n1+(i1+w+lag)]*data[i3*n2*n1+(i2+1)*n1+(i1+w+lag)];
				corr += data[i3*n2*n1+i2*n1+(i1+w)]*data[i3*n2*n1+(i2+1)*n1+(i1+w+lag)];
				
			    }
			}
			corr = fabsf(corr)/(sqrtf(done*dtwo)+FLT_EPSILON);
			if (corr > indat[i3*n2*n1+i2*n1+i1]) {
			    indat[i3*n2*n1+i2*n1+i1] = corr;
			}
		    }
		    if (1 != n3) {
			/* crossline m-lag crosscorrelation */
			crossdat[i3*n2*n1+i2*n1+i1] = 0.;
			for (lag=-lag2; lag <= lag2; lag++) {
			    done = 0.;
			    dtwo = 0.;
			    corr = 0.;
			    for (w=-win; w <= win; w++) {
				if (((i3+1) < n3) && ((i1+w) >=0) && ((i1+w) < n1) && ((i1+w+lag) >=0) && ((i1+w+lag) < n1)) {
				    done += data[i3*n2*n1+i2*n1+(i1+w)]*data[i3*n2*n1+i2*n1+(i1+w)];
				    dtwo += data[(i3+1)*n2*n1+i2*n1+(i1+w+lag)]*data[(i3+1)*n2*n1+i2*n1+(i1+w+lag)];
				    corr += data[i3*n2*n1+i2*n1+(i1+w)]*data[(i3+1)*n2*n1+i2*n1+(i1+w+lag)];
				    
				}
			    }
			    corr = fabsf(corr)/(sqrtf(done*dtwo)+FLT_EPSILON);
			    if (corr > crossdat[i3*n2*n1+i2*n1+i1]) {
				crossdat[i3*n2*n1+i2*n1+i1] = corr;
			    }
			}
		    }
		}
	    }
	}

	if (1 != n3) {
	    for (i3=0; i3 < n3; i3++) {
		for (i2=0; i2 < n2; i2++) {
		    for (i1=0; i1 < n1; i1++) {
			data[i3*n2*n1+i2*n1+i1] = sqrtf(indat[i3*n2*n1+i2*n1+i1]*crossdat[i3*n2*n1+i2*n1+i1]);
		    }
		}
	    }
	    sf_floatwrite(data,n1*n2*n3,out);
	} else {
	    sf_floatwrite(indat,n1*n2*n3,out);
	}
    }
    sf_warning(".");
    exit(0);
}
/* 	$Id$	 */
