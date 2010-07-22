/* Plot rays.

Takes: < rays.rsf > plot.vpl

Run "sfdoc stdplot" for more parameters.
*/
/*
  Copyright (C) 2004 University of Texas at Austin

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
#include <rsfplot.h>

int main(int argc, char* argv[])
{
    int n1, n2, ir, jr;
    float o1, d1, o2, d2, **traj;
    int nsr, it, nt;
    sf_file frame=NULL;

    sf_init (argc,argv);
    vp_init ();

    if (NULL != sf_getstring("frame")) frame = sf_input("frame");

    if (!sf_getint("n1",&n1) && 
	(NULL == frame || 
	 !sf_histint(frame,"n1",&n1))) sf_error("Need n1=");
    if (!sf_getint("n2",&n2) &&
	(NULL == frame || 
	 !sf_histint(frame,"n2",&n2))) sf_error("Need n2=");
    /* frame dimensions */
    if (!sf_getfloat("d1",&d1) &&
	(NULL == frame || 
	 !sf_histfloat(frame,"d1",&d1))) sf_error("Need d1=");
    if (!sf_getfloat("d2",&d2) &&
	(NULL == frame || 
	 !sf_histfloat(frame,"d2",&d2))) sf_error("Need d2=");
    /* frame sampling */
    if (!sf_getfloat("o1",&o1) &&
	(NULL == frame || 
	 !sf_histfloat(frame,"o1",&o1))) sf_error("Need o1=");
    if (!sf_getfloat("o2",&o2) &&
	(NULL == frame || 
	 !sf_histfloat(frame,"o2",&o2))) sf_error("Need o2=");
    /* frame origin */
    if (!sf_getint("nt",&nt)) nt=n1*n2;
    /* maximum ray length */
    if (!sf_getint("jr",&jr)) jr=1;
    /* skip rays */

    traj = sf_floatalloc2 (2,nt);

    /* transp and yreverse */
    vp_stdplot_init (o1,o1+(n1-1)*d1,o2,o2+(n2-1)*d2,
		     true,false,true,false);
    vp_plot_init(1);
    vp_plot_set(0);
 
    if (1 != fread(&nsr,sizeof(int),1,stdin)) sf_error("read error");

    for (ir=0; ir < nsr; ir++) {
	if (1 != fread(&it,sizeof(int),1,stdin) ||
	    (it+1)*2 != fread(traj[0],sizeof(float),(it+1)*2,stdin))
	    sf_error("read error");
	if (ir>0 && ir%jr) continue;  
	vp_umove(traj[it][1],traj[it][0]);
	while (--it >= 0) {
	    vp_udraw(traj[it][1],traj[it][0]);
	}
    }

    vp_frame();


    exit(0);
}
