/* Seislet Transform Denoising using Mask Operator*/
/*
  Copyright (C) 2013 University of Texas at Austin
   
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
#include <rsfpwd.h>

int main(int argc, char *argv[])
{
    int i, j, n1, n2, n3, n12; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/
    int order1, cut;
    float *data1, *data2, /* *adata, */ **dd;
    sf_file in, out, dip, slet, sletcut;
    float eps;
    char *type;
    bool unit=false, inv=true;

    sf_init(argc,argv);

    in = sf_input("in");
    dip=sf_input("dip");

    out = sf_output("out");
    slet=sf_output("slet");
    sletcut=sf_output("sletcut");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12=n1*n2;
    n3 = sf_leftsize(in,2);
    /* get the trace length (n1) and the number of traces (n2) and n3*/

    if(!sf_getint("cut",&cut)) cut=n2/4;
    /* cut threshold value */

    if(!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if(!sf_getint("order1",&order1))order1=1;
    /* accuracy order for seislet transform*/

    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */

    if(!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */
   

    data1 = sf_floatalloc(n1*n2);	/*allocate memory*/
    data2 = sf_floatalloc(n1*n2);	/*allocate memory*/
    /*   adata = sf_floatalloc(n1*n2);	allocate memory*/
    dd=sf_floatalloc2(n1,n2);

    seislet_init(n1,n2,inv,unit,eps,order1,type[0]);  /* unit=false inv=true */

    for(i=0;i<n3;i++)  {
	sf_floatread(data1,n1*n2,in);
	seislet_set(dd);
	sf_floatread(dd[0],n12,dip);

	seislet_lop(true,false,n12,n12,data2,data1); 
 	/*data1 is the input raw data and data2 is the unsorted seislet domain seismic data.*/

	if(NULL!=sf_getstring("slet"))	
	sf_floatwrite(data2,n1*n2,slet);
	/* seismic domain */

	/***********************************************************/
	/*applying mask operator*/
	for(i=cut; i<n2; i++)
		for(j=0;j<n1;j++)
			data2[i*n1+j]=0;
	/***********************************************************/

	if(NULL!=sf_getstring("sletcut"))	
	sf_floatwrite(data2,n1*n2,sletcut);
	/* cutted seislet domain */

	seislet_lop(false,false,n12,n12,data2,data1);	
	/*data2 is sorted seislet domain seismic data and data1 is denoised t-x domain seismic data*/

	sf_floatwrite(data1,n1*n2,out);
    }

    exit(0);
}

