/* Find max value and its sampled position along fast dimension */
/*
  Copyright (C) 2010 University of Texas at Austin

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

#include <math.h>
#include <rsf.h>

float find_max (const float* a, int n, int* max_ind_pt);

int main (int argc, char* argv[])
{
	bool verb;
	float o1,o2,d1,d2;	
	int n1,	n2,n3, shift;
	int i2,i3;	
	char *unit1,*label1;
	float *column=NULL,*max=NULL; 
	int *index=NULL;	
	sf_file in=NULL,out=NULL,max_val=NULL;
	
	sf_init (argc,argv);
	
	in=sf_input("in");
	out=sf_output("out");
	sf_settype (out, SF_INT);
	if (NULL != sf_getstring("max_val")) max_val=sf_output("max_val");
	
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
	
	label1=sf_charalloc(100);	
	unit1=sf_charalloc(100);	

	label1=sf_histstring(in,"label2");
	unit1=sf_histstring(in,"unit2");

	if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");

	if (!sf_getint("shift",&shift)) shift=0;
	/* shift */
	if (!sf_getbool("verb",&verb)) verb=false;
	

	column = sf_floatalloc(n1);
    max = sf_floatalloc(n2);
    index = sf_intalloc(n2);

	sf_putint(out,"n1",n2);
    sf_putfloat(out,"o1",o2);

    sf_putfloat(out,"d1",d2);

	if (!(label1==NULL)) sf_putstring(out,"label1",label1);
	if (!(unit1==NULL))  sf_putstring(out,"unit1",unit1);
	//sf_warning("Son qua d2=%s",label1);

	sf_putint(out,"n2",1);
	sf_putstring(out,"o2","");
    sf_putstring(out,"d2","");
	sf_putstring(out,"label2","");
	sf_putstring(out,"unit2","");
	
	if (!(max_val==NULL)) {
		sf_putint(max_val,"n1",n2);
    	sf_putfloat(max_val,"o1",o2);
		sf_putint(max_val,"n2",1);
		sf_putstring(max_val,"o2","");
    	sf_putstring(max_val,"d2","");
		sf_putstring(max_val,"label2","");
		sf_putstring(max_val,"unit2","");
	}
	/* reading the number of gahters in data*/
    n3 = sf_leftsize(in,2);	

	for (i3=0;i3<n3;i3++) { /*gahters loop */
	    sf_warning("Gather %d/%d",i3+1,n3);
		for (i2=0;i2<n2;i2++) {

			sf_floatread(column,n1,in);

			max[i2]=find_max (column, n1, index+i2);
			if (d1<0)
			index[i2]=n1-index[i2]+shift;
			else
			index[i2]-=shift;
	    	//sf_warning("Son qua, max=%f index=%d",max[i2],index[i2]);
		}
	
	if (!(max_val==NULL)) {
		sf_floatwrite(max,n2,max_val);
	}
		sf_intwrite(index,n2,out);	
	} /* END gahters loop */
    exit (0);
}

float find_max (const float* a, int n, int* max_ind_pt) {

        float max_value;
        int i;

        max_value = a[0];

        *max_ind_pt = 0;

        for (i = 1; i < n; i++)
		{
                if (a[i] > max_value)
				{
                        max_value = a[i];
                        *max_ind_pt  = i;
                }
        }
        return max_value;
}
