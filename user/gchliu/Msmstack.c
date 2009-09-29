/* Stack a dataset over the second dimensions by smart stacking. */

/*
  Copyright (C) 2009 China university of Petroleum, Beijing, China
                     and University of Texas at Austin
       
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

#include <stdio.h>
#include <rsf.h>
#include <math.h>


int main(int argc, char* argv[])
{
    int nt, nh, ncmp, it, icmp, ih, zero, s, l;
    bool ifwt;
    float *indata, *outdata, *stack; /*win1, *win2; */
    float dt, esp, ee, dh, dcmp,cmp0, t0, h0, weight, sumweight;
    sf_file in, out;
    
    /* initialization */
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    /* get rsf_file parameters */
    if (!sf_histint(in,"n1",&nt)) sf_error("Need n1=");
    if (!sf_histfloat(in,"d1",&dt)) dt=1.;
    if (!sf_histfloat(in,"o1",&t0)) t0=0.;
    
    if (!sf_histint(in,"n2",&nh)) sf_error("Need n2=");
    if (!sf_histfloat(in,"d2",&dh)) dh=1.;
    if (!sf_histfloat(in,"o2",&h0)) h0=0.;

    if (!sf_histint(in,"n3",&ncmp)) ncmp=1;
    if (!sf_histfloat(in,"d3",&dcmp)) dcmp=1.;
    if (!sf_histfloat(in,"o3",&cmp0)) cmp0=0.;

    /* get parameters */
    if (!sf_getint("s",&s)) s=1;
    /*exponent*/  
    if (!sf_getint("l",&l)) l=0;
    /* parameter for alpha-trimmed mean */
    if (!sf_getbool("ifwt",&ifwt)) ifwt= true;
    if (!sf_getfloat("esp",&esp)) esp=1e-10;
    if (!sf_getfloat("ee",&ee)) esp=1.;
    /*false: equal weight;true: smart stack*/
    if (l < 0 || l > nh)  sf_warning("l exceed the range of offset");

    /* change output data dimensions */
    sf_putint(out,"n1",nt);
    sf_putint(out,"n2",ncmp);

    sf_putfloat(out,"d1",dt);
    sf_putfloat(out,"o1",t0);
    sf_putfloat(out,"d2",dcmp);
    sf_putfloat(out,"o2",cmp0);
    sf_putint(out,"n3",1);
     
    indata = sf_floatalloc(nt*nh);
    outdata = sf_floatalloc(nt);
    stack = sf_floatalloc(nt);
 
    for (icmp=0; icmp < ncmp; icmp++){        
         sf_floatread(indata,nt*nh,in);
         for (it=0; it < nt; it++){
              stack[it] = 0;
              outdata[it] = 0;
         }

         /*compute the alpha-trimmed mean trace */
         for (it=0; it < nt; it++){
              zero=0;
              for(ih=l; ih < nh-l; ih++){
                 if (indata[ih*nt+it]!=0)
                     zero++;
                 stack[it]+= indata[ih*nt+it];
              }
              if (zero==0)
              stack[it]=0.;
              else 
              stack[it]=stack[it];
         }
            
         /* compute weights */
	 for (it=0; it < nt; it++){
              zero = 0;
	      sumweight = 0;
	      for (ih=0; ih < nh; ih++){
	            if (indata[ih*nt+it]*stack[it]>0){
                        if (indata[ih*nt+it]!=0)
                            zero++;   /* number of non-zero values */
		        if (ifwt)
		            weight=1.0/(pow(fabs(indata[ih*nt+it]-stack[it]),s)+esp);
			else {
			    weight=1.0;
                         }
		        sumweight += weight;
			outdata[it] +=indata[ih*nt+it]*weight;
	             }		        
              }
	      if (sumweight==0 || zero==0)
                  outdata[it]=0;
	      else
                  outdata[it]=outdata[it]/(zero*sumweight+ee);
	      
	 }     
         sf_floatwrite(outdata,nt,out);
         sf_warning("running cmp is = %d of %d",icmp, ncmp); 
     }

    exit(0);
}
