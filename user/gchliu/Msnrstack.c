/* Stack a dataset over the second dimensions by SNR weighted method. */

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
    int nt, nh, ncmp, it, icmp, ih, iw, w, zero, shift;
    float *indata, *outdata, *stack, *win1, *win2, *outweight;
    float dt,var, maxwin1, maxwin2, ee, esp, dh, dcmp,cmp0, t0, h0, /* sumab, sumwin1, sumwin2, */ sumweight, sft;

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

    /* get other parameters */
    if (!sf_getint("w",&w)) w=50;  
    /* sliding window size*/    
    if (!sf_getfloat("sft",&sft)) sft=1;  
    /*weight shift*/
    if (!sf_getfloat("ee",&ee)) ee=1.0;      
    if (!sf_getfloat("esp",&esp)) esp=1.0;       
    
    sf_putint(out,"n1",nt);
    sf_putint(out,"n2",ncmp);

    sf_putfloat(out,"d1",dt);
    sf_putfloat(out,"o1",t0);
    sf_putfloat(out,"d2",dcmp);
    sf_putfloat(out,"o2",cmp0);
    sf_putint(out,"n3",1);
    
    
                  
    indata = sf_floatalloc(nt*nh);
    outdata = sf_floatalloc(nt);
    outweight = sf_floatalloc(nt*nh);
    stack = sf_floatalloc(nt);
    win1 = sf_floatalloc(w);
    win2 = sf_floatalloc(w);
 
    for (icmp=0; icmp < ncmp; icmp++){ 
        
         sf_floatread(indata,nt*nh,in);

         for (it=0; it < nt; it++){
              stack[it] = 0;
              outdata[it] = 0;
	      for (ih=0;ih<nh;ih++){
	           outweight[nt*ih+it]=0;
	      }
         }
	 
	 
         /* computer the directly stack trace */
         for (it=0; it < nt; it++){
              zero=0;
              for(ih=0; ih < nh; ih++){
                 if (indata[ih*nt+it]!=0)
                     zero++;
                 
                 stack[it]+= indata[ih*nt+it];
              }
              if (zero==0)
              stack[it]=stack[it];
              else 
              stack[it]=stack[it]/zero;
         }         
         /* estimate the noise variances */
         for (it=0; it < nt; it++){ 
              zero = 0;
	       sumweight=0;
              for (ih=0; ih < nh; ih++){
		  /* sumwin1 = 0;
                    sumwin2 = 0;
                    sumab = 0; */
                   shift = SF_MAX(0,SF_MIN(nt-w, it-w/2-1));
		   /* weight=0; */
		   var=0;
		   maxwin1=0;
		   maxwin2=0;
                   for (iw=0; iw <w; iw++){
                        win1[iw] = indata[ih*nt+iw+shift];
                        win2[iw] = stack[iw+shift];
			if (fabs(win1[iw])>fabs(maxwin1))
			     maxwin1 = fabs(win1[iw]);
			if (fabs(win2[iw])>fabs(maxwin2))
			     maxwin2 = fabs(win2[iw]);
	           }
		   for (iw=0; iw <w; iw++){		     
			var += fabs(win1[iw]/maxwin1-win2[iw]/maxwin2)*fabs(win1[iw]/maxwin1-win2[iw]/maxwin2);    
                   }
                   outweight[nt*ih+it] =1./(var+ee); 
                   
                   if (indata[ih*nt+it]!=0){
                     zero++;
                     outdata[it] += pow(1./(var+ee),1)*indata[ih*nt+it];
                     sumweight +=  pow(1./(var+ee),1);
		   }
              }
              if (zero==0)
                  outdata[it]=0;              
              else     
                  outdata[it] = outdata[it]/(zero*sumweight+esp);
              
         }
	 
         sf_floatwrite(outdata,nt,out);
         sf_warning("running cmp is = %d of %d",icmp, ncmp); 
     }

    exit(0);
}
