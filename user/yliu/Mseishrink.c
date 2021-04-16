/* 2-D Seislet shrinkage denoising. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

int main(int argc, char* argv[])
{
    int nw, n1, n2, n12, i1, i2, i3, n3, j, iter, nperc, interp, order;
    float *mm, *dd, *dd2=NULL, **pp, *d1, *m1, eps, perc, *hcurv;
    float *data, *model, *ldata=NULL, min, max, dm, dperc;
    sf_complex *norm12;
    char *type;
    bool verb, dwt;
    sf_file in, out, dip, lcurve=NULL, hcurve=NULL, norm=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n3 = sf_leftsize(in,2);

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

    if (NULL == (type=sf_getstring("type"))) type="biorthogonal";
    /* [haar,linear,biorthogonal] wavelet type, the default is biorthogonal  */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (!sf_getbool("dwt",&dwt)) dwt = false;
    /* if y, dwt in vertical axis */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization parameter */

    if (!sf_getfloat("perc",&perc)) perc=90.;
    /* percentage for shrinkage */

    if (!sf_getint("nperc",&nperc)) nperc=1;
    /* number of percentage dimension */
    if (nperc < 1) nperc=1;

    if (NULL != sf_getstring ("lcurve")) {
	lcurve = sf_output("lcurve");
    }
    if (NULL != sf_getstring ("hcurve")) {
	hcurve = sf_output("hcurve");
    }
    if (NULL != sf_getstring ("norm")) {
	norm = sf_output("norm");
    }

    pp = sf_floatalloc2(n1,n2);
    mm = sf_floatalloc(n12);
    dd = sf_floatalloc(n12);
    dd2 = sf_floatalloc(n12);
    d1 = sf_floatalloc(n1);
    m1 = sf_floatalloc(n1);
    data = sf_floatalloc(nperc);
    model = sf_floatalloc(nperc);
    norm12 = sf_complexalloc(nperc);
    min = 0.;
    max = 0.;
    order = 1;
    dm = 1.e+20;
    dperc = 100.0/(nperc-1);

    seislet_init(n1,n2,true,true,eps,order,type[0]);
    sf_wavelet_init(n1,true,true,type[0]);

    seislet_set(pp);
    sf_floatread(pp[0],n12,dip);

    for (i3=0; i3 < n3; i3++) {
	sf_warning("slice %d of %d",i3+1,n3);

	sf_floatread(dd,n12,in);
	for (i1=0; i1 < n12; i1++) {
	    dd2[i1] = dd[i1];
	}
	if (1==nperc) {
	    sf_sharpen_init(n12,perc,0.5);
	    seislet_lop(true,false,n12,n12,mm,dd);
	    if (dwt) {
		for (i2=0; i2 < n2; i2++) {
		    for (j=0; j < n1; j++) {
			d1[j] = mm[i2*n1+j];
		    }
		    sf_wavelet_lop(false,false,n1,n1,d1,m1);
		    for (j=0; j < n1; j++) {
			mm[i2*n1+j] = m1[j];
		    }
		}
	    } /* Forward DWT */
	    sf_sharpen(mm);
	    sf_weight_apply(n12,mm);
	    if (dwt) {
		for (i2=0; i2 < n2; i2++) {
		    for (j=0; j < n1; j++) {
			m1[j] = mm[i2*n1+j];
		    }
		    sf_wavelet_lop(true,false,n1,n1,d1,m1);
		    for (j=0; j < n1; j++) {
			mm[i2*n1+j] = d1[j];
		    }
		}
	    }	 /* Inverse DWT */    
	    seislet_lop(false,false,n12,n12,mm,dd2);
	    sf_floatwrite (dd2,n12,out);
	} else {
	    dm = 1.e+20;
	    
	    /* L-curve */
	    for (iter=0; iter < nperc; iter++) {
		perc = dperc*iter;
		if (perc > 100.) perc=100.;
		if (verb)
		    sf_warning("Percentage %2.1f of 100",perc);
		sf_sharpen_init(n12,perc,0.5);
		seislet_lop(true,false,n12,n12,mm,dd);
		if (dwt) {
		    for (i2=0; i2 < n2; i2++) {
			for (j=0; j < n1; j++) {
			    d1[j] = mm[i2*n1+j];
			}
			sf_wavelet_lop(false,false,n1,n1,d1,m1);
			for (j=0; j < n1; j++) {
			    mm[i2*n1+j] = m1[j];
			}
		    }
		} /* Forward DWT */
		sf_sharpen(mm);
		sf_weight_apply(n12,mm);
		if (dwt) {
		    for (i2=0; i2 < n2; i2++) {
			for (j=0; j < n1; j++) {
			    m1[j] = mm[i2*n1+j];
			}
			sf_wavelet_lop(true,false,n1,n1,d1,m1);
			for (j=0; j < n1; j++) {
			    mm[i2*n1+j] = d1[j];
			}
		    }
		}	 /* Inverse DWT */    		
		seislet_lop(false,false,n12,n12,mm,dd2);

		data[iter] = 0.;
		model[iter] = 0.;
		for (j=0; j < n12; j++) {
		    data[iter] += (dd[j]-dd2[j])*(dd[j]-dd2[j]);
		    model[iter] += fabsf(mm[j]);
		}
		data[iter] = 10*logf(data[iter]);
		model[iter] = 10*logf(model[iter]);
		/* sf_warning("%d %f %f",iter,model[iter],data[iter]); */
		if (0!=iter && iter<nperc) {
		    if (dm > fabsf(model[iter]-model[iter-1]) 
			&& 0.!=fabsf(model[iter]-model[iter-1])) {
			dm = fabsf(model[iter]-model[iter-1]);
		    }
		}
		norm12[iter]=sf_cmplx(model[iter],data[iter]);
	    }

	    if (NULL != sf_getstring ("norm")) {
		sf_settype(norm,SF_COMPLEX);
		sf_putint(norm,"n1",nperc);
		sf_putfloat(norm,"d1",dperc);
		sf_putfloat(norm,"o1",0.);
		sf_putstring(norm,"label1","Percentage");
		sf_putstring(norm,"unit1","");
		sf_putint(norm,"n2",1);
		sf_putstring(norm,"label2","Norm1&2");
		sf_putstring(norm,"unit2","");
		sf_putint(norm,"n3",1);
		sf_complexwrite(norm12,nperc,norm);
	    }

	    max = model[0];
	    min = model[nperc-2];
	    /* interpolate model */
	    interp = (int)((max-min)/dm) +1;
	    ldata = sf_floatalloc(interp);
	    /* sddata = sf_floatalloc(interp); */
	    order = 1;
	    for (j=0; j < interp; j++) {
		if ((max-dm*j) == model[order] && order < (nperc-1)) {
		    ldata[j] = data[order];
		    order ++;
		} else {
		    ldata[j] = data[order-1]+(data[order]-data[order-1])/
			(model[order]-model[order-1])*
			(max-dm*j-model[order-1]);
		    if ((max-dm*j) < model[order] && order < (nperc-1)) {
			order ++;
		    }
		}
	    }
	    if (NULL != sf_getstring ("lcurve")) {
		sf_putint(lcurve,"n1",interp);
		sf_putfloat(lcurve,"d1",-dm);
		sf_putfloat(lcurve,"o1",max);
		sf_putstring(lcurve,"label1","||m||\\_\\s75 1\\^\\s100 ");
		sf_putstring(lcurve,"unit1","");
		sf_putint(lcurve,"n2",1);
		sf_putstring(lcurve,"label2",
			     "||d-d*||\\_\\s75 2\\s300 \\^\\s75 2 \\s100 ");
		sf_putstring(lcurve,"unit2","");
		sf_putint(lcurve,"n3",1);
		sf_floatwrite(ldata,interp,lcurve);
	    }

	    hcurv = sf_floatalloc(nperc);
	    /* H-curve */
	    for (j=1; j < nperc-1; j++) {
	        hcurv[j] = ((model[j+1]-model[j-1])*
			    (data[j+1]-2*data[j]+data[j-1])-
			    (data[j+1]-data[j-1])*
			    (model[j+1]-2*model[j]+model[j-1]))/
		            powf(((model[j+1]-model[j-1])*
				  (model[j+1]-model[j-1])+
				  (data[j+1]-data[j-1])*
				  (data[j+1]-data[j-1])),1.5);
	    }
	    hcurv[0] = 0.;
	    hcurv[nperc-1] = 0.;

	    if (NULL != sf_getstring ("hcurve")) {
		sf_putint(hcurve,"n1",nperc);
		sf_putfloat(hcurve,"d1",dperc);
		sf_putfloat(hcurve,"o1",0.);
		sf_putstring(hcurve,"label1","Percentage");
		sf_putstring(hcurve,"unit1","");
		sf_putint(hcurve,"n2",1);
		sf_putstring(hcurve,"label2","Curvature");
		sf_putstring(hcurve,"unit2","");
		sf_putint(hcurve,"n3",1);
		sf_floatwrite(hcurv,nperc,hcurve);
	    }
	    max = hcurv[nperc/2];
	    perc = 0.;
	    for (j=(int)(0.1*nperc); j < nperc; j++) {
		if (max < hcurv[j]) {
		    max = hcurv[j];
		    perc = j*dperc;
		}
	    }
	    sf_warning("Key thresholding is percentage=%2.1f",perc);

	    /* Shrinkage with key thresholding */
	    sf_sharpen_init(n12,perc,0.5);
	    seislet_lop(true,false,n12,n12,mm,dd);
	    if (dwt) {
		for (i2=0; i2 < n2; i2++) {
		    for (j=0; j < n1; j++) {
			d1[j] = mm[i2*n1+j];
		    }
		    sf_wavelet_lop(false,false,n1,n1,d1,m1);
		    for (j=0; j < n1; j++) {
			mm[i2*n1+j] = m1[j];
		    }
		}
	    } /* Forward DWT */
	    sf_sharpen(mm);
	    sf_weight_apply(n12,mm);
	    if (dwt) {
		for (i2=0; i2 < n2; i2++) {
		    for (j=0; j < n1; j++) {
			m1[j] = mm[i2*n1+j];
		    }
		    sf_wavelet_lop(true,false,n1,n1,d1,m1);
		    for (j=0; j < n1; j++) {
			mm[i2*n1+j] = d1[j];
		    }
		}
	    }	 /* Inverse DWT */    
	    seislet_lop(false,false,n12,n12,mm,dd2);
	    sf_floatwrite (dd2,n12,out);

	}
    }
    exit(0);
}

/* 	$Id$	 */
