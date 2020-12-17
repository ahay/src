/* Simple acoustic RTM code
#
# Originally written by Yangkang Chen
# Modified by Shan Qu
# 
# CODE for DTV regularized FWI
# Reference:
# Qu, S., E. Verschuur, and Y. Chen, 2019, FWI/JMI with an automatic directional total variation constraint, Geophysics, 84, R175-R183.
#
*/

#include <rsf.h>
#include <time.h>

static float c0=-30./12.,c1=+16./12.,c2=- 1./12.;

void oneway_abc(float **uo, float **um, float **vm, int nx, int nz, int nb, float dx, float dz, float dt, bool iffree)
/* oneway - absorbing condition*/
{
    int iz,ix,ib;
    float q;
    float qx,qt,qxt,b=0.5;/*b=0.5 is the best.*/

    for(ix=0;ix<nx;ix++) {
	for(ib=0;ib<nb;ib++) {

	    /* top BC */
	    if(!iffree) { /* not free surface, apply ABC */
		iz = nb-ib;
		q = vm[ix][           nb  ] *dt/dz; 
		qx=(b*(1+q)-q)/(1+q)/(1-b);
		qt=(b*(1+q)-1)/(1+q)/(1-b);
		qxt=b/(b-1);
				
 	    uo      [ix  ][iz] 
		= um[ix][iz+1]*(-qxt) 
		+ um[ix][iz]*(-qt)
		+ uo[ix][iz+1]*(-qx); 		    
	    }

	    /* bottom BC */
	    iz = nz-nb+ib-1;
	    q = vm[ix][nz-nb-1] *dt/dz; 

		qx=(b*(1+q)-q)/(1+q)/(1-b);
		qt=(b*(1+q)-1)/(1+q)/(1-b);
		qxt=b/(b-1);
				
 	    uo      [ix  ][iz] 
		= um[ix][iz-1]*(-qxt) 
		+ um[ix][iz]*(-qt)
		+ uo[ix][iz-1]*(-qx); 

	}
    }

    for(iz=0;iz<nz;iz++) {
	for(ib=0;ib<nb;ib++) {

	    /* left BC */
	    ix = nb-ib;
	    q = vm[           nb  ][iz] *dt/dx;
		qx=(b*(1+q)-q)/(1+q)/(1-b);
		qt=(b*(1+q)-1)/(1+q)/(1-b);
		qxt=b/(b-1);
		
	    uo      [ix  ][iz] 
		= um[ix+1][iz]*(-qxt) 
		+ um[ix  ][iz]*(-qt)
		+ uo[ix+1][iz]*(-qx);		


	    /* right BC */
	    ix = nx-nb+ib-1;
	    q = vm[nx-nb-1][iz] *dt/dx; 
		qx=(b*(1+q)-q)/(1+q)/(1-b);
		qt=(b*(1+q)-1)/(1+q)/(1-b);
		qxt=b/(b-1);
			    
	    uo      [ix  ][iz] 
		= um[ix-1][iz]*(-qxt) 
		+ um[ix  ][iz]*(-qt)
		+ uo[ix-1][iz]*(-qx);			
	}
    }

}

void step_forward(float **p0, float **p1, float **p2, float **vm, float dt, float dz, float dx, int nz, int nx, int nb, bool iffree)
/*< forward modeling step, Clayton-Enquist ABC incorporated >*/
{
    int ix,iz;
    float tmp;
    float idx,idz;
    idx=1./(dx*dx);idz=1./(dz*dz);
    
	//sf_warning("idx=%g,idz=%g",idx,idz);	    
	for (iz=nb; iz<nz-nb; iz++) {
	    for (ix=nb; ix<nx-nb; ix++) {
		tmp = 
		    c0* p1[ix  ][iz  ] * (idx+idz) + 
		    c1*(p1[ix-1][iz  ] + p1[ix+1][iz  ])*idx +
		    c2*(p1[ix-2][iz  ] + p1[ix+2][iz  ])*idx +
		    c1*(p1[ix  ][iz-1] + p1[ix  ][iz+1])*idz +
		    c2*(p1[ix  ][iz-2] + p1[ix  ][iz+2])*idz;

		    p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vm[ix][iz]*vm[ix][iz]*dt*dt*tmp;	 

	    }
	}
	oneway_abc(p2, p1, vm, nx, nz, nb, dx, dz, dt, iffree);
	
}

void step_backward(float **illum, float **lap, float **p0, float **p1, float **p2, float **vm, float dt, float dz, float dx, int nz, int nx, int nb)
/*< step backward >*/
{   
    int ix,iz;
    float tmp;
    float idx,idz;
    idx=1./(dx*dx);idz=1./(dz*dz);
    
	//sf_warning("idx=%g,idz=%g",idx,idz);	    
	for (iz=nb; iz<nz-nb; iz++)
	    for (ix=nb; ix<nx-nb; ix++) {
		tmp = 
		    c0* p1[ix  ][iz  ] * (idx+idz) + 
		    c1*(p1[ix-1][iz  ] + p1[ix+1][iz  ])*idx +
		    c2*(p1[ix-2][iz  ] + p1[ix+2][iz  ])*idx +
		    c1*(p1[ix  ][iz-1] + p1[ix  ][iz+1])*idz +
		    c2*(p1[ix  ][iz-2] + p1[ix  ][iz+2])*idz;
		    lap[ix][iz]=tmp;
		    p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vm[ix][iz]*vm[ix][iz]*lap[ix][iz]*dt*dt;
		    illum[ix][iz]+=p1[ix][iz]*p1[ix][iz];	 
	    }
}

void source(float **p, float *source, int *sxz, int ns, int nz, bool add)
/*< add/subtract seismic sources >*/
{
	int is, sx, sz;
	if(add){
		for(is=0;is<ns; is++){
			sx=sxz[is]/nz;
			sz=sxz[is]%nz;
			p[sx][sz]+=source[is];
		}
	}else{
		for(is=0;is<ns; is++){
			sx=sxz[is]/nz;
			sz=sxz[is]%nz;
			p[sx][sz]-=source[is];
		}
	}
}

void transpose(float *matrix, float *trans, int n1, int n2)
/*< matrix transpose: matrix tansposed to be trans >*/
{
	int i1, i2;

	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	    trans[i2+n2*i1]=matrix[i1+n1*i2];
}

void record(float *seis_it, int *gxz, float **p, int ng, int nz)
/*< record seismogram at time it into a vector length of ng >*/
{
	int ig, gx, gz;
	for(ig=0;ig<ng; ig++)
	{
		gx=gxz[ig]/nz;
		gz=gxz[ig]%nz;
		seis_it[ig]=p[gx][gz];
	}
}

void sg_position(int *sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns, int nz)
/*<  shot/geophone positions >*/
{
	int is, sz, sx;
	for(is=0; is<ns; is++)
	{
		sz=szbeg+is*jsz;
		sx=sxbeg+is*jsx;
		sxz[is]=sz+nz*sx;
	}
}

void rw_boundary(float *bd, float **p, int nz, int nx, bool write)
/*< if write==true, write/save boundaries out of variables;
 else  read boundaries into wavefields (for 4th order FD) >*/
{
	int ix,iz;
	if(write){
		        for(ix=0; ix<nx; ix++)
                for(iz=0; iz<2; iz++)
                {
                        bd[iz+4*ix]=p[ix][iz];
                        bd[iz+2+4*ix]=p[ix][nz-1-iz];
                }
                
		        for(iz=0; iz<nz; iz++)
                for(ix=0; ix<2; ix++)
                {
                        bd[4*nx+iz+nz*ix]=p[ix][iz];
                        bd[4*nx+iz+nz*(ix+2)]=p[nx-1-ix][iz];
                }	
	}else{
                for(ix=0; ix<nx; ix++)
                for(iz=0; iz<2; iz++)
                {
                        p[ix][iz]=bd[iz+4*ix];
                        p[ix][nz-1-iz]=bd[iz+2+4*ix];
                }
                
                for(iz=0; iz<nz; iz++)
                for(ix=0; ix<2; ix++)
                {
                        p[ix][iz]=bd[4*nx+iz+nz*ix];
                        p[nx-1-ix][iz]=bd[4*nx+iz+nz*(ix+2)];
                }
	}
}

void source_illum(float **imag2, float **imag1, float **imagtmp, float **illum, int nz, int nx, int nb)
/*< apply source illumination >*/
{
  int ix, iz;
  for(ix=nb; ix<nx-nb; ix++){
    for(iz=nb; iz<nz-nb; iz++){
	imag2[ix][iz]+=imagtmp[ix][iz]/(illum[ix][iz]+SF_EPS);
	imag1[ix][iz]+=imagtmp[ix][iz];
    }
  }
}

void cal_residual(float *dcal, float *dobs, float *dres, int ng)
/*< calculate residual >*/
{
  int ig;
  for(ig=0; ig<ng; ig++){
    dres[ig]=dcal[ig]-dobs[ig];
  }
}

void xcorr(float **imag1, float **sp, float **gp, int nz, int nx)
/*< cross-correlation >*/
{
  int ix, iz;
  for(ix=0; ix<nx; ix++){
    for(iz=0; iz<nz; iz++){
      imag1[ix][iz]+=sp[ix][iz]*gp[ix][iz];
    }
  }
}

int main(int argc, char *argv[])
{
	bool verb, roll, iffree;
	int is, it, distx, distz, ifr;
	int nz, nx, nt, ns, ng, nb;
	int sxbeg, szbeg, gxbeg, gzbeg, jsx, jsz, jgx, jgz;/*  parameters of acquisition geometery */
	float dx, dz, fm, dt, tmp, amp;
	float *dobs, *dcal, *dref, *wlt, *bd, *trans;
	int *sxz, *gxz;		
	float **vm, **vw, **imag1, **imag2, **illum, **imagtmp, **lap, **sp0, **sp1, **sp2, **gp0, **gp1, **gp2, **ptr=NULL;	
	clock_t start=0, stop;/* timer */
	sf_file vmig,vwater, shots, img1,img2,illums;/* I/O files */

    /* initialize Madagascar */
    sf_init(argc,argv);

    /* set up I/O files */
    vmig=sf_input ("in");   /* migration velocity model, unit=m/s */
    vwater=sf_input ("vwater");   /* water velocity model, unit=m/s */    
	shots=sf_input("shots"); /* recorded shots */
    img1=sf_output("out"); /* RTM image without source illumination*/ 
    img2=sf_output("img2"); /* RTM image with source illumination*/ 
        
	illums=sf_output("illums");/* source illumination */

    /* get parameters from velocity model and recorded shots */
	if (!sf_getbool("verb",&verb)) verb=true;/* vebosity */
    if (!sf_histint(vmig,"n1",&nz)) sf_error("no n1");/* nz */
    if (!sf_histint(vmig,"n2",&nx)) sf_error("no n2");/* nx */
    if (!sf_histfloat(vmig,"d1",&dz)) sf_error("no d1");/* dz */
   	if (!sf_histfloat(vmig,"d2",&dx)) sf_error("no d2");/* dx */

   	if (!sf_histint(shots,"n1",&nt)) sf_error("no nt");
	/* total modeling time steps */
   	if (!sf_histint(shots,"n2",&ng)) sf_error("no ng");
	/* total receivers in each shot */
   	if (!sf_histint(shots,"n3",&ns)) sf_error("no ns");
	/* number of shots */
   	if (!sf_histfloat(shots,"d1",&dt)) sf_error("no dt");
	/* time sampling interval */
   	if (!sf_histfloat(shots,"amp",&amp)) sf_error("no amp");
	/* maximum amplitude of ricker */
   	if (!sf_histfloat(shots,"fm",&fm)) sf_error("no fm");
	/* dominant freq of ricker */
   	if (!sf_histint(shots,"sxbeg",&sxbeg)) sf_error("no sxbeg");
	/* x-begining index of sources, starting from 0 */
   	if (!sf_histint(shots,"szbeg",&szbeg)) sf_error("no szbeg");
	/* x-begining index of sources, starting from 0 */
   	if (!sf_histint(shots,"gxbeg",&gxbeg)) sf_error("no gxbeg");
	/* x-begining index of receivers, starting from 0 */
   	if (!sf_histint(shots,"gzbeg",&gzbeg)) sf_error("no gzbeg");
	/* x-begining index of receivers, starting from 0 */
   	if (!sf_histint(shots,"jsx",&jsx)) sf_error("no jsx");
	/* source x-axis  jump interval  */
   	if (!sf_histint(shots,"jsz",&jsz)) sf_error("no jsz");
	/* source z-axis jump interval  */
   	if (!sf_histint(shots,"jgx",&jgx)) sf_error("no jgx");
	/* receiver x-axis jump interval  */
   	if (!sf_histint(shots,"jgz",&jgz)) sf_error("no jgz");
	/* receiver z-axis jump interval  */
   	if (!sf_histint(shots,"roll",&ifr)) sf_error("roll or not required");
	/* default: rolling acquisition*/

    if(! sf_getbool("free",&iffree)) iffree=false; /*free surface*/
    if(! sf_getint("nb",&nb)) nb=2; /*boundary condition width*/
    
	sf_putint(img1,"n1",nz);	
	sf_putint(img1,"n2",nx);
	sf_putfloat(img1,"d1",dz);
	sf_putfloat(img1,"d2",dx);
	sf_putstring(img1,"label1","Depth");
	sf_putstring(img1,"label2","Position");
	sf_putint(img2,"n1",nz);	
	sf_putint(img2,"n2",nx);
	sf_putfloat(img2,"d1",dz);
	sf_putfloat(img2,"d2",dx);
	sf_putstring(img2,"label1","Depth");
	sf_putstring(img2,"label2","Position");
	sf_putint(illums,"n1",nz);	
	sf_putint(illums,"n2",nx);
	sf_putfloat(illums,"d1",dz);
	sf_putfloat(illums,"d2",dx);
	sf_putint(illums,"n3",ns);
	sf_putint(illums,"d3",1);
	sf_putint(illums,"o3",1);

    /*below output the reflection data for comparison purpose*/
    sf_file refs;
    refs=sf_output("refs");
	sf_putint(refs,"n1",ng);	
	sf_putint(refs,"n2",nt);
	sf_putfloat(refs,"d1",jgx*dx);
	sf_putfloat(refs,"d2",dt);
	sf_putint(refs,"n3",ns);
	sf_putint(refs,"d3",1);
	sf_putint(refs,"o3",1);    
    /*above output the reflection data for comparison purpose*/
    
	roll=(ifr>0)?true:false;

	vm=sf_floatalloc2(nz, nx);/* migration velocity */
	vw=sf_floatalloc2(nz, nx);/* water velocity */
	imag1=sf_floatalloc2(nz, nx);/* RTM image without source illumination */
	imag2=sf_floatalloc2(nz, nx);/* RTM image with source illumination */
	
	sp0=sf_floatalloc2(nz, nx);/* source wavefield p0 */
	sp1=sf_floatalloc2(nz, nx);/* source wavefield p1 */
	sp2=sf_floatalloc2(nz, nx);/* source wavefield p2 */
	gp0=sf_floatalloc2(nz, nx);/* geophone/receiver wavefield p0 */
	gp1=sf_floatalloc2(nz, nx);/* geophone/receiver wavefield p1 */
	gp2=sf_floatalloc2(nz, nx);/* geophone/receiver wavefield p2 */
	lap=sf_floatalloc2(nz, nx);/* laplace of the source wavefield */
	illum=sf_floatalloc2(nz, nx);/* illumination of the source wavefield */
	imagtmp=sf_floatalloc2(nz, nx);/* the numerator (tmp image) for illumination */	
	wlt=(float*)malloc(nt*sizeof(float));/* ricker wavelet */
	sxz=(int*)malloc(ns*sizeof(int)); /* source positions */
	gxz=(int*)malloc(ng*sizeof(int)); /* geophone positions */
	bd=(float*)malloc(nt*(4*(nz+nx))*sizeof(float));/* boundaries for wavefield reconstruction */
	trans=(float*)malloc(ng*nt*sizeof(float));/* transposed one shot */
	dobs=(float*)malloc(ng*nt*sizeof(float));/* observed seismic data */
	dcal=(float*)malloc(ng*sizeof(float));	/* calculated/synthetic seismic data */
	dref=(float*)malloc(ns*ng*nt*sizeof(float));/* residual/error between synthetic and observation */

	/* initialize varibles */
	sf_floatread(vm[0], nz*nx, vmig);
	sf_floatread(vw[0], nz*nx, vwater);	
	memset(sp0[0], 0, nz*nx*sizeof(float));
	memset(sp1[0], 0, nz*nx*sizeof(float));
	memset(sp2[0], 0, nz*nx*sizeof(float));
	memset(gp0[0], 0, nz*nx*sizeof(float));
	memset(gp1[0], 0, nz*nx*sizeof(float));
	memset(gp2[0], 0, nz*nx*sizeof(float));
	memset(imag1[0], 0, nz*nx*sizeof(float));
	memset(imag2[0], 0, nz*nx*sizeof(float));
	memset(lap[0], 0, nz*nx*sizeof(float));
	
	for(it=0;it<nt;it++){
	    tmp=SF_PI*fm*(it*dt-1.0/fm);tmp*=tmp;
		/*tmp=SF_PI*fm*(it*dt-1.0/fm);tmp*=tmp;*//*kt=1.0/fm/dt*/
		wlt[it]=(1.0-2.0*tmp)*expf(-tmp);
	}
	
	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz))	
	{ sf_warning("sources exceeds the computing zone!\n"); exit(1);}
	sg_position(sxz, szbeg, sxbeg, jsz, jsx, ns, nz);
	distx=sxbeg-gxbeg;
	distz=szbeg-gzbeg;
	if (roll)	{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz &&
		(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz))	
		{ sf_warning("geophones exceeds the computing zone!\n"); exit(1); }
	}
	else{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
		{ sf_warning("geophones exceeds the computing zone!\n"); exit(1); }
	}
	sg_position(gxz, gzbeg, gxbeg, jgz, jgx, ng, nz);
	memset(bd, 0, nt*(2*nz+nx)*sizeof(float));
	memset(dobs, 0, ng*nt*sizeof(float));
	memset(dcal, 0, ng*sizeof(float));
	
		if(verb) start=clock();// record starting time
		sf_seek(shots, 0L, SEEK_SET);		
		for(is=0;is<ns;is++)
		{
			sf_floatread(trans, ng*nt, shots);
			transpose(trans, dobs, nt, ng);
			if (roll)	{
				gxbeg=sxbeg+is*jsx-distx;
				sg_position(gxz, gzbeg, gxbeg, jgz, jgx, ng, nz);
			}
			memset(sp0[0], 0, nz*nx*sizeof(float));
			memset(sp1[0], 0, nz*nx*sizeof(float));
			
			/*The following is used to subtract direct waves*/
			for(it=0; it<nt; it++)
			{
				source(sp1, &wlt[it], &sxz[is], 1, nz, true);			
				step_forward(sp0, sp1, sp2, vw, dt, dz, dx, nz, nx, nb, iffree);
				ptr=sp0; sp0=sp1; sp1=sp2; sp2=ptr; /*switch wavefield pointer at t-1,t,t+1*/
				record(dcal, gxz, sp0, ng, nz);     /*record the wavefield from water velocity*/
				cal_residual(dcal, &dobs[it*ng], &dref[is*ng*nt+it*ng], ng);/*get reflection data*/
			}

			
			/*The above is used to subtract direct waves*/
			
			/*The following is used to save source wavefield to boundary*/
			memset(sp0[0], 0, nz*nx*sizeof(float));/*first re-initialize the wavefield*/
	        memset(sp1[0], 0, nz*nx*sizeof(float));/*first re-initialize the wavefield*/
	        memset(sp2[0], 0, nz*nx*sizeof(float));/*first re-initialize the wavefield*/
	      //  memset(dcal, 0, ng*sizeof(float));/*first re-initialize the recorded data*/
	        
			for(it=0; it<nt; it++)
			{
				source(sp1, &wlt[it], &sxz[is], 1, nz, true);			
				step_forward(sp0, sp1, sp2, vm, dt, dz, dx, nz, nx, nb, iffree);
				ptr=sp0; sp0=sp1; sp1=sp2; sp2=ptr;
				
				/*save boundary*/
				rw_boundary(&bd[it*4*(nx+nz)], sp0, nz, nx, true);
	
			}
			/*The above is used to save source wavefield to boundary*/

			
			ptr=sp0; sp0=sp1; sp1=ptr;/*bd=u_n;sp0=u_n+1;sp1=bd=u_n*/
			memset(gp0[0], 0, nz*nx*sizeof(float));
			memset(gp1[0], 0, nz*nx*sizeof(float));
	        memset(illum[0], 0, nz*nx*sizeof(float));
	        memset(imagtmp[0], 0, nz*nx*sizeof(float));
			for(it=nt-1; it>-1; it--)
			{
				/* read boundary for time-reversal modeling */
				rw_boundary(&bd[it*4*(nx+nz)], sp1, nz, nx, false);
				step_backward(illum, lap, sp0, sp1, sp2, vm, dt, dz, dx, nz, nx, nb);
				
				source(gp1, &dref[is*ng*nt+it*ng], gxz, ng, nz, true);
				step_forward(gp0, gp1, gp2, vm, dt, dz, dx, nz, nx, nb, iffree);

//				xcorr(imagtmp, lap, gp1, nz, nx);/*temporary image for one shot*/
				xcorr(imagtmp, sp1, gp1, nz, nx);/*temporary image for one shot*/				
				ptr=sp0; sp0=sp1; sp1=sp2; sp2=ptr;
				ptr=gp0; gp0=gp1; gp1=gp2; gp2=ptr;
			}
			
			/*imagtmp is the numerator and illum is the denominator*/
			source_illum(imag2, imag1, imagtmp, illum, nz, nx, nb);	 /*apply source illumination*/
			
		    /* write the source illumination part */
		    sf_floatwrite(illum[0], nz*nx, illums);
		    if(verb) {// output important computation information 
			sf_warning("Shot %d is finished",is+1);
		    }
		}
			
		/*write the current velocity*/
		sf_floatwrite(imag1[0], nz*nx, img1);
		sf_floatwrite(imag2[0], nz*nx, img2);
		
		/*output reflection data*/
		sf_floatwrite(dref,ng*nt*ns,refs);		
		
		if(verb) {// output important computation information 
			stop=clock();// record ending time 
			sf_warning("RTM is finished in %f s",((float)(stop-start))/CLOCKS_PER_SEC);
		}

	free(*vm); free(vm);
	free(*sp0); free(sp0);
	free(*sp1); free(sp1);
	free(*sp2); free(sp2);
	free(*gp0); free(gp0);
	free(*gp1); free(gp1);
	free(*gp2); free(gp2);
	free(*imag1); free(imag1);
	free(*imag2); free(imag2);
	free(*lap); free(lap);
	free(*illum); free(illum);
	free(*imagtmp); free(imagtmp);	
	free(wlt);
	free(sxz);
	free(gxz);
	free(bd);
	free(trans);
	free(dobs);
	free(dcal);
	free(dref);


	exit(0);
}


