/* Simple acoustic wave propagation code
#
# Originally written by Yangkang Chen
# Modified by Shan Qu
# 
# CODE for DTV regularized FWI
# Reference:
# Qu, S., E. Verschuur, and Y. Chen, 2019, FWI/JMI with an automatic directional total variation constraint, Geophysics, 84, R175-R183.
#
*/

#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <rsf.h>

// #include <iostream>
// #include <complex>
// #include <sstream>
// #include <assert.h>

#define M_PI 3.1415926535897932385E0
// typedef struct { float Re; float Im; }  complex;         /* for complex number */

float **alloc2D(int size_x, int size_y);
void delete2D(float **arr, int size_x, int size_y);
void sg_position(int *sxz, int sxbeg, int jsx, int ns);
void step_forward(float **p0, float **p1, float **p2, float **vv, float dt, float dz, float dx, int nz, int nx, int nb, int nd, bool iffree);
void oneway_abc(float **uo, float **um, float **vv, int nx, int nz, int nb, int nd, float dx, float dz, float dt, bool iffree);
void save_2d_float(float **p, int n1, int n2, char *filename);
void save_1d_float(float *p, int n1, char *filename);
void save_1d_int(int *p, int n1, char *filename);
void rw_boundary(float *bd, float **p, int nz, int nx, int nb, int nd, bool write);
void step_backward(float **illum, float **lap, float **p0, float **p1, float **p2, float **vv, float dt, float dz, float dx, int nz, int nx, int nb, int nd);
void pre_gradient(float **grad, float **vv, float **illum, int nz, int nx, int nb, int nd, float stab, bool precon);
void smoothz(float **g, float **smg, int rb, int nz, int nx);
void smoothx(float **g, float **smg, int rb, int nz, int nx);
float max(float x, float y);
float min(float x, float y);
float cal_beta(float **g0, float **g1, float **cg, int nz, int nx);
float cal_eps(float **vv, float **cg, float stab, int nz, int nx);
void sum_alpha(float *alpha1, float *alpha2, float *dcaltmp, float *dobs, float *derr, int ng);
float cal_alpha(float *alpha1, float *alpha2, float eps, float stab, int ng);
// std::string IntToString (int a);
//void fft(complex x[], int n, int flag);
//void fft2D(complex x[], int n1, int n2, int flag);
//void fft3D(complex x[], int n1, int n2, int n3, int flag);
//void cooley_tukey(complex x[], int n, int flag, int n2);
//static void stockham(complex x[], int n, int flag, int n2, complex y[]);

void TVfilter(float **tvx, float **tvz, float **vv, float **theta, float alp, int nx, int nz, int ndb);
void TVTfilter(float **tvx, float **tvz, float **dtv, float **theta, float alp, int nx, int nz, int ndb);
float get_coef(float **dv, float ** dtv, float mu, int nx, int nz, float stab);
void shrink(float **tvx, float **b1, float **d1, float mu, int nx, int nz);
void DTV(float **vv, float **dv, float **dtv, float **theta, float **b1, float **b2, float **d1, float **d2, float alp, float mu, float lambda, int nx, int nz, int ndb, float stab);

int main (int argc, char* argv[]) 
{    

	sf_file f_vel, f_shots;
	FILE *fp_input;

	sf_init (argc, argv);
	
	f_vel=sf_input("in");
	f_shots=sf_output("out");
	
	int ix, iz, it, is, ig;
	int nxo, nzo;
	int nb, nd;
	int nx, nz;
	float dx, dz, dt;
	int nt;
	int ng, ns, jsx, jgx, sxbeg, szbeg, gxbeg, gzbeg;
	bool roll, iffree,is_direct;
	float amp, fm, tmp;
	int distx, distz;
	char fullname[20];

	float** vel = NULL;
	float** vv = NULL;
	float** sp0 = NULL;		//source wavefield p0
	float** sp1 = NULL;	    //source wavefield p1
	float** sp2 = NULL;	    //source wavefield p2
	float** ptr = NULL;

	float* wlt = NULL;		//ricker wavelet
	int* sxz = NULL;		//source position
	int* gxz = NULL;		//receiver position

	float* dobs = NULL;		//observed seismic data
	float* shots = NULL;	// observed data

	int jsz,jgz;
	
// 	nxo = 334;
// 	nzo = 92;
// 	nb = 2;     // boundary condition width
// 	nd = 2;     // 2rd order fd
// 	nx = nxo + 2*(nb+nd);
// 	nz = nzo + 2*(nb+nd);
// 	dx = 12;
// 	dz = 12;

    if (!sf_histint(f_vel,"n2",&nxo)) sf_error("no n2");/* nx */
    if (!sf_histint(f_vel,"n1",&nzo)) sf_error("no n1");/* nx */
    if (!sf_histfloat(f_vel,"d2",&dx)) sf_error("no d2");/* nx */
    if (!sf_histfloat(f_vel,"d1",&dz)) sf_error("no n1");/* nx */
     
// 	nt = 2000;
// 	dt = 0.001;
// 	ng = 334;        // total receivers in each shot
// 	ns = 23;         // number of shots
// 	jsx = 15;        // source x jump interval
// 	jgx = 1;         // receiver x jump interval
// 	sxbeg = 2;       // x begining index of sources
// 	szbeg = 0;       // z begining index of sources
// 	gxbeg = 0;       // x begining index of receivers
// 	gzbeg = 0;       // z begining index of receivers

// 	roll = false;    // rolling survey
// 	iffree = false;  // free surface
	is_direct = false;
// 	amp = 1;         // max amp of ricker
// 	fm = 14;         // dominant freq of ricker

   	if (!sf_getint("nt",&nt)) sf_error("no nt");
	/* total modeling time steps */
   	if (!sf_getint("ng",&ng)) sf_error("no ng");
	/* total receivers in each shot */
   	if (!sf_getint("ns",&ns)) sf_error("no ns");
	/* number of shots */
   	if (!sf_getfloat("dt",&dt)) sf_error("no dt");
	/* time sampling interval */
	
 	if (!sf_getfloat("amp",&amp)) amp=1;
	/* maximum amplitude of ricker */
    	if (!sf_getfloat("fm",&fm)) fm=14;	
	/* dominant freq of ricker */
    	if (!sf_getint("jsx",&jsx))   sf_error("no jsx");
	/* source x-axis  jump interval  */
    	if (!sf_getint("jsz",&jsz))   jsz=0; /*not used*/
	/* source z-axis jump interval  */
    	if (!sf_getint("jgx",&jgx))   jgx=1;
	/* receiver x-axis jump interval */
    	if (!sf_getint("jgz",&jgz))   jgz=0; /*not used*/
	/* receiver z-axis jump interval */
    	if (!sf_getint("sxbeg",&sxbeg))   sf_error("no sxbeg");
	/* x-begining index of sources, starting from 0 */
    	if (!sf_getint("szbeg",&szbeg))   sf_error("no szbeg");
	/* z-begining index of sources, starting from 0 */
    	if (!sf_getint("gxbeg",&gxbeg))   sf_error("no gxbeg");
	/* x-begining index of receivers, starting from 0 */
    	if (!sf_getint("gzbeg",&gzbeg))   sf_error("no gzbeg");
	/* z-begining index of receivers, starting from 0 */
	if (!sf_getbool("roll",&roll)) roll=false;
	/* default, rolling acquisition*/
    if(! sf_getbool("free",&iffree)) iffree=false; /*free surface*/
    if(! sf_getint("nb",&nb)) nb=2; /*boundary condition width*/
    
    nd=2;
	nx = nxo + 2*(nb+nd);
	nz = nzo + 2*(nb+nd);
    
	sf_putint(f_shots,"n1",ng);
	sf_putint(f_shots,"n2",nt);
	sf_putint(f_shots,"n3",ns);
	sf_putfloat(f_shots,"d1",jgx*dx);
	sf_putfloat(f_shots,"d2",dt);
	sf_putfloat(f_shots,"d3",jsx*dx);
	sf_putfloat(f_shots,"o1",gxbeg*dx);
	sf_putfloat(f_shots,"o2",0);
	sf_putfloat(f_shots,"o3",sxbeg*dx);
	sf_putstring(f_shots,"label2","Time");
	sf_putstring(f_shots,"label1","Position");
	sf_putstring(f_shots,"label3","Shot");
	sf_putstring(f_shots,"unit2","s");
	sf_putstring(f_shots,"unit1","m");
	sf_putfloat(f_shots,"amp",amp);
	sf_putfloat(f_shots,"fm",fm);
	sf_putint(f_shots,"ng",ng);
	sf_putint(f_shots,"szbeg",szbeg);
	sf_putint(f_shots,"sxbeg",sxbeg);
	sf_putint(f_shots,"gzbeg",gzbeg);
	sf_putint(f_shots,"gxbeg",gxbeg);
	sf_putint(f_shots,"jsx",jsx);
	sf_putint(f_shots,"jsz",jsz);
	sf_putint(f_shots,"jgx",jgx);
	sf_putint(f_shots,"jgz",jgz);
	sf_putint(f_shots,"roll",roll?1:0);
	
	szbeg=0;
	gzbeg=0;
	
    vel = alloc2D(nxo, nzo);
	vv = alloc2D(nx, nz);
	sp0 = alloc2D(nx, nz);
	sp1 = alloc2D(nx, nz);
	sp2 = alloc2D(nx, nz);


	wlt = (float *)malloc(nt * sizeof(float));
	sxz = (int *)malloc(ns * sizeof(int));
	gxz = (int *)malloc(ng * sizeof(int));
	dobs = (float *)malloc(ng*nt * sizeof(float));
	shots = (float *)malloc(ng*nt*ns * sizeof(float));


	/* read vv from vel.bin and extend it*/
 	if ((fp_input = fopen("vel.bin", "rb")) == NULL)
 	{
 	    printf("Cannot open inputfile.\n");
 	    exit(0);
 	}
 	
	
 	for (iz=0; iz<nzo; iz++)
 	    for (ix=0; ix<nxo; ix++)
 			fread(&(vel[ix][iz]), sizeof(float), 1, fp_input);
 	fclose(fp_input);
	
//	sf_floatread(vel[0],nzo*nxo,f_vel);

	for (ix=nb+nd; ix<nx-nb-nd; ix++)
		for (iz=nb+nd; iz<nz-nb-nd; iz++)
			vv[ix][iz] = vel[ix-(nb+nd)][iz-(nb+nd)];
	
	for (ix=nb+nd; ix<nx-nb-nd; ix++)
	{
		for (iz=0; iz<nb+nd; iz++)
			vv[ix][iz] = vv[ix][nb+nd];
		for (iz=nz-nb-nd; iz<nz; iz++)
			vv[ix][iz] = vv[ix][nz-nb-nd-1];
	}
	for (iz=nb+nd; iz<nz-nb-nd; iz++)
	{
		for (ix=0; ix<nb+nd; ix++)
			vv[ix][iz] = vv[nb+nd][iz];
		for (ix=nx-nb-nd; ix<nx; ix++)
			vv[ix][iz] = vv[nx-nb-nd-1][iz];
	}
	for (ix=0; ix<nb+nd; ix++)
	{
		for (iz=0; iz<nb+nd; iz++)
			vv[ix][iz] = vv[nb+nd][nb+nd];
		for (iz=nz-nb-nd; iz<nz; iz++)
			vv[ix][iz] = vv[nb+nd][nz-nb-nd-1];
	}
	for (ix=nx-nb-nd; ix<nx; ix++)
	{
		for (iz=0; iz<nb+nd; iz++)
			vv[ix][iz] = vv[nx-nb-nd-1][nb+nd];
		for (iz=nz-nb-nd; iz<nz; iz++)
			vv[ix][iz] = vv[nx-nb-nd-1][nz-nb-nd-1];
	}

	if (is_direct)
	for (ix=0; ix<nx; ix++)
		for (iz=0; iz<nz;iz++)
			vv[ix][iz] = 1500.0;


	/* initialize varibles */
	for (ix=0; ix<nx; ix++)
		for (iz=0; iz<nz;iz++)
		{
			sp0[ix][iz] = 0.0;
			sp1[ix][iz] = 0.0;
			sp2[ix][iz] = 0.0;
		}

	for (ix=0; ix<ng*nt; ix++)
		dobs[ix] = 0.0;

	/* define wavelet */
	for(it=0; it<nt; it++)
	{
		tmp = M_PI * fm * (it*dt-1.0/fm);
		tmp *= tmp;
		wlt[it] = (1.0-2.0*tmp) * expf(-tmp);
	}
	

	if (!(sxbeg>=0 && szbeg>=0 && 
		sxbeg+(ns-1)*jsx<nxo && szbeg<nzo))	
	{ 
		printf("sources exceeds the computing zone!\n"); 
		exit(1);
	}

	sg_position(sxz, sxbeg, jsx, ns);


	distx = sxbeg - gxbeg;

	if (roll)	
	{
		if (!(gxbeg>=0 && gzbeg>=0 && 
			gxbeg+(ng-1)*jgx<nxo && gzbeg<nzo && 
			(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nxo))
		{ 
			printf("geophones exceeds the computing zone!\n"); 
			exit(1); 
		}
	}
	else
	{
		if (!(gxbeg>=0 && gzbeg>=0 && 
			gxbeg+(ng-1)*jgx<nxo && gzbeg<nzo))	
		{ 
			printf("geophones exceeds the computing zone!\n"); 
			exit(1); 
		}
	}
	sg_position(gxz, gxbeg, jgx, ng);




	for (is=0; is<ns; is++)
	{
		sf_warning("is=%d \n", is);
		if (roll)
		{
			gxbeg=sxbeg+is*jsx-distx;
			sg_position(gxz, gxbeg, jgx, ng);
		}
		for (ix=0; ix<nx; ix++)
			for (iz=0; iz<nz;iz++)
			{
				sp0[ix][iz] = 0.0;
				sp1[ix][iz] = 0.0;
			}

		for (it=0; it<nt; it++)
		{
			/* add source*/
			sp1[sxz[is]+nb+nd][szbeg+nb+nd] += wlt[it];
			step_forward(sp0, sp1, sp2, vv, dt, dz, dx, nz, nx, nb, nd, iffree);
			ptr = sp0; 
			sp0 = sp1; 
			sp1 = sp2; 
			sp2 = ptr;

			for (ig=0; ig<ng; ig++)
			{
				dobs[it*ng+ig] = sp0[gxz[ig]+nb+nd][gzbeg+nb+nd];
				shots[is*ng*nt+it*ng+ig] =  dobs[it*ng+ig];
			}
			
		}
	}
	// save data
// 	if (is_direct)
// 	{		
// 		snprintf(fullname,20,"shots_direct.bin"); /*note I change here*/
// 		save_1d_float(shots, ng*nt*ns, fullname);
// 	}
// 	else
// 	{
// 		snprintf(fullname,20,"shots.bin"); /*note I change here*/
// 		save_1d_float(shots, ng*nt*ns, fullname);
// 	}



	sf_floatwrite(shots,ng*nt*ns,f_shots);


	// clear memory
	delete2D(vel, nxo, nzo);
	delete2D(vv, nx, nz);
	delete2D(sp0, nx, nz);
	delete2D(sp1, nx, nz);
	delete2D(sp2, nx, nz);

	free(wlt);
	free(sxz);
	free(gxz);

	free(dobs);
	free(shots);

	return 0;
}

float **alloc2D(int size_x, int size_y)
{
	int i;
// 	float **arr = NULL;

// 	arr = new float *[size_x];
// 
// 	for (i=0; i<size_x; i++)
// 		arr[i] = new float [size_y];


    float **arr = (float **)malloc(size_x * sizeof(float *));
	arr[0] = (float *)malloc(size_x * size_y * sizeof(float));
	for (i=1; i < size_x; i++) {
		arr[i] = arr[0] + i * size_y;
    }
	return arr;
}

void delete2D(float **arr, int size_x, int size_y)
{
	free(*arr);
	free(arr);
}

void sg_position(int *sxz, int sxbeg, int jsx, int ns)
/*<  shot/geophone positions >*/
{
	int is, sx;
	for (is=0; is<ns; is++)
	{
		sx = sxbeg+is*jsx;
		sxz[is] = sx;
	}
}

void oneway_abc(float **uo, float **um, float **vv, int nx, int nz, int nb, int nd, float dx, float dz, float dt, bool iffree)
/* oneway - absorbing condition */
/* Reference: Absorbing boundary conditions for elastic waves, R.L. Higdon, Geophysics, 1991, 56, 2, 231-241.
Reverse time migration using analytical wavefield and wavefield decomposition imaging condition, Enjiang Wang and Yang Liu, SEG abstract (86th), 2016, 4461-4465 */
{
    int iz,ix,ib;
    float q;
    float qx,qt,qxt,b;

	b = 0.5;   //b=0.5 is the best.

    for (ix=nd; ix<nx-nd; ix++)
	{
		for (ib=nd; ib<nb+nd; ib++)
		{
			/* top BC */
			if (!iffree) //not free surface, apply ABC
			{
				iz = nb + nd - ib;
				q = vv[ix][nb] * dt / dz;
				qx = (b*(1+q)-q) / (1+q) / (1-b);
				qt = (b*(1+q)-1) / (1+q) / (1-b);
				qxt = b / (b-1);

 				uo[ix][iz] = um[ix][iz+1] * (-qxt) +
							 um[ix][iz  ] * (-qt)  +
							 uo[ix][iz+1] * (-qx);
			}

			/* bottom BC */
			iz = nz - nb - nd + ib - 1;
			q = vv[ix][nz-nb-1] * dt / dz;

			qx = (b*(1+q)-q) / (1+q) / (1-b);
			qt = (b*(1+q)-1) / (1+q) / (1-b);
			qxt = b / (b-1);

 			uo[ix][iz] = um[ix][iz-1] * (-qxt) +
						 um[ix][iz  ] * (-qt)  +
						 uo[ix][iz-1] * (-qx);

		}
    }

    for (iz=nd; iz<nz-nd; iz++)
	{
		for (ib=nd; ib<nb+nd; ib++)
		{
			/* left BC */
			ix = nb + nd - ib;
			q = vv[nb][iz] * dt / dx;
			qx = (b*(1+q)-q) / (1+q) / (1-b);
			qt = (b*(1+q)-1) / (1+q) / (1-b);
			qxt = b / (b-1);

			uo[ix][iz] = um[ix+1][iz] * (-qxt) +
						 um[ix  ][iz] * (-qt)  +
						 uo[ix+1][iz] * (-qx);


			/* right BC */
			ix = nx - nb - nd + ib - 1;
			q = vv[nx-nb-1][iz] * dt / dx;
			qx = (b*(1+q)-q) / (1+q) / (1-b);
			qt = (b*(1+q)-1) / (1+q) / (1-b);
			qxt = b / (b-1);

			uo[ix][iz] = um[ix-1][iz] * (-qxt) +
						 um[ix  ][iz] * (-qt)  +
						 uo[ix-1][iz] * (-qx);
		}
    }

}


void step_forward(float **p0, float **p1, float **p2, float **vv, float dt, float dz, float dx, int nz, int nx, int nb, int nd, bool iffree)
/*< forward modeling step, Clayton-Enquist ABC incorporated >*/
{
    int ix,iz;
    float tmp;
    float idx,idz;
	float c0, c1, c2;
	c0 = -30./12.;
	c1 = +16./12.;
	c2 = - 1./12.;
    idx = 1.0 / (dx*dx);
	idz = 1.0 / (dz*dz);

	for (iz=nd; iz<nz-nd; iz++)
	{
	    for (ix=nd; ix<nx-nd; ix++)
		{
			tmp =
		    c0* p1[ix  ][iz  ] * (idx+idz) +
		    c1*(p1[ix-1][iz  ] + p1[ix+1][iz  ])*idx +
		    c2*(p1[ix-2][iz  ] + p1[ix+2][iz  ])*idx +
		    c1*(p1[ix  ][iz-1] + p1[ix  ][iz+1])*idz +
		    c2*(p1[ix  ][iz-2] + p1[ix  ][iz+2])*idz;

		    p2[ix][iz] = 2*p1[ix][iz] - p0[ix][iz] +
				vv[ix][iz]*vv[ix][iz]*dt*dt*tmp;

	    }
	}
	oneway_abc(p2, p1, vv, nx, nz, nb, nd, dx, dz, dt, iffree);

}

void save_2d_float(float **p, int n1, int n2, char *filename)
{
	int i,j;
	FILE *fp_write;
	fp_write = fopen(filename, "wb");
	for (i=0; i<n1; i++)
		for (j=0; j<n2; j++)
			fwrite(&(p[i][j]), sizeof(float), 1, fp_write);
	fclose(fp_write);
}
void save_1d_float(float *p, int n1, char *filename)
{
	int i;
	FILE *fp_write;
	fp_write = fopen(filename, "wb");
	for (i=0; i<n1; i++)
		fwrite(&(p[i]), sizeof(float), 1, fp_write);
	fclose(fp_write);
}

void save_1d_int(int *p, int n1, char *filename)
{
	int i;
	FILE *fp_write;
	fp_write = fopen(filename, "wb");
	for (i=0; i<n1; i++)
		fwrite(&(p[i]), sizeof(int), 1, fp_write);
	fclose(fp_write);
}


void rw_boundary(float *bd, float **p, int nz, int nx, int nb, int nd, bool write)
/*< if write==true, write/save boundaries out of variables;
 else  read boundaries into wavefields (for 4th order FD) >*/
{
	int ix,iz;
	if(write)
	{
		for (ix=0; ix<nx; ix++)
			for (iz=0; iz<nb; iz++)
                {
                        bd[iz   +2*nb*ix] = p[ix][iz];
                        bd[iz+nb+2*nb*ix] = p[ix][nz-1-iz];
                }
        for (iz=0; iz<nz; iz++)
            for (ix=0; ix<nb; ix++)
                {
                        bd[2*nb*nx+iz+nz*ix     ] = p[ix][iz];
                        bd[2*nb*nx+iz+nz*(ix+nb)] = p[nx-1-ix][iz];
                }

		/*for (ix=nd; ix<nx-nd; ix++)
		{
			for (iz=0; iz<nd; iz++)
				bd[2*nb*nx+iz+nz*ix] = 0.0;
			for (iz=nz-nd; iz<nz; iz++)
				bd[2*nb*nx+iz+nz*ix] = 0.0;
		}
		for (iz=nd; iz<nz-nd; iz++)
		{
			for (ix=0; ix<nd; ix++)
				bd[2*nb*nx+iz+nz*ix] = 0.0;
			for (ix=nx-nd; ix<nx; ix++)
				bd[2*nb*nx+iz+nz*ix] = 0.0;
		}
		for (ix=0; ix<nd; ix++)
		{
			for (iz=0; iz<nd; iz++)
				bd[2*nb*nx+iz+nz*ix] = 0.0;
			for (iz=nz-nd; iz<nz; iz++)
				bd[2*nb*nx+iz+nz*ix] = 0.0;
		}
		for (ix=nx-nd; ix<nx; ix++)
		{
			for (iz=0; iz<nd; iz++)
				bd[2*nb*nx+iz+nz*ix] = 0.0;
			for (iz=nz-nd; iz<nz; iz++)
				bd[2*nb*nx+iz+nz*ix] = 0.0;
		}*/
	}
	else
	{
		for (ix=0; ix<nx; ix++)
			for (iz=0; iz<nb; iz++)
                {
                        p[ix][iz]      = bd[iz   +2*nb*ix];
                        p[ix][nz-1-iz] = bd[iz+nb+2*nb*ix];
                }
        for (iz=0; iz<nz; iz++)
			for (ix=0; ix<nb; ix++)
                {
                        p[ix][iz]      = bd[2*nb*nx+iz+nz*ix];
                        p[nx-1-ix][iz] = bd[2*nb*nx+iz+nz*(ix+nb)];
                }
	}
}


void step_backward(float **illum, float **lap, float **p0, float **p1, float **p2, float **vv, float dt, float dz, float dx, int nz, int nx, int nb, int nd)
/*< backward modeling step, Clayton-Enquist ABC incorporated >*/
{
    int ix,iz;
    float tmp;
    float idx,idz;
	float c0, c1, c2;
	c0 = -30./12.;
	c1 = +16./12.;
	c2 = - 1./12.;
    idx = 1.0 / (dx*dx);
	idz = 1.0 / (dz*dz);

	for (iz=nd; iz<nz-nd; iz++)
		for (ix=nd; ix<nx-nd; ix++)
		{
			tmp =
		    c0* p1[ix  ][iz  ] * (idx+idz) +
		    c1*(p1[ix-1][iz  ] + p1[ix+1][iz  ])*idx +
		    c2*(p1[ix-2][iz  ] + p1[ix+2][iz  ])*idx +
		    c1*(p1[ix  ][iz-1] + p1[ix  ][iz+1])*idz +
		    c2*(p1[ix  ][iz-2] + p1[ix  ][iz+2])*idz;
			lap[ix][iz] = tmp;

		    p2[ix][iz] = 2*p1[ix][iz] - p0[ix][iz] +
				vv[ix][iz]*vv[ix][iz]*dt*dt*tmp;

			illum[ix][iz] += p1[ix][iz]*p1[ix][iz];
	    }
}


void pre_gradient(float **grad, float **vv, float **illum, int nz, int nx, int nb, int nd, float stab, bool precon)
/*< scale gradient >*/
{
	int ix, iz;
	float a;

	for (ix=nb+nd; ix<nx-nb-nd; ix++)
	{
		for (iz=nb+nd; iz<nz-nb-nd; iz++)
		{
			a=vv[ix][iz];
			if (precon)
				a = a * sqrtf(illum[ix][iz]+stab); /*precondition with wavefield illumination*/

			grad[ix][iz] = grad[ix][iz] * 2.0/a;
		}
	}

	for (ix=nb+nd; ix<nx-(nb+nd); ix++)
	{
		for (iz=0; iz<nb+nd; iz++)
			grad[ix][iz] = grad[ix][nb+nd];
		for (iz=nz-(nb+nd); iz<nz; iz++)
			grad[ix][iz] = grad[ix][nz-(nb+nd)-1];
	}
	for (iz=nb+nd; iz<nz-(nb+nd); iz++)
	{
		for (ix=0; ix<nb+nd; ix++)
			grad[ix][iz] = grad[nb+nd][iz];
		for (ix=nx-(nb+nd); ix<nx; ix++)
			grad[ix][iz] = grad[nx-(nb+nd)-1][iz];
	}
	for (ix=0; ix<nb+nd; ix++)
	{
		for (iz=0; iz<nb+nd; iz++)
			grad[ix][iz] = grad[nb+nd][nb+nd];
		for (iz=nz-(nb+nd); iz<nz; iz++)
			grad[ix][iz] = grad[nb+nd][nz-(nb+nd)-1];
	}
	for (ix=nx-(nb+nd); ix<nx; ix++)
	{
		for (iz=0; iz<nb+nd; iz++)
			grad[ix][iz] = grad[nx-(nb+nd)-1][nb+nd];
		for (iz=nz-(nb+nd); iz<nz; iz++)
			grad[ix][iz] = grad[nx-(nb+nd)-1][nz-(nb+nd)-1];
	}



	//for (ix=nd; ix<nx-nd; ix++)
	//{
	//	for (iz=nd; iz<nz-nd; iz++)
	//	{
	//		a=vv[ix][iz];
	//		if (precon)
	//			a *= sqrtf(illum[ix][iz] + stab); /*precondition with wavefield illumination*/
	//		grad[ix][iz] *= 2.0/a;
	//	}
	//}

	//for (ix=nd; ix<nx-nd; ix++)
	//{
	//	for (iz=0; iz<nd; iz++)
	//		grad[ix][iz] = grad[ix][nd];
	//	for (iz=nz-nd; iz<nz; iz++)
	//		grad[ix][iz] = grad[ix][nz-nd-1];
	//}
	//for (iz=nd; iz<nz-nd; iz++)
	//{
	//	for (ix=0; ix<nd; ix++)
	//		grad[ix][iz] = grad[nd][iz];
	//	for (ix=nx-nd; ix<nx; ix++)
	//		grad[ix][iz] = grad[nx-nd-1][iz];
	//}
	//for (ix=0; ix<nd; ix++)
	//{
	//	for (iz=0; iz<nd; iz++)
	//		grad[ix][iz] = grad[nd][nd];
	//	for (iz=nz-(nb+nd); iz<nz; iz++)
	//		grad[ix][iz] = grad[nd][nz-nd-1];
	//}
	//for (ix=nx-nd; ix<nx; ix++)
	//{
	//	for (iz=0; iz<nb+nd; iz++)
	//		grad[ix][iz] = grad[nx-nd-1][nb+nd];
	//	for (iz=nz-nd; iz<nz; iz++)
	//		grad[ix][iz] = grad[nx-nd-1][nz-nd-1];
	//}

}



void smoothz(float **g, float **smg, int rb, int nz, int nx)
/*< gaussian bell smoothing for z-axis >*/
{
	int ix, iz, i;
	float s;

	for (ix=0; ix<nx; ix++)
		for (iz=0; iz<nz; iz++)
		{
			s = 0.0;
			for (i=-rb; i<=rb; i++)
				if (iz+i>=0 && iz+i<nz)
					s += expf(-(2.0*i*i)/rb) * g[ix][iz+i];
			smg[ix][iz] = s;
		}
}

void smoothx(float **g, float **smg, int rb, int nz, int nx)
/*< gaussian bell smoothing for x-axis >*/
{
	int ix, iz, i;
	float s;

	for (ix=0; ix<nx; ix++)
		for (iz=0; iz<nz; iz++)
		{
			s = 0.0;
			for (i=-rb; i<=rb; i++)
				if(ix+i>=0 && ix+i<nx)
					s+=expf(-(2.0*i*i)/rb)*g[ix+i][iz];
			smg[ix][iz] = s;
		}
}


float cal_beta(float **g0, float **g1, float **cg, int nz, int nx)
/*< calculate beta >*/
{
	int ix, iz;
	float a,b,c,d;

	a = b = c = d = 0;

	for (ix=0; ix<nx; ix++)
		for (iz=0; iz<nz; iz++)
		{
			a += (-g1[ix][iz]) * (-g1[ix][iz]-(-g0[ix][iz]));// numerator of HS
			b += cg[ix][iz] * (-g1[ix][iz]-(-g0[ix][iz]));// denominator of HS,DY
			c += g1[ix][iz] * g1[ix][iz];		// numerator of DY
			d += g0[ix][iz] * g0[ix][iz];		// denominator of PR
		}

	float beta_HS = -a/b;
	beta_HS	= max(0.0, beta_HS);

	float beta_DY = -c/b;
	beta_DY = max(0.0, beta_DY);

	float beta_PR = a/d;
	beta_PR = max(0.0, beta_PR);

	return max(0.0, min(beta_HS, beta_DY));
	//return beta_PR;
}


float max(float x, float y)
{
    if (x >= y)
    {
        return x;
    }
    return y;
}

float min(float x, float y)
{
    if (x <= y)
    {
        return x;
    }
    return y;
}



float cal_eps(float **vv, float **cg, float stab, int nz, int nx)
/*< calculate eps >*/
{
	int ix, iz;
	float vvmax, cgmax;
	vvmax = cgmax = 0.0;

	for(ix=0; ix<nx; ix++)
		for(iz=0; iz<nz; iz++)
		{
			vvmax = max(vvmax, fabsf(vv[ix][iz]));
			cgmax = max(cgmax, fabsf(cg[ix][iz]));
		}

	return 0.01 * vvmax / (cgmax + stab);
}


void sum_alpha(float *alpha1, float *alpha2, float *dcaltmp, float *dobs, float *derr, int ng)
/*< calculate numerator and denominator of alpha >*/
{
	int ig;
	float a, b, c;
	for (ig=0; ig<ng; ig++)
	{
		c = derr[ig];
		a = dobs[ig] + c;/* since f(mk)-dobs[id]=derr[id], thus f(mk)=b+c; */
		b = dcaltmp[ig] - a;/* f(mk+eps*cg)-f(mk) */
		alpha1[ig] -= b*c;
		alpha2[ig] += b*b;
	}
}

float cal_alpha(float *alpha1, float *alpha2, float eps, float stab, int ng)
/*< calculate alpha >*/
{
	int ig;
	float a, b;

	a = b = 0;
	for (ig=0; ig<ng; ig++)
	{
		a += alpha1[ig];
		b += alpha2[ig];
	}

	return (a*eps/(b+stab));
}

// std::string IntToString (int a)
// {
//     std::ostringstream temp;
//     temp << a;
//     return temp.str();
// }



//void bandlimit_vel_gradient(float **grad, float fmax, float **vel, int nx, int nz, float dx, float dz)
//{
//	float dkx, dkz;
//	int nkx, nkz;
//	int ix, iz;
//	float vel_max;
//	float **Filterx, **Filterz;
//
//
//	nkx = nx/2;
//	nkz = nz/2;
//
//	dkx = 2 * M_PI / (nx*dx);
//	dkz = 2 * M_PI / (nz*dz);
//
//	vel_min = 0.0;
//
//	for(ix=0; ix<nx; ix++)
//		for(iz=0; iz<nz; iz++)
//		{
//			vel_min = min(vel_min, vel[ix][iz]);
//		}
//
//	kmax = 2 * M_PI * fmax / vel_min;
//
//}


// static void stockham(complex x[], int n, int flag, int n2, complex y[])
// {
//    complex  *y_orig, *tmp;
//    int  i, j, k, k2, Ls, r, jrs;
//    int  half, m, m2;
//    float  wr, wi, tr, ti;
// 
//    y_orig = y;
//    r = half = n >> 1;
//    Ls = 1;                                         /* Ls=L/2 is the L star */
// 
//    while(r >= n2) {                              /* loops log2(n/n2) times */
//       tmp = x;                           /* swap pointers, y is always old */
//       x = y;                                   /* x is always for new data */
//       y = tmp;
//       m = 0;                        /* m runs over first half of the array */
//       m2 = half;                             /* m2 for second half, n2=n/2 */
//       for(j = 0; j < Ls; ++j) {
//          wr = cos(M_PI*j/Ls);                   /* real and imaginary part */
//          wi = -flag * sin(M_PI*j/Ls);                      /* of the omega */
//          jrs = j*(r+r);
//          for(k = jrs; k < jrs+r; ++k) {           /* "butterfly" operation */
//             k2 = k + r;
//             tr =  wr*y[k2].Re - wi*y[k2].Im;      /* complex multiply, w*y */
//             ti =  wr*y[k2].Im + wi*y[k2].Re;
//             x[m].Re = y[k].Re + tr;
//             x[m].Im = y[k].Im + ti;
//             x[m2].Re = y[k].Re - tr;
//             x[m2].Im = y[k].Im - ti;
//             ++m;
//             ++m2;
//          }
//       }
//       r  >>= 1;
//       Ls <<= 1;
//    };
// 
//    if (y != y_orig) {                     /* copy back to permanent memory */
//       for(i = 0; i < n; ++i) {               /* if it is not already there */
//          y[i] = x[i];               /* performed only if log2(n/n2) is odd */
//       }
//    }
// 
//    assert(Ls == n/n2);                        /* ensure n is a power of 2  */
//    assert(1 == n || m2 == n);           /* check array index within bound  */
// }


/* The Cooley-Tukey multiple column algorithm, see page 124 of Loan.
   x[] is input data, overwritten by output, viewed as n/n2 by n2
   array. flag = 1 for forward and -1 for backward transform.
*/
// void cooley_tukey(complex x[], int n, int flag, int n2)
// {
//    complex c;
//    int i, j, k, m, p, n1;
//    int Ls, ks, ms, jm, dk;
//    float wr, wi, tr, ti;
// 
//    n1 = n/n2;                               /* do bit reversal permutation */
//    for(k = 0; k < n1; ++k) {        /* This is algorithms 1.5.1 and 1.5.2. */
//       j = 0;
//       m = k;
//       p = 1;                               /* p = 2^q,  q used in the book */
//       while(p < n1) {
//          j = 2*j + (m&1);
//          m >>= 1;
//          p <<= 1;
//       }
//       assert(p == n1);                   /* make sure n1 is a power of two */
//       if(j > k) {
//          for(i = 0; i < n2; ++i) {                     /* swap k <-> j row */
//             c = x[k*n2+i];                              /* for all columns */
//             x[k*n2+i] = x[j*n2+i];
//             x[j*n2+i] = c;
//          }
//       }
//    }
//                                               /* This is (3.1.7), page 124 */
//    p = 1;
//    while(p < n1) {
//       Ls = p;
//       p <<= 1;
//       jm = 0;                                                /* jm is j*n2 */
//       dk = p*n2;
//       for(j = 0; j < Ls; ++j) {
//          wr = cos(M_PI*j/Ls);                   /* real and imaginary part */
//          wi = -flag * sin(M_PI*j/Ls);                      /* of the omega */
//          for(k = jm; k < n; k += dk) {                      /* "butterfly" */
//             ks = k + Ls*n2;
//             for(i = 0; i < n2; ++i) {                      /* for each row */
//                m = k + i;
//                ms = ks + i;
//                tr =  wr*x[ms].Re - wi*x[ms].Im;
//                ti =  wr*x[ms].Im + wi*x[ms].Re;
//                x[ms].Re = x[m].Re - tr;
//                x[ms].Im = x[m].Im - ti;
//                x[m].Re += tr;
//                x[m].Im += ti;
//             }
//          }
//          jm += n2;
//       }
//    }
// }






void TVfilter(float **tvx, float **tvz, float **vv, float **theta, float alp, int nx, int nz, int ndb)
{
	float h, a, b;
	float ux, uz;
	int ix, iz;
	h = 0.5;


	for (ix=ndb; ix<nx-ndb-1; ix++)
		for (iz=ndb; iz<nz-ndb-1; iz++)
		{
			a = cos(theta[ix][iz]);
			b = sin(theta[ix][iz]);
			//printf("a=%f, b=%f\n ", a, b);
			ux = vv[ix  ][iz  ] * h +
				 vv[ix  ][iz+1] * h +
				 vv[ix+1][iz  ] * (-h) +
				 vv[ix+1][iz+1] * (-h);
			uz = vv[ix  ][iz  ] * h +
				 vv[ix  ][iz+1] * (-h) +
				 vv[ix+1][iz  ] * h +
				 vv[ix+1][iz+1] * (-h);
			tvx[ix][iz] = alp * (a*ux+b*uz);
			tvz[ix][iz] = (2.0-alp) * (a*uz-b*ux);
		}


	for (ix=ndb; ix<nx-ndb-1; ix++)
	{
		tvx[ix][nz-ndb-1] = tvx[ix][nz-ndb-1-1];
		tvz[ix][nz-ndb-1] = tvz[ix][nz-ndb-1-1];
	}
	for (iz=ndb; iz<nz-ndb-1; iz++)
	{
		tvx[nx-ndb-1][iz] = tvx[nx-ndb-1-1][iz];
		tvz[nx-ndb-1][iz] = tvz[nx-ndb-1-1][iz];
	}

	tvx[nx-ndb-1][nz-ndb-1] = tvx[nx-ndb-1-1][nz-ndb-1-1];
	tvz[nx-ndb-1][nz-ndb-1] = tvz[nx-ndb-1-1][nz-ndb-1-1];

}

void TVTfilter(float **tvx, float **tvz, float **dtv, float **theta, float alp, int nx, int nz, int ndb)
{
	float h, a, b;
	float **ux, **uz;
	int ix, iz;
	h = 0.5;

	ux = alloc2D(nx, nz);
	uz = alloc2D(nx, nz);

	for (ix=0; ix<nx; ix++)
		for (iz=0; iz<nz;iz++)
		{
			ux[ix][iz] = 0.0;
			uz[ix][iz] = 0.0;
		}

	for (ix=ndb; ix<nx-ndb; ix++)
		for (iz=ndb; iz<nz-ndb; iz++)
		{
			ux[ix][iz] = alp * tvx[ix][iz];
			uz[ix][iz] = (2.0-alp) * tvz[ix][iz];
			a = cos(theta[ix][iz]);
			b = sin(theta[ix][iz]);
			tvx[ix][iz] = a*ux[ix][iz]-b*uz[ix][iz];
			tvz[ix][iz] = a*uz[ix][iz]+b*ux[ix][iz];
		}

	for (ix=ndb+1; ix<nx-ndb; ix++)
		for (iz=ndb+1; iz<nz-ndb; iz++)
		{
			ux[ix][iz] = tvx[ix-1][iz-1] * (-h) +
						 tvx[ix-1][iz  ] * (-h) +
						 tvx[ix  ][iz-1] * h +
						 tvx[ix  ][iz  ] * h;
			uz[ix][iz] = tvz[ix-1][iz-1] * (-h) +
						 tvz[ix-1][iz  ] * h +
						 tvz[ix  ][iz-1] * (-h) +
						 tvz[ix  ][iz  ] * h;
		}

	for (ix=ndb+1; ix<nx-ndb; ix++)
	{
		ux[ix][ndb] = ux[ix][ndb+1];
		uz[ix][ndb] = uz[ix][ndb+1];
	}
	for (iz=ndb+1; iz<nz-ndb; iz++)
	{
		ux[ndb][iz] = ux[ndb+1][iz];
		uz[ndb][iz] = uz[ndb+1][iz];
	}
	ux[ndb][ndb] = ux[ndb+1][ndb+1];
	uz[ndb][ndb] = uz[ndb+1][ndb+1];

	for (ix=ndb; ix<nx-ndb; ix++)
		for (iz=ndb; iz<nz-ndb; iz++)
			dtv[ix][iz] = ux[ix][iz] + uz[ix][iz];

	delete2D(ux, nx, nz);
	delete2D(uz, nx, nz);
}



void DTV(float **vv, float **dv, float **dtv, float **theta, float **b1, float **b2, float **d1, float **d2, float alp, float mu, float lambda, int nx, int nz, int ndb, float stab)
{
	int ix, iz;
	float **tvx, **tvz;
	float **tmp1, **tmp2;

	tvx = alloc2D(nx, nz);
	tvz = alloc2D(nx, nz);
	tmp1 = alloc2D(nx, nz);
	tmp2 = alloc2D(nx, nz);
	for (ix=0; ix<nx; ix++)
		for (iz=0; iz<nz;iz++)
		{
			tvx[ix][iz] = 0.0;
			tvz[ix][iz] = 0.0;
			dtv[ix][iz] = 0.0;
			tmp1[ix][iz] = 0.0;
			tmp2[ix][iz] = 0.0;
		}

	// get tvx tvz
	TVfilter(tvx, tvz, vv, theta, alp, nx, nz, ndb);


	for (ix=0; ix<nx; ix++)
		for (iz=0; iz<nz;iz++)
		{
			tvx[ix][iz] = d1[ix][iz] - tvx[ix][iz] - b1[ix][iz];
			tvz[ix][iz] = d2[ix][iz] - tvz[ix][iz] - b2[ix][iz];
		}


	// get dtv
	TVTfilter(tvx, tvz, dtv, theta, alp, nx, nz, ndb);

	for (ix=0; ix<nx; ix++)
		for (iz=0; iz<nz;iz++)
		{
			dtv[ix][iz] *= lambda;
			vv[ix][iz] += dtv[ix][iz] + dv[ix][iz];
		}


	TVfilter(tvx, tvz, vv, theta, alp, nx, nz, ndb);


	shrink(tvx, b1, d1, mu, nx, nz);
	shrink(tvz, b2, d2, mu, nx, nz);

	for (ix=0; ix<nx; ix++)
		for (iz=0; iz<nz;iz++)
		{
			tmp1[ix][iz] = tvx[ix][iz] + b1[ix][iz];
			tmp2[ix][iz] = tvz[ix][iz] + b2[ix][iz];
		}

	save_2d_float(tmp1, nx, nz, "tvx.bin");
	save_2d_float(tmp2, nx, nz, "tvz.bin");
	save_2d_float(d1, nx, nz, "d1.bin");
	save_2d_float(d2, nx, nz, "d2.bin");

	for (ix=0; ix<nx; ix++)
		for (iz=0; iz<nz;iz++)
		{
			b1[ix][iz] += tvx[ix][iz] - d1[ix][iz];
			b2[ix][iz] += tvz[ix][iz] - d2[ix][iz];
		}
	delete2D(tvx, nx, nz);
	delete2D(tvz, nx, nz);
}

float get_coef(float **dv, float **dtv, float mu, int nx, int nz, float stab)
{
	int ix, iz;
	float a,b;
	a = b = 0;
	for (ix=0; ix<nx; ix++)
		for (iz=0; iz<nz;iz++)
		{
			a += fabsf(dv[ix][iz]);
			b += fabsf(dtv[ix][iz]);
		}

	return mu * a / (b + 0.005 * b);

}

void shrink(float **tvx, float **b1, float **d1, float mu, int nx, int nz)
{
	int ix, iz;
	float tvxb_max, lambda;
	tvxb_max = 0;
	lambda = 1.0 / mu;

	for(ix=0; ix<nx; ix++)
		for(iz=0; iz<nz; iz++)
			d1[ix][iz] = max(fabsf(tvx[ix][iz]+b1[ix][iz])-lambda, 0)
				/ (max(fabsf(tvx[ix][iz]+b1[ix][iz])-lambda, 0)+lambda)
				* (tvx[ix][iz]+b1[ix][iz]);


}

