/* Estimation of depth-delay of common-image gathers */
/*
 Copyright (C) 2016 University of Texas at Austin
 
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

int nx, nz, seed;
float idx, **p, eps2;

void sf_conv_lop1(bool adj, bool add, int nm, int nd, float *xx, float *yy)
/*< linear operator >*/
{
    int ix, iz;

    sf_adjnull(adj, add, nm, nd, xx, yy);

    if(adj){

        // first column
        xx[0] = xx[0] - yy[0]*(p[0][0]+idx);
        for(iz=1; iz<nz-2; iz++)
            xx[iz] = xx[iz] + yy[iz-1]*p[0][iz-1] - yy[iz]*(p[0][iz]+idx);
        xx[nz-2] = xx[nz-2] + yy[nz-3]*p[0][nz-3] - yy[nz-2]*(p[0][nz-2]+idx) - yy[nz-1]*p[0][nz-1];
        xx[nz-1] = xx[nz-1] + yy[nz-2]*p[0][nz-2] + yy[nz-1]*(p[0][nz-1]-idx);

        // second column
        for(ix=1; ix<nx-2; ix++)
            for(iz=1; iz<nz-2; iz++){
                xx[ix*nz+iz] = xx[ix*nz+iz] + yy[(ix-1)*nz+iz]*idx + yy[ix*nz+iz-1]*p[ix][iz-1] - yy[ix*nz+iz]*(p[ix][iz]+idx);
            }
        for(ix=1; ix<nx-2; ix++){
            xx[ix*nz] = xx[ix*nz] + yy[(ix-1)*nz]*idx - yy[ix*nz]*(p[ix][0]+idx);
            xx[ix*nz+nz-2] = xx[ix*nz+nz-2] + yy[(ix-1)*nz+nz-2]*idx + yy[ix*nz+nz-3]*p[ix][nz-3] - yy[ix*nz+nz-2]*(p[ix][nz-2]+idx) - yy[ix*nz+nz-1]*p[ix][nz-1];
            xx[ix*nz+nz-1] = xx[ix*nz+nz-1] + yy[(ix-1)*nz+nz-1]*idx + yy[ix*nz+nz-2]*p[ix][nz-2] + yy[ix*nz+nz-1]*(p[ix][nz-1]-idx);
        }

        // third column
        xx[(nx-2)*nz] = xx[(nx-2)*nz] + yy[(nx-3)*nz]*idx - yy[(nx-2)*nz]*(p[nx-2][0]+idx) - yy[(nx-1)*nz]*idx;
        for(iz=1; iz<nz-2; iz++)
            xx[(nx-2)*nz+iz] = xx[(nx-2)*nz+iz] + yy[(nx-3)*nz+iz]*idx + yy[(nx-2)*nz+iz-1]*p[nx-2][iz-1] - yy[(nx-2)*nz+iz]*(p[nx-2][iz]+idx) - yy[(nx-1)*nz+iz]*idx;
        xx[(nx-2)*nz+nz-2] = xx[(nx-2)*nz+nz-2] + yy[(nx-3)*nz+nz-2]*idx + yy[(nx-2)*nz+nz-3]*p[nx-2][nz-3] - yy[(nx-2)*nz+nz-2]*(p[nx-2][nz-2]+idx) - yy[(nx-2)*nz+nz-1]*p[nx-2][nz-1] - yy[(nx-1)*nz+nz-2]*idx;
        xx[(nx-2)*nz+nz-1] = xx[(nx-2)*nz+nz-1] + yy[(nx-3)*nz+nz-1]*idx + yy[(nx-2)*nz+nz-2]*p[nx-2][nz-2] + yy[(nx-2)*nz+nz-1]*(p[nx-2][nz-1]-idx) - yy[(nx-1)*nz+nz-1]*idx;

        // last column
        xx[(nx-1)*nz] = xx[(nx-1)*nz] + yy[(nx-2)*nz]*idx + yy[(nx-1)*nz]*(idx-p[nx-1][0]);
        for(iz=1; iz<nz-2; iz++)
            xx[(nx-1)*nz+iz] = xx[(nx-1)*nz+iz] + yy[(nx-2)*nz+iz]*idx + yy[(nx-1)*nz+iz-1]*p[nx-1][iz-1] + yy[(nx-1)*nz+iz]*(idx-p[nx-1][iz]);
        xx[(nx-1)*nz+nz-2] = xx[(nx-1)*nz+nz-2] + yy[(nx-2)*nz+nz-2]*idx + yy[(nx-1)*nz+nz-3]*p[nx-1][nz-3] + yy[(nx-1)*nz+nz-2]*(idx-p[nx-1][nz-2]) - yy[(nx-1)*nz+nz-1]*p[nx-1][nz-1];
        xx[(nx-1)*nz+nz-1] = xx[(nx-1)*nz+nz-1] + yy[(nx-2)*nz+nz-1]*idx + yy[(nx-1)*nz+nz-2]*p[nx-1][nz-2] + yy[(nx-1)*nz+nz-1]*(idx+p[nx-1][nz-1]);

    }else{

        // left-top part
        for(ix=0; ix<nx-1; ix++)
            for(iz=0; iz<nz-1; iz++){
                yy[iz+ix*nz] = yy[iz+ix*nz] - xx[iz+ix*nz]*(idx+p[ix][iz]) + xx[iz+1+ix*nz]*p[ix][iz] + xx[iz+(ix+1)*nz]*idx;
            }

        // bottom boundary
        for(ix=0; ix<nx-1; ix++)
            yy[nz-1+ix*nz] = yy[nz-1+ix*nz] - xx[nz-2+ix*nz]*p[ix][nz-1] + xx[nz-1+ix*nz]*(p[ix][nz-1]-idx) + xx[nz-1+(ix+1)*nz]*idx;

        // right boundary
        for(iz=0; iz<nz-1; iz++)
            yy[iz+(nx-1)*nz] = yy[iz+(nx-1)*nz] - xx[iz+(nx-2)*nz]*idx + xx[iz+(nx-1)*nz]*(idx-p[nx-1][iz]) + xx[iz+1+(nx-1)*nz]*p[nx-1][iz];

        // last element
        yy[nz-1+(nx-1)*nz] = yy[nz-1+(nx-1)*nz] - xx[nz-1+(nx-2)*nz]*idx - xx[nz-2+(nx-1)*nz]*p[nx-1][nz-1] + xx[nz-1+(nx-1)*nz]*(idx+p[nx-1][nz-1]);
    }
}

void sf_conv_lop2(bool adj, bool add, int nm, int nd, float *xx, float *yy)
/*< linear operator >*/
{
    int ix, iz;

    sf_adjnull(adj, add, nm, nd, xx, yy);

    if(adj){

        // first column
        xx[0] = xx[0] - yy[0]*(p[0][0]+idx);
        for(iz=1; iz<nz-2; iz++)
            xx[iz] = xx[iz] + yy[iz-1]*p[0][iz-1] - yy[iz]*(p[0][iz]+idx);
        xx[nz-2] = xx[nz-2] + yy[nz-3]*p[0][nz-3] - yy[nz-2]*(p[0][nz-2]+idx) - yy[nz-1]*p[0][nz-1];
        xx[nz-1] = xx[nz-1] + yy[nz-2]*p[0][nz-2] + yy[nz-1]*(p[0][nz-1]-idx);

        // second column
        for(ix=1; ix<nx-2; ix++)
            for(iz=1; iz<nz-2; iz++){
                xx[ix*nz+iz] = xx[ix*nz+iz] + yy[(ix-1)*nz+iz]*idx + yy[ix*nz+iz-1]*p[ix][iz-1] - yy[ix*nz+iz]*(p[ix][iz]+idx);
            }
        for(ix=1; ix<nx-2; ix++){
            xx[ix*nz] = xx[ix*nz] + yy[(ix-1)*nz]*idx - yy[ix*nz]*(p[ix][0]+idx);
            xx[ix*nz+nz-2] = xx[ix*nz+nz-2] + yy[(ix-1)*nz+nz-2]*idx + yy[ix*nz+nz-3]*p[ix][nz-3] - yy[ix*nz+nz-2]*(p[ix][nz-2]+idx) - yy[ix*nz+nz-1]*p[ix][nz-1];
            xx[ix*nz+nz-1] = xx[ix*nz+nz-1] + yy[(ix-1)*nz+nz-1]*idx + yy[ix*nz+nz-2]*p[ix][nz-2] + yy[ix*nz+nz-1]*(p[ix][nz-1]-idx);
        }

        // third column
        xx[(nx-2)*nz] = xx[(nx-2)*nz] + yy[(nx-3)*nz]*idx - yy[(nx-2)*nz]*(p[nx-2][0]+idx) - yy[(nx-1)*nz]*idx;
        for(iz=1; iz<nz-2; iz++)
            xx[(nx-2)*nz+iz] = xx[(nx-2)*nz+iz] + yy[(nx-3)*nz+iz]*idx + yy[(nx-2)*nz+iz-1]*p[nx-2][iz-1] - yy[(nx-2)*nz+iz]*(p[nx-2][iz]+idx) - yy[(nx-1)*nz+iz]*idx;
        xx[(nx-2)*nz+nz-2] = xx[(nx-2)*nz+nz-2] + yy[(nx-3)*nz+nz-2]*idx + yy[(nx-2)*nz+nz-3]*p[nx-2][nz-3] - yy[(nx-2)*nz+nz-2]*(p[nx-2][nz-2]+idx) - yy[(nx-2)*nz+nz-1]*p[nx-2][nz-1] - yy[(nx-1)*nz+nz-2]*idx;
        xx[(nx-2)*nz+nz-1] = xx[(nx-2)*nz+nz-1] + yy[(nx-3)*nz+nz-1]*idx + yy[(nx-2)*nz+nz-2]*p[nx-2][nz-2] + yy[(nx-2)*nz+nz-1]*(p[nx-2][nz-1]-idx) - yy[(nx-1)*nz+nz-1]*idx;

        // last column
        xx[(nx-1)*nz] = xx[(nx-1)*nz] + yy[(nx-2)*nz]*idx + yy[(nx-1)*nz]*(idx-p[nx-1][0]);
        for(iz=1; iz<nz-2; iz++)
            xx[(nx-1)*nz+iz] = xx[(nx-1)*nz+iz] + yy[(nx-2)*nz+iz]*idx + yy[(nx-1)*nz+iz-1]*p[nx-1][iz-1] + yy[(nx-1)*nz+iz]*(idx-p[nx-1][iz]);
        xx[(nx-1)*nz+nz-2] = xx[(nx-1)*nz+nz-2] + yy[(nx-2)*nz+nz-2]*idx + yy[(nx-1)*nz+nz-3]*p[nx-1][nz-3] + yy[(nx-1)*nz+nz-2]*(idx-p[nx-1][nz-2]) - yy[(nx-1)*nz+nz-1]*p[nx-1][nz-1];
        xx[(nx-1)*nz+nz-1] = xx[(nx-1)*nz+nz-1] + yy[(nx-2)*nz+nz-1]*idx + yy[(nx-1)*nz+nz-2]*p[nx-1][nz-2] + yy[(nx-1)*nz+nz-1]*(idx+p[nx-1][nz-1]);

        // padded part
        for(iz=0; iz<nz; iz++)
            xx[iz+(seed-1)*nz] = xx[iz+(seed-1)*nz] + eps2*yy[iz+nm+(seed-1)*nz];
    }else{

        // left-top part
        for(ix=0; ix<nx-1; ix++)
            for(iz=0; iz<nz-1; iz++){
                yy[iz+ix*nz] = yy[iz+ix*nz] - xx[iz+ix*nz]*(idx+p[ix][iz]) + xx[iz+1+ix*nz]*p[ix][iz] + xx[iz+(ix+1)*nz]*idx;
            }

        // bottom boundary
        for(ix=0; ix<nx-1; ix++)
            yy[nz-1+ix*nz] = yy[nz-1+ix*nz] - xx[nz-2+ix*nz]*p[ix][nz-1] + xx[nz-1+ix*nz]*(p[ix][nz-1]-idx) + xx[nz-1+(ix+1)*nz]*idx;

        // right boundary
        for(iz=0; iz<nz-1; iz++)
            yy[iz+(nx-1)*nz] = yy[iz+(nx-1)*nz] - xx[iz+(nx-2)*nz]*idx + xx[iz+(nx-1)*nz]*(idx-p[nx-1][iz]) + xx[iz+1+(nx-1)*nz]*p[nx-1][iz];

        // last element
        yy[nz-1+(nx-1)*nz] = yy[nz-1+(nx-1)*nz] - xx[nz-1+(nx-2)*nz]*idx - xx[nz-2+(nx-1)*nz]*p[nx-1][nz-1] + xx[nz-1+(nx-1)*nz]*(idx+p[nx-1][nz-1]);

        // padded part
        for(ix=nm; ix<nd; ix++)
            yy[ix] = yy[ix];
        for(iz=0; iz<nz; iz++)
            yy[iz+nm+(seed-1)*nz] = yy[iz+nm+(seed-1)*nz] + eps2*xx[iz+(seed-1)*nz]; 
    }
}

void sf_conv_lop3(bool adj, bool add, int nm, int nd, float *xx, float *yy)
/*< linear operator >*/
{
    float *temp;
    temp=sf_floatalloc(2*nm);

    if(adj){
        sf_conv_lop2(false,false,nm,2*nm,yy,temp);
        sf_conv_lop2(true,add,nm,2*nm,xx,temp);
    }else{
        sf_conv_lop2(false,false,nm,2*nm,xx,temp);
        sf_conv_lop2(true,add,nm,2*nm,yy,temp);
    }

    free(temp);
}

void sf_conv_lop4(bool adj, bool add, int nm, int nd, float *xx, float *yy)
/*< linear operator >*/
{
    int ix, iz;
	float *temp;
    temp=sf_floatalloc(nz);

    sf_adjnull(adj, add, nm, nd, xx, yy);

    if(adj){

		for(iz=0; iz<nz; iz++)
			temp[iz] = xx[iz+(seed-1)*nz];

        // first column
        xx[0] = xx[0] - yy[0]*(p[0][0]+idx);
        for(iz=1; iz<nz-2; iz++)
            xx[iz] = xx[iz] + yy[iz-1]*p[0][iz-1] - yy[iz]*(p[0][iz]+idx);
        xx[nz-2] = xx[nz-2] + yy[nz-3]*p[0][nz-3] - yy[nz-2]*(p[0][nz-2]+idx) - yy[nz-1]*p[0][nz-1];
        xx[nz-1] = xx[nz-1] + yy[nz-2]*p[0][nz-2] + yy[nz-1]*(p[0][nz-1]-idx);

        // second column
        for(ix=1; ix<nx-2; ix++)
            for(iz=1; iz<nz-2; iz++){
                xx[ix*nz+iz] = xx[ix*nz+iz] + yy[(ix-1)*nz+iz]*idx + yy[ix*nz+iz-1]*p[ix][iz-1] - yy[ix*nz+iz]*(p[ix][iz]+idx);
            }
        for(ix=1; ix<nx-2; ix++){
            xx[ix*nz] = xx[ix*nz] + yy[(ix-1)*nz]*idx - yy[ix*nz]*(p[ix][0]+idx);
            xx[ix*nz+nz-2] = xx[ix*nz+nz-2] + yy[(ix-1)*nz+nz-2]*idx + yy[ix*nz+nz-3]*p[ix][nz-3] - yy[ix*nz+nz-2]*(p[ix][nz-2]+idx) - yy[ix*nz+nz-1]*p[ix][nz-1];
            xx[ix*nz+nz-1] = xx[ix*nz+nz-1] + yy[(ix-1)*nz+nz-1]*idx + yy[ix*nz+nz-2]*p[ix][nz-2] + yy[ix*nz+nz-1]*(p[ix][nz-1]-idx);
        }

        // third column
        xx[(nx-2)*nz] = xx[(nx-2)*nz] + yy[(nx-3)*nz]*idx - yy[(nx-2)*nz]*(p[nx-2][0]+idx) - yy[(nx-1)*nz]*idx;
        for(iz=1; iz<nz-2; iz++)
            xx[(nx-2)*nz+iz] = xx[(nx-2)*nz+iz] + yy[(nx-3)*nz+iz]*idx + yy[(nx-2)*nz+iz-1]*p[nx-2][iz-1] - yy[(nx-2)*nz+iz]*(p[nx-2][iz]+idx) - yy[(nx-1)*nz+iz]*idx;
        xx[(nx-2)*nz+nz-2] = xx[(nx-2)*nz+nz-2] + yy[(nx-3)*nz+nz-2]*idx + yy[(nx-2)*nz+nz-3]*p[nx-2][nz-3] - yy[(nx-2)*nz+nz-2]*(p[nx-2][nz-2]+idx) - yy[(nx-2)*nz+nz-1]*p[nx-2][nz-1] - yy[(nx-1)*nz+nz-2]*idx;
        xx[(nx-2)*nz+nz-1] = xx[(nx-2)*nz+nz-1] + yy[(nx-3)*nz+nz-1]*idx + yy[(nx-2)*nz+nz-2]*p[nx-2][nz-2] + yy[(nx-2)*nz+nz-1]*(p[nx-2][nz-1]-idx) - yy[(nx-1)*nz+nz-1]*idx;

        // last column
        xx[(nx-1)*nz] = xx[(nx-1)*nz] + yy[(nx-2)*nz]*idx + yy[(nx-1)*nz]*(idx-p[nx-1][0]);
        for(iz=1; iz<nz-2; iz++)
            xx[(nx-1)*nz+iz] = xx[(nx-1)*nz+iz] + yy[(nx-2)*nz+iz]*idx + yy[(nx-1)*nz+iz-1]*p[nx-1][iz-1] + yy[(nx-1)*nz+iz]*(idx-p[nx-1][iz]);
        xx[(nx-1)*nz+nz-2] = xx[(nx-1)*nz+nz-2] + yy[(nx-2)*nz+nz-2]*idx + yy[(nx-1)*nz+nz-3]*p[nx-1][nz-3] + yy[(nx-1)*nz+nz-2]*(idx-p[nx-1][nz-2]) - yy[(nx-1)*nz+nz-1]*p[nx-1][nz-1];
        xx[(nx-1)*nz+nz-1] = xx[(nx-1)*nz+nz-1] + yy[(nx-2)*nz+nz-1]*idx + yy[(nx-1)*nz+nz-2]*p[nx-1][nz-2] + yy[(nx-1)*nz+nz-1]*(idx+p[nx-1][nz-1]);

		// projection operator
		for(iz=0; iz<nz; iz++)
			xx[iz+(seed-1)*nz] = temp[iz];

		free(temp);

    }else{

		// projection operator
		for(iz=0; iz<nz; iz++)
			xx[iz+(seed-1)*nz] = 0.;

        // left-top part
        for(ix=0; ix<nx-1; ix++)
            for(iz=0; iz<nz-1; iz++){
                yy[iz+ix*nz] = yy[iz+ix*nz] - xx[iz+ix*nz]*(idx+p[ix][iz]) + xx[iz+1+ix*nz]*p[ix][iz] + xx[iz+(ix+1)*nz]*idx;
            }

        // bottom boundary
        for(ix=0; ix<nx-1; ix++)
            yy[nz-1+ix*nz] = yy[nz-1+ix*nz] - xx[nz-2+ix*nz]*p[ix][nz-1] + xx[nz-1+ix*nz]*(p[ix][nz-1]-idx) + xx[nz-1+(ix+1)*nz]*idx;

        // right boundary
        for(iz=0; iz<nz-1; iz++)
            yy[iz+(nx-1)*nz] = yy[iz+(nx-1)*nz] - xx[iz+(nx-2)*nz]*idx + xx[iz+(nx-1)*nz]*(idx-p[nx-1][iz]) + xx[iz+1+(nx-1)*nz]*p[nx-1][iz];

        // last element
        yy[nz-1+(nx-1)*nz] = yy[nz-1+(nx-1)*nz] - xx[nz-1+(nx-2)*nz]*idx - xx[nz-2+(nx-1)*nz]*p[nx-1][nz-1] + xx[nz-1+(nx-1)*nz]*(idx+p[nx-1][nz-1]);
    }
}


int main(int argc, char *argv[])
{
    bool adj, inv, shape;
    int ix, iz, n, niter, rect1, rect2;
    float dx, dz, eps1;
    float *s, *np, *pp, *b;

    sf_file Fp, Fs, Fnp;

    sf_init(argc, argv);

    if(!sf_getbool("adj", &adj)) adj=false;
	/* if adj=y, adjoint operator */
    if(!sf_getbool("inv", &inv)) inv=false;
	/* if inv=y, perform inversion */
    if(!sf_getbool("shape", &shape)) shape=false;
	/* if shape=y, use projection method */

    Fp=sf_input("Fp");
    if(adj || inv){
        Fnp=sf_input("in");
        Fs=sf_output("out");
    }else{
        Fs=sf_input("in");
        Fnp=sf_output("out");
    }

    if(inv){
        if(!sf_getint("niter", &niter)) niter=100;
		/* number of iterations */

        if(!sf_getfloat("eps1", &eps1)) eps1=0.;
		/* shaping regularization parameter */
        if(!sf_getint("rect1", &rect1)) rect1=3;
		/* shaping smoothing parameter in 1st axis */
        if(!sf_getint("rect2", &rect2)) rect2=3;
		/* shaping smoothing parameter in 2nd axis */

        if(!sf_getint("seed", &seed)) seed=0;
		/* index of reference trace */
        if(!sf_getfloat("eps2", &eps2)) eps2=1e3;
		/* regularization parameter in model constraint */
    }

    if(!sf_histint(Fp, "n1", &nz)) sf_error("No n1= in Fp");
    if(!sf_histfloat(Fp, "d1", &dz)) sf_error("No d1= in Fp");
    if(!sf_histint(Fp, "n2", &nx)) sf_error("No n2= in Fp");
    if(!sf_histfloat(Fp, "d2", &dx)) sf_error("No d2= in Fp");

    n=nz*nx;
    p=sf_floatalloc2(nz, nx);
    s=sf_floatalloc(n);
    np=sf_floatalloc(n);
    if(seed>0 && !shape) b=sf_floatalloc(n*2);
    if(inv) pp=sf_floatalloc(n);

    sf_floatread(p[0], n, Fp);

    idx=-1./dx;
    for(ix=0; ix<nx; ix++){
        for(iz=0; iz<nz; iz++){
            p[ix][iz] = -p[ix][iz]/dz;
        }
    }

    if(adj || inv){
        if(seed>0 && !shape)
            sf_floatread(b, n, Fnp);
        else
            sf_floatread(np, n, Fnp);
        // pad data if seed>0 and shape=false
        if(seed>0 && !shape){
            for(ix=0; ix<n; ix++){
                b[ix+n]=0.;
            }
            sf_conv_lop2(true, false, n, 2*n, np, b);
        }
    }else{
        sf_floatread(s, n, Fs);
    }

    if(inv){
		/* shaping regularization with eps1 != 0 */
        sf_triangle2_init(rect1, rect2, nz, nx, 1);
        sf_conjgrad_init(n, n, n, n, eps1, 1.e-8, true, false);
        if(seed==0){ /* no model constraint */
            sf_conjgrad(NULL, sf_conv_lop1, sf_triangle2_lop, pp, s, np, niter);
        }else{ /* model constraint */
            if(shape) sf_conjgrad(NULL, sf_conv_lop4, sf_triangle2_lop, pp, s, np, niter);
			else sf_conjgrad(NULL, sf_conv_lop3, sf_triangle2_lop, pp, s, np, niter);
        }
        sf_conjgrad_close();
    }else{
		/* adj or not */
        sf_conv_lop1(adj, false, n, n, s, np);
    }

    if(adj || inv){
        sf_floatwrite(s, n, Fs);
    }else{
        sf_floatwrite(np, n, Fnp);
    }

    exit(0);
}
