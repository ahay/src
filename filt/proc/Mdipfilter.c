/*
#	pass=1	Pass dips in band
#	    =0	Reject dips in band
#
#	angle=0 Define band in terms of velocities
#	     =1 Define band in terms of angles
#
#	dim=2	2D dip-filter
#          =3	3D dip-filter
#
#	taper=1 linear taper
#            =2 cosine taper
#
#	Velocity/Dip band defined by:
#
#	    v2_____v3
#	     /     \
#	    /       \
#	 __/         \___
#	  v1         v4
#
# AUTHOR: James Rickett - Oct 97
*/
#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1,n2,n3, n1pad,n2pad,n3pad, nw, dim,nfr;
    float d1,d2,d3,kx0,dkx,kx,ky0,dky,ky,w0,dw,w,vel,degtan,v;
    float ang1, ang2, ang3, ang4, v1,v2,v3,v4,factor,taperfn;
    float ***data;
    float complex ***cdata;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_getint("dim",&dim)) dim=2;

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_histfloat(in,"d3",&d3)) d3=1.;

/* determine frequency sampling */
    n1pad = sf_npfar(n1);
    nw = n1pad/2+1;
    dw = 2.0*SF_PI/(n1pad*d1);
    w0 = 0.;

/* determine wavenumber sampling (for real to complex FFT) */
    n2pad = sf_npfa(n2);
    dkx = 2.0*SF_PI/(n2pad*d2);
    kx0 = -SF_PI/d2;

    if (dim==2) {
	nfr=n3; 
	n3=1; 
	n3pad=1;
    } else {
	nfr=1;
	n3pad = sf_npfa(n3);
	dky = 2.0*SF_PI/(n3pad*d3);
	ky0 = -SF_PI/d3;
    }

    data = sf_floatalloc3(n1,n2,n3);
    cdata = sf_complexalloc3(n1pad,n2pad,n3pad);

    exit(0);
}

#ifdef jkhgghdf

from par: real v1=0., v2=0.1, v3=99999., v4=999999., v=-1.
from par: real ang1=-50., ang2=-45., ang3=45., ang4=50.
from par: integer pass=1, angle=0, taper=1

if (v.lt.0.) v=(d2*(n2-1))/(d1*(n1-1))

if (angle==1) {
	v1=degtan(ang1)*v
	v2=degtan(ang2)*v
	v3=degtan(ang3)*v
	v4=degtan(ang4)*v
}
if ((v1>=v2)||(v2>=v3)||(v3>=v4)) call seperr('Need v1 < v2 < v3 < v4')

ifr=1	# Loop over panels
while (ifr <= nfr) { ifr=ifr+1

	# Read the data from stdin
	call sreed('in',data,n1*n2*n3*4)

	# Pad to power of 2
	do i3=1,n3pad {
	do i2=1,n2pad {
	do i1=1,n1pad {
		if ((i1 <= n1) && (i2 <= n2) && (i3 <= n3)) {
			cdata(i1,i2,i3)=cmplx(data(i1,i2,i3),0)
		} else {
			cdata(i1,i2,i3)=cmplx(0,0)
		}
	}}}

	# 3D FFT
	call ft3D1axis(0, 1.,n1pad,n2pad,n3pad,cdata)
	call ft3D2axis(0,-1.,n1pad,n2pad,n3pad,cdata)
	if (dim==3) call ft3D3axis(0,-1.,n1pad,n2pad,n3pad,cdata)

	dw =1/(n1pad*d1); w0 =-((n1pad-1)/2)*dw
	dkx=1/(n2pad*d2); kx0=-((n2pad-1)/2)*dkx
	dky=1/(n3pad*d3); ky0=-((n3pad-1)/2)*dky

	# Apply mute
	do i3=1, n3pad { ky=(i3-1)*dky+ky0
		if (dim==2) ky=0.
	do i2=1, n2pad { kx=(i2-1)*dkx+kx0
	do i1=1, n1pad { w =(i1-1)*dw +w0
		vel=w/(sqrt(kx**2+ky**2)+0.000001); factor=1.
		if ((vel>=v1) && (vel<=v2)) factor=1-taperfn(v1,v2,vel,taper)
		if ((vel>=v2) && (vel<=v3)) factor=0.
		if ((vel>=v3) && (vel<=v4)) factor=taperfn(v3,v4,vel,taper)
		if (pass == 1) factor=1.-factor
		cdata(i1,i2,i3)=cdata(i1,i2,i3)*factor
	}}}

	# Inverse 3D FFT
	call ft3D1axis(1, 1.,n1pad,n2pad,n3pad,cdata)
	call ft3D2axis(1,-1.,n1pad,n2pad,n3pad,cdata)
	if (dim==3) call ft3D3axis(1,-1.,n1pad,n2pad,n3pad,cdata)

	# Truncate to original size
	do i3=1,n3 { 
	do i2=1,n2 { 
	do i1=1,n1 { 
		data(i1,i2,i3)=real(cdata(i1,i2,i3))
	}}}

	# Write to stdout
	call srite('out',data,n1*n2*n3*4)
}
return;end

real function degtan(theta)
real theta,pi
pi=3.14159267
degtan=tan(pi*theta/180)
return;end

real function taperfn(start,end,value,type)
real start,end,value,pi
integer type
pi=3.14159267
taperfn=(value-start)/(end-start)
if (type.eq.2) {
	taperfn=sin(taperfn*pi/2)
}
return;end

#endif
