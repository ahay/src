function cc = thoms2stiff_vti(thoms)

%thoms[vp0,eps,del,gam]
vp0 = thoms(1);
vs0 = thoms(2);
eps = thoms(3);
del = thoms(4);
gam = thoms(5);

vp2 = vp0*vp0;
vs2 = vs0*vs0;
ep2 = 1+2*eps;
de2 = 1+2*del;
ga2 = 1+2*gam;

c11 = vp2*ep2;
c22 = c11;
c33 = vp2;
c44 = vs2;
c55 = c44;
c66 = vs2*ga2;

c13 = sqrt((de2*c33-c55)*(c33-c55)) - c44;
c23 = c13;
c12 = c11 - 2*c66;

cc = [c11 c12 c13 c22 c23 c33 c44 c55 c66];