function cc = thoms2stiff_ort(thoms)

%thoms[vp0,vs0,eps1,eps2,del1,del2,del3,gam1,gam2]
vp0 = thoms(1);
vs0 = thoms(2);
eps1 = thoms(3);
eps2 = thoms(4);
del1 = thoms(5);
del2 = thoms(6);
del3 = thoms(7);
gam1 = thoms(8);
gam2 = thoms(9);

vp2 = vp0*vp0;
vs2 = vs0*vs0;
ep1 = 1+2*eps1;
ep2 = 1+2*eps2;
de1 = 1+2*del1;
de2 = 1+2*del2;
de3 = 1+2*del3;
ga1 = 1+2*gam1;
ga2 = 1+2*gam2;

c11 = vp2*ep2;
c22 = vp2*ep1;
c33 = vp2;
c55 = vs2;
c66 = vs2*ga1;
c44 = c66/ga2;
c13 = sqrt((de2*c33-c55)*(c33-c55)) - c55;
c23 = sqrt((de1*c33-c44)*(c33-c44)) - c44;
c12 = sqrt((de3*c11-c66)*(c11-c66)) - c66;

cc = [c11 c12 c13 c22 c23 c33 c44 c55 c66];