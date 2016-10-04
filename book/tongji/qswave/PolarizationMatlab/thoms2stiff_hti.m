function cc = thoms2stiff_hti(thoms)

%thoms[vp0,vs0,eps1,eps2,del1,del2,del3,gam1,gam2]
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

c33 = vp2;
c55 = vs2;
c66 = c55;
c11 = vp2*ep2;
c22 = c33;
c44 = c55/ga2;

c13 = sqrt((de2*c33-c55)*(c33-c55)) - c55;
c23 = c33 - 2*c44;
c12 = c13;

cc = [c11 c12 c13 c22 c23 c33 c44 c55 c66];