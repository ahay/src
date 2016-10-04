function cc = thoms2stiff_iso(thoms)

%thoms[vp0,eps,del,gam]
vp0 = thoms(1);
vs0 = thoms(2);

c11 = vp2;
c22 = c11;
c33 = c11;
c44 = vs2;
c55 = c44;
c66 = c44;

c12 = c11 - 2*c44;
c13 = c12;
c23 = c12;

cc = [c11 c12 c13 c22 c23 c33 c44 c55 c66];