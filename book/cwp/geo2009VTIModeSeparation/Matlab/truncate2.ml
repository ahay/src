N=10;
vp=2;
vs=1;
epsilon=0.25;
delta=-0.29;

boxn=5;
PolarizationBox(N,vp,vs,epsilon,delta,boxn);

%boxn=15;
%PolariNzationBox(N,vp,vs,epsilon,delta,boxn);

%boxn=25;
%PolarizationBox(N,vp,vs,epsilon,delta,boxn);

%boxn=N/2;
%PolarizationBox(N,vp,vs,epsilon,delta,boxn);


%boxn=N;
%PolarizationBox(N,vp,vs,epsilon,delta,boxn);


print -depsc junk_ml.eps;

quit