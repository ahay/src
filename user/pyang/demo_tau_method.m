function demo_tau_method
%% Copyright (c) Xi'an Jiaotong University 2014
% % Originally written by Guowei Zhang, modified by Pengliang Yang
% %
% % Reference: Joakim O. Blanch, Johan O.A. Robertsson, and William W. Symmes,
% %  Modeling of a constant Q: methodology and algorithm for an efficient and 
% %  optimally inexpensive viscoelastic technique, Geophysics, vol 60 no.1
% % Note that there are some errors in the formulae appeared in the paper.

clc,clear

Qp=30;
Qs=20;
f0=25.0;

f=2:0.1:125;
w=2*pi*f;
wa=5*2*pi;
wb=125*2*pi;
tausig=1/2/f0/pi;

Qp0=28;
Qs0=18;

I0=fun1(wb,tausig)-fun1(wa,tausig);
I1=fun2(wb,tausig)-fun2(wa,tausig);

taup=1/Qp0*I0/I1;
taus=1/Qs0*I0/I1;
epsilop=(taup+1)*tausig;
epsilos=(taus+1)*tausig;

%% only use 1 SLS (standard linear solid )
Q1=findQ(w, tausig, taup);
plot(f,Q1)


Qp=30.0;
Qs=20.0;

f1=5.0;
f2=25.0;
f3=125.0;

f=2:0.1:125;
w=2*pi*f;
wa=5*2*pi;
wb=125*2*pi;

tau1=1/2/pi/f2/10;
tau2=1/2/f2/pi;
tau3=1/2/f2/pi*10;


Qp0=28.0;

I01=fun1(wb,tau1)-fun1(wa,tau1);
I02=fun1(wb,tau2)-fun1(wa,tau2);
I03=fun1(wb,tau3)-fun1(wa,tau3);

I11=fun2(wb,tau1)-fun2(wa,tau1);
I12=fun2(wb,tau2)-fun2(wa,tau2);
I13=fun2(wb,tau3)-fun2(wa,tau3);

I221=fun3(wb,tau1,tau2)-fun3(wa,tau1,tau2);
I231=fun3(wb,tau1,tau3)-fun3(wa,tau1,tau3);
I232=fun3(wb,tau2,tau3)-fun3(wa,tau2,tau3);

s0=I01+I02+I03;
s1=I11+I12+I13;
s3=I221+I231+I232;
tau=3/Qp0*s0/(s1+2*s3);


taueps1=(tau+1)*tau1;
taueps2=(tau+1)*tau2;
taueps3=(tau+1)*tau3;
s1=(1+w.^2*taueps1*tau1)./(1+w.^2*tau1.^2);
s2=(1+w.^2*taueps2*tau2)./(1+w.^2*tau2.^2);
s3=(1+w.^2*taueps3*tau3)./(1+w.^2*tau3.^2);
t1=w*(taueps1-tau1)./(1+w.^2*tau1.^2);
t2=w*(taueps2-tau2)./(1+w.^2*tau2.^2);
t3=w*(taueps3-tau3)./(1+w.^2*tau3.^2);
%% use 3 SLS for eaxct value 
Q=(s1+s2+s3)./(t1+t2+t3);
hold on
semilogx(f,Q)


function Q=findQ(w, tausig, tau)
Q=(1+w.^2*tausig.^2*(1+tau))./(w*tausig*tau);

function y=fun1(w, tau)
a=w*tau;
y=log(1+a.^2)./(2*tau);

function y=fun2(w, tau)
a=w*tau;
y=(atan(a)-a./(1+a.^2))./(2*tau);

function y=fun3(w, taul, tauk)
a=atan(w*taul)./taul-atan(w*tauk)/tauk;
b=taul.*tauk./(tauk.^2-taul.^2);
y=a*b;
