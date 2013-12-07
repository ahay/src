% by Mehdi Eftekhari Far
%last updated: 12/3/2010
function [weight]=cg_avo(n_unknown,cons,bkbm,tolerancer,itmax,I)


M=n_unknown;
w=zeros(1,M);
B=cons';
g=-B;
p=g;
A=bkbm;
for i=1:I
    alpha=sum(g.*g)/sum(p'.*(A*p'));
    z=g.*g;
    w=w-alpha.*p;
    g=g-alpha.*(A*p')';
    d=g.*g;
    beta=d/z;
    p=g+beta*p;
end


g=A*w'-cons;
g=g';
p=g;
for i=1:itmax
    alpha1=sum(g.*g)/sum(p'.*(A*p'));
    z=g.*g;
    w=w-alpha1.*p;
    g=g-alpha1.*(A*p')';
    d=g.*g;
    beta=d/z;
    p=g+beta*p;
    r=sum(g.*g)/sum(-B.*(-B));
%     if r<=tolerancer, break, end
end

weight=w';
end