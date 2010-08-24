function [a b c dd] = christofel3dtest(cc,k)

c11=cc(1);
c12=cc(2);
c13=cc(3);
c22=cc(4);
c23=cc(5);
c33=cc(6);
c44=cc(7);
c55=cc(8);
c66=cc(9);

kx=k(1);
ky=k(2);
kz=k(3);

G(1,1)=c11*kx^2+c66*ky^2+c55*kz^2;
G(2,2)=c66*kx^2+c22*ky^2+c44*kz^2;
G(3,3)=c55*kx^2+c44*ky^2+c33*kz^2;

G(1,2)=(c12+c66)*kx*ky; G(2,1)=G(1,2);
G(1,3)=(c13+c55)*kx*kz; G(3,1)=G(1,3);
G(2,3)=(c23+c44)*ky*kz; G(3,2)=G(2,3);

[v,d]=eig(G);

%----------


dd=[d(1,1) d(2,2) d(3,3)];
[dd,ind]=sort(dd,'descend');

a=v(:,ind(1));
b=v(:,ind(2));
c=v(:,ind(3));

if (k*a< 0)
    a=-a;
end





