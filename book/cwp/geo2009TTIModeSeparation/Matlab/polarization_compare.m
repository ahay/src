function [Nuxx,Nuzz]=polarization_compare(N,vp,vs,epsilon,delta,nu,color)

[kxx,kzz]=meshgrid(-N:N);
kk= sqrt(kxx.^2 + kzz.^2);

cc=zeros(6,6);


c33=vp*vp;
c55=vs*vs;
c11=(2*epsilon+1)*c33;
c13=sqrt(2*c33*(c33-c55)*delta+(c33-c55)^2)-c55;


cc(3,3)=c33;cc(5,5)=c55;cc(1,1)=c11;cc(1,3)=c13;
    
cc;
nu=nu*pi/180;
nu
m11=c11*cos(nu)^4+2*(c13+2*c55)*cos(nu)^2*sin(nu)^2+c33*sin(nu)^4;
m13=(c11+6*c13+c33-4*c55-(c11-2*c13+c33-4*c55)*cos(4*nu))/8;
m15=(c11-c33+(c11-2*c13+c33-4*c55)*cos(2*nu))*sin(2*nu)/4;
m33=c33*cos(nu)^4+2*(c13+2*c55)*cos(nu)^2*sin(nu)^2+c11*sin(nu)^4;
m35=-(-c11+c33+(c11-2*c13+c33-4*c55)*cos(2*nu))*sin(2*nu)/4;
m55=(c11-2*c13+c33+4*c55-(c11-2*c13+c33-4*c55)*cos(4*nu))/8;


cc(3,3)=m33;cc(5,5)=m55;cc(1,1)=m11;cc(1,3)=m13;cc(1,5)=m15;cc(3,5)=m35;
cc;



[m,n] = size(kxx);

for i=1:m
    for j=1:n
       G(1,1)=m11*kxx(i,j)^2+m55*kzz(i,j)^2+2*m15*kxx(i,j)*kzz(i,j);
       G(1,2)=m15*kxx(i,j)^2+m35*kzz(i,j)^2+(m13+m55)*kxx(i,j)*kzz(i,j);
       G(2,1)=G(1,2);
       G(2,2)=m55*kxx(i,j)^2+m33*kzz(i,j)^2+2*m35*kxx(i,j)*kzz(i,j);
        
       [v,d]= eig(G);
         
       sind=1;
       pind=2;
       if (d(1:1)>d(2,2))
           sind=2;
           pind=1;
       end

       uxx(i,j)=v(1,pind);
       uzz(i,j)=v(2,pind);
       if (uxx(i,j)*kxx(i,j) + uzz(i,j)*kzz(i,j)<=0)
            uxx(i,j)=-v(1,pind);
            uzz(i,j)=-v(2,pind);
       end
       
       vxx(i,j)=v(1,sind);
       vzz(i,j)=v(2,sind);
       if (vxx(i,j)*kzz(i,j) - vzz(i,j)*kxx(i,j)<=0)
            vxx(i,j)=-v(1,sind);
            vzz(i,j)=-v(2,sind);
       end
        
    end
end

uxx= uxx.*kk;
uzz= uzz.*kk;
vxx= vxx.*kk;
vzz= vzz.*kk;



h=quiver(kxx/N*3.14,kzz/N*3.14,uxx,uzz);
set(h, 'Color', color, 'LineWidth',1);

daspect([1 1 1]);
xlabel('k_x (radians) ');
ylabel('k_z (radians) ');
