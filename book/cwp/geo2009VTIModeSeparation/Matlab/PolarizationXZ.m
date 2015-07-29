function [Nuxx,Nuzz]=PolarizationBox(N,vp,vs,epsilon,delta)

[kxx,kzz]=meshgrid(-N:N);
kk= sqrt(kxx.^2 + kzz.^2);

c33=vp*vp;
c55=vs*vs;
c11=(2*epsilon+1)*c33;
c13=sqrt(2*c33*(c33-c55)*delta+(c33-c55)^2)-c55;

[m,n] = size(kxx);

for i=1:m
    for j=1:n
       G(1,1)=c11*kxx(i,j)^2+c55*kzz(i,j)^2;
       G(1,2)=(c13+c55)*kxx(i,j)*kzz(i,j);
       G(2,1)=G(1,2);
       G(2,2)=c55*kxx(i,j)^2+c33*kzz(i,j)^2;
        
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


figure;
subplot(121)
h=imagesc([-1:.1:1],[-1:.1:1],uxx/N);colorbar;daspect([1 1 1]);xlabel('k_x');ylabel('k_z');

subplot(122)
h=imagesc([-1:.1:1],[-1:.1:1],uzz/N);colorbar;daspect([1 1 1]);xlabel('k_x');ylabel('k_z');


