function [a b c] = cseed(cc,k)

c11=cc(1);
c12=cc(2);
c13=cc(3);
c22=cc(4);
c23=cc(5);
c33=cc(6);
c44=cc(7);
c55=cc(8);
c66=cc(9);

kx=k(1); ky=k(2); kz=k(3);

G(1,1)=c11*kx^2+c66*ky^2+c55*kz^2;
G(2,2)=c66*kx^2+c22*ky^2+c44*kz^2;
G(3,3)=c55*kx^2+c44*ky^2+c33*kz^2;

G(1,2)=(c12+c66)*kx*ky; G(2,1)=G(1,2);
G(1,3)=(c13+c55)*kx*kz; G(3,1)=G(1,3);
G(2,3)=(c23+c44)*ky*kz; G(3,2)=G(2,3);

[v,d]=eig(G);

    if     (d(1,1)>d(2,2) & d(1,1)>d(3,3))
                inda=1; 
                if(d(2,2)>d(3,3))
                    indb=2; 
                    indc=3;
                else
                    indb=3;
                    indc=2;
                end
    elseif (d(2,2)>d(1,1) & d(2,2)>d(3,3))
                inda=2; 
                if(d(1,1)>d(3,3))
                    indb=1; 
                    indc=3;
                else
                    indb=3;
                    indc=1;
                end
    else
                inda=3; 
                if(d(1,1)>d(2,2))
                    indb=1; 
                    indc=2;
                else
                    indb=2;
                    indc=1;
                end
    end
        
    %largest eigenvalue->eigenvector
    ax=v(1,inda);
    ay=v(2,inda);
    az=v(3,inda);
    
    if (ax*kx + ay*ky + az*kz <= 0)
        ax=-ax;
        ay=-ay;
        az=-az;
    end
    
    bx=v(1,indb);
    by=v(2,indb);
    bz=v(3,indb);
    
    cx=v(1,indc);
    cy=v(2,indc);
    cz=v(3,indc);
    
a=[ax ay az];
b=[bx by bz];
c=[cx cy cz];
